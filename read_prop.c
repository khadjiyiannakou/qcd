/* disconnected strange contribution to nucleon
 *
 * reads forward propagators
 * and creates the disconnected loop
 *
 * Kyriakos Hadjiyiannakou 2011
 ****************************************/
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>
#include "projectors.h"

// definitions
#define NX 32
#define NY 32
#define NZ 32 
#define NT 64
#define Nmu 4
#define Nic 3
#define N_noise_max 100
//functions

qcd_complex_16 conju(qcd_complex_16 a){
  qcd_complex_16 res;
  res.re=a.re;
  res.im=-a.im;
  return res;
}

qcd_complex_16 mul(qcd_complex_16 a,qcd_complex_16 b){
qcd_complex_16 res;
res.re=a.re*b.re-a.im*b.im;
res.im=a.re*b.im+a.im*b.re;
return res;
}

qcd_complex_16 add(qcd_complex_16 a,qcd_complex_16 b){
qcd_complex_16 res;
res.re=a.re+b.re;
res.im=a.im+b.im;
return res;
}

// important macros
#define ind(ix,iy,iz,it,mu,ic) (ix)*NY*NZ*NT*Nmu*Nic + (iy)*NZ*NT*Nmu*Nic + (iz)*NT*Nmu*Nic + (it)*Nmu*Nic + (mu)*Nic + ic 
//////
int main(int argc,char* argv[]){

   qcd_uint_2 mu,nu,ku,lu,c1,c2,c3,c1p,c2p,c3p;// various loop variables
   qcd_uint_2 id1,id2,id3,cc1,cc2,al,be;
   qcd_uint_4 i,j,k, v,lx,ly,lz,ip1,im1,v3; 
   qcd_int_4 x,y,z;
   qcd_uint_2 ic1,ic2,ic3;                    //
   qcd_uint_4 x_src[4];                       // source and sink coordinates
   qcd_uint_4 t_sink, t_start, t_stop, t,lt;
   qcd_real_8 tmp;                            // general purpuse
   FILE *fp_momlist;
   FILE *fp_prop;                           // output file
   FILE *fp_sollist;
   FILE *fp_loop_p;

   int params_len;                            // needed to read inputfiles
   char *params;                              // needed to read inputfiles
   char momlist_name[qcd_MAX_STRING_LENGTH];    // name of momenta-list file
   char gauge_name[qcd_MAX_STRING_LENGTH];      // name of gauge-configuration file
   char prop_name[qcd_MAX_STRING_LENGTH];     // name of output file proton 2pt function
   char loop_p_name[qcd_MAX_STRING_LENGTH];
   char param_name[qcd_MAX_STRING_LENGTH];      // name of parameter file  
//   char sollist_name[qcd_MAX_STRING_LENGTH];    // name of momenta-list file
   char noise_sol[N_noise_max][qcd_MAX_STRING_LENGTH];
   char list_sol[qcd_MAX_STRING_LENGTH];
   char uprop_name[qcd_MAX_STRING_LENGTH];      // file names of up and down quark propagators
   char dprop_name[qcd_MAX_STRING_LENGTH];      

   qcd_geometry geo;                            // geometry structure
   qcd_propagator uprop;                        // propagator
   qcd_propagator dprop;                        // propagator
   qcd_vector vec;                             // needed when smearing
   qcd_vector Phi;
   qcd_gaugeField u;                            // gauge field 
   qcd_gaugeField uAPE;                         // APE smeared gaugeField
   qcd_gaugeField *u_ptr, *uAPE_ptr, *utmp_ptr;

   qcd_uint_4 nsmear, nsmearAPE;       // gaussian and APE smearing: n
   qcd_real_8 alpha, alphaAPE;         // gaussian and APE smearing: alpha

   qcd_real_8 theta[4] = {M_PI,0.0,0.0,0.0};    // antiperiodic b.c. in time
   qcd_real_8 mu_value;
   qcd_uint_2 L[4];
   qcd_uint_2 P[4];
   qcd_complex_16 phase_factor;         
   qcd_complex_16 z1, z2;                       // temp variables
   qcd_complex_16 C, C2;   
   qcd_complex_16 corr, corr2;
   qcd_real_8 plaq;
   qcd_int_4 ctr, ctr2;				//count the none zero values of gamma matrices
   qcd_int_2 gamma_ind[16][2];		//store the indices for gamma matrices
   qcd_complex_16 gamma_val[16];		//store the values
   qcd_int_4 gamma_index,n_noise;
   qcd_complex_16 *block;                       // to store the block (2pt function before FT)

   qcd_int_4 (*mom)[3];                         // momenta-list
   qcd_uint_2 Nmom;
   int myid,numprocs, namelen;    
   char processor_name[MPI_MAX_PROCESSOR_NAME];
   FILE *pt_up,*pt_down;
             
             
   //set up MPI
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);         // num. of processes taking part in the calculation
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);             // each process gets its ID
   MPI_Get_processor_name(processor_name,&namelen); // 
               
   
   //////////////////// READ INPUT FILE /////////////////////////////////////////////
   if(myid ==0){
     pt_up = fopen("prop_up.dat","w");
     pt_down = fopen("prop_down.dat","w");
   }
   if(argc!=2)
   {
      if(myid==0) fprintf(stderr,"No input file specified\n");
      exit(EXIT_FAILURE);
   }

   strcpy(param_name,argv[1]);
   if(myid==0)
   {
      i=0;
      printf("opening input file %s\n",param_name);
      params=qcd_getParams(param_name,&params_len);
      if(params==NULL)
      {
         i=1;
      }
   }
   MPI_Bcast(&i,1,MPI_INT, 0, MPI_COMM_WORLD);
   if(i==1) exit(EXIT_FAILURE);
   MPI_Bcast(&params_len, 1, MPI_INT, 0, MPI_COMM_WORLD);
   if(myid!=0) params = (char*) malloc(params_len*sizeof(char));
   MPI_Bcast(params, params_len, MPI_CHAR, 0, MPI_COMM_WORLD);
   
   sscanf(qcd_getParam("<processors_txyz>",params,params_len),"%hd %hd %hd %hd",&P[0], &P[1], &P[2], &P[3]);
   sscanf(qcd_getParam("<lattice_txyz>",params,params_len),"%hd %hd %hd %hd",&L[0], &L[1], &L[2], &L[3]);
   if(P[0] != 1)
   {
      if(myid==0) fprintf(stderr,"Error! Number of processors in t-direction must be one.\n");
      exit(EXIT_FAILURE);
   }
   if(qcd_initGeometry(&geo,L,P, theta, myid, numprocs)) exit(EXIT_FAILURE);   
   
   if(myid==0) printf(" Local lattice: %i x %i x %i x %i\n",geo.lL[0],geo.lL[1],geo.lL[2],geo.lL[3]);
      
                   
   strcpy(uprop_name,qcd_getParam("<propagator_u>",params,params_len));
   if(myid==0) printf("Got propagator file name: %s\n",uprop_name);
   strcpy(dprop_name,qcd_getParam("<propagator_d>",params,params_len));
   if(myid==0) printf("Got propagator file name: %s\n",dprop_name);

 
           
//   strcpy(sollist_name,qcd_getParam("<solution_list>",params,params_len));
//   if(myid==0) printf("Got solution-list file name: %s\n",sollist_name);
    
   
   free(params);

   j = 0;
   j += qcd_initPropagator(&uprop, &geo);
   j += qcd_initPropagator(&dprop, &geo);
   //   j += qcd_initVector(&vec, &geo);
   // j += qcd_initVector(&Phi, &geo);

   

            
   MPI_Allreduce(&j, &k, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   if(k>0)
   {
      if(myid==0) printf("not enough memory\n");
      exit(EXIT_FAILURE);
   }
   if(myid==0) printf("memory for propagators and gauge-field allocated\n");

   // load propagators                                                                                                                                 

   if(qcd_getPropagator(uprop_name,qcd_PROP_LIME, &uprop)) exit(EXIT_FAILURE); // read u propagators from binary in ram                              
   if(myid==0) printf("up propagator loaded\n");
   if(qcd_getPropagator(dprop_name,qcd_PROP_LIME, &dprop)) exit(EXIT_FAILURE); // read d propagators from binary in ram                              
   if(myid==0) printf("down propagator loaded\n");

   if(myid == 0){
   for(nu=0; nu<4; nu++)
     for(c2=0; c2<3; c2++)
       for(x=0; x< geo.lL[1]; x++)
	 for(y=0; y< geo.lL[2]; y++)
	   for(z=0; z< geo.lL[3]; z++)
	     for(t=0; t< geo.lL[0]; t++)
	       for(mu=0; mu<4; mu++)
		 for(c1=0; c1<3; c1++)
		   {
		     fprintf(pt_up,"%+e %+e\n",uprop.D[qcd_LEXIC(t,x,y,z,geo.lL)][mu][nu][c1][c2].re,uprop.D[qcd_LEXIC(t,x,y,z,geo.lL)][mu][nu][c1][c2].im);
		     fprintf(pt_down,"%+e %+e\n",dprop.D[qcd_LEXIC(t,x,y,z,geo.lL)][mu][nu][c1][c2].re,dprop.D[qcd_LEXIC(t,x,y,z,geo.lL)][mu][nu][c1][c2].im);
		   }
   }
   
   //##############################################################################
   MPI_Finalize();
}
