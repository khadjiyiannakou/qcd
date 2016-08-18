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
#define N_noise 1
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

   int params_len;                            // needed to read inputfiles
   char *params;                              // needed to read inputfiles

   char gauge_name[qcd_MAX_STRING_LENGTH];      // name of gauge-configuration file
   char prop_name[qcd_MAX_STRING_LENGTH];     // name of output file proton 2pt function
   
   char param_name[qcd_MAX_STRING_LENGTH];      // name of parameter file  
//   char sollist_name[qcd_MAX_STRING_LENGTH];    // name of momenta-list file
   char noise_sol[N_noise][qcd_MAX_STRING_LENGTH];
   char list_sol[qcd_MAX_STRING_LENGTH];
//   char uprop_name[qcd_MAX_STRING_LENGTH];      // file names of up and down quark propagators
//   char dprop_name[qcd_MAX_STRING_LENGTH];      

   qcd_geometry geo;                            // geometry structure
//   qcd_propagator uprop;                        // propagator
//   qcd_propagator dprop;                        // propagator
   qcd_vector vec;                             // needed when smearing
   qcd_vector Phi;
   qcd_gaugeField u;                            // gauge field 
   qcd_gaugeField uAPE;                         // APE smeared gaugeField
   qcd_gaugeField *u_ptr, *uAPE_ptr, *utmp_ptr;

   qcd_uint_4 nsmear, nsmearAPE;       // gaussian and APE smearing: n
   qcd_real_8 alpha, alphaAPE;         // gaussian and APE smearing: alpha

   qcd_real_8 theta[4] = {M_PI,0.0,0.0,0.0};    // antiperiodic b.c. in time
   qcd_uint_2 L[4];
   qcd_uint_2 P[4];
   qcd_complex_16 phase_factor;         
   qcd_complex_16 z1, z2;                       // temp variables
   qcd_complex_16 C, C2;   
   qcd_complex_16 corr, corr2;
   qcd_real_8 plaq;
   qcd_int_4 ctr, ctr2;
   qcd_int_2 cg5cg5_ind[16*16][4];
   qcd_complex_16 cg5cg5_val[16*16];
   
   qcd_complex_16 *block;                       // to store the block (2pt function before FT)

   qcd_int_4 (*mom)[3];                         // momenta-list

   int myid,numprocs, namelen;    
   char processor_name[MPI_MAX_PROCESSOR_NAME];
   				 
//   complex *Phi_kr = (complex*)malloc(NX*NY*NZ*NT*Nmu*Nic*sizeof(complex));				 
   qcd_complex_16 loop[NT];
   qcd_int_4 gamma5[4];
   gamma5[0]=1;gamma5[1]=1;gamma5[2]=-1;gamma5[3]=-1;
   // initialize loop values to zero
   for(int i=0;i<NT;i++){
   loop[i].re=0.;
   loop[i].im=0.;
   }           
             
             
   //set up MPI
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);         // num. of processes taking part in the calculation
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);             // each process gets its ID
   MPI_Get_processor_name(processor_name,&namelen); // 
               
   
   //////////////////// READ INPUT FILE /////////////////////////////////////////////
      
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
      
   sscanf(qcd_getParam("<t>",params,params_len),"%d %d",&t_start, &t_stop);
   if(myid==0) printf("Got sink time slices: %d ... %d\n",t_start,t_stop);
                     	  
   sscanf(qcd_getParam("<source_pos_txyz>",params,params_len),"%d %d %d %d",&x_src[0],&x_src[1],&x_src[2],&x_src[3]);
   if(myid==0) printf("Got source coords: %d %d %d %d\n",x_src[0],x_src[1],x_src[2],x_src[3]);
     
   strcpy(list_sol,qcd_getParam("<list_sol>",params,params_len));
   if(myid==0) printf("Got list_sol file name: %s\n",list_sol);

   
   strcpy(gauge_name,qcd_getParam("<cfg_name>",params,params_len));
   if(myid==0) printf("Got conf name: %s\n",gauge_name);
          
          
          
   strcpy(prop_name,qcd_getParam("<prop>",params,params_len));
   if(myid==0) printf("Got output file name: %s\n",prop_name);
           
//   strcpy(sollist_name,qcd_getParam("<solution_list>",params,params_len));
//   if(myid==0) printf("Got solution-list file name: %s\n",sollist_name);
    
   sscanf(qcd_getParam("<alpha_gauss>",params,params_len),"%lf",&alpha);
   if(myid==0) printf(" Got alpha_gauss: %lf\n",alpha);
   sscanf(qcd_getParam("<nsmear_gauss>",params,params_len),"%d",&nsmear);
   if(myid==0) printf(" Got nsmear_gauss: %d\n",nsmear);
   sscanf(qcd_getParam("<alpha_APE>",params,params_len),"%lf",&alphaAPE);
   if(myid==0) printf(" Got alpha_APE: %lf\n",alphaAPE);
   sscanf(qcd_getParam("<nsmear_APE>",params,params_len),"%d",&nsmearAPE);
   if(myid==0) printf(" Got nsmear_APE: %d\n",nsmearAPE); 
   free(params);



         
   //#####################################################################   
   // allocate memory
   // load gauge-field and APE-smear it
   j=0;
   j += qcd_initGaugeField(&u, &geo);
   j += qcd_initGaugeField(&uAPE,&geo);
   MPI_Allreduce(&j, &k, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   if(k>0)
   {
      if(myid==0) fprintf(stderr,"Error, not enough memory\n");
      exit(EXIT_FAILURE);
   }
      
   if(qcd_getGaugeField(gauge_name,qcd_GF_LIME,&u)) 
   {
      if(myid==0) fprintf(stderr,"Error reading gauge field\n");
      exit(EXIT_FAILURE);
   }
   if(myid==0) printf("gauge-field loaded\n");   
   plaq = qcd_calculatePlaquette(&u);
   if(myid==0) printf("plaquette = %e\n",plaq);
   qcd_communicateGaugePM(&u);
   qcd_waitall(&geo);
   u_ptr = &u;
   uAPE_ptr = &uAPE;   
   for(i=0; i<nsmearAPE; i++)
   {
      qcd_apeSmear3d(uAPE_ptr, u_ptr, alphaAPE);
      utmp_ptr=u_ptr; u_ptr=uAPE_ptr; uAPE_ptr=utmp_ptr;   
   }
   utmp_ptr=u_ptr; u_ptr=uAPE_ptr; uAPE_ptr=utmp_ptr; //reverse the last swap. Also needed when nsmearAPE=0
   qcd_destroyGaugeField(u_ptr);
   uAPE = *uAPE_ptr;
  
   j = 0;
//   j += qcd_initPropagator(&uprop, &geo);
//   j += qcd_initPropagator(&dprop, &geo);
   j += qcd_initVector(&vec, &geo);
   j += qcd_initVector(&Phi, &geo);
   
   block = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
            
   MPI_Allreduce(&j, &k, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   if(k>0)
   {
      if(myid==0) printf("not enough memory\n");
      exit(EXIT_FAILURE);
   }
   if(myid==0) printf("memory for propagators and gauge-field allocated\n");
         
   
   //##############################################################################
   
   // load solution vector
   fp_sollist=fopen(list_sol,"r");
   for(int i=0;i<N_noise;i++){
   fscanf(fp_sollist,"%s",&noise_sol[i]);
   }
// the heart of the program //
   for(int ns=0;ns<N_noise;ns++){
    if(qcd_getVector(noise_sol[ns],qcd_PROP_LIME,0,0,&Phi)) exit(EXIT_FAILURE);
     for(int it=0;it<Phi.geo->lL[0];it++){
      for(int ix=0;ix<Phi.geo->lL[1];ix++)
       for(int iy=0;iy<Phi.geo->lL[2];iy++)
        for(int iz=0;iz<Phi.geo->lL[3];iz++)
	 for(int mu=0;mu<Nmu;mu++)
	  for(int ic=0;ic<Nic;ic++){
           loop[it].re += gamma5[mu]*(mul(conju(Phi.D[qcd_LEXIC(it,ix,iy,iz,geo.lL)][mu][ic]),Phi.D[qcd_LEXIC(it,ix,iy,iz,geo.lL)][mu][ic]).re);
	  // if(it==0)printf("%+e %+e \n",loop[it].re,gamma5[mu]*(mul(conju(Phi.D[qcd_LEXIC(it,ix,iy,iz,geo.lL)][mu][ic]),Phi.D[qcd_LEXIC(it,ix,iy,iz,geo.lL)][mu][ic]).re));
	  }
     }
   }

   for(int i=0;i<NT;i++){
//    loop[i].im=(loop[i].re*2*0.0035)/N_noise;
    printf("%+e\n",loop[i].re);
   }

//	printf("%s\n",noise_sol[1]);
//   if(qcd_getVector(noise_sol,qcd_PROP_LIME,0,0,&Phi)) exit(EXIT_FAILURE);
/*
   fp_prop=fopen(prop_name,"w");
   if(myid==0){   
   for(x=0; x< Phi.geo->lL[1]; x++)
   for(y=0; y< Phi.geo->lL[2]; y++)
   for(z=0; z< Phi.geo->lL[3]; z++)
   for(t=0; t< Phi.geo->lL[0]; t++)
   for(mu=0; mu<4; mu++)
   for(c1=0; c1<3; c1++)
   {
      fprintf(fp_prop,"%+e %+e\n",Phi.D[qcd_LEXIC(t,x,y,z,geo.lL)][mu][c1].re,Phi.D[qcd_LEXIC(t,x,y,z,geo.lL)][mu][c1].im);
   }

   }
*/
   // load propagators
/*   
   if(qcd_getPropagator(uprop_name,qcd_PROP_LIME, &uprop)) exit(EXIT_FAILURE);
   if(myid==0) printf("up propagator loaded\n");
   if(qcd_getPropagator(dprop_name,qcd_PROP_LIME, &dprop)) exit(EXIT_FAILURE);
   if(myid==0) printf("down propagator loaded\n");   
*/
   if(myid==0) printf("cleaning up...\n");   
   free(block);
   free(mom);
   qcd_destroyVector(&Phi);
   qcd_destroyGeometry(&geo);
   MPI_Finalize();
}
