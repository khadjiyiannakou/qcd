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
   
   sscanf(qcd_getParam("<mu_value>",params,params_len),"%lf",&mu_value);
   if(myid==0) printf("Got mu_value: %lf\n",mu_value);
   sscanf(qcd_getParam("<gamma_insertion>",params,params_len),"%d",&gamma_index);
   if(myid==0)printf("Got gamma matrix index: %d\n",gamma_index);
   sscanf(qcd_getParam("<number_noise>",params,params_len),"%d",&n_noise);
   if(myid==0)printf("Got number of noise vector: %d\n",n_noise);
   strcpy(loop_p_name,qcd_getParam("<loop_name>",params,params_len));
   if(myid==0) printf("Got output loop file name: %s\n",loop_p_name);
   strcpy(momlist_name,qcd_getParam("<momenta_list>",params,params_len));
   if(myid==0) printf("Got momenta-list file name: %s\n",momlist_name);

   free(params);

   //loop definition
   qcd_complex_16 loop,loop_test;

   //calculate non zero matrix elements
   ctr=0;
   for(mu=0;mu<Nmu;mu++)
   for(nu=0;nu<Nmu;nu++){
    C = qcd_GAMMAG5[gamma_index][mu][nu];
    if(qcd_NORM(C) > 1e-3){
     gamma_val[ctr].re = C.re;
     gamma_val[ctr].im = C.im;
     gamma_ind[ctr][0] = mu;
     gamma_ind[ctr][1] = nu;
     ctr++;
    }
   }
/*
if(myid==0){
for(int i=0;i<ctr;i++)printf("%+e %+e \n",gamma_val[i].re,gamma_val[i].im);
}
*/
//open output file to write in
    if(myid==0)
   {
      fp_loop_p = fopen(loop_p_name,"w");
      if(fp_loop_p==NULL)
      {
         printf("failed to open %s for writing\n",loop_p_name);
         k=1;
      }
      fp_momlist = fopen(momlist_name,"r");
      if(fp_momlist==NULL)
      {
        printf("failed to open %s for reading\n",momlist_name);
         k=1;
      }
   }
   MPI_Bcast(&k,1,MPI_INT, 0, MPI_COMM_WORLD);
      if(k==1) exit(EXIT_FAILURE);

   //load momenta-list   
   if(myid==0) fscanf(fp_momlist,"%i\n",&Nmom);
   MPI_Bcast(&Nmom,1,MPI_INT, 0, MPI_COMM_WORLD);
   if(myid==0) printf("will read %i momenta combinations\n",Nmom);

   mom = malloc(Nmom*3*sizeof(qcd_int_4));

   if(myid==0)
   {
      for(j=0; j<Nmom; j++)
      {
         fscanf(fp_momlist,"%i %i %i\n",&(mom[j][0]),&(mom[j][1]),&(mom[j][2]));
         //printf("got combination %i %i %i\n",mom[j][0],mom[j][1],mom[j][2]);  
      }
      fclose(fp_momlist);
   }
   MPI_Bcast(&(mom[0][0]),Nmom*3,MPI_INT,0, MPI_COMM_WORLD);
   if(myid==0) printf("momenta list read and broadcasted\n");
   
   // another loop with index for momenta
	qcd_complex_16 loop_final[Nmom][L[0]];
	for(int i=0;i<Nmom;i++)
	for(int j=0;j<L[0];j++){
		loop_final[i][j].re = 0;
		loop_final[i][j].im = 0;
	}
   qcd_complex_16 loop2[n_noise*Nmom*L[0]];
   // initialize loop values to zero
   for(int i=0;i<n_noise*Nmom*L[0];i++){
   loop2[i].re=0;
   loop2[i].im=0;
   }

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
   for(int i=0;i<n_noise;i++){
   fscanf(fp_sollist,"%s",&noise_sol[i]);
   }
// the heart of the program //

//calculate the traces//
   for(int ns=0;ns<n_noise;ns++){
    if(qcd_getVector(noise_sol[ns],qcd_PROP_LIME,0,0,&Phi)) exit(EXIT_FAILURE);
    for(int it=0;it<Phi.geo->lL[0];it++){
     for(int v3=0; v3 < geo.lV3; v3++) block[v3]= (qcd_complex_16) {0,0};
     for(int ctr2=0;ctr2<ctr;ctr2++){     
      for(int ic=0;ic<Nic;ic++){
       for(int lx=0;lx<Phi.geo->lL[1];lx++){
        for(int ly=0;ly<Phi.geo->lL[2];ly++){
         for(int lz=0;lz<Phi.geo->lL[3];lz++){
           v3 = qcd_LEXIC0(lx,ly,lz,geo.lL);
	   mu = gamma_ind[ctr2][0];
           nu = gamma_ind[ctr2][1];
           block[v3] =qcd_CADD(block[v3],qcd_CMUL(gamma_val[ctr2],qcd_CMUL(qcd_CONJ(Phi.D[qcd_LEXIC(it,lx,ly,lz,geo.lL)][mu][ic]),Phi.D[qcd_LEXIC(it,lx,ly,lz,geo.lL)][nu][ic])));
	 //  if(myid==0 && it==0)printf("%+e %+e %+e %+e %+e %+e \n",qcd_CMUL(qcd_CONJ(Phi.D[qcd_LEXIC(it,lx,ly,lz,geo.lL)][mu][ic]),Phi.D[qcd_LEXIC(it,lx,ly,lz,geo.lL)][nu][ic]).re,qcd_CMUL(qcd_CONJ(Phi.D[qcd_LEXIC(it,lx,ly,lz,geo.lL)][mu][ic]),Phi.D[qcd_LEXIC(it,lx,ly,lz,geo.lL)][nu][ic]).im,Phi.D[qcd_LEXIC(it,lx,ly,lz,geo.lL)][nu][ic].re,Phi.D[qcd_LEXIC(it,lx,ly,lz,geo.lL)][nu][ic].im,qcd_CONJ(Phi.D[qcd_LEXIC(it,lx,ly,lz,geo.lL)][nu][ic]).re,qcd_CONJ(Phi.D[qcd_LEXIC(it,lx,ly,lz,geo.lL)][nu][ic]).im);
	  } // close inner loop
         } // close y
        } // close x
       } //close color 
      } //close spin
// fourier transform inside every time-slice
      for(int j=0;j<Nmom;j++){
       loop = (qcd_complex_16) {0,0};
            for(lx=0; lx<geo.lL[1]; lx++)
            for(ly=0; ly<geo.lL[2]; ly++)
            for(lz=0; lz<geo.lL[3]; lz++)
            {
               v3 = qcd_LEXIC0(lx,ly,lz,geo.lL);
               x=lx+geo.Pos[1]*geo.lL[1];
               y=ly+geo.Pos[2]*geo.lL[2];
               z=lz+geo.Pos[3]*geo.lL[3];
               tmp = (((double) mom[j][0]*x)/geo.L[1] + ((double) mom[j][1]*y)/geo.L[2] + ((double) mom[j][2]*z)/geo.L[3])*2*M_PI;
               C2=(qcd_complex_16) {cos(tmp), -sin(tmp)}; //TABULATE FOR LARGE SPEEDUP!!!
               loop=qcd_CADD(loop, qcd_CMUL(block[v3],C2));
            }

MPI_Reduce(&(loop.re),&(loop2[(ns)*Nmom*L[0]+(j)*L[0]+it].re),2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      }// close j for all momenta
     } // close time
   } // close ns Number of noise vector

if(myid==0){
 for(int imom=0;imom<Nmom;imom++)
 for(int it=0;it< L[0];it++)
 for(int ns=0;ns<n_noise;ns++){
  loop_final[imom][it] = qcd_CADD(loop_final[imom][it],loop2[(ns)*Nmom*L[0]+(imom)*L[0]+it]);
 }
}


if(myid==0){
   for(int imom=0;imom<Nmom;imom++)
   for(int i=0;i<L[0];i++){
//    loop[i].im=(loop[i].re*2*0.0035)/N_noise;
    fprintf(fp_loop_p,"%+e %+e\n",loop_final[imom][i].re/n_noise,loop_final[imom][i].im/n_noise);
   }
}
/*
if(myid==0){
 //  for(int imom=0;imom<Nmom;imom++)
   for(int i=0;i<L[0];i++){
//    loop[i].im=(loop[i].re*2*0.0035)/N_noise;
    fprintf(fp_loop_p,"%+e %+e\n",loop2[i].re/n_noise,loop2[i].im/n_noise);
   }
}
*/

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
