/* source.c
 *
 * creates smeared sources
 * Gaussian smearing with 
 * APE-smeared gauge fields
 *
 * Tomasz Korzec 2009
 ****************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include <qcd.h>
 
enum {
  POINT_SOURCE,
  NOISE_SOURCE
} src_type;
 
int main(int argc,char* argv[])
{
   FILE*  pfile;
   char*  params;
   char   gauge_name[qcd_MAX_STRING_LENGTH];
   char   param_name[qcd_MAX_STRING_LENGTH];
   char   out_name[qcd_MAX_STRING_LENGTH];
   char   vec_names[1024][qcd_MAX_STRING_LENGTH];
   char   src_type_s[qcd_MAX_STRING_LENGTH];
   qcd_int_4   x_src[4],lx_src[4],i,j,t,isource,nsmear,nsmearAPE;
   qcd_uint_2   mu,nu,col,c1,c2,s;
   qcd_real_8   alpha,alphaAPE,plaq;
   qcd_uint_4   params_len;   


   qcd_geometry geo;
   qcd_gaugeField u, uAPE;
   qcd_gaugeField *u_ptr, *uAPE_ptr, *utmp_ptr;
   qcd_propagator source;
   qcd_vector vec;
   qcd_uint_2 P[4];
   qcd_uint_2 L[4];
   qcd_real_8 theta[4]={M_PI,0.,0.,0.}; // boundary conditions

   qcd_uint_4    ismear,nsources;
   qcd_uint_4    ape_ismear;
   
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
   if(qcd_initGeometry(&geo,L,P, theta, myid, numprocs)) exit(EXIT_FAILURE);
   
   if(myid==0) printf(" Local lattice: %i x %i x %i x %i\n",geo.lL[0],geo.lL[1],geo.lL[2],geo.lL[3]);  
   

   sscanf(qcd_getParam("<alpha_gauss>",params,params_len),"%lf",&alpha);
   if(myid==0) printf(" Got alpha_gauss: %lf\n",alpha);
   sscanf(qcd_getParam("<nsmear_gauss>",params,params_len),"%d",&nsmear);
   if(myid==0) printf(" Got nsmear_gauss: %d\n",nsmear);
   sscanf(qcd_getParam("<alpha_APE>",params,params_len),"%lf",&alphaAPE);
   if(myid==0) printf(" Got alpha_APE: %lf\n",alphaAPE);
   sscanf(qcd_getParam("<nsmear_APE>",params,params_len),"%d",&nsmearAPE);
   if(myid==0) printf(" Got nsmear_APE: %d\n",nsmearAPE);   
   strcpy(gauge_name,qcd_getParam("<cfg_name>",params,params_len));
   if(myid==0) printf(" Got conf name: %s\n",gauge_name);
   strcpy(out_name,qcd_getParam("<source>",params,params_len));
   if(myid==0) printf(" Got source name %s\n",out_name);


   free(params);      
   ///////////////////////////////////////////////////////////////////////////////////////////////////
  


       qcd_initGaugeField(&u,&geo);
       qcd_initGaugeField(&uAPE,&geo);

       if(qcd_getGaugeField(gauge_name,qcd_GF_LIME,&u))
	 {
	   fprintf(stderr,"process %i: Error reading gauge field!\n",myid);
	   exit(EXIT_FAILURE);
	 }
       
       if(myid==0) printf("gauge-field loaded\n");
       plaq = qcd_calculatePlaquette(&u);
       if(myid==0) printf("plaquette = %e\n",plaq);
       
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
    
       if(myid==0) printf("gauge-field APE-smeared\n");
       plaq = qcd_calculatePlaquette(&uAPE);
       if(myid==0) printf("plaquette = %e\n",plaq); 


       /*
       qcd_initVector(&vec, &geo);
       vec.D[0][0][0].re = 1.;
       for(i=0; i<nsmear; i++)
         {
           if(qcd_gaussIteration3d(&vec,&uAPE,alpha,0))
             {
               fprintf(stderr,"process %i: Error while smearing!\n",geo.myid);
               exit(EXIT_FAILURE);
             }
         }

       
       if(myid==0){
	 FILE *ptr_test;
	 ptr_test = fopen("kale.dat","w");
	 // for(int iv = 0 ; iv < 24*24*24*48 ; iv++)
	 for(int z = 0 ; z < 24 ; z++)
	   for(int y = 0 ; y < 24 ; y++)
	     for(int x =0 ; x < 24 ; x++){
	       int iv = qcd_LEXIC(0,x,y,z,geo.lL);
	       for(int mu = 0 ; mu < 4 ; mu++)
		 for(int ic = 0 ; ic < 3 ; ic++)
		   fprintf(ptr_test,"%d %d %d \t %e %e\n",x,y,z,vec.D[iv][mu][ic].re,vec.D[iv][mu][ic].im);
	     }
       }
       //qcd_initPropagator(&source,&geo);
       
       */

   ////////////////////////////////////// CLEAN UP AND EXIT ///////////////////////////////////////////
   //qcd_destroyPropagator(&source);
   qcd_destroyGeometry(&geo);
   MPI_Finalize();
   return(EXIT_SUCCESS);
}//end main
