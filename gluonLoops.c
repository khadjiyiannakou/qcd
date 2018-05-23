/* test the stout smearing
**
 *  * APE-smeared gauge fields
 *
 * Tomasz Korzec 2009
 ****************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include <qcd.h>
 
int main(int argc,char* argv[])
{
   FILE*  pfile;
   char*  params = NULL;
   char   gauge_name[qcd_MAX_STRING_LENGTH];
   char   param_name[qcd_MAX_STRING_LENGTH];
   char   output_name[qcd_MAX_STRING_LENGTH];

   qcd_int_4   i,nsmearstout;
   qcd_real_8   plaq,rho;
   int   params_len;   

   qcd_geometry geo;
   qcd_gaugeField u,u_out,utmp;
   qcd_gaugeField *u_ptr, *uAPE_ptr, *utmp_ptr;
   qcd_vector vec;
   qcd_uint_2 P[4];
   qcd_uint_2 L[4];
   qcd_real_8 theta[4]={M_PI,0.,0.,0.}; // boundary conditions
   
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
   
   sscanf(qcd_getParam("<rho_stout>",params,params_len),"%lf",&rho);
   if(myid==0) printf(" Got rho_stout: %lf\n",rho);
   sscanf(qcd_getParam("<nstout>",params,params_len),"%d",&nsmearstout);
   if(myid==0) printf(" Got nstout: %d\n",nsmearstout);   
   strcpy(gauge_name,qcd_getParam("<cfg_name>",params,params_len));
   if(myid==0) printf(" Got conf name: %s\n",gauge_name);

   strcpy(output_name,qcd_getParam("<output_name>",params,params_len));
   if(myid==0) printf(" Got output name: %s\n",output_name);

   if(myid==0)
   {
      pfile = fopen(output_name,"w");   
      if(pfile == NULL)
      {
         printf("failed to open %s for writing\n",output_name);
         i=1;
      }
   }
   MPI_Bcast(&i,1,MPI_INT, 0, MPI_COMM_WORLD);
   if(i==1) exit(EXIT_FAILURE);
  


   free(params);      
   ///////////////////////////////////////////////////////////////////////////////////////////////////

   if(P[0] != 1)
   {
      if(myid==0) fprintf(stderr,"Error! Number of processors in t-direction must be one.\n");
      exit(EXIT_FAILURE);
   }
   
   qcd_initGaugeField(&u,&geo);
   qcd_initGaugeField(&u_out,&geo);
   qcd_initGaugeField(&utmp,&geo);

   if(qcd_getGaugeField(gauge_name,qcd_GF_LIME,&u))
     {
       fprintf(stderr,"process %i: Error reading gauge field!\n",myid);
       exit(EXIT_FAILURE);
     }
       
   if(myid==0) printf("gauge-field loaded\n");

   qcd_communicateGaugePM(&u);
   qcd_waitall(&geo);
   plaq = qcd_calculatePlaquette(&u);
   if(myid==0) printf("plaquette not smeared= %e\n",plaq);

   qcd_copyGaugeField(&utmp,&u);
   /** test >=0 stout smearing **/
   for(i=0; i<nsmearstout; i++)
     {
       qcd_stoutsmearing(&u_out, &utmp, rho, 4);
       plaq=qcd_calculatePlaquette(&u_out);
       if(myid==0) printf("iter=%d, plaquette stout=%e\n",i,plaq);
       qcd_copyGaugeField(&utmp,&u_out);
     }   
   
   double *gLoops = malloc(geo.lL[0]*sizeof(double));
   qcd_calculateGluonLoops(&u_out,gLoops);
   for(int it =0 ; it < geo.lL[0] ; it++)
     if(myid==0)
       fprintf(pfile,"%d %+e\n",it,gLoops[it]);

   free(gLoops);
   qcd_destroyGaugeField(&utmp); 
   qcd_destroyGaugeField(&u_out); 
   qcd_destroyGaugeField(&u); 
   
   ////////////////////////////////////// CLEAN UP AND EXIT ///////////////////////////////////////////
   qcd_destroyGeometry(&geo);
   MPI_Finalize();
   return(EXIT_SUCCESS);
}//end main
