/*
  Kyriakos Hadjiyiannakou 2018
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
   FILE *fp_momlist;
   char*  params = NULL;
   char   gauge_name[qcd_MAX_STRING_LENGTH];
   char   param_name[qcd_MAX_STRING_LENGTH];
   char   output_name[qcd_MAX_STRING_LENGTH];
   char momlist_name[qcd_MAX_STRING_LENGTH];
   char sumT[qcd_MAX_STRING_LENGTH]; // if this string is yes then it sums the time slices
   char doSg[qcd_MAX_STRING_LENGTH]; // if this string is yes we will compute also Sg for trace anomaly
   char doTraceless[qcd_MAX_STRING_LENGTH]; // In case we want to make the traceless the field after 

   qcd_int_4   i,k,j,nsmearstout;
   qcd_int_4 printStep;
   qcd_int_4 Nprint;
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

   qcd_int_4 (*mom)[3];                         // momenta-list
   int Nmom;
   

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

   sscanf(qcd_getParam("<printStep>",params,params_len),"%d",&printStep);
   if(myid==0) printf(" Got printStep: %d\n",printStep);

   if( (nsmearstout%printStep) != 0){
     fprintf(stderr,"Error! printStep should divide the number of smearing steps");
     exit(EXIT_FAILURE);
   }
   Nprint = nsmearstout/printStep + 1; // +1 because we want also the unsmeared case

   strcpy(gauge_name,qcd_getParam("<cfg_name>",params,params_len));
   if(myid==0) printf(" Got conf name: %s\n",gauge_name);

   strcpy(output_name,qcd_getParam("<output_name>",params,params_len));
   if(myid==0) printf(" Got output name: %s\n",output_name);

   strcpy(momlist_name,qcd_getParam("<momenta_list>",params,params_len));
   if(myid==0) printf("Got momenta-list file name: %s\n",momlist_name);

   strcpy(sumT,qcd_getParam("<sumT>",params,params_len));
   if(myid==0) printf("Do you want to perform the sum over time-slices?: %s\n",sumT);

   strcpy(doSg,qcd_getParam("<doSg>",params,params_len));
   if(myid==0) printf("Do you want to compute also trace anomaly?: %s\n",doSg);

   strcpy(doTraceless,qcd_getParam("<doTraceless>",params,params_len));
   if(myid==0) printf("Do you want to use traceless in the gluon fields?: %s\n",doTraceless);




   if(myid==0)
   {
      pfile = fopen(output_name,"w");   
      if(pfile == NULL)
	{
	  printf("failed to open %s for writing\n",output_name);
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

   free(params);      

   //load momenta-list   
   if(myid==0) fscanf(fp_momlist,"%d\n",&Nmom);
   MPI_Bcast(&Nmom,1,MPI_INT, 0, MPI_COMM_WORLD);
   if(myid==0) printf("will read %d momenta combinations\n",Nmom);

   mom = malloc(Nmom*3*sizeof(qcd_int_4));

   if(myid==0)
   {
      for(j=0; j<Nmom; j++)
      {
         fscanf(fp_momlist,"%i %i %i\n",&(mom[j][0]),&(mom[j][1]),&(mom[j][2]));
      }
      fclose(fp_momlist);   
   }
   MPI_Bcast(&(mom[0][0]),Nmom*3,MPI_INT,0, MPI_COMM_WORLD);
   if(myid==0) printf("momenta list read and broadcasted\n");   

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

   if( strcmp(doTraceless,"yes") == 0  || strcmp(doTraceless,"YES") == 0){
     traceless_gluon_Field(&utmp,&u);
     qcd_copyGaugeField(&u,&utmp);
   }

   qcd_communicateGaugePM(&u);
   qcd_waitall(&geo);
   plaq = qcd_calculatePlaquette(&u);
   if(myid==0) printf("plaquette not smeared= %e\n",plaq);
     
   qcd_complex_16 *gLoops = malloc(Nprint*geo.lL[0]*Nmom*2*sizeof(double));
   qcd_complex_16 *gLoops_Sg;
   if( strcmp(doSg,"yes") == 0  || strcmp(doSg,"YES") == 0)
     gLoops_Sg = malloc(Nprint*geo.lL[0]*Nmom*2*sizeof(double));
   else
     gLoops_Sg = NULL;

   if( strcmp(doSg,"yes") == 0  || strcmp(doSg,"YES") == 0)
     qcd_calculateGluonLoops_T00g_Sg(&u,gLoops,gLoops_Sg,mom,Nmom);
   else
     qcd_calculateGluonLoops(&u,gLoops,mom,Nmom); // gluon loops without any smearing

   qcd_copyGaugeField(&utmp,&u);

   int ii=0;
   for(i=0; i<nsmearstout; i++)
     {
       qcd_stoutsmearing(&u_out, &utmp, rho, 4);
       plaq=qcd_calculatePlaquette(&u_out);
       if(myid==0) printf("iter=%d, plaquette stout=%e\n",i,plaq);
       qcd_copyGaugeField(&utmp,&u_out);
       if((i+1)%printStep == 0){
	 ii+=1;
	 if( strcmp(doSg,"yes") == 0  || strcmp(doSg,"YES") == 0)
	   qcd_calculateGluonLoops_T00g_Sg(&u_out,&(gLoops[ii*geo.lL[0]*Nmom]),&(gLoops_Sg[ii*geo.lL[0]*Nmom]),mom,Nmom);
	 else
	   qcd_calculateGluonLoops(&u_out,&(gLoops[ii*geo.lL[0]*Nmom]),mom,Nmom);
       }
     }   

   qcd_complex_16  sumLoop;
   qcd_complex_16 sumLoopSg;

   if( strcmp(sumT,"yes") == 0  || strcmp(sumT,"YES") == 0){
     if(myid == 0) printf("You have enable sum of the loops on time");
     for(int st = 0 ; st < Nprint ; st++)
       for(int imom = 0 ; imom < Nmom ; imom++){
	 sumLoop.re=0.; sumLoop.im=0.;
	 if(myid==0){
	   if( strcmp(doSg,"yes") == 0  || strcmp(doSg,"YES") == 0){
	     for(int it =0 ; it < geo.lL[0] ; it++){
	       sumLoop.re += gLoops[st*geo.lL[0]*Nmom + it*Nmom + imom].re;
	       sumLoop.im += gLoops[st*geo.lL[0]*Nmom + it*Nmom + imom].im;
	       sumLoopSg.re += gLoops_Sg[st*geo.lL[0]*Nmom + it*Nmom + imom].re;
	       sumLoopSg.im += gLoops_Sg[st*geo.lL[0]*Nmom + it*Nmom + imom].im;
	     } 
	     fprintf(pfile,"%d %+d %+d %+d %+e %+e %+e %+e\n",st*printStep,mom[imom][0],mom[imom][1],mom[imom][2],sumLoop.re,sumLoop.im,
		     sumLoopSg.re,sumLoopSg.im);
	   }
	   else{
	     for(int it =0 ; it < geo.lL[0] ; it++){
	       sumLoop.re += gLoops[st*geo.lL[0]*Nmom + it*Nmom + imom].re;
	       sumLoop.im += gLoops[st*geo.lL[0]*Nmom + it*Nmom + imom].im;
	     } 
	     fprintf(pfile,"%d %+d %+d %+d %+e %+e\n",st*printStep,mom[imom][0],mom[imom][1],mom[imom][2],sumLoop.re,sumLoop.im);
	   }
	 }
       }
   }
   else{
     for(int st = 0 ; st < Nprint ; st++)
       for(int it =0 ; it < geo.lL[0] ; it++)
	 for(int imom = 0 ; imom < Nmom ; imom++)
	   if(myid==0){
	     if( strcmp(doSg,"yes") == 0  || strcmp(doSg,"YES") == 0)
	       fprintf(pfile,"%d %d %+d %+d %+d %+e %+e %+e %+e\n",st*printStep,it,mom[imom][0],mom[imom][1],mom[imom][2],
		       gLoops[st*geo.lL[0]*Nmom + it*Nmom + imom].re, gLoops[st*geo.lL[0]*Nmom + it*Nmom + imom].im,
		       gLoops_Sg[st*geo.lL[0]*Nmom + it*Nmom + imom].re, gLoops_Sg[st*geo.lL[0]*Nmom + it*Nmom + imom].im);	       
	     else
	       fprintf(pfile,"%d %d %+d %+d %+d %+e %+e\n",st*printStep,it,mom[imom][0],mom[imom][1],mom[imom][2],
		       gLoops[st*geo.lL[0]*Nmom + it*Nmom + imom].re, gLoops[st*geo.lL[0]*Nmom + it*Nmom + imom].im);
	   }
   }

   free(gLoops);
   free(gLoops_Sg);
   free(mom);
   qcd_destroyGaugeField(&utmp); 
   qcd_destroyGaugeField(&u_out); 
   qcd_destroyGaugeField(&u); 
   
   ////////////////////////////////////// CLEAN UP AND EXIT ///////////////////////////////////////////
   qcd_destroyGeometry(&geo);
   MPI_Finalize();
   return(EXIT_SUCCESS);
}//end main
