// Kyriakos Hadjiyiannakou 
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>


int main(int argc,char* argv[])
{

   qcd_uint_2 mu,nu,rho,ku,lu,c1,c2;          // various loop variables
   qcd_uint_4 i,j,k,lt,x,y,z,t,ip1,im1; 
   qcd_uint_2 id1,id2,id3,id4;
   qcd_uint_2 ic1,ic2,ic3,ic4;                    
   qcd_real_8 tmp;                            // general purpuse
   qcd_uint_2 xx[4];

   FILE *fp_vfun;
   FILE *fp_gloops;
   FILE *fp_props;

   int params_len;               // needed to read inputfiles
   char *params;                 // needed to read inputfiles

   char vfun_name[qcd_MAX_STRING_LENGTH];     
   char param_name[qcd_MAX_STRING_LENGTH];      // name of parameter file
   char pprop_name[qcd_MAX_STRING_LENGTH];

   qcd_geometry geo;                            // geometry structure
   char gloops_name[qcd_MAX_STRING_LENGTH];

   qcd_real_8 theta[4] = {M_PI,0.0,0.0,0.0};    // antiperiodic b.c. in time
   qcd_uint_2 L[4];
   qcd_uint_2 P[4];
   qcd_complex_16 uprop[4][4][3][3];
   qcd_complex_16 dprop[4][4][3][3];
   qcd_real_8 *gloops;
     
   qcd_int_2  pn[4];                            //integer momentum
   qcd_real_8 px[4];                             //= n*2*pi/L

   int myid,numprocs, namelen;    
   char processor_name[MPI_MAX_PROCESSOR_NAME];     

   int nStout;
		      
   //////////////////////////////////////////////////////////////////////////////////////          
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

   strcpy(gloops_name,qcd_getParam("<gloops_name>",params,params_len));
   if(myid==0) printf("Got gluon loop file name: %s\n",gloops_name);
                                 
   strcpy(vfun_name,qcd_getParam("<vfd_name>",params,params_len));
   if(myid==0) printf("Got output file name for vertex function: %s\n",vfun_name);

   strcpy(pprop_name,qcd_getParam("<pprop_name>",params,params_len));
   if(myid==0) printf("Got input file name for up/down propagator: %s\n",pprop_name);
            
   sscanf(qcd_getParam("<nstout>",params,params_len),"%hd",&nStout);
   if(myid==0) printf("Got number of stouts including without smearing: %i \n",nStout);


   gloops = malloc(nStout*sizeof(qcd_real_8));
   int d_dummy;
   double lf_dummy;

   free(params);
   

   if(P[0]*P[1]*P[2]*P[3] != 1){
     fprintf(stderr,"Error only one process is needed for this code (very cheap calculation)\n");
     exit(EXIT_FAILURE);
   }
   //#####################################################################   

   if(myid==0)
   {
      j=0;
      if( (fp_vfun=fopen(vfun_name,"w"))==NULL) j++;
      if( (fp_gloops=fopen(gloops_name,"r"))==NULL) j++;
      if( (fp_props=fopen(pprop_name,"r"))==NULL ) j++;
      if(j>0) fprintf(stderr,"Error while opening input/output files!\n");
   }
   MPI_Bcast(&j,1,MPI_INT, 0, MPI_COMM_WORLD);
   if(j>0) exit(EXIT_FAILURE);
   
      

   ////////////////////////////////////////////////////////////////////////////////////////////////////////

   for(mu = 0 ; mu < 4 ; mu++)
     for(nu = 0 ; nu < 4 ; nu++)
       for(c1 = 0 ; c1 < 3 ; c1++)
	 for(c2 = 0 ; c2 < 3 ; c2++){
	   fscanf(fp_props,"%d %d %d %d %lf %lf %lf %lf",&d_dummy,&d_dummy,&d_dummy,&d_dummy,&(uprop[mu][nu][c1][c2].re), &(uprop[mu][nu][c1][c2].im),
		  &(dprop[mu][nu][c1][c2].re), &(dprop[mu][nu][c1][c2].im));
	 }


     // read the gluon loop
     for(int i = 0 ; i <= nStout ; i++) fscanf(fp_gloops,"%d %d %d %d %lf %lf",&d_dummy,&d_dummy,&d_dummy,&d_dummy,&(gloops[i]),&lf_dummy);

     for(int i = 0 ; i <= nStout ; i++)
       for(mu = 0 ; mu < 4 ; mu++)
	 for(nu = 0 ; nu < 4 ; nu++)
	   for(c1 = 0 ; c1 < 3 ; c1++)
	     for(c2 = 0 ; c2 < 3 ; c2++)
	       fprintf(fp_vfun, "%d %d %d %d %d \t %+e %+e \t %+e %+e\n",i, mu, nu, c1, c2, uprop[mu][nu][c1][c2].re * gloops[i], uprop[mu][nu][c1][c2].im * gloops[i], dprop[mu][nu][c1][c2].re * gloops[i], dprop[mu][nu][c1][c2].im * gloops[i]);
     
   free(gloops);
   qcd_destroyGeometry(&geo);
   MPI_Finalize();
}//end main 

