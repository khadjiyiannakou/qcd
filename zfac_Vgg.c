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
   char gauge_name[qcd_MAX_STRING_LENGTH];
   char param_name[qcd_MAX_STRING_LENGTH];      // name of parameter file
   char loops_name[qcd_MAX_STRING_LENGTH];
   char pprop_name[qcd_MAX_STRING_LENGTH];

   qcd_geometry geo;                            // geometry structure
   qcd_gaugeField u;                            // gauge field 
   char gloops_name[qcd_MAX_STRING_LENGTH];

   qcd_real_8 theta[4] = {M_PI,0.0,0.0,0.0};    // antiperiodic b.c. in time
   qcd_uint_2 L[4];
   qcd_uint_2 P[4];
   qcd_complex_16 phase_factor_x, phase_factor_y;
   qcd_complex_16 ux_lFT[4][3][3]; //local
   qcd_complex_16 uy_lFT[4][3][3];
   qcd_complex_16 ux_gFT[4][3][3]; //global
   qcd_complex_16 uy_gFT[4][3][3];
   qcd_complex_16 gprop[4][4][3][3];
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
   if(myid==0) printf("Got gluon loop file name: %s\n",loops_name);
                                 
   strcpy(vfun_name,qcd_getParam("<vfd_name>",params,params_len));
   if(myid==0) printf("Got output file name for vertex function: %s\n",vfun_name);

   strcpy(pprop_name,qcd_getParam("<pprop_name>",params,params_len));
   if(myid==0) printf("Got output file name for gluon propagator: %s\n",pprop_name);
            
   sscanf(qcd_getParam("<momentum>",params,params_len),"%hd %hd %hd %hd",&pn[0], &pn[1], &pn[2], &pn[3]);
   if(myid==0) printf("Got momentum: [(%i+0.5)*2*pi/%i, %i*2*pi/%i, %i*2*pi/%i, %i*2*pi/%i] \n",pn[0],L[0],pn[1],L[1],pn[2],L[2],pn[3],L[3]);

   sscanf(qcd_getParam("<nstout>",params,params_len),"%hd",&nStout);
   if(myid==0) printf("Got number of stouts including without smearing: %i \n",nStout);

   strcpy(gauge_name,qcd_getParam("<cfg_name>",params,params_len));
   if(myid==0) printf("Got conf name: %s\n",gauge_name);

   px[0]=(pn[0]*2*M_PI+theta[0])/L[0]; 
   px[1]=(pn[1]*2*M_PI+theta[1])/L[1];
   px[2]=(pn[2]*2*M_PI+theta[2])/L[2];
   px[3]=(pn[3]*2*M_PI+theta[3])/L[3];

   gloops = malloc(nStout*sizeof(qcd_real_8));
   int d_dummy;
   double lf_dummy;

   free(params);
   

   if(P[0] != 1){
     fprintf(stderr,"Error only one process in t direction is allowed\n");
     exit(EXIT_FAILURE);
   }
   //#####################################################################   
   // allocate memory
  
   j = 0;
   j += qcd_initGaugeField(&u, &geo);
   MPI_Allreduce(&j, &k, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   if(k>0)
   {
      if(myid==0) printf("not enough memory\n");
      exit(EXIT_FAILURE);
   }
   
   // load gauge-field
   if(qcd_getGaugeField(gauge_name,qcd_GF_LIME,&u)) exit(EXIT_FAILURE);
   if(myid==0) printf("gauge-field loaded\n");   
   

   if(myid==0)
   {
      j=0;
      if( (fp_vfun=fopen(vfun_name,"w"))==NULL) j++;
      if( (fp_gloops=fopen(gloops_name,"r"))==NULL) j++;
      if( (fp_props=fopen(pprop_name,"w"))==NULL ) j++;
      if(j>0) fprintf(stderr,"Error while opening input/output files!\n");
   }
   MPI_Bcast(&j,1,MPI_INT, 0, MPI_COMM_WORLD);
   if(j>0) exit(EXIT_FAILURE);
   
      

   ////////////////////////////////////////////////////////////////////////////////////////////////////////
   //create and print out the vertex functions
   memset(&(gprop[0][0][0][0].re),0,4*4*3*3*sizeof(qcd_complex_16));
   memset(&(ux_lFT[0][0][0].re),0,4*3*3*sizeof(qcd_complex_16));
   memset(&(uy_lFT[0][0][0].re),0,4*3*3*sizeof(qcd_complex_16));
   memset(&(ux_gFT[0][0][0].re),0,4*3*3*sizeof(qcd_complex_16));
   memset(&(uy_gFT[0][0][0].re),0,4*3*3*sizeof(qcd_complex_16));


   for(t=0; t<geo.lL[0]; t++)
     for(x=0; x<geo.lL[1]; x++)
       for(y=0; y<geo.lL[2]; y++)
	 for(z=0; z<geo.lL[3]; z++)
	   {
	     tmp = 0;
	     tmp = (t+geo.lL[0]*geo.Pos[0])*px[0];
	     tmp+= (x+geo.lL[1]*geo.Pos[1])*px[1];
	     tmp+= (y+geo.lL[2]*geo.Pos[2])*px[2];
	     tmp+= (z+geo.lL[3]*geo.Pos[3])*px[3];
	     phase_factor_x = (qcd_complex_16) {cos(tmp), -sin(tmp)};
	     phase_factor_y=qcd_CONJ(phase_factor_x);
	     
	     i = qcd_LEXIC(t,x,y,z,geo.lL);
	     for(mu=0; mu < 4 ; mu++)
	       for(c1=0; c1<3; c1++)
		 for(c2=0; c2<3; c2++)
		   {
		     ux_lFT[mu][c1][c2] = qcd_CADD(ux_lFT[mu][c1][c2], qcd_CMUL(phase_factor_x, u.D[i][mu][c1][c2]));
		     uy_lFT[mu][c1][c2] = qcd_CADD(uy_lFT[mu][c1][c2], qcd_CMUL(phase_factor_y, u.D[i][mu][c1][c2]));
		   }
	     
	   }
   


   MPI_Reduce(&(ux_lFT[0][0][0].re), &(ux_gFT[0][0][0].re), 4*3*3*2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(&(uy_lFT[0][0][0].re), &(uy_gFT[0][0][0].re), 4*3*3*2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

   // compute the gluon propagator in momentum space
   if(myid == 0){
     for(mu = 0 ; mu < 4 ; mu++)
       for(nu = 0 ; nu < 4 ; nu++)
	 qcd_MUL3x3(gprop[mu][nu], ux_gFT[mu], uy_gFT[nu]);

     // write gluon propagator
     for(mu = 0 ; mu < 4 ; mu++)
       for(nu = 0 ; nu < 4 ; nu++)
	 for(c1 = 0 ; c1 < 3 ; c1++)
	   for(c2 = 0 ; c2 < 3 ; c2++)
	     fprintf(fp_props, "%d %d %d %d \t %+e %+e\n", mu, nu, c1, c2, gprop[mu][nu][c1][c2].re/geo.V, gprop[mu][nu][c1][c2].im/geo.V);

     // read the gluon loop
     for(int i = 0 ; i <= nStout ; i++) fscanf(fp_gloops,"%d %d %d %d %lf %lf",&d_dummy,&d_dummy,&d_dummy,&d_dummy,&(gloops[i]),&lf_dummy);

     for(int i = 0 ; i <= nStout ; i++)
       for(mu = 0 ; mu < 4 ; mu++)
	 for(nu = 0 ; nu < 4 ; nu++)
	   for(c1 = 0 ; c1 < 3 ; c1++)
	     for(c2 = 0 ; c2 < 3 ; c2++)
	       fprintf(fp_vfun, "%d %d %d %d %d \t %+e %+e\n",i, mu, nu, c1, c2, gprop[mu][nu][c1][c2].re * gloops[i] / geo.V, gprop[mu][nu][c1][c2].im * gloops[i] / geo.V);
     
   }
   free(gloops);
   qcd_destroyGeometry(&geo);
   qcd_destroyGaugeField(&u);
   MPI_Finalize();
}//end main 

