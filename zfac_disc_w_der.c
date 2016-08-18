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
   FILE *fp_loops;
   FILE *fp_props;

   FILE *fp_vfun_vD, *fp_vfun_aD, *fp_vfun_tD;
   FILE *fp_loops_vD, *fp_loops_aD, *fp_loops_tD;

   int params_len;               // needed to read inputfiles
   char *params;                 // needed to read inputfiles

   char vfun_name[qcd_MAX_STRING_LENGTH];     
   char vfun_name_vD[qcd_MAX_STRING_LENGTH];
   char vfun_name_aD[qcd_MAX_STRING_LENGTH];
   char vfun_name_tD[qcd_MAX_STRING_LENGTH];

   char param_name[qcd_MAX_STRING_LENGTH];      // name of parameter file

   qcd_geometry geo;                            // geometry structure
   char uprop_name[qcd_MAX_STRING_LENGTH];       // name of inverted Fourier source
   char dprop_name[qcd_MAX_STRING_LENGTH];       // name of inverted Fourier source
   char loops_name[qcd_MAX_STRING_LENGTH];
   char loops_name_vD[qcd_MAX_STRING_LENGTH];
   char loops_name_aD[qcd_MAX_STRING_LENGTH];
   char loops_name_tD[qcd_MAX_STRING_LENGTH];
   char pprop_name[qcd_MAX_STRING_LENGTH];

   qcd_propagator uprop, dprop;           // propagator

   qcd_real_8 theta[4] = {M_PI,0.0,0.0,0.0};    // antiperiodic b.c. in time
   qcd_uint_2 L[4];
   qcd_uint_2 P[4];
   qcd_complex_16 phase_factor, C;
   
   qcd_complex_16 u_vfun[16][4][4][3][3];
   qcd_complex_16 d_vfun[16][4][4][3][3];

   qcd_complex_16 u_vfun_vD[10][4][4][3][3];
   qcd_complex_16 d_vfun_vD[10][4][4][3][3];

   qcd_complex_16 u_vfun_aD[10][4][4][3][3];
   qcd_complex_16 d_vfun_aD[10][4][4][3][3];

   qcd_complex_16 u_vfun_tD[64][4][4][3][3];
   qcd_complex_16 d_vfun_tD[64][4][4][3][3];

   qcd_complex_16 uprop_sum[4][4][3][3];
   qcd_complex_16 uprop_sum_tmp[4][4][3][3];

   qcd_complex_16 dprop_sum[4][4][3][3];
   qcd_complex_16 dprop_sum_tmp[4][4][3][3];


   qcd_real_8 loops[16];
   qcd_real_8 loops_vD[10];
   qcd_real_8 loops_aD[10];
   qcd_real_8 loops_tD[64];
     
   qcd_int_2  pn[4];                            //integer momentum
   qcd_real_8 p[4];                             //= n*2*pi/L


   int myid,numprocs, namelen;    
   char processor_name[MPI_MAX_PROCESSOR_NAME];     
   				 
		      
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
        
   strcpy(uprop_name,qcd_getParam("<upropagator>",params,params_len));
   if(myid==0) printf("Got propagator file name: %s\n",uprop_name);

   strcpy(dprop_name,qcd_getParam("<dpropagator>",params,params_len));
   if(myid==0) printf("Got propagator file name: %s\n",dprop_name);


   strcpy(loops_name,qcd_getParam("<loops_name>",params,params_len));
   if(myid==0) printf("Got loop file name: %s\n",loops_name);

   strcpy(loops_name_vD,qcd_getParam("<loops_name_vD>",params,params_len));
   if(myid==0) printf("Got loop file name: %s\n",loops_name_vD);

   strcpy(loops_name_aD,qcd_getParam("<loops_name_aD>",params,params_len));
   if(myid==0) printf("Got loop file name: %s\n",loops_name_aD);

   strcpy(loops_name_tD,qcd_getParam("<loops_name_tD>",params,params_len));
   if(myid==0) printf("Got loop file name: %s\n",loops_name_tD);

                                 
   strcpy(vfun_name,qcd_getParam("<vfd_ultra_local_name>",params,params_len));
   if(myid==0) printf("Got output file name: %s\n",vfun_name);

   strcpy(vfun_name_vD,qcd_getParam("<vfd_vD_name>",params,params_len));
   if(myid==0) printf("Got output file name: %s\n",vfun_name_vD);

   strcpy(vfun_name_aD,qcd_getParam("<vfd_aD_name>",params,params_len));
   if(myid==0) printf("Got output file name: %s\n",vfun_name_aD);

   strcpy(vfun_name_tD,qcd_getParam("<vfd_tD_name>",params,params_len));
   if(myid==0) printf("Got output file name: %s\n",vfun_name_tD);


   strcpy(pprop_name,qcd_getParam("<pprop_name>",params,params_len));
   if(myid==0) printf("Got output file name: %s\n",pprop_name);
            
   sscanf(qcd_getParam("<momentum>",params,params_len),"%hd %hd %hd %hd",&pn[0], &pn[1], &pn[2], &pn[3]);
   if(myid==0) printf("Got momentum: [(%i+0.5)*2*pi/%i, %i*2*pi/%i, %i*2*pi/%i, %i*2*pi/%i] \n",pn[0],L[0],pn[1],L[1],pn[2],L[2],pn[3],L[3]);

   p[0]=(pn[0]*2*M_PI+theta[0])/L[0]; //!!!!!!!!!! this I am not sure need to be checked if needed here
   p[1]=(pn[1]*2*M_PI+theta[1])/L[1];
   p[2]=(pn[2]*2*M_PI+theta[2])/L[2];
   p[3]=(pn[3]*2*M_PI+theta[3])/L[3];
    
   free(params);
   

   if(P[0] != 1){
     fprintf(stderr,"Error only one process in t direction is allowed\n");
     exit(EXIT_FAILURE);
   }
   //#####################################################################   
   // allocate memory
  
   j = 0;
   j += qcd_initPropagator(&uprop, &geo);
   j += qcd_initPropagator(&dprop, &geo);
   
   MPI_Allreduce(&j, &k, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   if(k>0)
   {
      if(myid==0) printf("not enough memory\n");
      exit(EXIT_FAILURE);
   }
   if(myid==0) printf("memory for propagators allocated\n");
   
   
   //##############################################################################   
   // load propagator
   if(qcd_getPropagator(uprop_name,qcd_PROP_LIME, &uprop)) exit(EXIT_FAILURE);
   if(myid==0) printf("propagator loaded\n");   

   if(qcd_getPropagator(dprop_name,qcd_PROP_LIME, &dprop)) exit(EXIT_FAILURE);
   if(myid==0) printf("propagator loaded\n");      
   //################################################################################
   // transform propagators to basis with theta-periodic boundaries in the temporal direction

   for(lt=0; lt<geo.lL[0]; lt++)
   {
      t = lt + geo.Pos[0] * geo.lL[0];
      phase_factor   = (qcd_complex_16) {cos(theta[0]*t/geo.L[0]),sin(theta[0]*t/geo.L[0])};
      qcd_mulPropagatorC3d(&uprop, phase_factor, lt);
      qcd_mulPropagatorC3d(&dprop, phase_factor, lt);
   }
   if(myid==0) printf("propagators transformed to basis with theta-periodic boundary conditions\n");

   ////////////////////
   qcd_tranformPropagatorPhysicalPlus(&uprop,&geo,0, 0, geo.L[0] );
   qcd_tranformPropagatorPhysicalMinus(&dprop,&geo,0, 0, geo.L[0] );
   
   if(myid==0) printf("propagators transformed to physical basis\n");

   if(myid==0)
   {
      j=0;
      if( (fp_vfun=fopen(vfun_name,"w"))==NULL) j++;
      if( (fp_vfun_vD=fopen(vfun_name_vD,"w"))==NULL) j++;
      if( (fp_vfun_aD=fopen(vfun_name_aD,"w"))==NULL) j++;
      if( (fp_vfun_tD=fopen(vfun_name_tD,"w"))==NULL) j++;

      if( (fp_loops=fopen(loops_name,"r"))==NULL) j++;
      if( (fp_loops_vD=fopen(loops_name_vD,"r"))==NULL) j++;
      if( (fp_loops_aD=fopen(loops_name_aD,"r"))==NULL) j++;
      if( (fp_loops_tD=fopen(loops_name_tD,"r"))==NULL) j++;

      if( (fp_props=fopen(pprop_name,"w"))==NULL ) j++;
      if(j>0) fprintf(stderr,"Error while opening input/output files!\n");
   }
   MPI_Bcast(&j,1,MPI_INT, 0, MPI_COMM_WORLD);
   if(j>0) exit(EXIT_FAILURE);
   
      

   ////////////////////////////////////////////////////////////////////////////////////////////////////////
   //create and print out the vertex functions
   memset(&(uprop_sum[0][0][0][0].re),0,4*4*3*3*sizeof(qcd_complex_16));
   memset(&(uprop_sum_tmp[0][0][0][0].re),0,4*4*3*3*sizeof(qcd_complex_16));

   memset(&(dprop_sum[0][0][0][0].re),0,4*4*3*3*sizeof(qcd_complex_16));
   memset(&(dprop_sum_tmp[0][0][0][0].re),0,4*4*3*3*sizeof(qcd_complex_16));

   memset(&(u_vfun[0][0][0][0][0].re),0,16*4*4*3*3*sizeof(qcd_complex_16));
   memset(&(d_vfun[0][0][0][0][0].re),0,16*4*4*3*3*sizeof(qcd_complex_16));

   memset(&(u_vfun_vD[0][0][0][0][0].re),0,10*4*4*3*3*sizeof(qcd_complex_16));
   memset(&(d_vfun_vD[0][0][0][0][0].re),0,10*4*4*3*3*sizeof(qcd_complex_16));

   memset(&(u_vfun_aD[0][0][0][0][0].re),0,10*4*4*3*3*sizeof(qcd_complex_16));
   memset(&(d_vfun_aD[0][0][0][0][0].re),0,10*4*4*3*3*sizeof(qcd_complex_16));
   
   memset(&(u_vfun_tD[0][0][0][0][0].re),0,64*4*4*3*3*sizeof(qcd_complex_16));
   memset(&(d_vfun_tD[0][0][0][0][0].re),0,64*4*4*3*3*sizeof(qcd_complex_16));

   for(t=0; t<geo.lL[0]; t++)
     for(x=0; x<geo.lL[1]; x++)
       for(y=0; y<geo.lL[2]; y++)
	 for(z=0; z<geo.lL[3]; z++)
	   {
	     tmp = (t+geo.lL[0]*geo.Pos[0])*p[0];
	     tmp+= (x+geo.lL[1]*geo.Pos[1])*p[1];
	     tmp+= (y+geo.lL[2]*geo.Pos[2])*p[2];
	     tmp+= (z+geo.lL[3]*geo.Pos[3])*p[3];
	     phase_factor = (qcd_complex_16) {cos(tmp), -sin(tmp)};
	     
	     i = qcd_LEXIC(t,x,y,z,geo.lL);
	     for(mu=0; mu<4; mu++)
	       for(nu=0; nu<4; nu++)
		 for(c1=0; c1<3; c1++)
		   for(c2=0; c2<3; c2++)
		     {
		       uprop_sum[mu][nu][c1][c2] = qcd_CADD(uprop_sum[mu][nu][c1][c2], qcd_CMUL(phase_factor, uprop.D[i][mu][nu][c1][c2]));
		       dprop_sum[mu][nu][c1][c2] = qcd_CADD(dprop_sum[mu][nu][c1][c2], qcd_CMUL(phase_factor, dprop.D[i][mu][nu][c1][c2]));
		     }
	     
	   }
   

   if(myid==0) printf("local vertex functions calculated\n");
   MPI_Reduce(&(uprop_sum[0][0][0][0].re), &(uprop_sum_tmp[0][0][0][0].re), 4*4*3*3*2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(&(dprop_sum[0][0][0][0].re), &(dprop_sum_tmp[0][0][0][0].re), 4*4*3*3*2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);



   //global sums and output into files      
   if(myid == 0){
     
     for(int i = 0 ; i < 16 ; i++) fscanf(fp_loops,"%lf",&(loops[i]));
     for(int i = 0 ; i < 10 ; i++) fscanf(fp_loops_vD,"%lf",&(loops_vD[i]));
     for(int i = 0 ; i < 10 ; i++) fscanf(fp_loops_aD,"%lf",&(loops_aD[i]));
     for(int i = 0 ; i < 64 ; i++) fscanf(fp_loops_tD,"%lf",&(loops_tD[i]));

     for(mu=0; mu<4; mu++)
       for(nu=0; nu<4; nu++)
	 for(c1=0; c1<3; c1++)
	   for(c2=0; c2<3; c2++)
	     {
	       ////////////////// up /////////////
	       u_vfun[0][mu][nu][c1][c2].re = loops[0] * uprop_sum_tmp[mu][nu][c1][c2].re; //scalar
	       u_vfun[0][mu][nu][c1][c2].im = loops[0] * uprop_sum_tmp[mu][nu][c1][c2].im;

	       u_vfun[1][mu][nu][c1][c2].re = loops[1] * uprop_sum_tmp[mu][nu][c1][c2].re; //pseudoscalar
	       u_vfun[1][mu][nu][c1][c2].im = loops[1] * uprop_sum_tmp[mu][nu][c1][c2].im;

	       for(int i = 0 ; i < 4 ; i++){
		 u_vfun[2+i][mu][nu][c1][c2].re = -loops[2+i] * uprop_sum_tmp[mu][nu][c1][c2].im; // for the vector loops are imaginary
		 u_vfun[2+i][mu][nu][c1][c2].im = loops[2+i] * uprop_sum_tmp[mu][nu][c1][c2].re;
	       }

	       for(int i = 0 ; i < 4 ; i++){
		 u_vfun[6+i][mu][nu][c1][c2].re = loops[6+i] * uprop_sum_tmp[mu][nu][c1][c2].re; // axial
		 u_vfun[6+i][mu][nu][c1][c2].im = loops[6+i] * uprop_sum_tmp[mu][nu][c1][c2].im;
	       }

	       for(int i = 0 ; i < 6 ; i++){
		 u_vfun[10+i][mu][nu][c1][c2].re = -loops[10+i] * uprop_sum_tmp[mu][nu][c1][c2].im; // for the tensor loops are imaginary
		 u_vfun[10+i][mu][nu][c1][c2].im = loops[10+i] * uprop_sum_tmp[mu][nu][c1][c2].re;		 
	       }
	       ////////////////// down /////////////
	       d_vfun[0][mu][nu][c1][c2].re = loops[0] * dprop_sum_tmp[mu][nu][c1][c2].re; //scalar
	       d_vfun[0][mu][nu][c1][c2].im = loops[0] * dprop_sum_tmp[mu][nu][c1][c2].im;

	       d_vfun[1][mu][nu][c1][c2].re = loops[1] * dprop_sum_tmp[mu][nu][c1][c2].re; //pseudoscalar
	       d_vfun[1][mu][nu][c1][c2].im = loops[1] * dprop_sum_tmp[mu][nu][c1][c2].im;

	       for(int i = 0 ; i < 4 ; i++){
		 d_vfun[2+i][mu][nu][c1][c2].re = -loops[2+i] * dprop_sum_tmp[mu][nu][c1][c2].im; // for the vector loops are imaginary
		 d_vfun[2+i][mu][nu][c1][c2].im = loops[2+i] * dprop_sum_tmp[mu][nu][c1][c2].re;
	       }

	       for(int i = 0 ; i < 4 ; i++){
		 d_vfun[6+i][mu][nu][c1][c2].re = loops[6+i] * dprop_sum_tmp[mu][nu][c1][c2].re; // axial
		 d_vfun[6+i][mu][nu][c1][c2].im = loops[6+i] * dprop_sum_tmp[mu][nu][c1][c2].im;
	       }
	       
	       for(int i = 0 ; i < 6 ; i++){
		 d_vfun[10+i][mu][nu][c1][c2].re = -loops[10+i] * dprop_sum_tmp[mu][nu][c1][c2].im; // for the tensor loops are imaginary
		 d_vfun[10+i][mu][nu][c1][c2].im = loops[10+i] * dprop_sum_tmp[mu][nu][c1][c2].re;		 
	       }

	       // derivatives vD are real
	       for(int i = 0 ; i < 10 ; i++){ // tha order the data as what Martha asked for (see create_localLoops_w_der.c)
		 u_vfun_vD[i][mu][nu][c1][c2].re = loops_vD[i] * uprop_sum_tmp[mu][nu][c1][c2].re;
		 u_vfun_vD[i][mu][nu][c1][c2].im = loops_vD[i] * uprop_sum_tmp[mu][nu][c1][c2].im;
		 d_vfun_vD[i][mu][nu][c1][c2].re = loops_vD[i] * dprop_sum_tmp[mu][nu][c1][c2].re;
		 d_vfun_vD[i][mu][nu][c1][c2].im = loops_vD[i] * dprop_sum_tmp[mu][nu][c1][c2].im;
	       }

	       // derivatives aD are imag
	       for(int i = 0 ; i < 10 ; i++){
		 u_vfun_aD[i][mu][nu][c1][c2].re = -loops_aD[i] * uprop_sum_tmp[mu][nu][c1][c2].im; 
		 u_vfun_aD[i][mu][nu][c1][c2].im = loops_aD[i] * uprop_sum_tmp[mu][nu][c1][c2].re; 
		 d_vfun_aD[i][mu][nu][c1][c2].re = -loops_aD[i] * dprop_sum_tmp[mu][nu][c1][c2].im; 
		 d_vfun_aD[i][mu][nu][c1][c2].im = loops_aD[i] * dprop_sum_tmp[mu][nu][c1][c2].re; 
	       }

	       // derivatives tD are real
	       for(int i = 0 ; i < 64 ; i++){
		 u_vfun_tD[i][mu][nu][c1][c2].re = loops_tD[i] * uprop_sum_tmp[mu][nu][c1][c2].re;
		 u_vfun_tD[i][mu][nu][c1][c2].im = loops_tD[i] * uprop_sum_tmp[mu][nu][c1][c2].im;
		 d_vfun_tD[i][mu][nu][c1][c2].re = loops_tD[i] * dprop_sum_tmp[mu][nu][c1][c2].re;
		 d_vfun_tD[i][mu][nu][c1][c2].im = loops_tD[i] * dprop_sum_tmp[mu][nu][c1][c2].im;
	       }

	     }
     // write data

       for(mu=0; mu<4; mu++)
	 for(nu=0; nu<4; nu++)
	   for(c1=0; c1<3; c1++)
	     for(c2=0; c2<3; c2++)
	       {
		 fprintf(fp_props,"%d %d %d %d \t %+e %+e \t %+e %+e\n",mu,nu,c1,c2,uprop_sum_tmp[mu][nu][c1][c2].re/geo.V, uprop_sum_tmp[mu][nu][c1][c2].im/geo.V, dprop_sum_tmp[mu][nu][c1][c2].re/geo.V, dprop_sum_tmp[mu][nu][c1][c2].im/geo.V);
	       }

     for(int i = 0 ; i < 16 ; i++)
       for(mu=0; mu<4; mu++)
	 for(nu=0; nu<4; nu++)
	   for(c1=0; c1<3; c1++)
	     for(c2=0; c2<3; c2++)
	       {
		 fprintf(fp_vfun,"%d \t %d %d %d %d \t %+e %+e \t %+e %+e\n",i,mu,nu,c1,c2,u_vfun[i][mu][nu][c1][c2].re/geo.V, u_vfun[i][mu][nu][c1][c2].im/geo.V, d_vfun[i][mu][nu][c1][c2].re/geo.V, d_vfun[i][mu][nu][c1][c2].im/geo.V);
	       }

     for(int i = 0 ; i < 10 ; i++)
       for(mu=0; mu<4; mu++)
	 for(nu=0; nu<4; nu++)
	   for(c1=0; c1<3; c1++)
	     for(c2=0; c2<3; c2++)
	       {
		 fprintf(fp_vfun_vD,"%d \t %d %d %d %d \t %+e %+e \t %+e %+e\n",i,mu,nu,c1,c2,u_vfun_vD[i][mu][nu][c1][c2].re/geo.V, u_vfun_vD[i][mu][nu][c1][c2].im/geo.V, d_vfun_vD[i][mu][nu][c1][c2].re/geo.V, d_vfun_vD[i][mu][nu][c1][c2].im/geo.V);
		 fprintf(fp_vfun_aD,"%d \t %d %d %d %d \t %+e %+e \t %+e %+e\n",i,mu,nu,c1,c2,u_vfun_aD[i][mu][nu][c1][c2].re/geo.V, u_vfun_aD[i][mu][nu][c1][c2].im/geo.V, d_vfun_aD[i][mu][nu][c1][c2].re/geo.V, d_vfun_aD[i][mu][nu][c1][c2].im/geo.V);
	       }

     for(int i = 0 ; i < 64 ; i++)
       for(mu=0; mu<4; mu++)
	 for(nu=0; nu<4; nu++)
	   for(c1=0; c1<3; c1++)
	     for(c2=0; c2<3; c2++)
	       {
		 fprintf(fp_vfun_tD,"%d \t %d %d %d %d \t %+e %+e \t %+e %+e\n",i,mu,nu,c1,c2,u_vfun_tD[i][mu][nu][c1][c2].re/geo.V, u_vfun_tD[i][mu][nu][c1][c2].im/geo.V, d_vfun_tD[i][mu][nu][c1][c2].re/geo.V, d_vfun_tD[i][mu][nu][c1][c2].im/geo.V);
	       }

   }
   
   qcd_destroyPropagator(&uprop);
   qcd_destroyPropagator(&dprop);
   qcd_destroyGeometry(&geo);
   MPI_Finalize();
}//end main 

   /*
   FILE *fp_uprop;
   FILE *fp_dprop;
   
   fp_uprop = fopen("uprop.0100.dat","w");
   fp_dprop = fopen("dprop.0100.dat","w");

   for(t=0; t<geo.lL[0]; t++)
     for(x=0; x<geo.lL[1]; x++)
       for(y=0; y<geo.lL[2]; y++)
	 for(z=0; z<geo.lL[3]; z++)
	   {
	     int iv = qcd_LEXIC(t,x,y,z,geo.lL);
	     for(mu=0; mu<4; mu++)
	       for(nu=0; nu<4; nu++)
		 for(c1=0; c1<3; c1++)
		   for(c2=0; c2<3; c2++)
		     {
		       fprintf(fp_uprop,"%+16.15e %+16.15e\n",uprop.D[iv][mu][nu][c1][c2].re, uprop.D[iv][mu][nu][c1][c2].im);
		       fprintf(fp_dprop,"%+16.15e %+16.15e\n",dprop.D[iv][mu][nu][c1][c2].re, dprop.D[iv][mu][nu][c1][c2].im);
		     }
	   }
   exit(-1);
   */
   ////////////////////
