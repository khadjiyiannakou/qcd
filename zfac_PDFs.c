/* zfac_PDFs.c
 *
 * reads half fourier transformed propagator and
 * writes out:
 *
 * 1) the propagator with momentum p (12x12 matrix)
 *    S(p) = sum_xy exp(-ip(x-y)) < d(x) \bar d(y) >
 *
 * 2) the bare vertex function with momentum p (12x12 matrix)
 *    for  the insertion \bar{u}(w) L(w->w+z) d(w+z)
 *    and \bar{u}(w) L(w->w-z) d(w-z) where L is a link line
 *
 * Kyriakos Hadjiyiannakou 2016
 ************************************************************************/
 
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>

typedef struct
{
  qcd_complex_16 (*D)[3][3]; //data
  qcd_geometry *geo;
  qcd_uint_2 initialized;
  int mu;                   // direction of the wilson line
  int Nz;
} qcd_WilsonLine;

int qcd_initWilsonLine(qcd_WilsonLine *Wl, qcd_geometry *geo, int mu)
{

  if(geo->L[mu] != geo->lL[mu]){
    fprintf(stderr,"error in qcd_initWilsonLine! The direction where you put the Wilson line cannot be partitioned.\n");
    return(1);
  }

  qcd_uint_8 size;
  if(!geo->initialized)
    {
      fprintf(stderr,"error in qcd_initWilsonLine! Geometry must be properly initialized.\n");
      return(1);
    }
  size = geo->lV;
  Wl->D = malloc(3*3*size*(geo->L[mu]/2)*sizeof(qcd_complex_16));
  if(Wl->D == NULL)
    {
      fprintf(stderr,"process %i: Error in qcd_initWilsonLine! Out of memory\n",geo->myid);
      return(1);
    }
  Wl->geo = geo;
  Wl->mu = mu;
  Wl->Nz = geo->L[mu]/2;
  Wl->initialized = 1;
  return(0);
}

void qcd_destroyWilsonLine(qcd_WilsonLine *Wl){
  if(Wl->initialized != 1){
    fprintf(stderr,"Error in qcd_destroyWilsonLine! qcd_WilsonLine not initialized\n");
    return;
  }
  free(Wl->D);
}

int qcd_createWilsonLineFwd(qcd_WilsonLine *Wl, qcd_gaugeField *u){
  if(Wl->initialized != 1){
    fprintf(stderr,"error in qcd_createWilsonLineFwd! initialize first before use\n");
    return(1);
  }

  qcd_uint_8 pointPlus;
  qcd_uint_8 i;
  for(int z=0;z<Wl->geo->lL[3];z++)
    for(int y=0;y<Wl->geo->lL[2];y++)
      for(int x=0;x<Wl->geo->lL[1];x++)
	for(int t=0;t<Wl->geo->lL[0];t++)
	  {
	    i = qcd_LEXIC(t,x,y,z,Wl->geo->lL);
	    for(int iz = 0 ; iz < Wl->Nz ; iz++){
	      if(iz==0){
		for(int c1 = 0 ; c1 < 3 ; c1++)
		  for(int c2 = 0 ; c2 < 3 ; c2++){
		    Wl->D[i*Wl->Nz+iz][c1][c2].re = (c1==c2) ? 1. : 0.;
		    Wl->D[i*Wl->Nz+iz][c1][c2].im = 0.;
		  }
	      }
	      else if(iz==1){
		for(int c1 = 0 ; c1 < 3 ; c1++)
		  for(int c2 = 0 ; c2 < 3 ; c2++)
		    Wl->D[i*Wl->Nz+iz][c1][c2] = u->D[i][Wl->mu][c1][c2];
	      }
	      else{
		if(Wl->mu == 0 )
		  pointPlus = qcd_LEXIC((t+iz-1)%Wl->geo->lL[0],x,y,z,Wl->geo->lL);
		if(Wl->mu == 1 )
		  pointPlus = qcd_LEXIC(t,(x+iz-1)%Wl->geo->lL[1],y,z,Wl->geo->lL);
		if(Wl->mu == 2 )
		  pointPlus = qcd_LEXIC(t,x,(y+iz-1)%Wl->geo->lL[2],z,Wl->geo->lL);
		if(Wl->mu == 3 )
		  pointPlus = qcd_LEXIC(t,x,y,(z+iz-1)%Wl->geo->lL[3],Wl->geo->lL);

		int stepback = iz - 1;
		for(int c1 = 0 ; c1 < 3 ; c1++)
                  for(int c2 = 0 ; c2 < 3 ; c2++){
		    Wl->D[i*Wl->Nz+iz][c1][c2] = (qcd_complex_16) {0.,0.};
		    for(int c3 = 0 ; c3 < 3 ; c3++)
		      Wl->D[i*Wl->Nz+iz][c1][c2] = qcd_CADD(Wl->D[i*Wl->Nz+iz][c1][c2], qcd_CMUL(Wl->D[i*Wl->Nz+stepback][c1][c3],
												 u->D[pointPlus][Wl->mu][c3][c2]));
		  }
		    
	      }
	    }

	  }
  return(0);
}


int qcd_createWilsonLineBwd(qcd_WilsonLine *Wl, qcd_gaugeField *u){
  if(Wl->initialized != 1){
    fprintf(stderr,"error in qcd_createWilsonLineBwd! initialize first before use\n");
    return(1);
  }

  qcd_uint_8 pointMinus;
  qcd_uint_8 i;
  for(int z=0;z<Wl->geo->lL[3];z++)
    for(int y=0;y<Wl->geo->lL[2];y++)
      for(int x=0;x<Wl->geo->lL[1];x++)
	for(int t=0;t<Wl->geo->lL[0];t++)
	  {
	    i = qcd_LEXIC(t,x,y,z,Wl->geo->lL);
	    for(int iz = 0 ; iz < Wl->Nz ; iz++){
	      if(iz==0){
		for(int c1 = 0 ; c1 < 3 ; c1++)
		  for(int c2 = 0 ; c2 < 3 ; c2++){
		    Wl->D[i*Wl->Nz+iz][c1][c2].re = (c1==c2) ? 1. : 0.;
		    Wl->D[i*Wl->Nz+iz][c1][c2].im = 0.;
		  }
	      }
	      else if(iz==1){

		if(Wl->mu == 0 )
		  pointMinus = qcd_LEXIC((t-iz+Wl->geo->lL[0])%Wl->geo->lL[0],x,y,z,Wl->geo->lL);
		if(Wl->mu == 1 )
		  pointMinus = qcd_LEXIC(t,(x-iz+Wl->geo->lL[1])%Wl->geo->lL[1],y,z,Wl->geo->lL);
		if(Wl->mu == 2 )
		  pointMinus = qcd_LEXIC(t,x,(y-iz+Wl->geo->lL[2])%Wl->geo->lL[2],z,Wl->geo->lL);
		if(Wl->mu == 3 )
		  pointMinus = qcd_LEXIC(t,x,y,(z-iz+Wl->geo->lL[3])%Wl->geo->lL[3],Wl->geo->lL);

		for(int c1 = 0 ; c1 < 3 ; c1++)
		  for(int c2 = 0 ; c2 < 3 ; c2++)
		    Wl->D[i*Wl->Nz+iz][c1][c2] = qcd_CONJ(u->D[pointMinus][Wl->mu][c2][c1]);
	      }
	      else{

		if(Wl->mu == 0 )
		  pointMinus = qcd_LEXIC((t-iz+Wl->geo->lL[0])%Wl->geo->lL[0],x,y,z,Wl->geo->lL);
		if(Wl->mu == 1 )
		  pointMinus = qcd_LEXIC(t,(x-iz+Wl->geo->lL[1])%Wl->geo->lL[1],y,z,Wl->geo->lL);
		if(Wl->mu == 2 )
		  pointMinus = qcd_LEXIC(t,x,(y-iz+Wl->geo->lL[2])%Wl->geo->lL[2],z,Wl->geo->lL);
		if(Wl->mu == 3 )
		  pointMinus = qcd_LEXIC(t,x,y,(z-iz+Wl->geo->lL[3])%Wl->geo->lL[3],Wl->geo->lL);

		int stepback = iz - 1;
		for(int c1 = 0 ; c1 < 3 ; c1++)
                  for(int c2 = 0 ; c2 < 3 ; c2++){
		    Wl->D[i*Wl->Nz+iz][c1][c2] = (qcd_complex_16) {0.,0.};
		    for(int c3 = 0 ; c3 < 3 ; c3++)
		      Wl->D[i*Wl->Nz+iz][c1][c2] = qcd_CADD(Wl->D[i*Wl->Nz+iz][c1][c2], qcd_CMUL(Wl->D[i*Wl->Nz+stepback][c1][c3],
												 qcd_CONJ(u->D[pointMinus][Wl->mu][c2][c3])));
		  }
		    
	      }
	    }

	  }
  return(0);
}


int main(int argc,char* argv[])
{
   qcd_uint_2 mu,nu,rho,ku,lu,c1,c2;          // various loop variables
   qcd_uint_4 i,j,k,lt,x,y,z,t,ip1,im1; 
   qcd_uint_2 id1,id2,id3,id4;
   qcd_uint_2 ic1,ic2,ic3,ic4;                    
   qcd_real_8 tmp;                            // general purpuse
   qcd_uint_2 xx[4];
   int faciip1, faciim1;                      // factors related to b.c.
  
   /* FILE *fp_vfun_v, *fp_vfun_a;   // output files */
   /* FILE *fp_vfun_s, *fp_vfun_p; */
   /* FILE *fp_vfun_t, *fp_vfun_vD, *fp_vfun_aD; */
   /* FILE *fp_vfun_tD, *fp_vfun_d1; */

   FILE *fp_vfun_v_WL;
   FILE *fp_vfun_a_WL;
   FILE *fp_vfun_t_WL;

   FILE *fp_pprop;      
  
   int params_len;               // needed to read inputfiles
   char *params;                 // needed to read inputfiles

   char gauge_name[qcd_MAX_STRING_LENGTH];      // name of gauge-configuration file
   /* char vfun_s_name[qcd_MAX_STRING_LENGTH];     // output file name, local scalar density vertex function */
   /* char vfun_p_name[qcd_MAX_STRING_LENGTH];     // output file name, local pseudoscalar density vertex function */
   /* char vfun_v_name[qcd_MAX_STRING_LENGTH];     // output file name, local vector current vertex function */
   /* char vfun_a_name[qcd_MAX_STRING_LENGTH];     // output file name, local axial current vertex function */
   /* char vfun_t_name[qcd_MAX_STRING_LENGTH];     // output file name, local tensor current vertex function */
   /* char vfun_vD_name[qcd_MAX_STRING_LENGTH];    // output file name, 1 derivative vector operator vertex function */
   /* char vfun_aD_name[qcd_MAX_STRING_LENGTH];    // output file name, 1 derivative axial operator vertex function */
   /* char vfun_tD_name[qcd_MAX_STRING_LENGTH];    // output file name, 1 derivative tensor operator vertex function */
   /* char vfun_d1_name[qcd_MAX_STRING_LENGTH];    // output file name, 1 derivative d1 operator vertex function */

   char vfun_v_WL_name[qcd_MAX_STRING_LENGTH];
   char vfun_a_WL_name[qcd_MAX_STRING_LENGTH];
   char vfun_t_WL_name[qcd_MAX_STRING_LENGTH];

   char pprop_name[qcd_MAX_STRING_LENGTH];      // name of output file, momentum propagator
   
   char param_name[qcd_MAX_STRING_LENGTH];      // name of parameter file

   qcd_geometry geo;                            // geometry structure
   qcd_propagator prop, lprop, rprop;           // propagator
   qcd_gaugeField u;                            // gauge field 
   qcd_WilsonLine Wl_fwd, Wl_bwd;

   char prop_name[qcd_MAX_STRING_LENGTH];       // name of inverted Fourier source

   qcd_real_8 theta[4] = {M_PI,0.0,0.0,0.0};    // antiperiodic b.c. in time
   qcd_uint_2 L[4];
   qcd_uint_2 P[4];
   qcd_complex_16 phase_factor, C, Cphase_factor;
   
   qcd_complex_16 Sup[4][4][3][3];              // Fourier transformed up-quark & down-quark
   qcd_complex_16 Sdown[4][4][3][3];            // Propagators

   /*
   qcd_complex_16 vfun_s[4][4][3][3];           // vertex function noether current
   qcd_complex_16 vfun_p[4][4][3][3];           // vertex function noether current
   qcd_complex_16 vfun_v[4][4][4][3][3];        // vertex function local vector current
   qcd_complex_16 vfun_a[4][4][4][3][3];        // vertex function local axial current
   qcd_complex_16 vfun_t[16][4][4][3][3];        // vertex function local tensor current
   qcd_complex_16 vfun_vD[16][4][4][3][3];      // vertex function 1 derivative vector operator
   qcd_complex_16 vfun_aD[16][4][4][3][3];      // vertex function 1 derivative axial operator
   qcd_complex_16 vfun_tD[64][4][4][3][3];      // vertex function 1 derivative tensor operator
   qcd_complex_16 vfun_d1[16][4][4][3][3];      // vertex function 1 derivative d1 operator
   qcd_complex_16 vfun_tmp[64][4][4][3][3];     // temp variable
   */

   qcd_complex_16 (*vfun_v_WL)[4][4][3][3];
   qcd_complex_16 (*vfun_a_WL)[4][4][3][3];
   qcd_complex_16 (*vfun_t_WL)[4][4][3][3];
   qcd_complex_16 (*vfun_tmp_WL)[4][4][3][3]; // temporal space for global reduction

   qcd_complex_16 *lxr;
   
   qcd_complex_16 g5sig[5][5][4][4];            // gamma_5 * [gamma_mu, gamma_nu] *1/2
   
   qcd_int_2  pn[4];                            //integer momentum
   qcd_real_8 p[4];                             //= n*2*pi/L


   int myid,numprocs, namelen;    
   char processor_name[MPI_MAX_PROCESSOR_NAME];     
   				 
   int dirLine;
             
             
   //////////////////////////////////////////////////////////////////////////////////////          
   //set up MPI
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);         // num. of processes taking part in the calculation
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);             // each process gets its ID
   MPI_Get_processor_name(processor_name,&namelen); // 
   
   for(i=0; i<5; i++)
   for(mu=0; mu<4; mu++)
   for(nu=0; nu<4; nu++)
   {  
      for(j=0; j<5; j++)
      {         
         g5sig[i][j][mu][nu]= (qcd_complex_16){0,0};

         for(ku=0; ku<4; ku++)
         for(lu=0; lu<4; lu++)
         {
            g5sig[i][j][mu][nu] = qcd_CADD(g5sig[i][j][mu][nu],qcd_CMUL(qcd_CMUL(qcd_GAMMA[5][mu][ku],
                                                                                 qcd_GAMMA[i][ku][lu]),
                                                                        qcd_GAMMA[j][lu][nu]));
            g5sig[i][j][mu][nu] = qcd_CSUB(g5sig[i][j][mu][nu],qcd_CMUL(qcd_CMUL(qcd_GAMMA[5][mu][ku],
                                                                                 qcd_GAMMA[j][ku][lu]),
                                                                        qcd_GAMMA[i][lu][nu]));
         }
         g5sig[i][j][mu][nu] = qcd_CSCALE(g5sig[i][j][mu][nu],0.5);
      }
   }   
   
   
   
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
        
   strcpy(prop_name,qcd_getParam("<propagator>",params,params_len));
   if(myid==0) printf("Got propagator file name: %s\n",prop_name);

   strcpy(gauge_name,qcd_getParam("<cfg_name>",params,params_len));
   if(myid==0) printf("Got conf name: %s\n",gauge_name);
               
   strcpy(vfun_v_WL_name,qcd_getParam("<vertex_function_v_WL_name>",params,params_len));
   if(myid==0) printf("Got output file name: %s\n",vfun_v_WL_name);

   strcpy(vfun_a_WL_name,qcd_getParam("<vertex_function_a_WL_name>",params,params_len));
   if(myid==0) printf("Got output file name: %s\n",vfun_a_WL_name);

   strcpy(vfun_t_WL_name,qcd_getParam("<vertex_function_t_WL_name>",params,params_len));
   if(myid==0) printf("Got output file name: %s\n",vfun_t_WL_name);

   strcpy(pprop_name,qcd_getParam("<pprop_name>",params,params_len));
   if(myid==0) printf("Got output file name: %s\n",pprop_name);
         
   sscanf(qcd_getParam("<momentum>",params,params_len),"%hd %hd %hd %hd",&pn[0], &pn[1], &pn[2], &pn[3]);
   if(myid==0) printf("Got momentum: [(%i+0.5)*2*pi/%i, %i*2*pi/%i, %i*2*pi/%i, %i*2*pi/%i] \n",pn[0],L[0],pn[1],L[1],pn[2],L[2],pn[3],L[3]);

   sscanf(qcd_getParam("<dirLine>",params,params_len),"%d",&dirLine);
   if(myid==0) printf("Got direction of the Wilson line:%i\n",dirLine);

   if(dirLine <= 0 || dirLine > 3){ // there is no meaning of Wilson line in the temporal direction
     fprintf(stderr,"error with the direction of the Wilson line\n");
     exit(EXIT_FAILURE);
   }


   p[0]=(pn[0]*2*M_PI+theta[0])/L[0];
   p[1]=(pn[1]*2*M_PI+theta[1])/L[1];
   p[2]=(pn[2]*2*M_PI+theta[2])/L[2];
   p[3]=(pn[3]*2*M_PI+theta[3])/L[3];
    
   free(params);
   
   
   //#####################################################################   
   // allocate memory
  
   j = 0;
   j += qcd_initPropagator(&prop, &geo);
   j += qcd_initPropagator(&lprop, &geo);
   j += qcd_initPropagator(&rprop, &geo);
   j += qcd_initGaugeField(&u, &geo);
   
   MPI_Allreduce(&j, &k, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   if(k>0)
   {
      if(myid==0) printf("not enough memory\n");
      exit(EXIT_FAILURE);
   }
   if(myid==0) printf("memory for propagators and gauge-field allocated\n");

   j = 0;
   j += qcd_initWilsonLine(&Wl_fwd,&geo,dirLine);
   j += qcd_initWilsonLine(&Wl_bwd,&geo,dirLine);

   MPI_Allreduce(&j, &k, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   if(k>0)
   {
      if(myid==0) printf("problem with qcd_initWilsonLine\n");
      exit(EXIT_FAILURE);
   }
   if(myid==0) printf("memory for Wilson lines allocated\n");

   //##############################################################################
   // load gauge-field
   if(qcd_getGaugeField(gauge_name,qcd_GF_LIME,&u)) exit(EXIT_FAILURE);
   if(myid==0) printf("gauge-field loaded\n");   
   
   qcd_communicateGaugePM(&u);
   
   // load propagator
   if(qcd_getPropagator(prop_name,qcd_PROP_LIME, &prop)) exit(EXIT_FAILURE);
   if(myid==0) printf("propagator loaded\n");   

   //#################################################################################
   if(qcd_createWilsonLineFwd(&Wl_fwd,&u)) exit(EXIT_FAILURE);
   if(qcd_createWilsonLineBwd(&Wl_bwd,&u)) exit(EXIT_FAILURE);

   // allocate memory
   j=0;
   vfun_v_WL = malloc(2*Wl_fwd.Nz*4*4*3*3*sizeof(qcd_complex_16));
   vfun_a_WL = malloc(2*Wl_fwd.Nz*4*4*3*3*sizeof(qcd_complex_16));
   vfun_t_WL = malloc(2*Wl_fwd.Nz*4*4*3*3*sizeof(qcd_complex_16));
   vfun_tmp_WL = malloc(2*Wl_fwd.Nz*4*4*3*3*sizeof(qcd_complex_16));
   lxr = malloc(2*Wl_fwd.Nz*sizeof(qcd_complex_16));
   if(vfun_v_WL == NULL)j++;
   if(vfun_a_WL == NULL)j++;
   if(vfun_t_WL == NULL)j++;
   if(lxr == NULL)j++;

   MPI_Allreduce(&j, &k, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   if(k>0)
   {
      if(myid==0) printf("problem allocate memory for vertex functions\n");
      exit(EXIT_FAILURE);
   }
   
   //################################################################################
   // transform propagators to basis with theta-periodic boundaries in the temporal direction
   for(lt=0; lt<geo.lL[0]; lt++)
   {
      t = lt + geo.Pos[0] * geo.lL[0];
      phase_factor   = (qcd_complex_16) {cos(theta[0]*t/geo.L[0]),sin(theta[0]*t/geo.L[0])};
      qcd_mulPropagatorC3d(&prop, phase_factor, lt);
   }
   if(myid==0) printf("propagators transformed to basis with theta-periodic boundary conditions\n");
   
   //transform propagators to physical basis
   for(i=0; i<geo.lV; i++)
   for(c1=0; c1<3; c1++)
   for(c2=0; c2<3; c2++)
   {
      //this works only with g5=diag(1 1 -1 -1)
      lprop.D[i][0][0][c1][c2] = (qcd_complex_16){-prop.D[i][0][0][c2][c1].im, -prop.D[i][0][0][c2][c1].re};
      rprop.D[i][0][0][c1][c2] = (qcd_complex_16){ prop.D[i][0][0][c1][c2].im, -prop.D[i][0][0][c1][c2].re};
      lprop.D[i][0][1][c1][c2] = (qcd_complex_16){-prop.D[i][1][0][c2][c1].im, -prop.D[i][1][0][c2][c1].re};
      rprop.D[i][0][1][c1][c2] = (qcd_complex_16){ prop.D[i][0][1][c1][c2].im, -prop.D[i][0][1][c1][c2].re};
      lprop.D[i][0][2][c1][c2] = (qcd_complex_16){ prop.D[i][2][0][c2][c1].re, -prop.D[i][2][0][c2][c1].im};
      rprop.D[i][0][2][c1][c2] = (qcd_complex_16){ prop.D[i][0][2][c1][c2].re,  prop.D[i][0][2][c1][c2].im};
      lprop.D[i][0][3][c1][c2] = (qcd_complex_16){ prop.D[i][3][0][c2][c1].re, -prop.D[i][3][0][c2][c1].im};
      rprop.D[i][0][3][c1][c2] = (qcd_complex_16){ prop.D[i][0][3][c1][c2].re,  prop.D[i][0][3][c1][c2].im};
      lprop.D[i][1][0][c1][c2] = (qcd_complex_16){-prop.D[i][0][1][c2][c1].im, -prop.D[i][0][1][c2][c1].re};
      rprop.D[i][1][0][c1][c2] = (qcd_complex_16){ prop.D[i][1][0][c1][c2].im, -prop.D[i][1][0][c1][c2].re};
      lprop.D[i][1][1][c1][c2] = (qcd_complex_16){-prop.D[i][1][1][c2][c1].im, -prop.D[i][1][1][c2][c1].re};
      rprop.D[i][1][1][c1][c2] = (qcd_complex_16){ prop.D[i][1][1][c1][c2].im, -prop.D[i][1][1][c1][c2].re};
      lprop.D[i][1][2][c1][c2] = (qcd_complex_16){ prop.D[i][2][1][c2][c1].re, -prop.D[i][2][1][c2][c1].im};
      rprop.D[i][1][2][c1][c2] = (qcd_complex_16){ prop.D[i][1][2][c1][c2].re,  prop.D[i][1][2][c1][c2].im};
      lprop.D[i][1][3][c1][c2] = (qcd_complex_16){ prop.D[i][3][1][c2][c1].re, -prop.D[i][3][1][c2][c1].im};
      rprop.D[i][1][3][c1][c2] = (qcd_complex_16){ prop.D[i][1][3][c1][c2].re,  prop.D[i][1][3][c1][c2].im};
      lprop.D[i][2][0][c1][c2] = (qcd_complex_16){ prop.D[i][0][2][c2][c1].re, -prop.D[i][0][2][c2][c1].im};
      rprop.D[i][2][0][c1][c2] = (qcd_complex_16){ prop.D[i][2][0][c1][c2].re,  prop.D[i][2][0][c1][c2].im};
      lprop.D[i][2][1][c1][c2] = (qcd_complex_16){ prop.D[i][1][2][c2][c1].re, -prop.D[i][1][2][c2][c1].im};
      rprop.D[i][2][1][c1][c2] = (qcd_complex_16){ prop.D[i][2][1][c1][c2].re,  prop.D[i][2][1][c1][c2].im};
      lprop.D[i][2][2][c1][c2] = (qcd_complex_16){ prop.D[i][2][2][c2][c1].im,  prop.D[i][2][2][c2][c1].re};
      rprop.D[i][2][2][c1][c2] = (qcd_complex_16){-prop.D[i][2][2][c1][c2].im,  prop.D[i][2][2][c1][c2].re};
      lprop.D[i][2][3][c1][c2] = (qcd_complex_16){ prop.D[i][3][2][c2][c1].im,  prop.D[i][3][2][c2][c1].re};
      rprop.D[i][2][3][c1][c2] = (qcd_complex_16){-prop.D[i][2][3][c1][c2].im,  prop.D[i][2][3][c1][c2].re};
      lprop.D[i][3][0][c1][c2] = (qcd_complex_16){ prop.D[i][0][3][c2][c1].re, -prop.D[i][0][3][c2][c1].im};
      rprop.D[i][3][0][c1][c2] = (qcd_complex_16){ prop.D[i][3][0][c1][c2].re,  prop.D[i][3][0][c1][c2].im};
      lprop.D[i][3][1][c1][c2] = (qcd_complex_16){ prop.D[i][1][3][c2][c1].re, -prop.D[i][1][3][c2][c1].im};
      rprop.D[i][3][1][c1][c2] = (qcd_complex_16){ prop.D[i][3][1][c1][c2].re,  prop.D[i][3][1][c1][c2].im};
      lprop.D[i][3][2][c1][c2] = (qcd_complex_16){ prop.D[i][2][3][c2][c1].im,  prop.D[i][2][3][c2][c1].re};
      rprop.D[i][3][2][c1][c2] = (qcd_complex_16){-prop.D[i][3][2][c1][c2].im,  prop.D[i][3][2][c1][c2].re};
      lprop.D[i][3][3][c1][c2] = (qcd_complex_16){ prop.D[i][3][3][c2][c1].im,  prop.D[i][3][3][c2][c1].re};
      rprop.D[i][3][3][c1][c2] = (qcd_complex_16){-prop.D[i][3][3][c1][c2].im,  prop.D[i][3][3][c1][c2].re};
   }
   if(myid==0) printf("propagators transformed to physical basis\n");
   qcd_waitall(&geo);

   
   ////////////////////////////////////////////////////////////////////////////////////////////////////////
   //create and print out the fourier transformed propagators
   memset(&(Sdown[0][0][0][0].re),0,12*12*sizeof(qcd_complex_16));
   memset(&(Sup[0][0][0][0].re),0,12*12*sizeof(qcd_complex_16));
   for(t=0; t<geo.lL[0]; t++)
   for(x=0; x<geo.lL[1]; x++)
   for(y=0; y<geo.lL[2]; y++)
   for(z=0; z<geo.lL[3]; z++)
   {
      tmp = (t+geo.lL[0]*geo.Pos[0])*p[0];
      tmp+= (x+geo.lL[1]*geo.Pos[1])*p[1];
      tmp+= (y+geo.lL[2]*geo.Pos[2])*p[2];
      tmp+= (z+geo.lL[3]*geo.Pos[3])*p[3];
      phase_factor = (qcd_complex_16) {cos(tmp), sin(tmp)};
      Cphase_factor= qcd_CONJ(phase_factor);
      i = qcd_LEXIC(t,x,y,z,geo.lL);
      for(mu=0; mu<4; mu++)
      for(nu=0; nu<4; nu++)
      for(c1=0; c1<3; c1++)
      for(c2=0; c2<3; c2++)
      {
         Sdown[mu][nu][c1][c2] = qcd_CADD(Sdown[mu][nu][c1][c2], qcd_CMUL(Cphase_factor,rprop.D[i][mu][nu][c1][c2]));
         Sup[mu][nu][c1][c2]   = qcd_CSUB(Sup[mu][nu][c1][c2],   qcd_CMUL( phase_factor,lprop.D[i][mu][nu][c1][c2]));
      }
   }   
   
   if(myid==0) printf("local momentum propagators calculated\n");

   if(myid==0)
   {
      j=0;
      if( (fp_vfun_v_WL=fopen(vfun_v_WL_name,"w"))==NULL) j++;
      if( (fp_vfun_a_WL=fopen(vfun_a_WL_name,"w"))==NULL) j++;
      if( (fp_vfun_t_WL=fopen(vfun_t_WL_name,"w"))==NULL) j++;
      if( (fp_pprop=fopen(pprop_name,"w"))==NULL) j++;
      if(j>0) fprintf(stderr,"Error while opening output files for writing!\n");
   }
   MPI_Bcast(&j,1,MPI_INT, 0, MPI_COMM_WORLD);
   if(j>0) exit(EXIT_FAILURE);
   
   if(myid==0) printf("all output files ready for writing\n");   
      
   for(mu=0; mu<4; mu++)
   for(nu=0; nu<4; nu++)
   for(c1=0; c1<3; c1++)
   for(c2=0; c2<3; c2++)
   {
      MPI_Reduce(&(Sdown[mu][nu][c1][c2].re), &(C.re), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&(Sdown[mu][nu][c1][c2].im), &(C.im), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      if(myid==0)
         fprintf(fp_pprop,"d %i %i %i %i %e %e\n",mu, nu, c1, c2, C.re/geo.V, C.im/geo.V);
      MPI_Reduce(&(Sup[mu][nu][c1][c2].re), &(C.re), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&(Sup[mu][nu][c1][c2].im), &(C.im), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      if(myid==0)
         fprintf(fp_pprop,"u %i %i %i %i %e %e\n",mu, nu, c1, c2, C.re/geo.V, C.im/geo.V);
   }

   if(myid==0) printf("global momentum propagators calculated and written\n");

   //================ create vertex functions ==================//
   memset(&(vfun_v_WL[0][0][0][0][0].re),0,2*Wl_fwd.Nz*4*4*3*3*sizeof(qcd_complex_16));
   memset(&(vfun_a_WL[0][0][0][0][0].re),0,2*Wl_fwd.Nz*4*4*3*3*sizeof(qcd_complex_16));
   memset(&(vfun_t_WL[0][0][0][0][0].re),0,2*Wl_fwd.Nz*4*4*3*3*sizeof(qcd_complex_16));
   int Nz = Wl_fwd.Nz;
   qcd_uint_8 pointPlus,pointMinus;
   for(int z=0;z<geo.lL[3];z++)
     for(int y=0;y<geo.lL[2];y++)
       for(int x=0;x<geo.lL[1];x++)
	 for(int t=0;t<geo.lL[0];t++)
	   {
	     qcd_uint_8 i = qcd_LEXIC(t,x,y,z,geo.lL);
	     for(id1=0; id1<4; id1++)
	       for(id2=0; id2<4; id2++)
		 for(id3=0; id3<4; id3++)
		   for(id4=0; id4<4; id4++)
		     for(ic1=0; ic1<3; ic1++)
		       for(ic4=0; ic4<3; ic4++)
			 { 
			   for(int iz = 0 ; iz < Nz; iz++){
			     lxr[iz] = (qcd_complex_16) {0,0};
			     lxr[Nz+iz] = (qcd_complex_16) {0,0};
			    
			     if(dirLine == 0 )
			       pointPlus = qcd_LEXIC((t+iz)%geo.lL[0],x,y,z,geo.lL);
			     if(dirLine == 1 )
			       pointPlus = qcd_LEXIC(t,(x+iz)%geo.lL[1],y,z,geo.lL);
			     if(dirLine == 2 )
			       pointPlus = qcd_LEXIC(t,x,(y+iz)%geo.lL[2],z,geo.lL);
			     if(dirLine == 3 )
			       pointPlus = qcd_LEXIC(t,x,y,(z+iz)%geo.lL[3],geo.lL);

			     if(dirLine == 0 )
			       pointMinus = qcd_LEXIC((t-iz+geo.lL[0])%geo.lL[0],x,y,z,geo.lL);
			     if(dirLine == 1 )
			       pointMinus = qcd_LEXIC(t,(x-iz+geo.lL[1])%geo.lL[1],y,z,geo.lL);
			     if(dirLine == 2 )
			       pointMinus = qcd_LEXIC(t,x,(y-iz+geo.lL[2])%geo.lL[2],z,geo.lL);
			     if(dirLine == 3 )
			       pointMinus = qcd_LEXIC(t,x,y,(z-iz+geo.lL[3])%geo.lL[3],geo.lL);

			     for(ic2=0; ic2<3; ic2++)
			       for(ic3=0; ic3<3; ic3++)
				 {
				   lxr[iz] = qcd_CADD(lxr[iz],
						      qcd_CMUL(lprop.D[i][id1][id2][ic1][ic2],
							       qcd_CMUL(Wl_fwd.D[i*Nz+iz][ic2][ic3],rprop.D[pointPlus][id3][id4][ic3][ic4])));
				   lxr[Nz+iz] = qcd_CADD(lxr[Nz+iz],
						      qcd_CMUL(lprop.D[i][id1][id2][ic1][ic2],
							       qcd_CMUL(Wl_bwd.D[i*Nz+iz][ic2][ic3],rprop.D[pointMinus][id3][id4][ic3][ic4])));

				 }
			     //============= Vector ==============//
			     if(qcd_NORM(qcd_GAMMA[dirLine][id2][id3])>1e-4){
			       vfun_v_WL[iz][id1][id4][ic1][ic4] = qcd_CADD(vfun_v_WL[iz][id1][id4][ic1][ic4],
									    qcd_CMUL(lxr[iz],qcd_GAMMA[dirLine][id2][id3]));
			       vfun_v_WL[Nz+iz][id1][id4][ic1][ic4] = qcd_CADD(vfun_v_WL[Nz+iz][id1][id4][ic1][ic4],
									    qcd_CMUL(lxr[Nz+iz],qcd_GAMMA[dirLine][id2][id3]));
			     }
			     //============ Axial ================//
			     if(qcd_NORM(qcd_G5GAMMA[dirLine][id2][id3])>1e-4){
			       vfun_a_WL[iz][id1][id4][ic1][ic4] = qcd_CADD(vfun_a_WL[iz][id1][id4][ic1][ic4],
									    qcd_CMUL(lxr[iz],qcd_G5GAMMA[dirLine][id2][id3]));
			       vfun_a_WL[Nz+iz][id1][id4][ic1][ic4] = qcd_CADD(vfun_a_WL[Nz+iz][id1][id4][ic1][ic4],
									    qcd_CMUL(lxr[Nz+iz],qcd_G5GAMMA[dirLine][id2][id3]));
			     }
			     //=========== Tensor ===============//
			     int ismn;
			     if(dirLine==1)
			       ismn=3;
			     if(dirLine==2)
			       ismn=3;
			     if(dirLine==3)
			       ismn=2;
			     if(qcd_NORM(g5sig[dirLine][ismn][id2][id3])>1e-4){
			       vfun_t_WL[iz][id1][id4][ic1][ic4] = qcd_CADD(vfun_t_WL[iz][id1][id4][ic1][ic4],
									    qcd_CMUL(lxr[iz],g5sig[dirLine][ismn][id2][id3]));
			       vfun_t_WL[Nz+iz][id1][id4][ic1][ic4] = qcd_CADD(vfun_t_WL[Nz+iz][id1][id4][ic1][ic4],
									       qcd_CMUL(lxr[Nz+iz],g5sig[dirLine][ismn][id2][id3]));
			     }
			     //=====================================//
			   }
			 }
	   }

   if(myid==0) printf("local vertex functions calculated\n");

   // !!!!!!!!!!!! Here I reinterpriter +z<->-z when I print the results
   // This is a convention and the same must be used in h matrix element

   //global sums and output into files

   //============ print vector ==========//
   MPI_Reduce(&(vfun_v_WL[0][0][0][0][0].re), &(vfun_tmp_WL[0][0][0][0][0].re), 2*Nz*4*4*3*3*2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   if(myid==0)
   for(int iz = 0 ; iz < Nz; iz++)  
   for(id1=0; id1<4; id1++)
   for(id4=0; id4<4; id4++)
   for(ic1=0; ic1<3; ic1++)
   for(ic4=0; ic4<3; ic4++)
     fprintf(fp_vfun_v_WL,"%i \t %i %i %i %i \t %+e %+e %+e %+e\n",iz,id1,id4,ic1,ic4,
	      vfun_tmp_WL[Nz+iz][id1][id4][ic1][ic4].re/geo.V,vfun_tmp_WL[Nz+iz][id1][id4][ic1][ic4].im/geo.V,
	      vfun_tmp_WL[iz][id1][id4][ic1][ic4].re/geo.V,vfun_tmp_WL[iz][id1][id4][ic1][ic4].im/geo.V);

   //============ print axial ==========//
   MPI_Reduce(&(vfun_a_WL[0][0][0][0][0].re), &(vfun_tmp_WL[0][0][0][0][0].re), 2*Nz*4*4*3*3*2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   if(myid==0)
   for(int iz = 0 ; iz < Nz; iz++)  
   for(id1=0; id1<4; id1++)
   for(id4=0; id4<4; id4++)
   for(ic1=0; ic1<3; ic1++)
   for(ic4=0; ic4<3; ic4++)
     fprintf(fp_vfun_a_WL,"%i \t %i %i %i %i \t %+e %+e %+e %+e\n",iz,id1,id4,ic1,ic4,
	      vfun_tmp_WL[Nz+iz][id1][id4][ic1][ic4].re/geo.V,vfun_tmp_WL[Nz+iz][id1][id4][ic1][ic4].im/geo.V,
	      vfun_tmp_WL[iz][id1][id4][ic1][ic4].re/geo.V,vfun_tmp_WL[iz][id1][id4][ic1][ic4].im/geo.V);

   //============ print tensor ==========//
   MPI_Reduce(&(vfun_t_WL[0][0][0][0][0].re), &(vfun_tmp_WL[0][0][0][0][0].re), 2*Nz*4*4*3*3*2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   if(myid==0)
   for(int iz = 0 ; iz < Nz; iz++)  
   for(id1=0; id1<4; id1++)
   for(id4=0; id4<4; id4++)
   for(ic1=0; ic1<3; ic1++)
   for(ic4=0; ic4<3; ic4++)
     fprintf(fp_vfun_t_WL,"%i \t %i %i %i %i \t %+e %+e %+e %+e\n",iz,id1,id4,ic1,ic4,
	      vfun_tmp_WL[Nz+iz][id1][id4][ic1][ic4].re/geo.V,vfun_tmp_WL[Nz+iz][id1][id4][ic1][ic4].im/geo.V,
	      vfun_tmp_WL[iz][id1][id4][ic1][ic4].re/geo.V,vfun_tmp_WL[iz][id1][id4][ic1][ic4].im/geo.V);

   if(myid==0){
     fclose(fp_vfun_v_WL);
     fclose(fp_vfun_a_WL);
     fclose(fp_vfun_t_WL);
   }


   free(vfun_v_WL);
   free(vfun_a_WL);
   free(vfun_t_WL);
   free(vfun_tmp_WL);
   free(lxr);
   

   qcd_destroyWilsonLine(&Wl_fwd);
   qcd_destroyWilsonLine(&Wl_bwd);
   qcd_destroyPropagator(&lprop);
   qcd_destroyPropagator(&rprop);
   qcd_destroyPropagator(&prop);
   qcd_destroyGaugeField(&u);
   qcd_destroyGeometry(&geo);
   MPI_Finalize();
}//end main 
