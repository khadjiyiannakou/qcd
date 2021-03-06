/* Proton nucleon three point function using all to all Propagator using TSM
 * This program takes the forward propagators                                                    
 * And the stochastic sources with the smeared stochastic solution vector
 * Then we calculate the three point function in momentum space
 * Kyriakos Hadjiyiannakou 2011                                                                 
 ****************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>
#include "projectors.h"

#define Nmu 4
#define Nic 3
#define NTERM 8
#define MAX_NUM_PHI 100000
#define MAX_NUM_CHI 1000


int main(int argc, char* argv[]){

  qcd_uint_2 mu,nu,ku,lu,lambda,xita,beta,rho,phita, delta, sigma; // dirac indices
  qcd_uint_2 c1,c2,c3,c4,c5,c6,tt;      // color indices
  qcd_uint_2 cc1,cc2;                   // two more color indices
  qcd_uint_2 i,j,k,l,m,it;                  // for temporal use only
  qcd_uint_2 num_momenta;                // store the number of total momenta
  qcd_int_4 gamma_index;               // this index takes the gamma matrix fo insertion operator
  FILE *fp_momlist;                     // file pointer to the list of momenta
  FILE *fp_threep;                      // file pointer to the output file where i will store the thee point functions
  int params_len;                       // the string legth of the parameters
  char *params;                         // pointer to the params
  qcd_uint_2 aa,bb,cc,ff,gg,hh;         // also color indices
  qcd_int_4 imom;                       // index for momenta

  char gauge_name[qcd_MAX_STRING_LENGTH];       // store the path of gauge configuration
  char param_name[qcd_MAX_STRING_LENGTH];      // name of parameter file                        
  char momlist_name[qcd_MAX_STRING_LENGTH];    // name of momenta-list file                     
  char uprop_name[qcd_MAX_STRING_LENGTH];      // file with the filenames of up propagators   
  char dprop_name[qcd_MAX_STRING_LENGTH];      // file with the filenames of down propagators
  char output_name[qcd_MAX_STRING_LENGTH];     // name of output file

  char phi_name[qcd_MAX_STRING_LENGTH];        // file with the filenames of phi low precision solution vector
  char xi_1_name[qcd_MAX_STRING_LENGTH];       // file list of the noise vector that gives the low precision solution vector
  char chi_name[qcd_MAX_STRING_LENGTH];        // file with list of the chi solution vector (phi_high - phi_low)
  char xi_2_name[qcd_MAX_STRING_LENGTH];       // file list of the noise vector that gives high and low precision solution vector chi

  FILE *fp_phi_name,*fp_xi_1_name,*fp_chi_name,*fp_xi_2_name; // file pointer to the different input files

  char phi_files[MAX_NUM_PHI][qcd_MAX_STRING_LENGTH];        // for each number of phi solution vector we store a string for the path
  char xi_1_files[MAX_NUM_PHI][qcd_MAX_STRING_LENGTH];       // same as xi_1
  char chi_files[MAX_NUM_CHI][qcd_MAX_STRING_LENGTH];        // same as chi
  char xi_2_files[MAX_NUM_CHI][qcd_MAX_STRING_LENGTH];       // same as xi_2

  qcd_uint_4 num_phi, num_xi_1, num_chi, num_xi_2; // the num_phi == num_xi_1, num_chi == num_xi_2 check this must be equal 
  qcd_uint_4 t_sink, t_start, t_stop, t,lt;        // i will use only t for global time and lt for local time
  qcd_complex_16 phase_factor;                     // complex number where we store the phase for fourier transform
  qcd_uint_4 v3,v;                                 // variables for v3 = lx*ly*lz , v = lx*ly*lz*lt
  qcd_uint_2 L[4];                                 // store global lattice dimensions
  qcd_uint_2 P[4];                                 // store the number of proccesors in each direction
  qcd_uint_4 x_src[4],t_cur;                       // x_src the coordinates of source and t_c time coordinate for insertion
  qcd_propagator uprop;                            // for up propagator
  qcd_propagator dprop;                           // for down propagator
  qcd_int_4 x,y,z;                                 // indices for spatial coordinates
  qcd_real_8 tmp;                                  // temporary complex number
  qcd_uint_4 lx,ly,lz;                             // local lattice variables
  qcd_vector phi,xi_1,chi,xi_2;                    // solution and noice vectors

  qcd_geometry geo;                                // a variable to start geometry

  qcd_int_4 (*mom)[3];                             // a list where we will store the momenta
  qcd_real_8 theta[4] = {M_PI,0.0,0.0,0.0};    // antiperiodic b.c. in time take back normal fermion fields                                                                                    
  qcd_real_8 plaq;                             // calculate the plaquette                    
  int myid,numprocs, namelen;                  // for mpi use
  char processor_name[MPI_MAX_PROCESSOR_NAME]; //mpi use

  // smearing variables
  qcd_vector vec;                             // temporary vector using for smearing
  qcd_gaugeField u, uAPE;                     // store gauge field before and after smearing
  qcd_gaugeField *u_ptr, *uAPE_ptr, *utmp_ptr; // some pointer to gauge fields
  qcd_uint_4 nsmear, nsmearAPE;               // smearing parameters
  qcd_real_8 alpha, alphaAPE;                  // smearing parameters

  qcd_int_4 ctr, ctr2, ctr11, ctr22;           // special spin indices trick for gammas multiplications  
  qcd_int_2 cg5cg5_ind[16*16][4], gamma_ind[16][2]; // store the non zero indices and gamma values
  qcd_complex_16 cg5cg5_val[16*16], gamma_val[16];  // %

  qcd_complex_16 z1, z2;                       // temp complex variables                                                                                                                                     
  qcd_complex_16 C, C2;                        // temp complex variables                                                                           
  qcd_complex_16 *W1[Nmu][Nic], *W2[Nmu][Nic], *W3[Nmu][Nmu][Nmu][Nic], *W5[Nmu][Nmu][Nmu][Nic], *W7[Nmu][Nic]; //W8=W6=W4=W2         W terms for contractions    
  qcd_complex_16 *W1p[Nmu][Nic], *W2p[Nmu][Nic], *W3p[Nmu][Nmu][Nmu][Nic], *W5p[Nmu][Nmu][Nmu][Nic], *W7p[Nmu][Nic]; //W8=W6=W4=W2    W terms after fourier
  qcd_complex_16 *W1p_r[Nmu][Nic], *W2p_r[Nmu][Nic], *W3p_r[Nmu][Nmu][Nmu][Nic], *W5p_r[Nmu][Nmu][Nmu][Nic], *W7p_r[Nmu][Nic];        // W after reduce to root
  qcd_complex_16 *W12[Nmu][Nmu], *W34[Nmu][Nmu], *W56[Nmu][Nmu], *W78[Nmu][Nmu];                                                    // after do all the color-spin sums
  qcd_complex_16 *threep[Nmu][Nmu];                                                                                                 // store three point function
 
 //set up MPI                                                                                 
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);           
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);            
  MPI_Get_processor_name(processor_name,&namelen);  
  
  if(argc != 2)
    {
      if(myid == 0) fprintf(stderr,"No input file specified \n"); // check for input file
      exit(EXIT_FAILURE);
    }
  
  strcpy(param_name, argv[1]); // copy the name of input file to a character variable

  if(myid == 0 )
    {
      i=0; // boolean check
      printf("Opening input file for reading parameters %s \n", param_name);
      params = qcd_getParams( param_name, &params_len); // retrieve params in string and params_len only for root
      if(params == NULL) i=1; // if fail to take params return 1
    }

  MPI_Bcast(&i, 1, MPI_INT, 0, MPI_COMM_WORLD); // broadcast success or fail to all
  if(i == 1) exit(EXIT_FAILURE); 

  MPI_Bcast(&params_len, 1, MPI_INT, 0, MPI_COMM_WORLD); // broadcast params_len to all

  if(myid != 0) params = (char*)malloc(params_len*sizeof(char)); // allocate memory for params for all exect root beacause a previous routine do that for root

  MPI_Bcast(params, params_len, MPI_CHAR, 0, MPI_COMM_WORLD); // Broadcast the params to all processors

  /////////// reading from strings for all processors //////////
  
  sscanf(qcd_getParam("<processors_txyz>",params,params_len),"%hd %hd %hd %hd",&P[0], &P[1], &P[2], &P[3]); // read number of processors in each direction
  sscanf(qcd_getParam("<lattice_txyz>",params,params_len),"%hd %hd %hd %hd",&L[0], &L[1], &L[2], &L[3]);    // read the dimensions of global lattice

  if(P[0] != 1)
    {
      if( myid == 0) fprintf(stderr, " Error! , Number of processors in t-direction must be always 1"); // only use one processor in time direction
      exit(EXIT_FAILURE);
    }

  if(qcd_initGeometry(&geo, L, P, theta, myid, numprocs)) exit(EXIT_FAILURE);                           // start the geometry each processor takes a segment of lattice
  if(myid == 0) printf( " Local Lattice : %i x %i x %i x %i \n", geo.lL[0], geo.lL[1], geo.lL[2], geo.lL[3]); // print the local lattice

  strcpy(output_name,qcd_getParam("<output_name>",params,params_len)); //get  output filename                                                                
  if(myid==0) printf("Got output file name: %s\n",output_name);

  strcpy(uprop_name,qcd_getParam("<propagator_u>",params,params_len)); //get  u propagator file list
  if(myid==0) printf("Got propagator file name: %s\n",uprop_name);

  strcpy(dprop_name,qcd_getParam("<propagator_d>",params,params_len)); //get d propagator file list
  if(myid==0) printf("Got propagator file name: %s\n",dprop_name);

  strcpy(momlist_name,qcd_getParam("<momenta_list>",params,params_len)); //get the filename for the momenta list
  if(myid==0) printf("Got momenta-list file name: %s\n",momlist_name);

  strcpy(phi_name,qcd_getParam("<phi_list>",params,params_len)); //get the filename of phi list                                                       
  if(myid==0) printf("Got phi-list file name: %s\n",phi_name);
  sscanf(qcd_getParam("<num_phi>",params,params_len),"%d",&num_phi); // read number of phi vectors 
  if(myid == 0) printf("Number of phi vectors is: %d\n",num_phi);

  strcpy(xi_1_name,qcd_getParam("<xi_1_list>",params,params_len)); //get the filename of xi_1 list                                   
  if(myid==0) printf("Got xi_1-list file name: %s\n",xi_1_name);
  sscanf(qcd_getParam("<num_xi_1>",params,params_len),"%d",&num_xi_1); // read number of xi_1 vectors                                                                  
  if(myid == 0)printf("Number of xi_1 vectors is: %d\n",num_xi_1);

  if( num_phi != num_xi_1 )
    {
      if(myid == 0)fprintf(stderr,"Error, num_phi must be the same with num_xi_1");
      exit(EXIT_FAILURE);
    }
  
  strcpy(chi_name,qcd_getParam("<chi_list>",params,params_len)); //get the filename of chi list                                                                        
  if(myid==0) printf("Got chi-list file name: %s\n",chi_name);
  sscanf(qcd_getParam("<num_chi>",params,params_len),"%d",&num_chi); // read number of chi vectors                                                                     
  if(myid == 0) printf("Number of chi vectors is: %d\n",num_chi);

  strcpy(xi_2_name,qcd_getParam("<xi_2_list>",params,params_len)); //get the filename of xi_2 list                                                                     
  if(myid==0) printf("Got xi_2-list file name: %s\n",xi_2_name);
  sscanf(qcd_getParam("<num_xi_2>",params,params_len),"%d",&num_xi_2); // read number of xi_2 vectors                                                                  
  if(myid == 0)printf("Number of xi_2 vectors is: %d\n",num_xi_2);

  if( num_chi != num_xi_2 )
    {
      if(myid == 0)fprintf(stderr,"Error, num_chi must be the same with num_xi_2");
      exit(EXIT_FAILURE);
    }

  strcpy(gauge_name,qcd_getParam("<cfg_name>",params,params_len));                // read gauge filename
  if(myid==0) printf("Got conf name: %s\n",gauge_name);
  
  sscanf(qcd_getParam("<alpha_gauss>",params,params_len),"%lf",&alpha);           // read alpha parameter for wuppertal smearing
  if(myid==0) printf(" Got alpha_gauss: %lf\n",alpha);
  
  sscanf(qcd_getParam("<nsmear_gauss>",params,params_len),"%d",&nsmear);          // number of iterations for wuppertal smearing
  if(myid==0) printf(" Got nsmear_gauss: %d\n",nsmear);
  
  sscanf(qcd_getParam("<alpha_APE>",params,params_len),"%lf",&alphaAPE);         // read alpha parameter for APE smearing
  if(myid==0) printf(" Got alpha_APE: %lf\n",alphaAPE);
  
  sscanf(qcd_getParam("<nsmear_APE>",params,params_len),"%d",&nsmearAPE);         // number of iterations for APE smearing
  if(myid==0) printf(" Got nsmear_APE: %d\n",nsmearAPE);

  sscanf(qcd_getParam("<source_pos_txyz>",params,params_len),"%d %d %d %d",&x_src[0],&x_src[1],&x_src[2],&x_src[3]);  // read the position where we put the source
  if(myid==0) printf("Got source coords: %d %d %d %d\n",x_src[0],x_src[1],x_src[2],x_src[3]);

  sscanf(qcd_getParam("<current_pos_t>",params,params_len),"%d",&t_cur);                                              // read the time position of insertion
  if(myid==0) printf("Got noise vector time-slice dilution: %d \n",t_cur);

  sscanf(qcd_getParam("<gamma_index>",params,params_len),"%d",&gamma_index);                                         // index to choose the gamma insertion matrix
  if(myid==0) printf("Got index for gamma matrix : %d \n",gamma_index);


  free(params); 

  //////////////////////// finish with parameters //////////////////////////////////////////////////

  //allocate memory for gauge-fields////////////////////////////////////////// no comments i just copied it from twop.c
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
  if(myid == 0 ) printf("start APE smearing\n");
  for(i=0; i<nsmearAPE; i++)
    {
      qcd_apeSmear3d(uAPE_ptr, u_ptr, alphaAPE);
      utmp_ptr=u_ptr; u_ptr=uAPE_ptr; uAPE_ptr=utmp_ptr;
    }
  utmp_ptr=u_ptr; u_ptr=uAPE_ptr; uAPE_ptr=utmp_ptr; //reverse the last swap. Also needed when nsmearAPE=0                                                                                           
  qcd_destroyGaugeField(u_ptr);
  uAPE = *uAPE_ptr;
  if(myid == 0) printf("printf finish APE smearing\n");
  //////////////////////////////////////////////////////////////////////////////////

  ///////////////////// start initialize the resources ///////////////////// only for propagator and vectors

  j = 0; // flag for ok allocation
  j += qcd_initPropagator(&uprop, &geo); // initialize memory for u propagator
  j += qcd_initPropagator(&dprop, &geo); // initialize memory for d propagator
  j += qcd_initVector(&vec, &geo);

  j += qcd_initVector(&phi, &geo);       // initialize phi vector
  j += qcd_initVector(&xi_1, &geo);      // initialize xi_1 vector
  j += qcd_initVector(&chi, &geo);       // initialize chi vector
  j += qcd_initVector(&xi_2, &geo);      // initialize xi_2 vector

  MPI_Allreduce(&j, &k, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); // readuce k for all to see if the allocation is ok!

  if(k>0) // check for memory allocation flag
    {
      if(myid==0) printf("not enough memory\n");
      exit(EXIT_FAILURE);
    }
  if(myid==0) printf("memory for propagators and tsm vectors allocated\n");

  //############################################################///////////////////

  // load propagators

  if(qcd_getPropagator(uprop_name,qcd_PROP_LIME, &uprop)) exit(EXIT_FAILURE); // read u propagators from binary in ram
  if(myid==0) printf("up propagator loaded\n");
  if(qcd_getPropagator(dprop_name,qcd_PROP_LIME, &dprop)) exit(EXIT_FAILURE); // read d propagators from binary in ram
  if(myid==0) printf("down propagator loaded\n");

  //################################################################################                                                                                                                 
  // transform propagators to basis with theta-periodic boundaries in the temporal direction  !!! remember to do the same fo solution vectors
  for(int lt=0; lt<geo.lL[0]; lt++)
    {
      t = lt + geo.Pos[0] * geo.lL[0];
      phase_factor   = (qcd_complex_16) {cos(theta[0]*t/geo.L[0]),sin(theta[0]*t/geo.L[0])};
      qcd_mulPropagatorC3d(&uprop, phase_factor, (t+x_src[0]) % geo.L[0]);                         // !!!!!!!!!!! ask giannis why he use  (t+x_src[0]) % geo.L[0])
      qcd_mulPropagatorC3d(&dprop, phase_factor, (t+x_src[0]) % geo.L[0]);
    }
  if(myid==0) printf("propagators transformed to basis with theta-periodic boundary conditions\n");

  //################################################################################                                                                                                                 
  // gaussian smearing of propagators (only the time-slices that will be used)                                                                                                                       
    
  for(mu=0;mu<4;mu++)
    for(c1=0;c1<3;c1++)
      {

			qcd_copyVectorPropagator(&vec,&uprop,mu,c1);             // copy each of 12 vectors for propagator to do smearing
		
	for(i=0; i<nsmear; i++)
	  {
	    
	    if(qcd_gaussIteration3dAll(&vec,&uAPE,alpha,i==0))   // smearing for up propagator
	      {
		fprintf(stderr,"process %i: Error while smearing!\n",geo.myid);
		exit(EXIT_FAILURE);
	      }
	    
	  }
		qcd_copyPropagatorVector(&uprop,&vec,mu,c1);
		qcd_copyVectorPropagator(&vec,&dprop,mu,c1);
	for(i=0; i<nsmear; i++)
	  {
	   
	    if(qcd_gaussIteration3dAll(&vec,&uAPE,alpha,i==0))        // smearing for down propagator
	      {
		fprintf(stderr,"process %i: Error while smearing!\n",geo.myid);
		exit(EXIT_FAILURE);
	      }
	    
	  }
		qcd_copyPropagatorVector(&dprop,&vec,mu,c1);
		
      }
  
  qcd_destroyVector(&vec);
                                               // destroy only the temporary vector i dont nead any more
  if(myid==0) printf("propagators smeared\n");
  /////////////////////////////////////////////////////////////////////////////////////////
      
  //calculate the dirac non zero value and indices
  ctr = 0;
  for(beta = 0; beta < 4 ; beta++)                // this are only for C * Gamma5 * \bar{C * Gamma5}
    for(nu = 0; nu < 4 ; nu++)
      for(lambda = 0 ; lambda < 4 ; lambda++)
	for(xita = 0 ; xita < 4 ; xita++)
	  {
	    C = qcd_CMUL(qcd_CGAMMA[5][beta][nu],qcd_BAR_CGAMMA[5][lambda][xita]);
	    if(qcd_NORM(C) > 1e-3)
	      {
		cg5cg5_val[ctr].re = C.re;                      // store values
		cg5cg5_val[ctr].im = C.im;
		cg5cg5_ind[ctr][0] = beta;                      // store indices
		cg5cg5_ind[ctr][1] = nu;
		cg5cg5_ind[ctr][2]= lambda;
		cg5cg5_ind[ctr][3]= xita;
		ctr++;                                         // important dont forget to increase this index
	      }
	  }
  
  ctr11 = 0 ;
  for(rho =0; rho < 4 ; rho++)                               // now for gamma inserion matrix only
    for(phita=0; phita < 4 ; phita++){
      C = qcd_GAMMA[gamma_index][rho][phita] ;
      if(qcd_NORM(C) > 1e-3)
	{
	  gamma_val[ctr11].re = C.re;
	  gamma_val[ctr11].im = C.im;
	  gamma_ind[ctr11][0] = rho;
	  gamma_ind[ctr11][1] = phita;
	  ctr11++;
	}
    }
 
  //open output file to write in and file for momenta                                                                                                            
  if(myid==0)
    {
      fp_threep = fopen(output_name,"w");
      if(fp_threep==NULL)
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
  if(k==1) exit(EXIT_FAILURE); // check fo success
  //////////////////////////////

  // load the momenta list///////////////////////////////
  if(myid==0) fscanf(fp_momlist,"%i\n",&num_momenta);
  MPI_Bcast(&num_momenta,1,MPI_INT, 0, MPI_COMM_WORLD);
  if(myid==0) printf("will read %i momenta combinations\n",num_momenta);
  mom = malloc(num_momenta*3*sizeof(qcd_int_4));
  if(myid==0)
    {
      for(j=0; j<num_momenta; j++)
	{
	  fscanf(fp_momlist,"%i %i %i\n",&(mom[j][0]),&(mom[j][1]),&(mom[j][2]));
	  //printf("got combination %i %i %i\n",mom[j][0],mom[j][1],mom[j][2]);                                                                                                                        
	}
      fclose(fp_momlist);
    }
  MPI_Bcast(&(mom[0][0]),num_momenta*3,MPI_INT,0, MPI_COMM_WORLD);
  if(myid==0) printf("momenta list read and broadcasted\n");
  ////////////////////// succesfuly pass the momenta to all cpu

  //allocate memory for all the W terms ///////////////////////////////////////////////
  for(i=0; i < Nmu ; i++)
    for(j=0; j < Nic ; j++){
      W1[i][j] = (qcd_complex_16*)malloc(geo.lV3 * sizeof(qcd_complex_16));
      W2[i][j] = (qcd_complex_16*)malloc(geo.lV3 * sizeof(qcd_complex_16));
      W7[i][j] = (qcd_complex_16*)malloc(geo.lV3 * sizeof(qcd_complex_16));
      W1p[i][j] = (qcd_complex_16*)malloc(num_momenta * sizeof(qcd_complex_16));
      W2p[i][j] =(qcd_complex_16*)malloc(num_momenta *sizeof(qcd_complex_16));
      W7p[i][j] =(qcd_complex_16*)malloc(num_momenta *sizeof(qcd_complex_16));
      W1p_r[i][j] =(qcd_complex_16*)malloc(num_momenta *sizeof(qcd_complex_16));
      W2p_r[i][j] =(qcd_complex_16*)malloc(num_momenta *sizeof(qcd_complex_16));
      W7p_r[i][j] =(qcd_complex_16*)malloc(num_momenta *sizeof(qcd_complex_16));
    }

  for( i=0; i < Nmu ; i++)
    for( j=0; j<Nmu ; j++)
      for(k=0; k<Nmu ; k++)
	for(l=0; l<Nic ; l++){
	  W3[i][j][k][l] = (qcd_complex_16*)malloc(geo.lV3 * sizeof(qcd_complex_16));
	  W5[i][j][k][l] = (qcd_complex_16*)malloc(geo.lV3 * sizeof(qcd_complex_16));
	  W3p[i][j][k][l] = (qcd_complex_16*)malloc(num_momenta * sizeof(qcd_complex_16));
	  W5p[i][j][k][l] =(qcd_complex_16*)malloc(num_momenta * sizeof(qcd_complex_16));
	  W3p_r[i][j][k][l] = (qcd_complex_16*)malloc(num_momenta * sizeof(qcd_complex_16));
	  W5p_r[i][j][k][l] =(qcd_complex_16*)malloc(num_momenta * sizeof(qcd_complex_16));
	}

  for( i = 0 ; i < Nmu ; i++)
    for(j = 0 ; j <Nmu ; j++)
      {
	W12[i][j] = (qcd_complex_16*)malloc(num_momenta*num_momenta*geo.L[0]);
	W34[i][j] = (qcd_complex_16*)malloc(num_momenta*num_momenta*geo.L[0]);
	W56[i][j] = (qcd_complex_16*)malloc(num_momenta*num_momenta*geo.L[0]);
	W78[i][j] = (qcd_complex_16*)malloc(num_momenta*num_momenta*geo.L[0]);
	threep[i][j] = (qcd_complex_16*)malloc(num_momenta*num_momenta*geo.L[0]);
      }

  for( i = 0 ; i < Nmu ; i++)
    for(j =0 ; j < Nmu ; j++)
      for(k=0; k < num_momenta*num_momenta*geo.L[0] ; k++)
	{
	  threep[i][j][k] = (qcd_complex_16) {0,0}; // initialize to zero the three point function
	}
  ///////////////////////////////// finish with the allocation for W terms

  // read the names of files from a list////////////////////////////////////////////////////////////
  fp_phi_name = fopen(phi_name,"r");
  fp_xi_1_name = fopen(xi_1_name,"r");
  fp_chi_name = fopen(chi_name, "r");
  fp_xi_2_name = fopen(xi_2_name,"r");

  if(fp_phi_name == NULL){
    if(myid == 0)fprintf(stderr,"Error opening list with phi vectors\n");
    exit(EXIT_FAILURE);
  }
  if(fp_xi_2_name == NULL){
    if(myid == 0)fprintf(stderr,"Error opening list with xi_2 vectors\n");
    exit(EXIT_FAILURE);
  }
  if(fp_xi_1_name == NULL){
    if(myid == 0)fprintf(stderr,"Error opening list with xi_1 vectors\n");
    exit(EXIT_FAILURE);
  }
  if(fp_chi_name == NULL){
    if(myid == 0)fprintf(stderr,"Error opening list with chi vectors\n");
    exit(EXIT_FAILURE);
  }

  imom = -1;
  while(!feof(fp_phi_name))
    {
      imom++;
      fscanf(fp_phi_name,"%s",&(phi_files[imom][0]));
    }
  fclose(fp_phi_name);
  if(imom!=num_phi){
    if(myid == 0)fprintf(stderr,"Error read right number of phi vector names ");
    exit(EXIT_FAILURE);
  }

  imom = -1;
  while(!feof(fp_xi_1_name))
    {
      imom++;
      fscanf(fp_xi_1_name,"%s",&(xi_1_files[imom][0]));
    }
  fclose(fp_xi_1_name);
  if(imom!=num_xi_1){
    if(myid == 0)fprintf(stderr,"Error read right number of xi_1 vector names ");
    exit(EXIT_FAILURE);
  }

  imom = -1;
  while(!feof(fp_chi_name))
    {
      imom++;
      fscanf(fp_chi_name,"%s",&(chi_files[imom][0]));
    }
  fclose(fp_chi_name);
  if(imom!=num_chi){
    if(myid == 0)fprintf(stderr,"Error read right number of chi vector names ");
    exit(EXIT_FAILURE);
  }

  imom = -1;
  while(!feof(fp_xi_2_name))
    {
      imom++;
      fscanf(fp_xi_2_name,"%s",&(xi_2_files[imom][0]));
    }
  fclose(fp_xi_2_name);
  if(imom!=num_xi_2){
    if(myid == 0)fprintf(stderr,"Error read right number of xi_2 vector names ");
    exit(EXIT_FAILURE);
  }





  //###################################### MAIN CALCULATIONS ####################################################
  for(int r = 0 ; r < num_phi ; r++){        // a summation over the ensemble of noise vectors

    for( i = 0 ; i < Nmu ; i++)
      for(j = 0 ; j <Nmu ; j++)
	for(k = 0 ; k < num_momenta*num_momenta*geo.L[0] ; k++)
	{
	  W12[i][j][k] = (qcd_complex_16) {0,0};         // for each noise vector we have to set the values of them to zero
	  W34[i][j][k] = (qcd_complex_16) {0,0};
	  W56[i][j][k] = (qcd_complex_16) {0,0};
	  W78[i][j][k] = (qcd_complex_16) {0,0};
	}

    //load the current phi and xi_1
    qcd_getVector(&(phi_files[r][0]),qcd_PROP_LIME,0,0,&phi);
    qcd_getVector(&(xi_1_files[r][0]),qcd_PROP_LIME,0,0,&xi_1);

    //we have to multiply phi by the exponential with theta
    for(int llt=0; llt<geo.lL[0]; llt++)
      {
	t = llt + geo.Pos[0] * geo.lL[0];
	phase_factor   = (qcd_complex_16) {cos(theta[0]*t/geo.L[0]),sin(theta[0]*t/geo.L[0])};
	qcd_mulVectorC3d(&phi, phase_factor, (t+t_cur) % geo.L[0]);                               // !!!!!!!!!!!!!!!!!!!!!! again ask giannis about that
      }
    
    //we have to do smearing on the phi
    for(i=0; i<nsmear; i++)
      {
	if(qcd_gaussIteration3dAll(&phi,&uAPE,alpha,i==0))
	  {
	    fprintf(stderr,"process %i: Error while smearing!\n",geo.myid);
	    exit(EXIT_FAILURE);
	  }
      }
    
    for( it = t_cur + 1 ; it < (t_cur + 1 + 10) ; it++)                                          // !!!!!!!!!!!!!!!!!!!!!!!!!! ask giannis about lattice close ring
      {
	if( it > L[0] ){
	  t = it - geo.L[0];
	}
	else {
	  t = it;
	}

	
	if(myid == 0)printf("t = %i \n",t);
	    
	// set to zero the W // this variables must set to zero for each time slice
	for(i=0; i < Nmu ; i++)
	  for(j=0; j < Nic ; j++)
	    for(v=0; v<geo.lV3; v++){
	      W1[i][j][v] = (qcd_complex_16) {0,0};
	      W2[i][j][v] = (qcd_complex_16) {0,0};
	      W7[i][j][v] = (qcd_complex_16) {0,0};
	    }
	for( i=0; i < Nmu ; i++)
	  for( j=0; j<Nmu ; j++)
	    for(k=0; k<Nmu ; k++)
	      for(l=0; l<Nic ; l++)
		for(v=0; v<geo.lV3 ; v++){
		  W3[i][j][k][l][v] = (qcd_complex_16) {0,0};
		  W5[i][j][k][l][v] = (qcd_complex_16) {0,0};
		}

	    ///////////////////////////////////////////////////////////////////////

	for(lx = 0 ; lx < geo.lL[1] ; lx++)                  // each processors do his work for its local lattice
	  for(ly = 0 ; ly< geo.lL[2] ; ly++)
	    for(lz = 0 ; lz < geo.lL[3] ; lz++){
	      v3 = qcd_LEXIC0(lx,ly,lz,geo.lL);
	      v =  qcd_LEXIC(lt,lx,ly,lz,geo.lL);
	    

	      for(ctr2 = 0; ctr2 < ctr ; ctr2++){
		beta = cg5cg5_ind[ctr2][0];
		nu = cg5cg5_ind[ctr2][1];
		lambda = cg5cg5_ind[ctr2][2];
		xita = cg5cg5_ind[ctr2][3];
	   
	      
		for(cc1 = 0 ; cc1 < 6 ; cc1++){
		  aa = qcd_EPS[cc1][0];
		  bb = qcd_EPS[cc1][1];
		  cc = qcd_EPS[cc1][2];
		      
		  for(cc2 = 0 ; cc2 < 6 ; cc2++){
		    ff = qcd_EPS[cc2][0];
		    gg = qcd_EPS[cc2][1];
		    hh = qcd_EPS[cc2][2];
		    
			// for W1 W7 
		    for(delta = 0 ; delta < Nmu ; delta++){
		      W1[delta][hh][v3] = qcd_CADD(W1[delta][hh][v3],qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( cg5cg5_val[ctr2] , 
													    qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),
												phi.D[v][delta][cc]),uprop.D[v][beta][lambda][aa][ff]), 
									      dprop.D[v][nu][xita][bb][gg]));
		      W7[delta][hh][v3] = qcd_CADD(W7[delta][hh][v3],qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( cg5cg5_val[ctr2] ,
													    qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),
												phi.D[v][beta][aa]),uprop.D[v][delta][lambda][cc][ff]),
									      dprop.D[v][nu][xita][bb][gg]));

		    }
			// for W3 W5
		    for(delta = 0 ; delta<Nmu ; delta++)
		      for(sigma = 0 ;sigma<Nmu; sigma++){
			W3[delta][sigma][lambda][ff][v3] = qcd_CADD(W3[delta][sigma][lambda][ff][v3],qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( cg5cg5_val[ctr2] ,
																	    qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),
																phi.D[v][delta][cc]),uprop.D[v][beta][sigma][aa][hh]),
													      dprop.D[v][nu][xita][bb][gg]));
			W5[delta][sigma][lambda][ff][v3] = qcd_CADD(W5[delta][sigma][lambda][ff][v3],qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( cg5cg5_val[ctr2] ,
																	    qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),
																phi.D[v][beta][aa]),uprop.D[v][delta][sigma][cc][hh]),
													      dprop.D[v][nu][xita][bb][gg]));
		      }

			//
			
		      
		  } // close color
		} //close color
	      } //close dirac indices
		  
	      // for W2 term
	      for(tt=0 ; tt< Nic ;tt++)
		for( ctr22 =0 ; ctr22 < ctr11 ; ctr22++){
		  rho = gamma_ind[ctr22][0];
		  phita = gamma_ind[ctr22][1];
		  
		  for(hh=0; hh<Nic ;hh++)
		    for(sigma =0 ;sigma<Nmu ; sigma++)
		      W2[sigma][hh][v3] = qcd_CADD(W2[sigma][hh][v3], qcd_CMUL(gamma_val[ctr22],qcd_CMUL(qcd_CONJ(xi_1.D[v][rho][tt]),uprop.D[v][phita][sigma][tt][hh])));

		}

		} //close spatial

	    // Fourier Transform //////////////////////////////////

	for(i=0; i < Nmu ; i++) // set to zero
	  for(j=0; j < Nic ; j++)
	    for(v=0; v< num_momenta; v++){
	      W1p[i][j][v] = (qcd_complex_16) {0,0};
	      W2p[i][j][v] = (qcd_complex_16) {0,0};
	      W7p[i][j][v] = (qcd_complex_16) {0,0};
	      W1p_r[i][j][v] = (qcd_complex_16) {0,0};
	      W2p_r[i][j][v] = (qcd_complex_16) {0,0};
	      W7p_r[i][j][v] = (qcd_complex_16) {0,0};
	    }

	for( i=0; i < Nmu ; i++) // set to zero
	  for( j=0; j<Nmu ; j++)
	    for(k=0; k<Nmu ; k++)
	      for(l=0; l<Nic ; l++)
		for(v=0; v<num_momenta ; v++){
		  W3p[i][j][k][l][v] = (qcd_complex_16) {0,0};
		  W5p[i][j][k][l][v] = (qcd_complex_16) {0,0};
		  W3p_r[i][j][k][l][v] = (qcd_complex_16) {0,0};
		  W5p_r[i][j][k][l][v] = (qcd_complex_16) {0,0};
		}

	for( j=0 ; j< num_momenta ; j++){ // begin fourier


	  for(lx=0; lx<geo.lL[1]; lx++)
	    for(ly=0; ly<geo.lL[2]; ly++)
	      for(lz=0; lz<geo.lL[3]; lz++)
		{

		  v3 = qcd_LEXIC0(lx,ly,lz,geo.lL);
		  x=lx+geo.Pos[1]*geo.lL[1] - x_src[1];                   /// !!!!!!!!!! ask giannis why substract x_src
		  y=ly+geo.Pos[2]*geo.lL[2] - x_src[2];
		  z=lz+geo.Pos[3]*geo.lL[3] - x_src[3];
		  tmp = (((double) mom[j][0]*x)/geo.L[1] + ((double) mom[j][1]*y)/geo.L[2] + ((double) mom[j][2]*z)/geo.L[3])*2*M_PI;
		  C2=(qcd_complex_16) {cos(tmp),sin(tmp)}; //TABULATE FOR LARGE SPEEDUP!!!                                                                                                   
			  
		  for( hh = 0 ; hh < Nic ; hh++)
		    for( delta = 0 ; delta < Nmu ; delta++){
		      W1p[delta][hh][j]=qcd_CADD(W1p[delta][hh][j], qcd_CMUL(W1[delta][hh][v3],C2));
		      W2p[delta][hh][j]=qcd_CADD(W2p[delta][hh][j], qcd_CMUL(W2[delta][hh][v3],C2));
		      W7p[delta][hh][j]=qcd_CADD(W7p[delta][hh][j], qcd_CMUL(W7[delta][hh][v3],C2));
		      
		      for( mu = 0 ; mu< Nmu ; mu++)
			for( nu = 0 ; nu <Nmu ; nu++){
			  W3p[delta][mu][nu][hh][j] = qcd_CADD(W3p[delta][mu][nu][hh][j], qcd_CMUL(W3[delta][mu][nu][hh][v3],C2));
			  W5p[delta][mu][nu][hh][j] = qcd_CADD(W5p[delta][mu][nu][hh][j], qcd_CMUL(W5[delta][mu][nu][hh][v3],C2));
			}
		    }
		}
	  
	  for( hh = 0 ; hh < Nic ; hh++) // reduce the values to root
	    for(delta = 0 ;delta <Nmu ; delta++){
	      MPI_Reduce(&(W1p[delta][hh][j].re), &(W1p_r[delta][hh][j].re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	      MPI_Reduce(&(W2p[delta][hh][j].re), &(W2p_r[delta][hh][j].re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	      MPI_Reduce(&(W7p[delta][hh][j].re), &(W7p_r[delta][hh][j].re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	      for( mu = 0 ; mu < Nmu ; mu++)
		for(nu = 0 ; nu < Nmu ; nu++){
		  MPI_Reduce(&(W3p[delta][mu][nu][hh][j].re), &(W3p_r[delta][mu][nu][hh][j].re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		  MPI_Reduce(&(W5p[delta][mu][nu][hh][j].re), &(W5p_r[delta][mu][nu][hh][j].re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		}
	    }
	  


	} // close loop for momentum for fourier
	    ////////////////////////////////////////////////////////////////////////
	    
	for( int p = 0 ; p < num_momenta ; p++)
	  for( int q = 0 ; q < num_momenta; q++)
	    for( delta = 0 ; delta < Nmu ; delta++)
	      for( sigma = 0 ; sigma < Nmu ; sigma++){
		
		for( hh = 0 ; hh < Nic ; hh++){
		  W12[delta][sigma][p*num_momenta*geo.L[0]+q*geo.L[0]+t] = qcd_CADD(W12[delta][sigma][p*num_momenta*geo.L[0]+q*geo.L[0]+t],qcd_CMUL(W1p_r[delta][hh][p] , W2p_r[sigma][hh][q]));
		  W78[delta][sigma][p*num_momenta*geo.L[0]+q*geo.L[0]+t] = qcd_CADD(W78[delta][sigma][p*num_momenta*geo.L[0]+q*geo.L[0]+t],qcd_CMUL(W7p_r[delta][hh][p] , W2p_r[sigma][hh][q]));
		  
		  for(lambda = 0 ; lambda < Nmu ; lambda++){
		    W34[delta][sigma][p*num_momenta*geo.L[0]+q*geo.L[0]+t] = qcd_CADD(W3p_r[delta][sigma][lambda][hh][p*num_momenta*geo.L[0]+q*geo.L[0]+t],
										      qcd_CMUL(W3p_r[delta][sigma][lambda][hh][p*num_momenta*geo.L[0]+q*geo.L[0]+t],
											       W2p_r[lambda][hh][q]));
		    W56[delta][sigma][p*num_momenta*geo.L[0]+q*geo.L[0]+t] = qcd_CADD(W5p_r[delta][sigma][lambda][hh][p*num_momenta*geo.L[0]+q*geo.L[0]+t],
										      qcd_CMUL(W5p_r[delta][sigma][lambda][hh][p*num_momenta*geo.L[0]+q*geo.L[0]+t],
											       W2p_r[lambda][hh][q]));
		  }
		  
		}
	      }

	for( int p = 0 ; p < num_momenta ; p++)
	  for( int q = 0 ; q < num_momenta; q++)
	    for( delta = 0 ; delta < Nmu ; delta++)
	      for( sigma = 0 ; sigma < Nmu ; sigma++){
		threep[delta][sigma][p*num_momenta*geo.L[0]+q*geo.L[0]+t] = qcd_CADD(threep[delta][sigma][p*num_momenta*geo.L[0]+q*geo.L[0]+t],qcd_CADD(qcd_CSUB(W12[delta][sigma][p*num_momenta*geo.L[0]+q*geo.L[0]+t],W34[delta][sigma][p*num_momenta*geo.L[0]+q*geo.L[0]+t]),qcd_CSUB(W56[delta][sigma][p*num_momenta*geo.L[0]+q*geo.L[0]+t], W78[delta][sigma][p*num_momenta*geo.L[0]+q*geo.L[0]+t])));
	      }
	 
	
      } //close time
    
  } // close r


	    




	  


	  





  //finish MPI library
  MPI_Finalize();


  return 0;
}
