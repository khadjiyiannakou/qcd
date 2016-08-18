/*
This program construct two point function for all particles for octet and decuplet
for decuplet we project on spin 3/2
 strarting from one time slice after the insertion you fixed in three point function
Editor : Kyriakos Hadjiyiannakou 2012
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>
#include "projectors_pb.h"
#include <gamma_J12.h>

#ifndef MAX_PARTICLES_2
#define MAX_PARTICLES_2 40
#endif

void* allocate(size_t size){
  void *ptr = malloc(size);
  if(ptr == NULL){
    fprintf(stderr,"not enought memory for malloc \n");
  }
  return ptr;
}

int main(int argc, char* argv[]){

  int malloc_flag, malloc_flag_check, input_flag;  // flags to check 
  int indices, gamma_index;                       
  int params_len;                       // the string legth of the parameters
  char *params;                         // pointer to the params
 
  char gauge_name[qcd_MAX_STRING_LENGTH];       // store the path of gauge configuration
  char param_name[qcd_MAX_STRING_LENGTH];      // name of parameter file                        
  char output_name[qcd_MAX_STRING_LENGTH];
  char uprop_name[qcd_MAX_STRING_LENGTH];      // file with the filenames of up propagators   
  char dprop_name[qcd_MAX_STRING_LENGTH];      // file with the filenames of down propagators
  char sprop_name[qcd_MAX_STRING_LENGTH];
  char cprop_name[qcd_MAX_STRING_LENGTH];
  
  char names[qcd_MAX_STRING_LENGTH];                           // temporary string for names of particles
  char momlist_name[qcd_MAX_STRING_LENGTH];    // name of momenta-list file

  qcd_uint_4 v3,v;                                 // variables for v3 = lx*ly*lz , v = lx*ly*lz*lt
  qcd_uint_2 L[4];                                 // store global lattice dimensions
  qcd_uint_2 P[4];                                 // store the number of proccesors in each direction
  qcd_uint_4 x_src[4],cur_time, number_tsinks, after_tcur;                       // x_src the coordinates of source and cur_time time coordinate for insertion

  qcd_uint_4  t,lt;        // i will use only t for global time and lt for local time
  qcd_complex_16 phase_factor, C2, C;                     // complex number where we store the phase for fourier transformed
  
  qcd_propagator uprop;                            // for up propagator
  qcd_propagator dprop;                           // for down propagator
  qcd_propagator sprop;				  // for strange propagator
  qcd_propagator cprop;
  
  qcd_int_4 x,y,z;                                 // indices for spatial coordinates
  qcd_real_8 tmp;                                  // fourier phase real number
  qcd_uint_4 lx,ly,lz;                             // local lattice variables
  qcd_vector vec_up[2], vec_down[2], vec_strange[2], vec_charm[2];        // temporary vectors to store for smearing
  qcd_geometry geo;                                // a variable to start geometry
  qcd_int_4 (*mom)[3];                         // momenta-list
  
   qcd_gaugeField uAPE, u;                         // APE smeared gaugeField
   qcd_gaugeField *u_ptr, *uAPE_ptr, *utmp_ptr;
   qcd_uint_4 nsmear, nsmearAPE;       // gaussian and APE smearing: n
   qcd_real_8 alphaGauss, alphaAPE;         // gaussian and APE smearing: alpha
   
   qcd_real_8 plaq;
   qcd_uint_4 nz_counter[10];
   qcd_int_2 cg5cg5b_ind[10][16*16][4];
   qcd_complex_16 cg5cg5b_val[10][16*16];
   qcd_complex_16 expon;
   
   int dummy;
   qcd_uint_4 num_momenta;
   qcd_complex_16 corr_temp[MAX_PARTICLES_2][16][10], corr[MAX_PARTICLES_2][16][10];
   qcd_complex_16 *block[MAX_PARTICLES_2][16][10];
   qcd_complex_16 TermJ12[3][4][4];
   qcd_complex_16 corr_J12[4][4], Term3[4][4], Term6[4][4];
   qcd_complex_16 corr_proj;
   // we have 9 elements that correspond to Cg_i * \bar{Cg_j} , the final element correspond to Cg_5 * \bar{Cg_5}, 
   //the final element will be used only for the octec, the other 9 elements for octet will be dummy
   
   char *particles[MAX_PARTICLES_2] = {"PROTON", "NEUTRON", "LAMBDA", "SIGMA_PLUS", "SIGMA_ZERO", "SIGMA_MINUS", "KSI_ZERO", "KSI_MINUS", "DELTA_PLUS_PLUS", "DELTA_PLUS", "DELTA_ZERO", "DELTA_MINUS", "SIGMA_STAR_PLUS", "SIGMA_STAR_ZERO", "SIGMA_STAR_MINUS", "KSI_STAR_ZERO", "KSI_STAR_MINUS", "OMEGA", "LAMBDA_C", "SIGMA_PLUS_PLUS_C", "SIGMA_PLUS_C", "SIGMA_ZERO_C", "KSI_PLUS_PLUS_C_C", "KSI_PLUS_C_C", "SIGMA_STAR_PLUS_PLUS_C", "SIGMA_STAR_PLUS_C", "SIGMA_STAR_ZERO_C", "KSI_STAR_PLUS_PLUS_C_C", "KSI_STAR_PLUS_C_C", "OMEGA_PLUS_PLUS_C_C_C", "OMEGA_ZERO_C", "KSI_PRIME_ZERO_C", "KSI_ZERO_C", "KSI_PRIME_PLUS_C", "KSI_PLUS_C", "OMEGA_PLUS_C_C", "OMEGA_STAR_ZERO_C", "KSI_STAR_ZERO_C", "KSI_STAR_PLUS_C", "OMEGA_STAR_PLUS_C_C"};
   
  FILE *fp_momlist, *fp_twop[MAX_PARTICLES_2];
  qcd_real_8 theta[4] = {M_PI,0.0,0.0,0.0};    // antiperiodic b.c. in time take back normal fermion fields           
  int myid,numprocs, namelen;                  // for mpi use
  char processor_name[MPI_MAX_PROCESSOR_NAME]; //mpi use
  
  qcd_uint_2 list_cgamma[4];
  
  list_cgamma[0] = 1 ;  // .> Cg1
  list_cgamma[1] = 2 ;  // .> Cg2
  list_cgamma[2] = 3 ;  // .> Cg3
  list_cgamma[3] = 5 ;  // .> Cg5
  
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);         // num. of processes taking part in the calculation
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);             // each process gets its ID
   MPI_Get_processor_name(processor_name,&namelen); // 

      ////////////////////// check for input file ////////////////////////////////////
    if(argc != 2)
    {
      if(myid == 0) fprintf(stderr,"No input file specified \n"); // check for input file
      exit(EXIT_FAILURE);
    }
  strcpy(param_name, argv[1]); // copy the name of input file to a character variable
  if(myid == 0 )
    {
      input_flag=0; // boolean check
      printf("Opening input file for reading parameters %s \n", param_name);
      params = qcd_getParams( param_name, &params_len); // retrieve params in string and params_len only for root
      if(params == NULL) input_flag=1; // if fail to take params return 1
    }
  MPI_Bcast(&input_flag, 1, MPI_INT, 0, MPI_COMM_WORLD); // broadcast success or fail to all
  if(input_flag == 1) exit(EXIT_FAILURE); 
  MPI_Bcast(&params_len, 1, MPI_INT, 0, MPI_COMM_WORLD); // broadcast params_len to all
  if(myid != 0) params = (char*)malloc(params_len*sizeof(char)); // allocate memory for params for all exect root beacause a previous routine do that for root
  MPI_Bcast(params, params_len, MPI_CHAR, 0, MPI_COMM_WORLD); // Broadcast the params to all processors
   //////////////////////////////////////////////////////////////////////////////
  
  sscanf(qcd_getParam("<processors_txyz>",params,params_len),"%hd %hd %hd %hd",&P[0], &P[1], &P[2], &P[3]); // read number of processors in each direction
  sscanf(qcd_getParam("<lattice_txyz>",params,params_len),"%hd %hd %hd %hd",&L[0], &L[1], &L[2], &L[3]);    // read the dimensions of global lattice
  
    if(P[0] != 1)
    {
      if( myid == 0) fprintf(stderr, " Error! , Number of processors in t-direction must be always 1"); // only use one processor in time direction
      exit(EXIT_FAILURE);
    }
    
  if(qcd_initGeometry(&geo, L, P, theta, myid, numprocs)) exit(EXIT_FAILURE);     // start the geometry each processor takes a segment of lattice
  if(myid == 0) printf( " Local Lattice : %d x %d x %d x %d \n", geo.lL[0], geo.lL[1], geo.lL[2], geo.lL[3]); // print the local lattice

  strcpy(gauge_name,qcd_getParam("<cfg_name>",params,params_len));                // read gauge filename
  if(myid==0) printf("Got conf name: %s\n",gauge_name);
  
   strcpy(momlist_name,qcd_getParam("<momenta_list>",params,params_len));
   if(myid==0) printf("Got momenta-list file name: %s\n",momlist_name);
  
      sscanf(qcd_getParam("<alpha_gauss>",params,params_len),"%lf",&alphaGauss);           // read alpha parameter for wuppertal smearing
    if(myid==0) printf(" Got alpha_gauss: %lf\n",alphaGauss);
    
    sscanf(qcd_getParam("<nsmear_gauss>",params,params_len),"%d",&nsmear);          // number of iterations for wuppertal smearing
    if(myid==0) printf(" Got nsmear_gauss: %d\n",nsmear);
    
    sscanf(qcd_getParam("<alpha_APE>",params,params_len),"%lf",&alphaAPE);         // read alpha parameter for APE smearing
    if(myid==0) printf(" Got alpha_APE: %lf\n",alphaAPE);
    
    sscanf(qcd_getParam("<nsmear_APE>",params,params_len),"%d",&nsmearAPE);         // number of iterations for APE smearing
    if(myid==0) printf(" Got nsmear_APE: %d\n",nsmearAPE);
  
  strcpy(output_name,qcd_getParam("<output_name>",params,params_len)); //get  output filename                                                                
  if(myid==0) printf("Got output file name: %s\n",output_name);

  strcpy(uprop_name,qcd_getParam("<propagator_u>",params,params_len)); //get  u propagator file list
  if(myid==0) printf("Got propagator file name: %s\n",uprop_name);

  strcpy(dprop_name,qcd_getParam("<propagator_d>",params,params_len)); //get d propagator file list
  if(myid==0) printf("Got propagator file name: %s\n",dprop_name);

  strcpy(sprop_name,qcd_getParam("<propagator_s>",params,params_len)); //get s propagator file list
  if(myid==0) printf("Got propagator file name: %s\n",sprop_name);
  
  strcpy(cprop_name,qcd_getParam("<propagator_c>",params,params_len)); //get c propagator file list
  if(myid==0) printf("Got propagator file name: %s\n",cprop_name);
  
  sscanf(qcd_getParam("<source_pos_txyz>",params,params_len),"%d %d %d %d",&x_src[0],&x_src[1],&x_src[2],&x_src[3]);  // read the position where we put the source
  if(myid==0) printf("Got source coords: %d %d %d %d\n",x_src[0],x_src[1],x_src[2],x_src[3]);

  sscanf(qcd_getParam("<current_fixed_t>",params,params_len),"%d",&cur_time);                 // read the time position of insertion
  if(myid==0) printf("Got time to fixed the current time slice: %d \n",cur_time);
  
  sscanf(qcd_getParam("<number_tsinks>",params,params_len),"%d",&number_tsinks);          // read number of sinks 
  if(myid==0) printf("Got the number of tsinks you want: %d \n",number_tsinks);
    
  sscanf(qcd_getParam("<after_tcur>",params,params_len),"%d",&after_tcur);          // read number after insertion current 
  if(myid==0) printf("Got the number of timeslices after current : %d \n", after_tcur);
  
  free(params);
   
  
      if(myid==0)
      {
	for(int iparticle =0 ; iparticle < MAX_PARTICLES_2 ; iparticle++){
	      sprintf(names,"%s_%s.dat",output_name,particles[iparticle]);
	      fp_twop[iparticle] = fopen(names,"w");
	      if(fp_twop[iparticle]==NULL)
		{
		  printf("failed to open %s for writing\n",names);
		  input_flag=1;
		}
	    }

	fp_momlist = fopen(momlist_name,"r");
	if(fp_momlist==NULL)
	  {
	    printf("failed to open %s for reading\n",momlist_name);
	    input_flag=1;
	  }
      }
    MPI_Bcast(&input_flag,1,MPI_INT, 0, MPI_COMM_WORLD);
    if(input_flag==1) exit(EXIT_FAILURE); // check fo success
    //////////////////////////////

    // load the momenta list///////////////////////////////
    if(myid==0) dummy = fscanf(fp_momlist,"%d\n",&num_momenta);
    MPI_Bcast(&num_momenta,1,MPI_INT, 0, MPI_COMM_WORLD);
    if(myid==0) printf("will read %d momenta combinations\n",num_momenta);
    mom = malloc(num_momenta*3*sizeof(qcd_int_4));
    if(myid==0)
      {
	for(int j=0; j<num_momenta; j++)
	  {
	    dummy = fscanf(fp_momlist,"%d %d %d\n",&(mom[j][0]),&(mom[j][1]),&(mom[j][2]));
	    //printf("got combination %i %i %i\n",mom[j][0],mom[j][1],mom[j][2]);                                                                                                                        
	  }
	fclose(fp_momlist);
      }
    MPI_Bcast(&(mom[0][0]),num_momenta*3,MPI_INT,0, MPI_COMM_WORLD);
    if(myid==0) printf("momenta list read and broadcasted\n");
    ////////////////////// succesfuly pass the momenta to all cpu
    
    malloc_flag=0;
    malloc_flag += qcd_initGaugeField(&u, &geo);
    malloc_flag += qcd_initGaugeField(&uAPE,&geo);
    MPI_Allreduce(&malloc_flag, &malloc_flag_check, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if(malloc_flag_check>0)
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
    for(int i=0; i<nsmearAPE; i++)
      {
	qcd_apeSmear3d(uAPE_ptr, u_ptr, alphaAPE);
	utmp_ptr=u_ptr; u_ptr=uAPE_ptr; uAPE_ptr=utmp_ptr;
      }
    utmp_ptr=u_ptr; u_ptr=uAPE_ptr; uAPE_ptr=utmp_ptr; //reverse the last swap. Also needed when nsmearAPE=0                                                                                           
    qcd_destroyGaugeField(u_ptr);
    uAPE = *uAPE_ptr;
    if(myid == 0) printf("printf finish APE smearing\n");
  ///////////////////////////////////////////////////////////////////////////////////////
    
    malloc_flag = 0; // flag for ok allocation
    malloc_flag += qcd_initPropagator(&uprop, &geo); // initialize memory for u propagator
    malloc_flag += qcd_initPropagator(&dprop, &geo); // initialize memory for d propagator
    malloc_flag += qcd_initPropagator(&sprop, &geo); // initialize memory for s propagator
    malloc_flag += qcd_initPropagator(&cprop, &geo); // initialize memory for s propagator
    
    malloc_flag += qcd_initVector(&(vec_up[0]), &geo);
    malloc_flag += qcd_initVector(&(vec_down[0]), &geo);
    malloc_flag += qcd_initVector(&(vec_strange[0]), &geo);
    malloc_flag += qcd_initVector(&(vec_charm[0]), &geo);
    
    malloc_flag += qcd_initVector(&(vec_up[1]), &geo);
    malloc_flag += qcd_initVector(&(vec_down[1]), &geo);
    malloc_flag += qcd_initVector(&(vec_strange[1]), &geo);
    malloc_flag += qcd_initVector(&(vec_charm[1]), &geo);

    MPI_Allreduce(&malloc_flag, &malloc_flag_check, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); // readuce k for all to see if the allocation is ok!

    if(malloc_flag_check>0) // check for memory allocation flag
      {
	if(myid==0) printf("not enough memory\n");
	exit(EXIT_FAILURE);
      }
    if(myid==0) printf("memory for propagators and vectors allocated\n");
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(qcd_getPropagator(uprop_name,qcd_PROP_LIME, &uprop)) exit(EXIT_FAILURE); // read u propagators from binary in ram
    if(myid==0) printf("up propagator loaded\n");
    if(qcd_getPropagator(dprop_name,qcd_PROP_LIME, &dprop)) exit(EXIT_FAILURE); // read d propagators from binary in ram
    if(myid==0) printf("down propagator loaded\n");
    if(qcd_getPropagator(sprop_name,qcd_PROP_LIME, &sprop)) exit(EXIT_FAILURE); // read s propagators from binary in ram
    if(myid==0) printf("strage propagator loaded\n");
    if(qcd_getPropagator(cprop_name,qcd_PROP_LIME, &cprop)) exit(EXIT_FAILURE); // read c propagators from binary in ram
    if(myid==0) printf("charm propagator loaded\n");
    
    ////////////////////////////////////////////////////////////////////////////////////////////////
    for(int lt=0; lt<geo.lL[0]; lt++)
      {
	t = lt + geo.Pos[0] * geo.lL[0];
	phase_factor = (qcd_complex_16) {cos(theta[0]*t/geo.L[0]),sin(theta[0]*t/geo.L[0])};
	qcd_mulPropagatorC3d(&uprop, phase_factor, (t+x_src[0]) % geo.L[0]);                      
	qcd_mulPropagatorC3d(&dprop, phase_factor, (t+x_src[0]) % geo.L[0]);
	qcd_mulPropagatorC3d(&sprop, phase_factor, (t+x_src[0]) % geo.L[0]);
        qcd_mulPropagatorC3d(&cprop, phase_factor, (t+x_src[0]) % geo.L[0]);
      }
    if(myid==0) printf("propagators transformed to basis with theta-periodic boundary conditions\n");
  ///////////////////////////////////////////////////////////////////////////////////
    ///////////////////// transform propagators from twisted to physical ///////////////
    qcd_tranformPropagatorPhysicalPlus(&uprop,&geo,cur_time, after_tcur, number_tsinks );
    qcd_tranformPropagatorPhysicalMinus(&dprop,&geo,cur_time, after_tcur, number_tsinks );
    qcd_tranformPropagatorPhysicalPlus(&sprop,&geo,cur_time, after_tcur, number_tsinks );
    qcd_tranformPropagatorPhysicalPlus(&cprop,&geo,cur_time, after_tcur, number_tsinks );

        ////////////////////////////////////////////////////////////////////////////////////////////////////
                                                                                                                                                                                                  
    // gaussian smearing of propagators (only the time-slices that will be used)                                                                                                                    
     qcd_communicateGaugeP(&uAPE);
     qcd_waitall(&geo);
        
    for(int mu=0;mu<4;mu++)
      for(int c1=0;c1<3;c1++)
	{
	  if(myid == 0)printf("smearing mu = %d , c1 = %d \n",mu,c1);
	  
         qcd_copyVectorPropagator(&(vec_up[0]),&uprop,mu,c1);             // copy each of 12 vectors for propagator to do smearing
          qcd_copyVectorPropagator(&(vec_down[0]),&dprop,mu,c1);             // copy each of 12 vectors for propagator to do smearing
          qcd_copyVectorPropagator(&(vec_strange[0]),&sprop,mu,c1);             // copy each of 12 vectors for propagator to do smearing
          qcd_copyVectorPropagator(&(vec_charm[0]),&cprop,mu,c1);             // copy each of 12 vectors for propagator to do smearing
	  
           qcd_gaussIteration3d_special_2pf(vec_up, vec_down, vec_strange, vec_charm, &geo, &uAPE, alphaGauss, nsmear, cur_time, number_tsinks, after_tcur,0);
	
	
	   qcd_copyPropagatorVector(&uprop,&(vec_up[0]),mu,c1);
          qcd_copyPropagatorVector(&dprop,&(vec_down[0]),mu,c1);
          qcd_copyPropagatorVector(&sprop,&(vec_strange[0]),mu,c1);
          qcd_copyPropagatorVector(&cprop,&(vec_charm[0]),mu,c1);
	} 
    
    for(int i=0;i<2;i++){
    qcd_destroyVector(&(vec_up[i]));                   // free memory 
    qcd_destroyVector(&(vec_down[i]));
    qcd_destroyVector(&(vec_strange[i]));
    qcd_destroyVector(&(vec_charm[i]));
    }
    qcd_destroyGaugeField(&uAPE);

                                               
    qcd_waitall(&geo);
    if(myid==0) printf("propagators smeared\n");
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////// allocate memory for block //////////////////////////////////////////
    malloc_flag = 0;

    for(int mu = 0 ; mu < 4 ; mu++)
      for(int nu = 0 ; nu < 4 ; nu++)
	for(int ij = 0 ; ij < 10 ; ij++)
	  for( int iparticle = 0 ; iparticle < MAX_PARTICLES_2 ; iparticle++){
	    block[iparticle][mu*4+nu][ij] = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
	    if(block[iparticle][mu*4+nu][ij] == NULL ) malloc_flag ++;
	  }
    
    MPI_Allreduce(&malloc_flag, &malloc_flag_check, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); // readuce k for all to see if the allocation is ok!
    if(malloc_flag_check>0) // check for memory allocation flag
      {
	if(myid==0) printf("not enough memory\n");
	exit(EXIT_FAILURE);
      }
    if(myid==0) printf("memory for blocks allocated\n");
  
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //////////////////// find non zero elements for gamma matrices compinations //////////////////
    
    for(int i=0; i < 10 ; i++) nz_counter[i] = 0;
    
    for(int icg1 = 0 ; icg1 < 3 ; icg1++)         // Cg_i * \bar{Cg_j}
      for(int icg2 = 0 ; icg2 < 3 ; icg2++)
	for(int alpha = 0 ; alpha < 4 ; alpha++)
	  for(int beta = 0 ; beta < 4 ; beta++)
	    for(int beta1 = 0 ; beta1 < 4 ; beta1++)
	      for(int alpha1 = 0 ; alpha1 < 4 ; alpha1++){
		C = qcd_CMUL(qcd_CGAMMA[list_cgamma[icg1]][alpha][beta],qcd_BAR_CGAMMA[list_cgamma[icg2]][beta1][alpha1]);
	      
		if(qcd_NORM(C) > 1e-3){
		  cg5cg5b_val[icg1*3+icg2][nz_counter[icg1*3+icg2]].re = C.re;  
		  cg5cg5b_val[icg1*3+icg2][nz_counter[icg1*3+icg2]].im = C.im;
		  cg5cg5b_ind[icg1*3+icg2][nz_counter[icg1*3+icg2]][0] = alpha;
		  cg5cg5b_ind[icg1*3+icg2][nz_counter[icg1*3+icg2]][1] = beta;
		  cg5cg5b_ind[icg1*3+icg2][nz_counter[icg1*3+icg2]][2] = beta1;
		  cg5cg5b_ind[icg1*3+icg2][nz_counter[icg1*3+icg2]][3] = alpha1;
		  nz_counter[icg1*3+icg2]++;  
		}   
	      }

    for(int alpha = 0 ; alpha < 4 ; alpha++)     // Cg_5 * \bar{Cg_5}
      for(int beta = 0 ; beta < 4 ; beta++)
	for(int beta1 = 0 ; beta1 < 4 ; beta1++)
	  for(int alpha1 = 0 ; alpha1 < 4 ; alpha1++){
	    C = qcd_CMUL(qcd_CGAMMA[list_cgamma[3]][alpha][beta],qcd_BAR_CGAMMA[list_cgamma[3]][beta1][alpha1]);

	    if(qcd_NORM(C) > 1e-3){
	      cg5cg5b_val[9][nz_counter[9]].re = C.re;
	      cg5cg5b_val[9][nz_counter[9]].im = C.im;
	      cg5cg5b_ind[9][nz_counter[9]][0] = alpha;
	      cg5cg5b_ind[9][nz_counter[9]][1] = beta;
	      cg5cg5b_ind[9][nz_counter[9]][2] = beta1;
	      cg5cg5b_ind[9][nz_counter[9]][3] = alpha1;
	      nz_counter[9]++;
	    }
	  }

 
        
    
    ///////////////////////////////////////////// start contractions /////////////////////////////////////////////////////////////
        for(int it = 0; it < number_tsinks ; it++){
	  t = ((it + cur_time + after_tcur)%geo.L[0]);
	  if(geo.myid == 0)printf("Time step %d\n",it);
	  
	  for(int mu = 0 ; mu < 4 ; mu++)
	    for(int nu = 0 ; nu < 4 ; nu++)
	      for(int ij = 0 ; ij < 10 ; ij++)
		for(int iparticle = 0 ; iparticle < MAX_PARTICLES_2 ; iparticle++)
		  for(int v3 =0 ; v3 < geo.lV3 ; v3++)
		    block[iparticle][mu*4+nu][ij][v3] = (qcd_complex_16) {0,0};
	   
	               
	  //	  octet_J12(&uprop, &dprop, &sprop, &cprop, &geo, block, t, nz_counter , cg5cg5b_ind , cg5cg5b_val);
	  // if(myid == 0) printf("finish octet\n");
	  
	  decuplet_J12(&uprop, &dprop, &sprop, &cprop, &geo, block, t, nz_counter , cg5cg5b_ind , cg5cg5b_val);	     
	  if(myid == 0)printf("finish decuplet\n");
	 
	  // octet_c_J12(&uprop, &dprop, &sprop, &cprop, &geo, block, t, nz_counter , cg5cg5b_ind , cg5cg5b_val);
	  //if(myid == 0) printf("finish octet_c\n");

	  decuplet_c_J12(&uprop, &dprop, &sprop, &cprop, &geo, block, t, nz_counter , cg5cg5b_ind , cg5cg5b_val);            
	  if(myid == 0)printf("finish decuplet_c\n");

	  // octet_final6_J12(&uprop, &dprop, &sprop, &cprop, &geo, block, t, nz_counter , cg5cg5b_ind , cg5cg5b_val);
          //if(myid == 0) printf("finish octet_final6\n");

	  decuplet_final4_J12(&uprop, &dprop, &sprop, &cprop, &geo, block, t, nz_counter , cg5cg5b_ind , cg5cg5b_val);
          if(myid == 0)printf("finish decuplet_final4\n");

	  
	 //////////////////////////////// start fourier /////////////////////////////////////////////
	 
	  int imom = 0;
	    //	   if(myid == 0)printf("momentum %d\n",imom);

	   for(int mu = 0 ; mu < 4 ; mu++)
	     for(int nu = 0 ; nu < 4 ; nu++)
	       for(int ij = 0 ; ij < 10 ; ij++)
		 for(int iparticle = 0 ; iparticle < MAX_PARTICLES_2 ; iparticle++)
		   corr_temp[iparticle][mu*4+nu][ij] = (qcd_complex_16) {0,0} ;
	    
	    for(int lx=0; lx < geo.lL[1]; lx++)
            for(int ly=0; ly < geo.lL[2]; ly++)
            for(int lz=0; lz < geo.lL[3]; lz++)
            {
               v3 = qcd_LEXIC0(lx,ly,lz,geo.lL);
               x=lx+geo.Pos[1]*geo.lL[1] - x_src[1];
               y=ly+geo.Pos[2]*geo.lL[2] - x_src[2];
               z=lz+geo.Pos[3]*geo.lL[3] - x_src[3];
               tmp = (((double) mom[imom][0]*x)/geo.L[1] + ((double) mom[imom][1]*y)/geo.L[2] + ((double) mom[imom][2]*z)/geo.L[3])*2*M_PI;
               expon=(qcd_complex_16) {cos(tmp), -sin(tmp)}; //TABULATE FOR LARGE SPEEDUP!!!
//               corr=qcd_CADD(corr, qcd_CMUL(block[v3],C2));
	      
	       for(int mu = 0 ; mu < 4 ; mu++)
		 for(int nu = 0 ; nu < 4 ; nu++)
		   for(int ij = 0 ; ij < 10 ; ij++)
		     for(int iparticle =0 ;iparticle < MAX_PARTICLES_2 ; iparticle++){
		       corr_temp[iparticle][mu*4+nu][ij] = qcd_CADD(corr_temp[iparticle][mu*4+nu][ij],qcd_CMUL(block[iparticle][mu*4+nu][ij][v3],expon));
		     } // close spin indices
		 
            } // close space sum
            
	    for(int mu = 0 ; mu < 4 ; mu++)
	      for(int nu =0 ; nu < 4 ; nu++)
		for(int ij = 0 ; ij < 10 ; ij++)
		  for(int iparticle = 0 ; iparticle < MAX_PARTICLES_2 ; iparticle++){
		    if( ((iparticle >= 8) && (iparticle <= 17)) || ((iparticle >= 24) && (iparticle <= 29)) || ((iparticle >= 36) && (iparticle <= 39)) )
		      MPI_Reduce(&(corr_temp[iparticle][mu*4+nu][ij].re),&(corr[iparticle][mu*4+nu][ij].re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		  }

	    if(myid == 0){
	      
	      for(int iparticle = 0 ; iparticle < MAX_PARTICLES_2 ; iparticle++){
		if( ((iparticle >= 8) && (iparticle <= 17)) || ((iparticle >= 24) && (iparticle <= 29)) || ((iparticle >= 36) && (iparticle <= 39)) ){		
		  for(int i = 0 ; i < 3 ; i++)
		    for(int mu = 0 ; mu < 4 ; mu++)
		      for(int nu = 0 ; nu < 4 ; nu++)
			TermJ12[i][mu][nu] = (qcd_complex_16) {0.,0.} ;

		  for(int mu = 0 ; mu < 4 ; mu++)
		    for(int nu = 0 ; nu < 4 ; nu++)
		      for(int lu = 0 ; lu < 4 ; lu++){
			TermJ12[0][mu][nu] = qcd_CADD(TermJ12[0][mu][nu] , qcd_CMUL( qcd_GAMMA_J12[0][mu][lu] , qcd_CSUB(corr[iparticle][lu*4+nu][1*3+0] , corr[iparticle][lu*4+nu][0*3+1] ) ) );
			TermJ12[1][mu][nu] = qcd_CADD(TermJ12[1][mu][nu] , qcd_CMUL( qcd_GAMMA_J12[1][mu][lu] , qcd_CSUB(corr[iparticle][lu*4+nu][2*3+0] , corr[iparticle][lu*4+nu][0*3+2] ) ) );
			TermJ12[2][mu][nu] = qcd_CADD(TermJ12[2][mu][nu] , qcd_CMUL( qcd_GAMMA_J12[2][mu][lu] , qcd_CSUB(corr[iparticle][lu*4+nu][2*3+1] , corr[iparticle][lu*4+nu][1*3+2] ) ) );
		      }

		  for(int mu = 0 ; mu < 4 ; mu++)
		    for(int nu = 0 ; nu < 4 ; nu++){
		      Term3[mu][nu] = qcd_CSCALE( qcd_CADD(qcd_CADD( corr[iparticle][mu*4+nu][0*3+0], corr[iparticle][mu*4+nu][1*3+1] ) , corr[iparticle][mu*4+nu][2*3+2] ) , 1./3. );
		      Term6[mu][nu] = qcd_CSCALE( qcd_CADD(qcd_CADD(TermJ12[0][mu][nu] , TermJ12[1][mu][nu] ) , TermJ12[2][mu][nu] ) , 1./3.) ;
		      corr_J12[mu][nu] = qcd_CADD( qcd_CADD( Term3[mu][nu] , Term6[mu][nu] ) , corr[iparticle][mu*4+nu][9] ) ;
		      // IMPORTANT : In the latter expresion the first term only decuplet contributes and the second term only octet contributes !!
		    }
		
				
		  corr_proj = (qcd_complex_16) {0.,0.} ;

		  for(int mu = 0 ; mu < 4 ; mu++)
		    for(int nu = 0 ; nu < 4 ; nu++)
		      corr_proj = qcd_CADD(corr_proj , qcd_CMUL( PROJECTOR[13][mu][nu] , corr_J12[nu][mu] ) );

		  fprintf(fp_twop[iparticle], "%d %d %+e %+e \n", it+after_tcur, imom, corr_proj.re, corr_proj.im);

	      
		} // if particles
	      } // close all particles
		

	    } // close only root if statement
	      

	 if(myid == 0) printf("finish time %d\n\n",it);
	} // close time
    
  
    
    
     if(myid == 0)printf("finish program \n");
  //finish MPI library
  MPI_Finalize();
  return 0;
}
