/*
This program take as input the sequential propagator for (up,down,strange) quark for fix current method and 
computes three point function for the octet and decuplet

Editor : Kyriakos Hadjiyiannakou 2012
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>
#include "projectors_pb.h"


void* allocate(size_t size){
  void *ptr = malloc(size);
  if(ptr == NULL){
    fprintf(stderr,"not enought memory for malloc \n");
  }
  return ptr;
}

#define MAX_PARTICLES 40
#define Nmu 4
#define NPROJECTORS 5

int main(int argc, char* argv[]){

#include "seq_fix_current_input/particles_part.h"          // this file includes the insertion parts for all particles
  int malloc_flag = 0, malloc_flag_check = 0, input_flag = 0, input_flag_check = 0, function_flag=0; // flags for check
  int  dummy, number_tsinks, particles_counter = 0, after_tcur;     // number_tsinks is for how many time slices move the sink, particle_counter count particles
  int params_len, string_length;                       // the string legth of the parameters
  char *params;                         // pointer to the params
  FILE *ptr_particles,*fp_threep[MAX_PARTICLES][3][NPROJECTORS]; // file pointers for particle list and particles three point function
  qcd_int_4 num_momenta;                // store the number of total momenta
  FILE *fp_momlist;                     // file pointer to the list of momenta

  char gauge_name[qcd_MAX_STRING_LENGTH];       // store the path of gauge configuration
  char param_name[qcd_MAX_STRING_LENGTH];      // name of parameter file                        
  char particle_list_name[qcd_MAX_STRING_LENGTH];              // name of the file where we have all particles 
  char particle_list[MAX_PARTICLES][qcd_MAX_STRING_LENGTH];    // list which has all particles names to use 
  char output_name[qcd_MAX_STRING_LENGTH];                     // output name for three point function
  char names[qcd_MAX_STRING_LENGTH];                           // temporary string for names of particles
  char momlist_name[qcd_MAX_STRING_LENGTH];    // name of momenta-list file                     


  char uprop_name[qcd_MAX_STRING_LENGTH], seq_uprop_name[qcd_MAX_STRING_LENGTH];      // file with the filenames of up propagators   
  char dprop_name[qcd_MAX_STRING_LENGTH], seq_dprop_name[qcd_MAX_STRING_LENGTH];      // file with the filenames of down propagators
  char sprop_name[qcd_MAX_STRING_LENGTH], seq_sprop_name[qcd_MAX_STRING_LENGTH] ;     // file with the filenames of strange propagators
  char cprop_name[qcd_MAX_STRING_LENGTH], seq_cprop_name[qcd_MAX_STRING_LENGTH] ;     // file with the filenames of charm propagators
  
  qcd_uint_2 particles_flag[MAX_PARTICLES] = {0};
  qcd_uint_4 v3,v;                                 // variables for v3 = lx*ly*lz , v = lx*ly*lz*lt
  qcd_uint_2 L[4];                                 // store global lattice dimensions
  qcd_uint_2 P[4];                                 // store the number of proccesors in each direction
  qcd_uint_4 x_src[4],cur_time, itmom;             // x_src the coordinates of source and cur_time time coordinate for insertion
                                                   // itmom is it*Nmom + imom
  qcd_uint_4  t,lt;        // i will use only t for global time and lt for local time
  qcd_complex_16 phase_factor,sum, sum_tpf[3];   // complex number where we store the phase for fourier transformed
  
  qcd_propagator uprop, seq_uprop;                            // for up propagator
  qcd_propagator dprop, seq_dprop;                           // for down propagator
  qcd_propagator sprop, seq_sprop;			     // for strange propagator
  qcd_propagator cprop, seq_cprop;                           // for charm propagator   
  
  qcd_int_4 x,y,z;                                 // indices for spatial coordinates
  qcd_real_8 tmp;                                  // temporary complex number
  qcd_uint_4 lx,ly,lz;                             // local lattice variables
  qcd_int_4 (*mom)[3];                             // a list where we will store the momenta

  qcd_gaugeField u, uAPE;                     // store gauge field before and after smearing
  qcd_gaugeField *u_ptr, *uAPE_ptr, *utmp_ptr; // some pointer to gauge fields
  qcd_uint_4 nsmear, nsmearAPE;               // smearing parameters
  qcd_real_8 alpha, alphaAPE;                  // smearing parameters
  
  qcd_vector vec_up[2],vec_down[2],vec_strange[2], vec_charm[2];                    // temporary vectors for smearing
  qcd_vector seq_vec_up[2],seq_vec_down[2],seq_vec_strange[2], seq_vec_charm[2];              
    
  qcd_geometry geo;                                // a variable to start geometry

  qcd_real_8 theta[4] = {M_PI,0.0,0.0,0.0};    // antiperiodic b.c. in time take back normal fermion fields                                                                                    
  qcd_real_8 plaq;                             // calculate the plaquette                    
  int myid,numprocs, namelen;                  // for mpi use
  char processor_name[MPI_MAX_PROCESSOR_NAME]; //mpi use
  qcd_real_8 start_time,end_time, all_start_time, all_end_time;
  
  qcd_complex_16 *part0[MAX_PARTICLES][Nmu][Nmu], *part1[MAX_PARTICLES][Nmu][Nmu], *part2[MAX_PARTICLES][Nmu][Nmu]; // part for each insertion type
  qcd_int_2 list_projectors[NPROJECTORS];       // list to store projectors that i will use


  char *particle_list_full[MAX_PARTICLES] = {"PROTON", "NEUTRON", "LAMBDA", "SIGMA_PLUS", "SIGMA_ZERO", "SIGMA_MINUS", "KSI_ZERO", "KSI_MINUS", "DELTA_PLUS_PLUS", "DELTA_PLUS", "DELTA_ZERO", "DELTA_MINUS", "SIGMA_STAR_PLUS", "SIGMA_STAR_ZERO", "SIGMA_STAR_MINUS", "KSI_STAR_ZERO", "KSI_STAR_MINUS", "OMEGA", "LAMBDA_C", "SIGMA_PLUS_PLUS_C", "SIGMA_PLUS_C", "SIGMA_ZERO_C", "KSI_PLUS_PLUS_C_C", "KSI_PLUS_C_C", "SIGMA_STAR_PLUS_PLUS_C", "SIGMA_STAR_PLUS_C", "SIGMA_STAR_ZERO_C", "KSI_STAR_PLUS_PLUS_C_C", "KSI_STAR_PLUS_C_C", "OMEGA_PLUS_PLUS_C_C_C", "OMEGA_ZERO_C", "KSI_PRIME_ZERO_C", "KSI_ZERO_C","KSI_PRIME_PLUS_C", "KSI_PLUS_C", "OMEGA_PLUS_C_C", "OMEGA_STAR_ZERO_C", "KSI_STAR_ZERO_C", "KSI_STAR_PLUS_C", "OMEGA_STAR_PLUS_C_C" };

  list_projectors[0] = 16; // -> g5*g1                                                                               
  list_projectors[1] = 15; // -> g5*g2                                                                                            
  list_projectors[2] = 4;  // -> g5*g3                                                                   
  list_projectors[3] = 3;  // -> g5*(g1 + g2 + g3)
  list_projectors[4] = 13; // -> (1+g0)
  
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
    
    if(qcd_initGeometry(&geo, L, P, theta, myid, numprocs)) exit(EXIT_FAILURE);                           // start the geometry each processor takes a segment of lattice
    if(myid == 0) printf( " Local Lattice : %d x %d x %d x %d \n", geo.lL[0], geo.lL[1], geo.lL[2], geo.lL[3]); // print the local lattice
    
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
  
    strcpy(seq_uprop_name,qcd_getParam("<seq_propagator_u>",params,params_len));    //get sequential  u propagator file list
    if(myid==0) printf("Got sequential propagator file name: %s\n",seq_uprop_name);
    
    strcpy(seq_dprop_name,qcd_getParam("<seq_propagator_d>",params,params_len)); //get sequential d propagator file list
    if(myid==0) printf("Got sequential propagator file name: %s\n",seq_dprop_name);
    
    strcpy(seq_sprop_name,qcd_getParam("<seq_propagator_s>",params,params_len)); //get sequential s propagator file list
    if(myid==0) printf("Got sequential propagator file name: %s\n",seq_sprop_name);
    
    strcpy(seq_cprop_name,qcd_getParam("<seq_propagator_c>",params,params_len)); //get sequential s propagator file list
    if(myid==0) printf("Got sequential propagator file name: %s\n",seq_cprop_name);
    
    strcpy(momlist_name,qcd_getParam("<momenta_list>",params,params_len)); //get the filename for the momenta list
    if(myid==0) printf("Got momenta-list file name: %s\n",momlist_name);
    
    
    sscanf(qcd_getParam("<source_pos_txyz>",params,params_len),"%d %d %d %d",&x_src[0],&x_src[1],&x_src[2],&x_src[3]);  // read the position where we put the source
    if(myid==0) printf("Got source coords: %d %d %d %d\n",x_src[0],x_src[1],x_src[2],x_src[3]);
    
    sscanf(qcd_getParam("<current_fixed_t>",params,params_len),"%d",&cur_time);              // read the time position of insertion
    if(myid==0) printf("Got time to fixed the current time slice: %d \n",cur_time);
  
    sscanf(qcd_getParam("<number_tsinks>",params,params_len),"%d",&number_tsinks);          // read number of sinks 
    if(myid==0) printf("Got the number of tsinks you want: %d \n",number_tsinks);
    
    sscanf(qcd_getParam("<after_tcur>",params,params_len),"%d",&after_tcur);          // read number after insertion current 
    if(myid==0) printf("Got the number of timeslices after current : %d \n", after_tcur);
    
    strcpy(particle_list_name,qcd_getParam("<particle_list_name>",params,params_len));     // read particles filename 
    if(myid==0) printf("Got particle list name: %s\n",particle_list_name);
  
    free(params);
    ///////////////////////////////////////////////////////////////////// 
    if (myid == 0){
      ptr_particles = fopen(particle_list_name,"r");
      if ( ptr_particles == NULL){
	fprintf(stderr,"fail to open list with particles for reading\n");
	input_flag = 1;
      }
    }
    MPI_Bcast(&input_flag,1,MPI_INT, 0, MPI_COMM_WORLD);
    if(input_flag==1) exit(EXIT_FAILURE); // check for success

    if(myid == 0){
      while(!feof(ptr_particles)){
	dummy = fscanf(ptr_particles,"%s",&(particle_list[particles_counter][0]));
	particles_counter++;
      }
      for(int i =0 ; i < particles_counter -1; i++)
	printf("Got particle : %s\n",particle_list[i]);
    }
  
    for(int i=0 ; i < MAX_PARTICLES ; i++){
      //    string_length = strlen(particle_list[i]);
      MPI_Bcast(particle_list[i],qcd_MAX_STRING_LENGTH, MPI_CHAR, 0, MPI_COMM_WORLD );  
    }
    
    MPI_Bcast(&particles_counter, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if( numprocs < MAX_PARTICLES ){
      if(myid == 0)fprintf(stderr,"Error number of processors must be at least equal to the number of particles\n");
      exit(EXIT_FAILURE);
    }

    if( (particles_counter-1) > MAX_PARTICLES ){
      if(myid == 0)fprintf(stderr,"Error, maximum number of particles in particles_list is 40\n");
      exit(EXIT_FAILURE);
    } 

  //////////////////////////////////// Start APE smearing ///////////////////////////
  
    //allocate memory for gauge-fields//////////////////////////////////////////
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
    //////////////////////////////////////////////////////////////////////////////////
    malloc_flag = 0; // flag for ok allocation
    malloc_flag += qcd_initPropagator(&uprop, &geo); // initialize memory for u propagator
    malloc_flag += qcd_initPropagator(&dprop, &geo); // initialize memory for d propagator
    malloc_flag += qcd_initPropagator(&sprop, &geo); // initialize memory for s propagator
    malloc_flag += qcd_initPropagator(&cprop, &geo); // initialize memory for c propagator
    

    malloc_flag += qcd_initPropagator(&seq_uprop, &geo); // initialize memory for u propagator
    malloc_flag += qcd_initPropagator(&seq_dprop, &geo); // initialize memory for d propagator
    malloc_flag += qcd_initPropagator(&seq_sprop, &geo); // initialize memory for s propagator
    malloc_flag += qcd_initPropagator(&seq_cprop, &geo); // initialize memory for c propagator
 
    malloc_flag += qcd_initVector(&(vec_up[0]), &geo);
    malloc_flag += qcd_initVector(&(vec_down[0]), &geo);
    malloc_flag += qcd_initVector(&(vec_strange[0]), &geo);
    malloc_flag += qcd_initVector(&(vec_charm[0]), &geo);

    malloc_flag += qcd_initVector(&(seq_vec_up[0]), &geo);
    malloc_flag += qcd_initVector(&(seq_vec_down[0]), &geo);
    malloc_flag += qcd_initVector(&(seq_vec_strange[0]), &geo);
    malloc_flag += qcd_initVector(&(seq_vec_charm[0]), &geo);
    
    malloc_flag += qcd_initVector(&(vec_up[1]), &geo);
    malloc_flag += qcd_initVector(&(vec_down[1]), &geo);
    malloc_flag += qcd_initVector(&(vec_strange[1]), &geo);
    malloc_flag += qcd_initVector(&(vec_charm[1]), &geo);

    malloc_flag += qcd_initVector(&(seq_vec_up[1]), &geo);
    malloc_flag += qcd_initVector(&(seq_vec_down[1]), &geo);
    malloc_flag += qcd_initVector(&(seq_vec_strange[1]), &geo);
    malloc_flag += qcd_initVector(&(seq_vec_charm[1]), &geo);

    
    MPI_Allreduce(&malloc_flag, &malloc_flag_check, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); // readuce k for all to see if the allocation is ok!

    if(malloc_flag_check>0) // check for memory allocation flag
      {
	if(myid==0) printf("not enough memory\n");
	exit(EXIT_FAILURE);
      }
    if(myid==0) printf("memory for propagators and vectors allocated\n");
    /////////////////////////////read propagators //////////////////////////////////////
    // load propagators
    if(myid == 0) start_time = MPI_Wtime();

    if(qcd_getPropagator(uprop_name,qcd_PROP_LIME, &uprop)) exit(EXIT_FAILURE); // read u propagators from binary in ram
    if(myid==0) printf("up propagator loaded\n");
    if(qcd_getPropagator(dprop_name,qcd_PROP_LIME, &dprop)) exit(EXIT_FAILURE); // read d propagators from binary in ram
    if(myid==0) printf("down propagator loaded\n");
    if(qcd_getPropagator(sprop_name,qcd_PROP_LIME, &sprop)) exit(EXIT_FAILURE); // read s propagators from binary in ram
    if(myid==0) printf("strage propagator loaded\n");
    if(qcd_getPropagator(cprop_name,qcd_PROP_LIME, &cprop)) exit(EXIT_FAILURE); // read c propagators from binary in ram
    if(myid==0) printf("charm propagator loaded\n");
    
    
    if(qcd_getPropagator(seq_uprop_name,qcd_PROP_LIME, &seq_uprop)) exit(EXIT_FAILURE); // read u propagators from binary in ram
    if(myid==0) printf("sequential up propagator loaded\n");
    if(qcd_getPropagator(seq_dprop_name,qcd_PROP_LIME, &seq_dprop)) exit(EXIT_FAILURE); // read d propagators from binary in ram
    if(myid==0) printf("sequential down propagator loaded\n");
    if(qcd_getPropagator(seq_sprop_name,qcd_PROP_LIME, &seq_sprop)) exit(EXIT_FAILURE); // read s propagators from binary in ram
    if(myid==0) printf("sequential strange propagator loaded\n");
    if(qcd_getPropagator(seq_cprop_name,qcd_PROP_LIME, &seq_cprop)) exit(EXIT_FAILURE); // read c propagators from binary in ram
    if(myid==0) printf("sequential charm propagator loaded\n");
    
    if(myid == 0) end_time = MPI_Wtime();

    if(myid == 0)printf("Timing for loading props is %e in minutes \n",(end_time-start_time)/60.);

  
    ///////////////////////////////////////////////////////////////////////////////////
    // transform propagators to basis with theta-periodic boundaries in the temporal direction  
    for(int lt=0; lt<geo.lL[0]; lt++)
      {
	t = lt + geo.Pos[0] * geo.lL[0];
	phase_factor = (qcd_complex_16) {cos(theta[0]*t/geo.L[0]),sin(theta[0]*t/geo.L[0])};
	qcd_mulPropagatorC3d(&uprop, phase_factor, (t+x_src[0]) % geo.L[0]);                      
	qcd_mulPropagatorC3d(&dprop, phase_factor, (t+x_src[0]) % geo.L[0]);
	qcd_mulPropagatorC3d(&sprop, phase_factor, (t+x_src[0]) % geo.L[0]);
        qcd_mulPropagatorC3d(&cprop, phase_factor, (t+x_src[0]) % geo.L[0]);
	qcd_mulPropagatorC3d(&seq_uprop, phase_factor, (t+x_src[0]) % geo.L[0]);
	qcd_mulPropagatorC3d(&seq_dprop, phase_factor, (t+x_src[0]) % geo.L[0]);
	qcd_mulPropagatorC3d(&seq_sprop, phase_factor, (t+x_src[0]) % geo.L[0]);
        qcd_mulPropagatorC3d(&seq_cprop, phase_factor, (t+x_src[0]) % geo.L[0]);
        
      }
    if(myid==0) printf("propagators transformed to basis with theta-periodic boundary conditions\n");
  ///////////////////////////////////////////////////////////////////////////////////
    ////////////////////////// transform propagators from twisted to physical base //////////////////
    qcd_tranformPropagatorPhysicalPlus(&uprop,&geo,cur_time, after_tcur, number_tsinks );
    qcd_tranformPropagatorPhysicalMinus(&dprop,&geo,cur_time, after_tcur, number_tsinks );
    qcd_tranformPropagatorPhysicalPlus(&sprop,&geo,cur_time, after_tcur, number_tsinks );
    qcd_tranformPropagatorPhysicalPlus(&cprop,&geo,cur_time, after_tcur, number_tsinks );
   
    qcd_tranformPropagatorPhysicalPlus(&seq_uprop,&geo,cur_time, after_tcur, number_tsinks );
    qcd_tranformPropagatorPhysicalMinus(&seq_dprop,&geo,cur_time, after_tcur, number_tsinks );
    qcd_tranformPropagatorPhysicalPlus(&seq_sprop,&geo,cur_time, after_tcur, number_tsinks );
    qcd_tranformPropagatorPhysicalPlus(&seq_cprop,&geo,cur_time, after_tcur, number_tsinks );
    
    if(myid == 0) printf("propagators transformed to physical base\n");
  //////////////////////////// read momenta ///////////////////////////////////
  
    if(myid == 0){
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
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    // gaussian smearing of propagators (only the time-slices that will be used)                                                                                                                       
        qcd_communicateGaugeP(&uAPE);
        qcd_waitall(&geo);
        if(myid == 0) start_time = MPI_Wtime();
    for(int mu=0;mu<4;mu++)
      for(int c1=0;c1<3;c1++)
	{
	  if(myid == 0)printf("smearing mu = %d , c1 = %d \n",mu,c1);
	  
          qcd_copyVectorPropagator_timerange(&(vec_up[0]),&uprop,mu,c1, cur_time, after_tcur, number_tsinks);             // copy each of 12 vectors for propagator to do smearing
          qcd_copyVectorPropagator_timerange(&(vec_down[0]),&dprop,mu,c1, cur_time, after_tcur, number_tsinks);             // copy each of 12 vectors for propagator to do smearing
          qcd_copyVectorPropagator_timerange(&(vec_strange[0]),&sprop,mu,c1, cur_time, after_tcur, number_tsinks);             // copy each of 12 vectors for propagator to do smearing
          qcd_copyVectorPropagator_timerange(&(vec_charm[0]),&cprop,mu,c1, cur_time, after_tcur, number_tsinks);             // copy each of 12 vectors for propagator to do smearing
          
          qcd_copyVectorPropagator_timerange(&(seq_vec_up[0]),&seq_uprop,mu,c1, cur_time, after_tcur, number_tsinks);             // copy each of 12 vectors for propagator to do smearing
          qcd_copyVectorPropagator_timerange(&(seq_vec_down[0]),&seq_dprop,mu,c1, cur_time, after_tcur, number_tsinks);             // copy each of 12 vectors for propagator to do smearing
          qcd_copyVectorPropagator_timerange(&(seq_vec_strange[0]),&seq_sprop,mu,c1, cur_time, after_tcur, number_tsinks);             // copy each of 12 vectors for propagator to do smearing
          qcd_copyVectorPropagator_timerange(&(seq_vec_charm[0]),&seq_cprop,mu,c1, cur_time, after_tcur, number_tsinks);             // copy each of 12 vectors for propagator to do smearing
	
	  
          qcd_gaussIteration3d_special_3pf(vec_up, vec_down, vec_strange, vec_charm, seq_vec_up
          , seq_vec_down, seq_vec_strange, seq_vec_charm, &geo, &uAPE, alpha, nsmear, cur_time, number_tsinks, after_tcur,0);
	
	
        
          qcd_copyPropagatorVector_timerange(&uprop,&(vec_up[0]),mu,c1, cur_time, after_tcur, number_tsinks);
          qcd_copyPropagatorVector_timerange(&dprop,&(vec_down[0]),mu,c1, cur_time, after_tcur, number_tsinks);
          qcd_copyPropagatorVector_timerange(&sprop,&(vec_strange[0]),mu,c1, cur_time, after_tcur, number_tsinks);
          qcd_copyPropagatorVector_timerange(&cprop,&(vec_charm[0]),mu,c1, cur_time, after_tcur, number_tsinks);
          qcd_copyPropagatorVector_timerange(&seq_uprop,&(seq_vec_up[0]),mu,c1, cur_time, after_tcur, number_tsinks);
          qcd_copyPropagatorVector_timerange(&seq_dprop,&(seq_vec_down[0]),mu,c1, cur_time, after_tcur, number_tsinks);
          qcd_copyPropagatorVector_timerange(&seq_sprop,&(seq_vec_strange[0]),mu,c1, cur_time, after_tcur, number_tsinks);
          qcd_copyPropagatorVector_timerange(&seq_cprop,&(seq_vec_charm[0]),mu,c1, cur_time, after_tcur, number_tsinks);
	  
	} // for \bar{d} \gamma^\mu d we dont need smearing on down propagators
    if(myid == 0) end_time = MPI_Wtime();
    if(myid == 0)printf("Timing for smearing is %e in minutes \n",(end_time-start_time)/60.);
    
    for(int i=0;i<2;i++){
    qcd_destroyVector(&(vec_up[i]));                   // free memory 
    qcd_destroyVector(&(vec_down[i]));
    qcd_destroyVector(&(vec_strange[i]));
    qcd_destroyVector(&(vec_charm[i]));
    qcd_destroyVector(&(seq_vec_up[i]));
    qcd_destroyVector(&(seq_vec_down[i]));
    qcd_destroyVector(&(seq_vec_strange[i]));
    qcd_destroyVector(&(seq_vec_charm[i]));
    }
    qcd_destroyGaugeField(&uAPE);

                                               
    if(myid==0) printf("propagators smeared\n");
    ///////////////////////////////ALLOCATE MEMORY FOR 3PF//////////////////////////////////////////////////////////
    for(int iparticle = 0 ; iparticle < MAX_PARTICLES ; iparticle++)
      for(int mu = 0 ; mu < Nmu ; mu++)
	for(int nu = 0 ; nu < Nmu ; nu++){
	  part0[iparticle][mu][nu] = (qcd_complex_16*)allocate(number_tsinks*num_momenta*sizeof(qcd_complex_16));
	  part1[iparticle][mu][nu] = (qcd_complex_16*)allocate(number_tsinks*num_momenta*sizeof(qcd_complex_16));
	  part2[iparticle][mu][nu] = (qcd_complex_16*)allocate(number_tsinks*num_momenta*sizeof(qcd_complex_16));
	}
    
    for(int iparticle = 0 ; iparticle < MAX_PARTICLES ; iparticle++)
      for(int mu = 0 ; mu < Nmu ; mu++)
	for(int nu = 0 ; nu < Nmu ; nu++)
	  for(int itmom =0 ; itmom < number_tsinks * num_momenta ; itmom++){
	    part0[iparticle][mu][nu][itmom].re = 0; part0[iparticle][mu][nu][itmom].im = 0;
	    part1[iparticle][mu][nu][itmom].re = 0; part1[iparticle][mu][nu][itmom].im = 0;
	    part2[iparticle][mu][nu][itmom].re = 0; part2[iparticle][mu][nu][itmom].im = 0;	  
	  }
    if(myid == 0 )printf("initialize 3pf ok\n");

    if(myid == 0) all_start_time = MPI_Wtime();
    //////////////////////// START PARTICLES /////////////////////////////////////////////
    for(int iparticle = 0 ; iparticle < MAX_PARTICLES ; iparticle++){
  
   ////////// from 0-7 octet particles ///////////////////////////////////////
      if( !(strcmp(particle_list[iparticle],"PROTON")) ){  // check if this is the right particle
    
	if(myid == 0) printf("Start Proton \n");
	if(myid == 0) start_time = MPI_Wtime();
	particles_flag[0] = 1 ;
	function_flag = tpf_proton(&uprop, &dprop, &sprop, &seq_uprop, &seq_dprop,
				   &seq_sprop, &geo, cur_time, number_tsinks, part0[0], 
				   part1[0], part2[0], &(mom[0]), num_momenta, x_src, after_tcur);
	if(myid == 0) end_time = MPI_Wtime();
	if(myid == 0)printf("Timing for PROTON is %e in minutes \n",(end_time-start_time)/60.);
      }
      
      if( !(strcmp(particle_list[iparticle],"NEUTRON")) ){  // check if this is the right particle
    
	if(myid == 0) printf("Start NEUTRON \n");
	if(myid == 0) start_time = MPI_Wtime();
	particles_flag[1] = 1 ;
	function_flag = tpf_neutron(&uprop, &dprop, &sprop, &seq_uprop, &seq_dprop,
				   &seq_sprop, &geo, cur_time, number_tsinks, part0[1], 
				   part1[1], part2[1], &(mom[0]), num_momenta, x_src, after_tcur);
	if(myid == 0) end_time = MPI_Wtime();
	if(myid == 0)printf("Timing for NEUTRON is %e in minutes \n",(end_time-start_time)/60.);
      }
      
      if( !(strcmp(particle_list[iparticle],"LAMBDA")) ){  // check if this is the right particle
    
	if(myid == 0) printf("Start LAMBDA \n");
	if(myid == 0) start_time = MPI_Wtime();
	particles_flag[2] = 1 ;
	function_flag = tpf_lambda(&uprop, &dprop, &sprop, &seq_uprop, &seq_dprop,
				   &seq_sprop, &geo, cur_time, number_tsinks, part0[2], 
				   part1[2], part2[2], &(mom[0]), num_momenta, x_src, after_tcur);
	if(myid == 0) end_time = MPI_Wtime();
	if(myid == 0)printf("Timing for LAMBDA is %e in minutes \n",(end_time-start_time)/60.);
      }
      
      if( !(strcmp(particle_list[iparticle],"SIGMA_PLUS")) ){  // check if this is the right particle
    
	if(myid == 0) printf("Start SIGMA_PLUS \n");
	if(myid == 0) start_time = MPI_Wtime();
	particles_flag[3] = 1 ;
	function_flag = tpf_sigma_plus(&uprop, &dprop, &sprop, &seq_uprop, &seq_dprop,
				   &seq_sprop, &geo, cur_time, number_tsinks, part0[3], 
				   part1[3], part2[3], &(mom[0]), num_momenta, x_src, after_tcur);
	if(myid == 0) end_time = MPI_Wtime();
	if(myid == 0)printf("Timing for SIGMA_PLUS is %e in minutes \n",(end_time-start_time)/60.);
      }
      
      if( !(strcmp(particle_list[iparticle],"SIGMA_ZERO")) ){  // check if this is the right particle
    
	if(myid == 0) printf("Start SIGMA_ZERO \n");
	if(myid == 0) start_time = MPI_Wtime();
	particles_flag[4] = 1 ;
	function_flag = tpf_sigma_zero_test(&uprop, &dprop, &sprop, &seq_uprop, &seq_dprop,
				   &seq_sprop, &geo, cur_time, number_tsinks, part0[4], 
				   part1[4], part2[4], &(mom[0]), num_momenta, x_src, after_tcur);
	if(myid == 0) end_time = MPI_Wtime();
	if(myid == 0)printf("Timing for SIGMA_ZERO is %e in minutes \n",(end_time-start_time)/60.);
      }
      
      if( !(strcmp(particle_list[iparticle],"SIGMA_MINUS")) ){  // check if this is the right particle
    
	if(myid == 0) printf("Start SIGMA_MINUS \n");
	if(myid == 0) start_time = MPI_Wtime();
	particles_flag[5] = 1 ;
	function_flag = tpf_sigma_minus(&uprop, &dprop, &sprop, &seq_uprop, &seq_dprop,
				   &seq_sprop, &geo, cur_time, number_tsinks, part0[5], 
				   part1[5], part2[5], &(mom[0]), num_momenta, x_src, after_tcur);
	if(myid == 0) end_time = MPI_Wtime();
	if(myid == 0)printf("Timing for SIGMA_MINUS is %e in minutes \n",(end_time-start_time)/60.);
      }
   
      if( !(strcmp(particle_list[iparticle],"KSI_ZERO")) ){  // check if this is the right particle
    
	if(myid == 0) printf("Start KSI_ZERO \n");
	if(myid == 0) start_time = MPI_Wtime();
	particles_flag[6] = 1 ;
	function_flag = tpf_ksi_zero(&uprop, &dprop, &sprop, &seq_uprop, &seq_dprop,
				   &seq_sprop, &geo, cur_time, number_tsinks, part0[6], 
				   part1[6], part2[6], &(mom[0]), num_momenta, x_src, after_tcur);
	if(myid == 0) end_time = MPI_Wtime();
	if(myid == 0)printf("Timing for KSI_ZERO is %e in minutes \n",(end_time-start_time)/60.);
      }
      
      if( !(strcmp(particle_list[iparticle],"KSI_MINUS")) ){  // check if this is the right particle
    
	if(myid == 0) printf("Start KSI_MINUS \n");
	if(myid == 0) start_time = MPI_Wtime();
	particles_flag[7] = 1 ;
	function_flag = tpf_ksi_minus(&uprop, &dprop, &sprop, &seq_uprop, &seq_dprop,
				   &seq_sprop, &geo, cur_time, number_tsinks, part0[7], 
				   part1[7], part2[7], &(mom[0]), num_momenta, x_src, after_tcur);
	if(myid == 0) end_time = MPI_Wtime();
	if(myid == 0)printf("Timing for KSI_MINUS is %e in minutes \n",(end_time-start_time)/60.);
      }
      /////////// finish octet particles ////////////////////////////////
      
      //////////// from 8-17 decuplet particles ////////////////////////
      if( !(strcmp(particle_list[iparticle],"DELTA_PLUS_PLUS")) ){  // check if this is the right particle
    
	if(myid == 0) printf("Start DELTA_PLUS_PLUS \n");
	if(myid == 0) start_time = MPI_Wtime();
	particles_flag[8] = 1 ;
	function_flag = tpf_delta_plus_plus(&uprop, &dprop, &sprop, &seq_uprop, &seq_dprop,
				   &seq_sprop, &geo, cur_time, number_tsinks, part0[8], 
				   part1[8], part2[8], &(mom[0]), num_momenta, x_src, after_tcur);
	if(myid == 0) end_time = MPI_Wtime();
	if(myid == 0)printf("Timing for DELTA_PLUS_PLUS is %e in minutes \n",(end_time-start_time)/60.);
      }
      
      if( !(strcmp(particle_list[iparticle],"DELTA_PLUS")) ){  // check if this is the right particle
    
	if(myid == 0) printf("Start DELTA_PLUS \n");
	if(myid == 0) start_time = MPI_Wtime();
	particles_flag[9] = 1 ;
	function_flag = tpf_delta_plus(&uprop, &dprop, &sprop, &seq_uprop, &seq_dprop,
				   &seq_sprop, &geo, cur_time, number_tsinks, part0[9], 
				   part1[9], part2[9], &(mom[0]), num_momenta, x_src, after_tcur);
	if(myid == 0) end_time = MPI_Wtime();
	if(myid == 0)printf("Timing for DELTA_PLUS is %e in minutes \n",(end_time-start_time)/60.);
      }
      
      if( !(strcmp(particle_list[iparticle],"DELTA_ZERO")) ){  // check if this is the right particle
    
	if(myid == 0) printf("Start DELTA_ZERO \n");
	if(myid == 0) start_time = MPI_Wtime();
	particles_flag[10] = 1 ;
	function_flag = tpf_delta_zero(&uprop, &dprop, &sprop, &seq_uprop, &seq_dprop,
				   &seq_sprop, &geo, cur_time, number_tsinks, part0[10], 
				   part1[10], part2[10], &(mom[0]), num_momenta, x_src, after_tcur);
	if(myid == 0) end_time = MPI_Wtime();
	if(myid == 0)printf("Timing for DELTA_ZERO is %e in minutes \n",(end_time-start_time)/60.);
      }
      
      if( !(strcmp(particle_list[iparticle],"DELTA_MINUS")) ){  // check if this is the right particle
    
	if(myid == 0) printf("Start DELTA_MINUS \n");
	if(myid == 0) start_time = MPI_Wtime();
	particles_flag[11] = 1 ;
	function_flag = tpf_delta_minus(&uprop, &dprop, &sprop, &seq_uprop, &seq_dprop,
				   &seq_sprop, &geo, cur_time, number_tsinks, part0[11], 
				   part1[11], part2[11], &(mom[0]), num_momenta, x_src, after_tcur);
	if(myid == 0) end_time = MPI_Wtime();
	if(myid == 0)printf("Timing for DELTA_MINUS is %e in minutes \n",(end_time-start_time)/60.);
      }
      
      if( !(strcmp(particle_list[iparticle],"SIGMA_STAR_PLUS")) ){  // check if this is the right particle
    
	if(myid == 0) printf("Start SIGMA_STAR_PLUS \n");
	if(myid == 0) start_time = MPI_Wtime();
	particles_flag[12] = 1 ;
	function_flag = tpf_sigma_star_plus(&uprop, &dprop, &sprop, &seq_uprop, &seq_dprop,
				   &seq_sprop, &geo, cur_time, number_tsinks, part0[12], 
				   part1[12], part2[12], &(mom[0]), num_momenta, x_src, after_tcur);
	if(myid == 0) end_time = MPI_Wtime();
	if(myid == 0)printf("Timing for SIGMA_STAR_PLUS is %e in minutes \n",(end_time-start_time)/60.);
      }
      
      if( !(strcmp(particle_list[iparticle],"SIGMA_STAR_ZERO")) ){  // check if this is the right particle
    
	if(myid == 0) printf("Start SIGMA_STAR_ZERO \n");
	if(myid == 0) start_time = MPI_Wtime();
	particles_flag[13] = 1 ;
	function_flag = tpf_sigma_star_zero(&uprop, &dprop, &sprop, &seq_uprop, &seq_dprop,
				   &seq_sprop, &geo, cur_time, number_tsinks, part0[13], 
				   part1[13], part2[13], &(mom[0]), num_momenta, x_src, after_tcur);
	if(myid == 0) end_time = MPI_Wtime();
	if(myid == 0)printf("Timing for SIGMA_STAR_ZERO is %e in minutes \n",(end_time-start_time)/60.);
      }   
      
      if( !(strcmp(particle_list[iparticle],"SIGMA_STAR_MINUS")) ){  // check if this is the right particle
    
	if(myid == 0) printf("Start SIGMA_STAR_MINUS \n");
	if(myid == 0) start_time = MPI_Wtime();
	particles_flag[14] = 1 ;
	function_flag = tpf_sigma_star_minus(&uprop, &dprop, &sprop, &seq_uprop, &seq_dprop,
				   &seq_sprop, &geo, cur_time, number_tsinks, part0[14], 
				   part1[14], part2[14], &(mom[0]), num_momenta, x_src, after_tcur);
	if(myid == 0) end_time = MPI_Wtime();
	if(myid == 0)printf("Timing for SIGMA_STAR_MINUS is %e in minutes \n",(end_time-start_time)/60.);
      } 
      
      if( !(strcmp(particle_list[iparticle],"KSI_STAR_ZERO")) ){  // check if this is the right particle
    
	if(myid == 0) printf("Start KSI_STAR_ZERO \n");
	if(myid == 0) start_time = MPI_Wtime();
	particles_flag[15] = 1 ;
	function_flag = tpf_ksi_star_zero(&uprop, &dprop, &sprop, &seq_uprop, &seq_dprop,
				   &seq_sprop, &geo, cur_time, number_tsinks, part0[15], 
				   part1[15], part2[15], &(mom[0]), num_momenta, x_src, after_tcur);
	if(myid == 0) end_time = MPI_Wtime();
	if(myid == 0)printf("Timing for KSI_STAR_ZERO is %e in minutes \n",(end_time-start_time)/60.);
      } 
      
           
      if( !(strcmp(particle_list[iparticle],"KSI_STAR_MINUS")) ){  // check if this is the right particle
    
	if(myid == 0) printf("Start KSI_STAR_MINUS \n");
	if(myid == 0) start_time = MPI_Wtime();
	particles_flag[16] = 1 ;
	function_flag = tpf_ksi_star_minus(&uprop, &dprop, &sprop, &seq_uprop, &seq_dprop,
				   &seq_sprop, &geo, cur_time, number_tsinks, part0[16], 
				   part1[16], part2[16], &(mom[0]), num_momenta, x_src, after_tcur);
	if(myid == 0) end_time = MPI_Wtime();
	if(myid == 0)printf("Timing for KSI_STAR_MINUS is %e in minutes \n",(end_time-start_time)/60.);
      } 
      
      if( !(strcmp(particle_list[iparticle],"OMEGA")) ){  // check if this is the right particle
    
	if(myid == 0) printf("Start OMEGA \n");
	if(myid == 0) start_time = MPI_Wtime();
	particles_flag[17] = 1 ;
	function_flag = tpf_omega(&uprop, &dprop, &sprop, &seq_uprop, &seq_dprop,
				   &seq_sprop, &geo, cur_time, number_tsinks, part0[17], 
				   part1[17], part2[17], &(mom[0]), num_momenta, x_src, after_tcur);
	if(myid == 0) end_time = MPI_Wtime();
	if(myid == 0)printf("Timing for OMEGA is %e in minutes \n",(end_time-start_time)/60.);
      } 
      
      //////////// start particles with charm propagator /////////////////////////////
      
      /* Here i will use same routine instead of strange i will use charm and with spart i will mean cpart */
      
      if( !(strcmp(particle_list[iparticle],"LAMBDA_C")) ){  // check if this is the right particle
    
        if(myid == 0) printf("Start LAMBDA_C \n");
        if(myid == 0) start_time = MPI_Wtime();
	particles_flag[18] = 1 ;
        function_flag = tpf_lambda(&uprop, &dprop, &cprop, &seq_uprop, &seq_dprop,
                                   &seq_cprop, &geo, cur_time, number_tsinks, part0[18], 
                                   part1[18], part2[18], &(mom[0]), num_momenta, x_src, after_tcur);
        if(myid == 0) end_time = MPI_Wtime();
        if(myid == 0)printf("Timing for LAMBDA_C is %e in minutes \n",(end_time-start_time)/60.);
      }
      
      if( !(strcmp(particle_list[iparticle],"SIGMA_PLUS_PLUS_C")) ){  // check if this is the right particle
    
        if(myid == 0) printf("Start SIGMA_PLUS_PLUS_C \n");
        if(myid == 0) start_time = MPI_Wtime();
	particles_flag[19] = 1 ;
        function_flag = tpf_sigma_plus(&uprop, &dprop, &cprop, &seq_uprop, &seq_dprop,
                                   &seq_cprop, &geo, cur_time, number_tsinks, part0[19], 
                                   part1[19], part2[19], &(mom[0]), num_momenta, x_src, after_tcur);
        if(myid == 0) end_time = MPI_Wtime();
        if(myid == 0)printf("Timing for SIGMA_PLUS_PLUS_C is %e in minutes \n",(end_time-start_time)/60.);
      }
      
      if( !(strcmp(particle_list[iparticle],"SIGMA_PLUS_C")) ){  // check if this is the right particle
    
        if(myid == 0) printf("Start SIGMA_PLUS_C \n");
        if(myid == 0) start_time = MPI_Wtime();
	particles_flag[20] = 1 ;
        function_flag = tpf_sigma_zero(&uprop, &dprop, &cprop, &seq_uprop, &seq_dprop,
                                   &seq_cprop, &geo, cur_time, number_tsinks, part0[20], 
                                   part1[20], part2[20], &(mom[0]), num_momenta, x_src, after_tcur);
        if(myid == 0) end_time = MPI_Wtime();
        if(myid == 0)printf("Timing for SIGMA_PLUS_C is %e in minutes \n",(end_time-start_time)/60.);
      }
      
      if( !(strcmp(particle_list[iparticle],"SIGMA_ZERO_C")) ){  // check if this is the right particle
    
        if(myid == 0) printf("Start SIGMA_ZERO_C \n");
        if(myid == 0) start_time = MPI_Wtime();
	particles_flag[21] = 1 ;
        function_flag = tpf_sigma_minus(&uprop, &dprop, &cprop, &seq_uprop, &seq_dprop,
                                   &seq_cprop, &geo, cur_time, number_tsinks, part0[21], 
                                   part1[21], part2[21], &(mom[0]), num_momenta, x_src, after_tcur);
        if(myid == 0) end_time = MPI_Wtime();
        if(myid == 0)printf("Timing for SIGMA_ZERO_C is %e in minutes \n",(end_time-start_time)/60.);
      }
   
      if( !(strcmp(particle_list[iparticle],"KSI_PLUS_PLUS_C_C")) ){  // check if this is the right particle
    
        if(myid == 0) printf("Start KSI_PLUS_PLUS_C_C \n");
        if(myid == 0) start_time = MPI_Wtime();
	particles_flag[22] = 1 ;
        function_flag = tpf_ksi_zero(&uprop, &dprop, &cprop, &seq_uprop, &seq_dprop,
                                   &seq_cprop, &geo, cur_time, number_tsinks, part0[22], 
                                   part1[22], part2[22], &(mom[0]), num_momenta, x_src, after_tcur);
        if(myid == 0) end_time = MPI_Wtime();
        if(myid == 0)printf("Timing for KSI_PLUS_PLUS_C_C is %e in minutes \n",(end_time-start_time)/60.);
      }
      
      if( !(strcmp(particle_list[iparticle],"KSI_PLUS_C_C")) ){  // check if this is the right particle
    
        if(myid == 0) printf("Start KSI_PLUS_C_C \n");
        if(myid == 0) start_time = MPI_Wtime();
	particles_flag[23] = 1 ;
        function_flag = tpf_ksi_minus(&uprop, &dprop, &cprop, &seq_uprop, &seq_dprop,
                                   &seq_cprop, &geo, cur_time, number_tsinks, part0[23], 
                                   part1[23], part2[23], &(mom[0]), num_momenta, x_src, after_tcur);
        if(myid == 0) end_time = MPI_Wtime();
        if(myid == 0)printf("Timing for KSI_PLUS_C_C is %e in minutes \n",(end_time-start_time)/60.);
      }
      
      //////// decaplet for charm ///////////////
      
      if( !(strcmp(particle_list[iparticle],"SIGMA_STAR_PLUS_PLUS_C")) ){  // check if this is the right particle
    
        if(myid == 0) printf("Start SIGMA_STAR_PLUS_PLUS_C \n");
        if(myid == 0) start_time = MPI_Wtime();
	particles_flag[24] = 1 ;
        function_flag = tpf_sigma_star_plus(&uprop, &dprop, &cprop, &seq_uprop, &seq_dprop,
                                   &seq_cprop, &geo, cur_time, number_tsinks, part0[24], 
                                   part1[24], part2[24], &(mom[0]), num_momenta, x_src, after_tcur);
        if(myid == 0) end_time = MPI_Wtime();
        if(myid == 0)printf("Timing for SIGMA_STAR_PLUS_PLUS_C is %e in minutes \n",(end_time-start_time)/60.);
      }
      
      if( !(strcmp(particle_list[iparticle],"SIGMA_STAR_PLUS_C")) ){  // check if this is the right particle
    
        if(myid == 0) printf("Start SIGMA_STAR_PLUS_C \n");
        if(myid == 0) start_time = MPI_Wtime();
	particles_flag[25] = 1 ;
        function_flag = tpf_sigma_star_zero(&uprop, &dprop, &cprop, &seq_uprop, &seq_dprop,
                                   &seq_cprop, &geo, cur_time, number_tsinks, part0[25], 
                                   part1[25], part2[25], &(mom[0]), num_momenta, x_src, after_tcur);
        if(myid == 0) end_time = MPI_Wtime();
        if(myid == 0)printf("Timing for SIGMA_STAR_PLUS_C is %e in minutes \n",(end_time-start_time)/60.);
      }   
      
      if( !(strcmp(particle_list[iparticle],"SIGMA_STAR_ZERO_C")) ){  // check if this is the right particle
    
        if(myid == 0) printf("Start SIGMA_STAR_ZERO_C \n");
        if(myid == 0) start_time = MPI_Wtime();
	particles_flag[26] = 1 ;
        function_flag = tpf_sigma_star_minus(&uprop, &dprop, &cprop, &seq_uprop, &seq_dprop,
                                   &seq_cprop, &geo, cur_time, number_tsinks, part0[26], 
                                   part1[26], part2[26], &(mom[0]), num_momenta, x_src, after_tcur);
        if(myid == 0) end_time = MPI_Wtime();
        if(myid == 0)printf("Timing for SIGMA_STAR_ZERO_C is %e in minutes \n",(end_time-start_time)/60.);
      } 
      
      if( !(strcmp(particle_list[iparticle],"KSI_STAR_PLUS_PLUS_C_C")) ){  // check if this is the right particle
    
        if(myid == 0) printf("Start KSI_STAR_PLUS_PLUS_C_C \n");
        if(myid == 0) start_time = MPI_Wtime();
	particles_flag[27] = 1 ;
        function_flag = tpf_ksi_star_zero(&uprop, &dprop, &cprop, &seq_uprop, &seq_dprop,
                                   &seq_cprop, &geo, cur_time, number_tsinks, part0[27], 
                                   part1[27], part2[27], &(mom[0]), num_momenta, x_src, after_tcur);
        if(myid == 0) end_time = MPI_Wtime();
        if(myid == 0)printf("Timing for KSI_STAR_PLUS_PLUS_C_C is %e in minutes \n",(end_time-start_time)/60.);
      } 
      
           
      if( !(strcmp(particle_list[iparticle],"KSI_STAR_PLUS_C_C")) ){  // check if this is the right particle
    
        if(myid == 0) printf("Start KSI_STAR_PLUS_C_C \n");
        if(myid == 0) start_time = MPI_Wtime();
	particles_flag[28] = 1 ;
        function_flag = tpf_ksi_star_minus(&uprop, &dprop, &cprop, &seq_uprop, &seq_dprop,
                                   &seq_cprop, &geo, cur_time, number_tsinks, part0[28], 
                                   part1[28], part2[28], &(mom[0]), num_momenta, x_src, after_tcur);
        if(myid == 0) end_time = MPI_Wtime();
        if(myid == 0)printf("Timing for KSI_STAR_PLUS_C_C is %e in minutes \n",(end_time-start_time)/60.);
      } 
      
      if( !(strcmp(particle_list[iparticle],"OMEGA_PLUS_PLUS_C_C_C")) ){  // check if this is the right particle
    
        if(myid == 0) printf("Start OMEGA_PLUS_PLUS_C_C_C \n");
        if(myid == 0) start_time = MPI_Wtime();
	particles_flag[29] = 1 ;
        function_flag = tpf_omega(&uprop, &dprop, &cprop, &seq_uprop, &seq_dprop,
                                   &seq_cprop, &geo, cur_time, number_tsinks, part0[29], 
                                   part1[29], part2[29], &(mom[0]), num_momenta, x_src, after_tcur);
        if(myid == 0) end_time = MPI_Wtime();
        if(myid == 0)printf("Timing for OMEGA_PLUS_PLUS_C_C_C is %e in minutes \n",(end_time-start_time)/60.);
      }

      /////////////// new 6 particles for octet style //////////////////////////////////
      if( !(strcmp(particle_list[iparticle],"OMEGA_ZERO_C")) ){  // check if this is the right particle
    
        if(myid == 0) printf("Start OMEGA_ZERO_C \n");
        if(myid == 0) start_time = MPI_Wtime();
	particles_flag[30] = 1 ;
        function_flag = tpf_proton(&sprop, &cprop, &cprop, &seq_sprop, &seq_cprop,              // here the third propagator will not be used
                                   &seq_cprop, &geo, cur_time, number_tsinks, part0[30], 
                                   part1[30], part2[30], &(mom[0]), num_momenta, x_src, after_tcur);
        if(myid == 0) end_time = MPI_Wtime();
        if(myid == 0)printf("Timing for OMEGA_ZERO_C is %e in minutes \n",(end_time-start_time)/60.);
      } 

      if( !(strcmp(particle_list[iparticle],"KSI_PRIME_ZERO_C")) ){  // check if this is the right particle
    
        if(myid == 0) printf("Start KSI_PRIME_ZERO_C \n");
        if(myid == 0) start_time = MPI_Wtime();
	particles_flag[31] = 1 ;
        function_flag = tpf_sigma_zero(&dprop, &sprop, &cprop, &seq_dprop, &seq_sprop,
                                   &seq_cprop, &geo, cur_time, number_tsinks, part0[31], 
                                   part1[31], part2[31], &(mom[0]), num_momenta, x_src, after_tcur);
        if(myid == 0) end_time = MPI_Wtime();
        if(myid == 0)printf("Timing for KSI_PRIME_ZERO_C is %e in minutes \n",(end_time-start_time)/60.);
      } 

      if( !(strcmp(particle_list[iparticle],"KSI_ZERO_C")) ){  // check if this is the right particle
    
        if(myid == 0) printf("Start KSI_ZERO_C \n");
        if(myid == 0) start_time = MPI_Wtime();
	particles_flag[32] = 1 ;
        function_flag = tpf_octet_uds(&dprop, &sprop, &cprop, &seq_dprop, &seq_sprop, 
                                   &seq_cprop, &geo, cur_time, number_tsinks, part0[32], 
                                   part1[32], part2[32], &(mom[0]), num_momenta, x_src, after_tcur);
        if(myid == 0) end_time = MPI_Wtime();
        if(myid == 0)printf("Timing for KSI_ZERO_C is %e in minutes \n",(end_time-start_time)/60.);
      } 

      if( !(strcmp(particle_list[iparticle],"KSI_PRIME_PLUS_C")) ){  // check if this is the right particle
    
        if(myid == 0) printf("Start KSI_PRIME_PLUS_C \n");
        if(myid == 0) start_time = MPI_Wtime();
	particles_flag[33] = 1 ;
        function_flag = tpf_sigma_zero(&uprop, &sprop, &cprop, &seq_uprop, &seq_sprop,
                                   &seq_cprop, &geo, cur_time, number_tsinks, part0[33], 
                                   part1[33], part2[33], &(mom[0]), num_momenta, x_src, after_tcur);
        if(myid == 0) end_time = MPI_Wtime();
        if(myid == 0)printf("Timing for KSI_PRIME_PLUS_C is %e in minutes \n",(end_time-start_time)/60.);
      } 

      if( !(strcmp(particle_list[iparticle],"KSI_PLUS_C")) ){  // check if this is the right particle
    
        if(myid == 0) printf("Start KSI_PLUS_C \n");
        if(myid == 0) start_time = MPI_Wtime();
	particles_flag[34] = 1 ;
        function_flag = tpf_octet_uds(&uprop, &sprop, &cprop, &seq_uprop, &seq_sprop, 
                                   &seq_cprop, &geo, cur_time, number_tsinks, part0[34], 
                                   part1[34], part2[34], &(mom[0]), num_momenta, x_src, after_tcur);
        if(myid == 0) end_time = MPI_Wtime();
        if(myid == 0)printf("Timing for KSI_PLUS_C is %e in minutes \n",(end_time-start_time)/60.);
      } 

      if( !(strcmp(particle_list[iparticle],"OMEGA_PLUS_C_C")) ){  // check if this is the right particle
    
        if(myid == 0) printf("Start OMEGA_PLUS_C_C \n");
        if(myid == 0) start_time = MPI_Wtime();
	particles_flag[35] = 1 ;
        function_flag = tpf_proton(&cprop, &sprop, &sprop, &seq_cprop, &seq_sprop,              // here the third propagator will not be used
                                   &seq_sprop, &geo, cur_time, number_tsinks, part0[35], 
                                   part1[35], part2[35], &(mom[0]), num_momenta, x_src, after_tcur);
        if(myid == 0) end_time = MPI_Wtime();
        if(myid == 0)printf("Timing for OMEGA_PLUS_C_C is %e in minutes \n",(end_time-start_time)/60.);
      } 

      /////////////////////// new 4 particles for decuplet style //////////////

      if( !(strcmp(particle_list[iparticle],"OMEGA_STAR_ZERO_C")) ){  // check if this is the right particle
    
	if(myid == 0) printf("Start OMEGA_STAR_ZERO_C \n");
	if(myid == 0) start_time = MPI_Wtime();
	particles_flag[36] = 1 ;
	function_flag = tpf_ksi_star_zero(&cprop, &dprop, &sprop, &seq_cprop, &seq_dprop,     // here the second propagator will not be used
				   &seq_sprop, &geo, cur_time, number_tsinks, part0[36], 
				   part1[36], part2[36], &(mom[0]), num_momenta, x_src, after_tcur);
	if(myid == 0) end_time = MPI_Wtime();
	if(myid == 0)printf("Timing for OMEGA_STAR_ZERO_C is %e in minutes \n",(end_time-start_time)/60.);
      } 
 
      if( !(strcmp(particle_list[iparticle],"KSI_STAR_ZERO_C")) ){  // check if this is the right particle
    
	if(myid == 0) printf("Start KSI_STAR_ZERO_C \n");
	if(myid == 0) start_time = MPI_Wtime();
	particles_flag[37] = 1 ;
	function_flag = tpf_decuplet_uds(&sprop, &dprop, &cprop, &seq_sprop, &seq_dprop,
				   &seq_cprop, &geo, cur_time, number_tsinks, part0[37], 
				   part1[37], part2[37], &(mom[0]), num_momenta, x_src, after_tcur);
	if(myid == 0) end_time = MPI_Wtime();
	if(myid == 0)printf("Timing for KSI_STAR_ZERO_C is %e in minutes \n",(end_time-start_time)/60.);
      } 

      if( !(strcmp(particle_list[iparticle],"KSI_STAR_PLUS_C")) ){  // check if this is the right particle
    
	if(myid == 0) printf("Start KSI_STAR_PLUS_C \n");
	if(myid == 0) start_time = MPI_Wtime();
	particles_flag[38] = 1 ;
	function_flag = tpf_decuplet_uds(&sprop, &uprop, &cprop, &seq_sprop, &seq_uprop,
				   &seq_cprop, &geo, cur_time, number_tsinks, part0[38], 
				   part1[38], part2[38], &(mom[0]), num_momenta, x_src, after_tcur);
	if(myid == 0) end_time = MPI_Wtime();
	if(myid == 0)printf("Timing for KSI_STAR_PLUS_C is %e in minutes \n",(end_time-start_time)/60.);
      } 

     
      if( !(strcmp(particle_list[iparticle],"OMEGA_STAR_PLUS_C_C")) ){  // check if this is the right particle
    
	if(myid == 0) printf("Start OMEGA_STAR_PLUS_C_C \n");
	if(myid == 0) start_time = MPI_Wtime();
	particles_flag[39] = 1 ;
	function_flag = tpf_ksi_star_zero(&sprop, &dprop, &cprop, &seq_sprop, &seq_dprop,     // here the second propagator will not be used
				   &seq_cprop, &geo, cur_time, number_tsinks, part0[39], 
				   part1[39], part2[39], &(mom[0]), num_momenta, x_src, after_tcur);
	if(myid == 0) end_time = MPI_Wtime();
	if(myid == 0)printf("Timing for OMEGA_STAR_PLUS_C_C is %e in minutes \n",(end_time-start_time)/60.);
      } 

      
    }
    if( myid == 0 ) all_end_time = MPI_Wtime();
    if(myid == 0)printf("Timing for all particles is %e in minutes \n",(all_end_time-all_start_time)/60.);
  ////////////////////////////////open files for writting & act with projector and write 3pf //////////////////////
    if(myid < MAX_PARTICLES ){
      if(particles_flag[myid] == 1){
	for(int ipart =0 ; ipart < 3 ; ipart++)
	  for(int iproj=0 ; iproj < NPROJECTORS ; iproj++){
	    if( (particles_parts[myid][ipart] != "z") && (particles_parts[myid][ipart] != "Z") ){
	      sprintf(names,"%s_%s_part_%s_proj_%d.dat",output_name,particle_list_full[myid],particles_parts[myid][ipart],list_projectors[iproj]);
	      fp_threep[myid][ipart][iproj] = fopen(names,"w");
	      if(fp_threep[myid][ipart][iproj]==NULL)
		{
		  printf("failed to open %s for writing\n",names);
		  input_flag=1;
		}
	    }
	  }
      }
    }
    
    MPI_Allreduce(&input_flag, &input_flag_check, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if(input_flag_check > 0){
      if(myid == 0)printf("failed open files for writting \n");
      exit(EXIT_FAILURE);
    }


    if(myid < MAX_PARTICLES ){
      if(particles_flag[myid] == 1){
	for(int iproj =0 ; iproj < NPROJECTORS ; iproj++)
	  for(int imom =0 ; imom < num_momenta ; imom++)
	    for(int it = 0 ; it < number_tsinks ; it++){
	 
	      for(int ipart=0; ipart < 3; ipart++)sum_tpf[ipart] = (qcd_complex_16) {0,0};
	      itmom = it * num_momenta + imom ;
	      for(int alpha=0 ; alpha < 4 ; alpha++)
		for(int beta=0; beta < 4 ; beta++){
		  sum_tpf[0] = qcd_CADD(sum_tpf[0],qcd_CMUL(PROJECTOR[list_projectors[iproj]][alpha][beta],part0[myid][beta][alpha][itmom]));
		  sum_tpf[1] = qcd_CADD(sum_tpf[1],qcd_CMUL(PROJECTOR[list_projectors[iproj]][alpha][beta],part1[myid][beta][alpha][itmom]));
		  sum_tpf[2] = qcd_CADD(sum_tpf[2],qcd_CMUL(PROJECTOR[list_projectors[iproj]][alpha][beta],part2[myid][beta][alpha][itmom]));

		}

	      for(int ipart =0 ; ipart < 3 ; ipart++)
		if( (particles_parts[myid][ipart] != "z") && (particles_parts[myid][ipart] != "Z") ){
		  fprintf(fp_threep[myid][ipart][iproj], "%d %d %+e %+e \n", imom, it + after_tcur, sum_tpf[ipart].re, sum_tpf[ipart].im);
		}

	    }
      }
    }
    ///////////////////////////////////////////////////////////////////////////////////////
    if(myid == 0)printf("finish program \n");
    //finish MPI library
    MPI_Finalize();
    return 0;
}
