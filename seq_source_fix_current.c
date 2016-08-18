/*
This program take a forward propagator for (up,down,strange) quark and multiply it with a local
current at a fix time with fix momentum, and write it in vector form

Editor : Kyriakos Hadjiyiannakou 2012
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>
#include "seq_fix_current_input/qcd_gamma_down.h"
#include "seq_fix_current_input/qcd_gamma_strange.h"
#include "seq_fix_current_input/qcd_gamma_up.h"
#include "seq_fix_current_input/qcd_gamma_charm.h"

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
  char vec_names_up[12][qcd_MAX_STRING_LENGTH], vec_names_down[12][qcd_MAX_STRING_LENGTH],
  vec_names_strange[12][qcd_MAX_STRING_LENGTH],vec_names_charm[12][qcd_MAX_STRING_LENGTH] ; // names with the source vectors
  
  qcd_uint_4 v3,v;                                 // variables for v3 = lx*ly*lz , v = lx*ly*lz*lt
  qcd_uint_2 L[4];                                 // store global lattice dimensions
  qcd_uint_2 P[4];                                 // store the number of proccesors in each direction
  qcd_uint_4 x_src[4],cur_time;                       // x_src the coordinates of source and cur_time time coordinate for insertion

  qcd_uint_4  t,lt;        // i will use only t for global time and lt for local time
  qcd_complex_16 phase_factor, C2;                     // complex number where we store the phase for fourier transformed
  
  qcd_propagator uprop;                            // for up propagator
  qcd_propagator dprop;                           // for down propagator
  qcd_propagator sprop;				  // for strange propagator
  qcd_propagator cprop;                           // for charm propagator
  
  qcd_int_4 x,y,z;                                 // indices for spatial coordinates
  qcd_real_8 tmp;                                  // fourier phase real number
  qcd_uint_4 lx,ly,lz;                             // local lattice variables
  qcd_vector vec_up, vec_down, vec_strange, vec_charm;        // temporary vectors to store insertion * propagator
  qcd_int_4 current_mom[3];                        // only one momentum for current because we use fix current method
  qcd_geometry geo;                                // a variable to start geometry

  qcd_real_8 theta[4] = {M_PI,0.0,0.0,0.0};    // antiperiodic b.c. in time take back normal fermion fields           
  int myid,numprocs, namelen;                  // for mpi use
  char processor_name[MPI_MAX_PROCESSOR_NAME]; //mpi use
  
  
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

  strcpy(gauge_name,qcd_getParam("<cfg_number>",params,params_len));                // read gauge filename
  if(myid==0) printf("Got conf name: %s\n",gauge_name);
  
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
  
  sscanf(qcd_getParam("<gamma_index>",params,params_len),"%d",&gamma_index);              // read the time position of insertion
  if(myid==0) printf("Got the index of insertion operator: %d \n",gamma_index);

  sscanf(qcd_getParam("<current_mom>",params,params_len),"%d %d %d",&current_mom[0], &current_mom[1], &current_mom[2]); // read momentum
  if(myid==0) printf("Got momentum for insertion current: %d %d %d\n",current_mom[0], current_mom[1], current_mom[2]);

  
  free(params);
  
  malloc_flag = 0;
  malloc_flag += qcd_initPropagator(&uprop, &geo); // initialize memory for u propagator
  malloc_flag += qcd_initPropagator(&dprop, &geo); // initialize memory for d propagator
  malloc_flag += qcd_initPropagator(&sprop, &geo); // initialize memory for s propagator
  malloc_flag += qcd_initPropagator(&cprop, &geo); // initialize memory for c propagator

  malloc_flag += qcd_initVector(&vec_up, &geo);    // initialize temporary vectors 
  malloc_flag += qcd_initVector(&vec_down, &geo);
  malloc_flag += qcd_initVector(&vec_strange,&geo);
  malloc_flag += qcd_initVector(&vec_charm,&geo);
  
  MPI_Allreduce(&malloc_flag, &malloc_flag_check, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); // readuce k for all to see if the allocation is ok!
  
    if(malloc_flag_check>0) // check for memory allocation flag
    {
      if(myid==0) printf("not enough memory\n");
      exit(EXIT_FAILURE);
    }
  if(myid==0) printf("memory for propagators allocated\n");
  
    //############################################################///////////////////

  // load propagators

  if(qcd_getPropagator(uprop_name,qcd_PROP_LIME, &uprop)) exit(EXIT_FAILURE); // read u propagators from binary in ram
  if(myid==0) printf("up propagator loaded\n");
  if(qcd_getPropagator(dprop_name,qcd_PROP_LIME, &dprop)) exit(EXIT_FAILURE); // read d propagators from binary in ram
  if(myid==0) printf("down propagator loaded\n");
  if(qcd_getPropagator(sprop_name,qcd_PROP_LIME, &sprop)) exit(EXIT_FAILURE); // read s propagators from binary in ram
  if(myid==0) printf("strange propagator loaded\n");
  if(qcd_getPropagator(cprop_name,qcd_PROP_LIME, &cprop)) exit(EXIT_FAILURE); // read c propagators from binary in ram
  if(myid==0) printf("charm propagator loaded\n");
  

  //################################################################################                                                                
 

  for(int sc =0 ; sc <12 ; sc++){
    sprintf(vec_names_up[sc],"%s_%s.%s.%05d",output_name,"up",gauge_name,sc);
    sprintf(vec_names_down[sc],"%s_%s.%s.%05d",output_name,"down",gauge_name,sc);
    sprintf(vec_names_strange[sc],"%s_%s.%s.%05d",output_name,"strange",gauge_name,sc);
    sprintf(vec_names_charm[sc],"%s_%s.%s.%05d",output_name,"charm",gauge_name,sc);
  }
  if(myid == 0)printf("finish give name to files\n");
 
  t = cur_time;             // fix the time to put insertion current
  
  for(int gamma =0 ; gamma < 4 ; gamma++) // slow indices
    for(int b =0 ; b < 3 ;b++){
      indices = gamma*3 + b;

      qcd_zeroVector(&vec_up);                    // zero vector for each slow spin-color index
      qcd_zeroVector(&vec_down);
      qcd_zeroVector(&vec_strange);
      qcd_zeroVector(&vec_charm);
      
      for(lx = 0 ; lx < geo.lL[1] ; lx++)                  // each processors do his work for its local lattice
	for(ly = 0 ; ly< geo.lL[2] ; ly++)
	  for(lz = 0 ; lz < geo.lL[3] ; lz++){
	    v =  qcd_LEXIC(t,lx,ly,lz,geo.lL);
	    
	    for(int alpha =0 ; alpha < 4 ; alpha++)
	      for(int a = 0 ; a < 3 ; a++){
		for(int beta = 0 ; beta < 4 ;beta++){ // this index will be sumed
		  vec_up.D[v][alpha][a] = qcd_CADD(vec_up.D[v][alpha][a],qcd_CMUL(qcd_GAMMA_OPERATORS_UP[gamma_index][alpha][beta], uprop.D[v][beta][gamma][a][b]));
		  vec_down.D[v][alpha][a] = qcd_CADD(vec_down.D[v][alpha][a],qcd_CMUL(qcd_GAMMA_OPERATORS_DOWN[gamma_index][alpha][beta], dprop.D[v][beta][gamma][a][b]));
		  vec_strange.D[v][alpha][a] = qcd_CADD(vec_strange.D[v][alpha][a],qcd_CMUL(qcd_GAMMA_OPERATORS_STRANGE[gamma_index][alpha][beta], sprop.D[v][beta][gamma][a][b]));
                  vec_charm.D[v][alpha][a] = qcd_CADD(vec_charm.D[v][alpha][a],qcd_CMUL(qcd_GAMMA_OPERATORS_CHARM[gamma_index][alpha][beta], cprop.D[v][beta][gamma][a][b]));
		}
		x=lx+geo.Pos[1]*geo.lL[1] - x_src[1];  // each processor calculate the global x,y,z from local lx,ly,lz
		y=ly+geo.Pos[2]*geo.lL[2] - x_src[2];
		z=lz+geo.Pos[3]*geo.lL[3] - x_src[3];
		tmp = (((double) current_mom[0]*x)/geo.L[1] + ((double) current_mom[1]*y)/geo.L[2] + ((double) current_mom[2]*z)/geo.L[3])*2*M_PI;
		C2=(qcd_complex_16) {cos(tmp), sin(tmp)}; // for current I put PLUS sign , PLEASE CHECK IT                              
		vec_up.D[v][alpha][a] = qcd_CMUL(vec_up.D[v][alpha][a], C2);   // multiply by the exponential
		vec_down.D[v][alpha][a] =qcd_CMUL(vec_down.D[v][alpha][a], C2);
		vec_strange.D[v][alpha][a] =qcd_CMUL(vec_strange.D[v][alpha][a], C2);
                vec_charm.D[v][alpha][a] =qcd_CMUL(vec_charm.D[v][alpha][a], C2);
	      } // close color a
	  } // close local spatial

      if(myid == 0)printf("finish spin=%d,color=%d\n",gamma,b);
      qcd_writeVectorLime(vec_names_up[indices], qcd_SOURCE_LIME, &vec_up);  // write sequential source in vector form 
      qcd_writeVectorLime(vec_names_down[indices], qcd_SOURCE_LIME, &vec_down);
      qcd_writeVectorLime(vec_names_strange[indices], qcd_SOURCE_LIME, &vec_strange);
      qcd_writeVectorLime(vec_names_charm[indices], qcd_SOURCE_LIME, &vec_charm);
      
    }       
  qcd_destroyGeometry(&geo);                       // free memory 
  qcd_destroyPropagator(&uprop);
  qcd_destroyPropagator(&dprop);
  qcd_destroyPropagator(&sprop);
  qcd_destroyPropagator(&cprop);
  qcd_destroyVector(&vec_up);
  qcd_destroyVector(&vec_down);
  qcd_destroyVector(&vec_strange);
  qcd_destroyVector(&vec_charm);
  
  if(myid == 0)printf("finish program \n");
  //finish MPI library
  MPI_Finalize();
  return 0;
}
