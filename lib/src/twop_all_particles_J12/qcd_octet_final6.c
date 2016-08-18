/* KYRIAKOS HADJIYIANNAKOU 
 twop for octet */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>

#ifndef MAX_PARTICLES_2
#define MAX_PARTICLES_2 40
#endif

int octet_final6_J12(qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator *sprop, qcd_propagator *cprop, qcd_geometry *geo,
          qcd_complex_16 *block[MAX_PARTICLES_2][16][10], qcd_uint_4 t, qcd_uint_4 nz_counter[10],qcd_int_2 cg5cg5b_ind[10][16*16][4],qcd_complex_16 cg5cg5b_val[10][16*16]){
 
   qcd_uint_2 gamma1, gamma, alpha1, alpha, beta1, beta;
   qcd_uint_2 a, b, c, a1, b1, c1;
   qcd_uint_4 v3,v;
   
             for(int inz = 0 ; inz < nz_counter[9] ; inz++){
               

             
                  alpha = cg5cg5b_ind[9][inz][0];
                  beta = cg5cg5b_ind[9][inz][1];
                  beta1 = cg5cg5b_ind[9][inz][2];
                  alpha1 = cg5cg5b_ind[9][inz][3];
                  
		  for(int gamma = 0 ; gamma < 4 ; gamma++)
                    for(int gamma1 = 0 ; gamma1 < 4 ; gamma1++)
               for(int cc1=0;cc1<6;cc1++)
               {
                  a = qcd_EPS[cc1][0];
                  b = qcd_EPS[cc1][1];
                  c = qcd_EPS[cc1][2];
                  for(int cc2=0;cc2<6;cc2++)
                  {          
                     a1 = qcd_EPS[cc2][0];
                     b1 = qcd_EPS[cc2][1];
                     c1 = qcd_EPS[cc2][2];
                     
                     for(int lx=0; lx < geo->lL[1]; lx++)
                     for(int ly=0; ly < geo->lL[2]; ly++)
                     for(int lz=0; lz < geo->lL[3]; lz++)
                     {
                        v3 = qcd_LEXIC0(lx,ly,lz,geo->lL);
                        v =  qcd_LEXIC(t,lx,ly,lz,geo->lL);
                        
                        ////////// omega_zero_c ////////////////////
                        block[30][gamma*4+gamma1][9][v3] = qcd_CSUB( block[30][gamma*4+gamma1][9][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[9][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),sprop->D[v][alpha][gamma1][a][c1] )
                                                ,cprop->D[v][beta][beta1][b][b1] ), sprop->D[v][gamma][alpha1][c][a1] ));
                        
                        block[30][gamma*4+gamma1][9][v3] = qcd_CADD( block[30][gamma*4+gamma1][9][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[9][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),sprop->D[v][alpha][alpha1][a][a1] )
                                                ,cprop->D[v][beta][beta1][b][b1] ), sprop->D[v][gamma][gamma1][c][c1] ));
                        /////////////////////////////////////
                        
                        ////////// ksi_prime_zero_c ////////////////////
                           
                        block[31][gamma*4+gamma1][9][v3] = qcd_CSUB( block[31][gamma*4+gamma1][9][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[9][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),dprop->D[v][alpha][gamma1][a][c1] )
                                                ,cprop->D[v][beta][beta1][b][b1] ), sprop->D[v][gamma][alpha1][c][a1] ), 0.5));
                        
                        block[31][gamma*4+gamma1][9][v3] = qcd_CSUB( block[31][gamma*4+gamma1][9][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[9][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),sprop->D[v][alpha][gamma1][a][c1] )
                                                ,cprop->D[v][beta][beta1][b][b1] ), dprop->D[v][gamma][alpha1][c][a1] ), 0.5));
                        
                        block[31][gamma*4+gamma1][9][v3] = qcd_CADD( block[31][gamma*4+gamma1][9][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[9][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),dprop->D[v][alpha][alpha1][a][a1] )
                                                ,cprop->D[v][beta][beta1][b][b1] ), sprop->D[v][gamma][gamma1][c][c1] ), 0.5));
                        
                        block[31][gamma*4+gamma1][9][v3] = qcd_CADD( block[31][gamma*4+gamma1][9][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[9][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),sprop->D[v][alpha][alpha1][a][a1] )
                                                ,cprop->D[v][beta][beta1][b][b1] ), dprop->D[v][gamma][gamma1][c][c1] ), 0.5));
                        
  
                        /////////////////////////////////////
                        
                        ///////// ksi_zero_c ////////////////////

			block[32][gamma*4+gamma1][9][v3] = qcd_CADD( block[32][gamma*4+gamma1][9][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE(
					         cg5cg5b_val[9][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),dprop->D[v][alpha][alpha1][a][a1] )
						 ,sprop->D[v][beta][beta1][b][b1] ), cprop->D[v][gamma][gamma1][c][c1] ));
                        
                        //////////////////////////////////////////
                        
                        
                        /////////////// ksi_prime_plus_c /////////////////
                        
                        block[33][gamma*4+gamma1][9][v3] = qcd_CSUB( block[33][gamma*4+gamma1][9][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[9][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][gamma1][a][c1] )
                                                ,cprop->D[v][beta][beta1][b][b1] ), sprop->D[v][gamma][alpha1][c][a1] ), 0.5));
                        
                        block[33][gamma*4+gamma1][9][v3] = qcd_CSUB( block[33][gamma*4+gamma1][9][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[9][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),sprop->D[v][alpha][gamma1][a][c1] )
                                                ,cprop->D[v][beta][beta1][b][b1] ), uprop->D[v][gamma][alpha1][c][a1] ), 0.5));
                        
                        block[33][gamma*4+gamma1][9][v3] = qcd_CADD( block[33][gamma*4+gamma1][9][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[9][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][alpha1][a][a1] )
                                                ,cprop->D[v][beta][beta1][b][b1] ), sprop->D[v][gamma][gamma1][c][c1] ), 0.5));
                        
                        block[33][gamma*4+gamma1][9][v3] = qcd_CADD( block[33][gamma*4+gamma1][9][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[9][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),sprop->D[v][alpha][alpha1][a][a1] )
                                                ,cprop->D[v][beta][beta1][b][b1] ), uprop->D[v][gamma][gamma1][c][c1] ), 0.5));
                        
                        //////////////////////////////////////////////

			///////////// ksi_plus_c //////////////////////

			block[34][gamma*4+gamma1][9][v3] = qcd_CADD( block[34][gamma*4+gamma1][9][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE(
					         cg5cg5b_val[9][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][alpha1][a][a1] )
						 ,sprop->D[v][beta][beta1][b][b1] ), cprop->D[v][gamma][gamma1][c][c1] ));

			//////////////////////////////////////////////

			//////////// omega_plus_c_c ////////////////

			block[35][gamma*4+gamma1][9][v3] = qcd_CSUB( block[35][gamma*4+gamma1][9][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE(
					        cg5cg5b_val[9][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),cprop->D[v][alpha][gamma1][a][c1] )
					        ,sprop->D[v][beta][beta1][b][b1] ), cprop->D[v][gamma][alpha1][c][a1] ));

                        block[35][gamma*4+gamma1][9][v3] = qcd_CADD( block[35][gamma*4+gamma1][9][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE(
						cg5cg5b_val[9][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),cprop->D[v][alpha][alpha1][a][a1] )
						,sprop->D[v][beta][beta1][b][b1] ), cprop->D[v][gamma][gamma1][c][c1] ));


			/////////////////////////////////////////////
                        
                        
                     } // close spatial 
                  } // close color 1
               } // close color 2
               
             }  // close spin indices
             
 return 0;
}
