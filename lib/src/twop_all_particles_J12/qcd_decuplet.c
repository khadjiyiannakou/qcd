/* KYRIAKOS HADJIYIANNAKOU 
 twop for decuplet */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>

#ifndef MAX_PARTICLES_2
#define MAX_PARTICLES_2 40
#endif

int decuplet_J12(qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator *sprop, qcd_propagator *cprop, qcd_geometry *geo,
          qcd_complex_16 *block[MAX_PARTICLES_2][16][10], qcd_uint_4 t, qcd_uint_4 nz_counter[10],qcd_int_2 cg5cg5b_ind[10][16*16][4],qcd_complex_16 cg5cg5b_val[10][16*16]){
 
   qcd_uint_2 gamma1, gamma, alpha1, alpha, beta1, beta;
   qcd_uint_2 a, b, c, a1, b1, c1;
   qcd_uint_4 v3,v;
   
   
        
        for( int icg1 = 0 ; icg1 < 3 ; icg1++)
	  for(int icg2 = 0 ; icg2 < 3 ; icg2++)
        for(int inz = 0 ; inz < nz_counter[icg1*3+icg2] ; inz++){
               
                  
                  alpha = cg5cg5b_ind[icg1*3+icg2][inz][0];
                  beta = cg5cg5b_ind[icg1*3+icg2][inz][1];
                  beta1 = cg5cg5b_ind[icg1*3+icg2][inz][2];
                  alpha1 = cg5cg5b_ind[icg1*3+icg2][inz][3];

		  for(int gamma = 0 ; gamma < 4 ; gamma++)
		    for(int gamma1 = 0 ;gamma1 < 4 ;gamma1++)
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

        
                  ////////// delta_plus_plus ///////////////////////////////////////////////////
                        
                        block[8][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CSUB( block[8][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][gamma1][a][c1] )
                                                ,uprop->D[v][beta][beta1][b][b1] ), uprop->D[v][gamma][alpha1][c][a1] ));
                        
                        block[8][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CADD( block[8][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][beta1][a][b1] )
                                                ,uprop->D[v][beta][gamma1][b][c1] ), uprop->D[v][gamma][alpha1][c][a1] ));
                        
                        block[8][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CADD( block[8][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][gamma1][a][c1] )
                                                ,uprop->D[v][beta][alpha1][b][a1] ), uprop->D[v][gamma][beta1][c][b1] ));
                        
                        block[8][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CSUB( block[8][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][alpha1][a][a1] )
                                                ,uprop->D[v][beta][gamma1][b][c1] ), uprop->D[v][gamma][beta1][c][b1] ));
                        
                        block[8][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CSUB( block[8][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][beta1][a][b1] )
                                                ,uprop->D[v][beta][alpha1][b][a1] ), uprop->D[v][gamma][gamma1][c][c1] ));
                        
                        block[8][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CADD( block[8][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][alpha1][a][a1] )
                                                ,uprop->D[v][beta][beta1][b][b1] ), uprop->D[v][gamma][gamma1][c][c1] ));
                  /////////////////////////////////////////////////////////////////////////////
                        
                //////////////// delta_plus //////////////////////////////////////////////////////
                        
                        block[9][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CSUB( block[9][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][beta1][a][b1] )
                                                ,uprop->D[v][beta][alpha1][b][a1] ), dprop->D[v][gamma][gamma1][c][c1] ), 1./3.));
                        
                        block[9][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CADD( block[9][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][alpha1][a][a1] )
                                                ,uprop->D[v][beta][beta1][b][b1] ), dprop->D[v][gamma][gamma1][c][c1] ), 1./3.));
                        // *2
                        block[9][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CADD( block[9][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][beta1][a][b1] )
                                                ,dprop->D[v][beta][gamma1][b][c1] ), uprop->D[v][gamma][alpha1][c][a1] ), 2./3.));
                        
                        block[9][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CADD( block[9][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][gamma1][a][c1] )
                                                ,uprop->D[v][beta][alpha1][b][a1] ), dprop->D[v][gamma][beta1][c][b1] ), 2./3.));
                        
                        block[9][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CSUB( block[9][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][alpha1][a][a1] )
                                                ,uprop->D[v][beta][gamma1][b][c1] ), dprop->D[v][gamma][beta1][c][b1] ), 2./3.));
                        
                        block[9][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CSUB( block[9][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][alpha1][a][a1] )
                                                ,dprop->D[v][beta][gamma1][b][c1] ), uprop->D[v][gamma][beta1][c][b1] ), 2./3.));
                        // * 4 
                        
                        block[9][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CSUB( block[9][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][gamma1][a][c1] )
                                                ,dprop->D[v][beta][beta1][b][b1] ), uprop->D[v][gamma][alpha1][c][a1] ), 4./3.));
                        
                        block[9][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CADD( block[9][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][alpha1][a][a1] )
                                                ,dprop->D[v][beta][beta1][b][b1] ), uprop->D[v][gamma][gamma1][c][c1] ), 4./3.));
                        
                //////////////////////////////////////////////
                        
                //////////////// delta_zero ////////////////////////
                        block[10][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CSUB( block[10][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),dprop->D[v][alpha][beta1][a][b1] )
                                                ,dprop->D[v][beta][alpha1][b][a1] ), uprop->D[v][gamma][gamma1][c][c1] ), 1./3.));
                        
                        block[10][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CADD( block[10][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),dprop->D[v][alpha][alpha1][a][a1] )
                                                ,dprop->D[v][beta][beta1][b][b1] ), uprop->D[v][gamma][gamma1][c][c1] ), 1./3.));
                        // *2
                        block[10][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CADD( block[10][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),dprop->D[v][alpha][beta1][a][b1] )
                                                ,uprop->D[v][beta][gamma1][b][c1] ), dprop->D[v][gamma][alpha1][c][a1] ), 2./3.));
                        
                        block[10][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CADD( block[10][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),dprop->D[v][alpha][gamma1][a][c1] )
                                                ,dprop->D[v][beta][alpha1][b][a1] ), uprop->D[v][gamma][beta1][c][b1] ), 2./3.));
                        
                        block[10][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CSUB( block[10][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),dprop->D[v][alpha][alpha1][a][a1] )
                                                ,dprop->D[v][beta][gamma1][b][c1] ), uprop->D[v][gamma][beta1][c][b1] ), 2./3.));
                        
                        block[10][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CSUB( block[10][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),dprop->D[v][alpha][alpha1][a][a1] )
                                                ,uprop->D[v][beta][gamma1][b][c1] ), dprop->D[v][gamma][beta1][c][b1] ), 2./3.));
                        // * 4 
                        
                        block[10][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CSUB( block[10][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),dprop->D[v][alpha][gamma1][a][c1] )
                                                ,uprop->D[v][beta][beta1][b][b1] ), dprop->D[v][gamma][alpha1][c][a1] ), 4./3.));
                        
                        block[10][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CADD( block[10][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),dprop->D[v][alpha][alpha1][a][a1] )
                                                ,uprop->D[v][beta][beta1][b][b1] ), dprop->D[v][gamma][gamma1][c][c1] ), 4./3.));
                /////////////////////////////////////////////////////   
                        
                ////////// delta_minus ///////////////////////////////
                        block[11][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CSUB( block[11][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),dprop->D[v][alpha][gamma1][a][c1] )
                                                ,dprop->D[v][beta][beta1][b][b1] ), dprop->D[v][gamma][alpha1][c][a1] ));
                        
                        block[11][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CADD( block[11][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),dprop->D[v][alpha][beta1][a][b1] )
                                                ,dprop->D[v][beta][gamma1][b][c1] ), dprop->D[v][gamma][alpha1][c][a1] ));
                        
                        block[11][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CADD( block[11][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),dprop->D[v][alpha][gamma1][a][c1] )
                                                ,dprop->D[v][beta][alpha1][b][a1] ), dprop->D[v][gamma][beta1][c][b1] ));
                        
                        block[11][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CSUB( block[11][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),dprop->D[v][alpha][alpha1][a][a1] )
                                                ,dprop->D[v][beta][gamma1][b][c1] ), dprop->D[v][gamma][beta1][c][b1] ));
                        
                        block[11][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CSUB( block[11][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),dprop->D[v][alpha][beta1][a][b1] )
                                                ,dprop->D[v][beta][alpha1][b][a1] ), dprop->D[v][gamma][gamma1][c][c1] ));
                        
                        block[11][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CADD( block[11][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),dprop->D[v][alpha][alpha1][a][a1] )
                                                ,dprop->D[v][beta][beta1][b][b1] ), dprop->D[v][gamma][gamma1][c][c1] ));
                /////////////////////////////////////////////////////
                        
                ////// sigma_star_plus ////////////////////////////
                        block[12][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CSUB( block[12][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][beta1][a][b1] )
                                                ,uprop->D[v][beta][alpha1][b][a1] ), sprop->D[v][gamma][gamma1][c][c1] ), 1./3.));
                        
                        block[12][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CADD( block[12][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][alpha1][a][a1] )
                                                ,uprop->D[v][beta][beta1][b][b1] ), sprop->D[v][gamma][gamma1][c][c1] ), 1./3.));
                        // *2
                        block[12][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CSUB( block[12][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][gamma1][a][c1] )
                                                ,uprop->D[v][beta][beta1][b][b1] ), sprop->D[v][gamma][alpha1][c][a1] ), 2./3.));
                        
                        block[12][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CADD( block[12][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][beta1][a][b1] )
                                                ,uprop->D[v][beta][gamma1][b][c1] ), sprop->D[v][gamma][alpha1][c][a1] ), 2./3.));
                        
                        block[12][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CSUB( block[12][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),sprop->D[v][alpha][gamma1][a][c1] )
                                                ,uprop->D[v][beta][beta1][b][b1] ), uprop->D[v][gamma][alpha1][c][a1] ), 2./3.));
                        
                        block[12][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CADD( block[12][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),sprop->D[v][alpha][gamma1][a][c1] )
                                                ,uprop->D[v][beta][alpha1][b][a1] ), uprop->D[v][gamma][beta1][c][b1] ), 2./3.));
                        // * 4 
                        
                        block[12][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CSUB( block[12][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),sprop->D[v][alpha][alpha1][a][a1] )
                                                ,uprop->D[v][beta][gamma1][b][c1] ), uprop->D[v][gamma][beta1][c][b1] ), 4./3.));
                        
                        block[12][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CADD( block[12][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),sprop->D[v][alpha][alpha1][a][a1] )
                                                ,uprop->D[v][beta][beta1][b][b1] ), uprop->D[v][gamma][gamma1][c][c1] ), 4./3.));
                //////////////////////////////////////////////////////
                        
                ///////////////////// sigma_star_zero ///////////////
                        block[13][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CADD( block[13][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),sprop->D[v][alpha][beta1][a][b1] )
                                                ,uprop->D[v][beta][gamma1][b][c1] ), dprop->D[v][gamma][alpha1][c][a1] ), 2./3.));
                        
                        block[13][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CADD( block[13][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][beta1][a][b1] )
                                                ,dprop->D[v][beta][gamma1][b][c1] ), sprop->D[v][gamma][alpha1][c][a1] ), 2./3.));
                        
                        block[13][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CADD( block[13][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),dprop->D[v][alpha][beta1][a][b1] )
                                                ,sprop->D[v][beta][gamma1][b][c1] ), uprop->D[v][gamma][alpha1][c][a1] ), 2./3.));
                        
                        block[13][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CADD( block[13][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),sprop->D[v][alpha][gamma1][a][c1] )
                                                ,uprop->D[v][beta][alpha1][b][a1] ), dprop->D[v][gamma][beta1][c][b1] ), 2./3.));
                        
                        block[13][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CADD( block[13][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][gamma1][a][c1] )
                                                ,dprop->D[v][beta][alpha1][b][a1] ), sprop->D[v][gamma][beta1][c][b1] ), 2./3.));
                        
                        block[13][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CADD( block[13][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),dprop->D[v][alpha][gamma1][a][c1] )
                                                ,sprop->D[v][beta][alpha1][b][a1] ), uprop->D[v][gamma][beta1][c][b1] ), 2./3.));
                        
                        block[13][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CADD( block[13][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),sprop->D[v][alpha][alpha1][a][a1] )
                                                ,uprop->D[v][beta][beta1][b][b1] ), dprop->D[v][gamma][gamma1][c][c1] ), 2./3.));
                        
                        block[13][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CADD( block[13][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][alpha1][a][a1] )
                                                ,dprop->D[v][beta][beta1][b][b1] ), sprop->D[v][gamma][gamma1][c][c1] ), 2./3.));
                        
                        block[13][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CADD( block[13][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),dprop->D[v][alpha][alpha1][a][a1] )
                                                ,sprop->D[v][beta][beta1][b][b1] ), uprop->D[v][gamma][gamma1][c][c1] ), 2./3.));
                        
                ////////////////////////////////////////////////////////
                        
                /////////////////////// sigma_star_minus /////////////
                        block[14][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CSUB( block[14][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),dprop->D[v][alpha][beta1][a][b1] )
                                                ,dprop->D[v][beta][alpha1][b][a1] ), sprop->D[v][gamma][gamma1][c][c1] ), 1./3.));
                        
                        block[14][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CADD( block[14][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),dprop->D[v][alpha][alpha1][a][a1] )
                                                ,dprop->D[v][beta][beta1][b][b1] ), sprop->D[v][gamma][gamma1][c][c1] ), 1./3.));
                        // *2
                        block[14][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CSUB( block[14][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),dprop->D[v][alpha][gamma1][a][c1] )
                                                ,dprop->D[v][beta][beta1][b][b1] ), sprop->D[v][gamma][alpha1][c][a1] ), 2./3.));
                        
                        block[14][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CADD( block[14][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),dprop->D[v][alpha][beta1][a][b1] )
                                                ,dprop->D[v][beta][gamma1][b][c1] ), sprop->D[v][gamma][alpha1][c][a1] ), 2./3.));
                        
                        block[14][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CSUB( block[14][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),sprop->D[v][alpha][gamma1][a][c1] )
                                                ,dprop->D[v][beta][beta1][b][b1] ), dprop->D[v][gamma][alpha1][c][a1] ), 2./3.));
                        
                        block[14][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CADD( block[14][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),sprop->D[v][alpha][gamma1][a][c1] )
                                                ,dprop->D[v][beta][alpha1][b][a1] ), dprop->D[v][gamma][beta1][c][b1] ), 2./3.));
                        // * 4 
                        
                        block[14][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CSUB( block[14][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),sprop->D[v][alpha][alpha1][a][a1] )
                                                ,dprop->D[v][beta][gamma1][b][c1] ), dprop->D[v][gamma][beta1][c][b1] ), 4./3.));
                        
                        block[14][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CADD( block[14][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),sprop->D[v][alpha][alpha1][a][a1] )
                                                ,dprop->D[v][beta][beta1][b][b1] ), dprop->D[v][gamma][gamma1][c][c1] ), 4./3.));
                        
                /////////////////////////////////////////////////////
                        
                //////////////////////// ksi_star_zero ////////////////////////
                        block[15][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CSUB( block[15][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),sprop->D[v][alpha][gamma1][a][c1] )
                                                ,uprop->D[v][beta][beta1][b][b1] ), sprop->D[v][gamma][alpha1][c][a1] ));
                        
                        block[15][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CADD( block[15][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),sprop->D[v][alpha][alpha1][a][a1] )
                                                ,uprop->D[v][beta][beta1][b][b1] ), sprop->D[v][gamma][gamma1][c][c1] ));
                //////////////////////////////////////////////////////////////
                        
                //////////////////////// ksi_star_minus ////////////////////////
                        block[16][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CSUB( block[16][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),sprop->D[v][alpha][gamma1][a][c1] )
                                                ,dprop->D[v][beta][beta1][b][b1] ), sprop->D[v][gamma][alpha1][c][a1] ));
                        
                        block[16][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CADD( block[16][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),sprop->D[v][alpha][alpha1][a][a1] )
                                                ,dprop->D[v][beta][beta1][b][b1] ), sprop->D[v][gamma][gamma1][c][c1] ));
                //////////////////////////////////////////////////////////////
                        
                ////////// omega ///////////////////////////////////////////////////
                        
                        block[17][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CSUB( block[17][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),sprop->D[v][alpha][gamma1][a][c1] )
                                                ,sprop->D[v][beta][beta1][b][b1] ), sprop->D[v][gamma][alpha1][c][a1] ));
                        
                        block[17][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CADD( block[17][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),sprop->D[v][alpha][beta1][a][b1] )
                                                ,sprop->D[v][beta][gamma1][b][c1] ), sprop->D[v][gamma][alpha1][c][a1] ));
                        
                        block[17][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CADD( block[17][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),sprop->D[v][alpha][gamma1][a][c1] )
                                                ,sprop->D[v][beta][alpha1][b][a1] ), sprop->D[v][gamma][beta1][c][b1] ));
                        
                        block[17][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CSUB( block[17][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),sprop->D[v][alpha][alpha1][a][a1] )
                                                ,sprop->D[v][beta][gamma1][b][c1] ), sprop->D[v][gamma][beta1][c][b1] ));
                        
                        block[17][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CSUB( block[17][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),sprop->D[v][alpha][beta1][a][b1] )
                                                ,sprop->D[v][beta][alpha1][b][a1] ), sprop->D[v][gamma][gamma1][c][c1] ));
                        
                        block[17][gamma*4+gamma1][icg1*3+icg2][v3] = qcd_CADD( block[17][gamma*4+gamma1][icg1*3+icg2][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                cg5cg5b_val[icg1*3+icg2][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),sprop->D[v][alpha][alpha1][a][a1] )
                                                ,sprop->D[v][beta][beta1][b][b1] ), sprop->D[v][gamma][gamma1][c][c1] ));
                  /////////////////////////////////////////////////////////////////////////////
                    } // close spatial 
                  } // close color 1
               } // close color 2
               
             }  // close spin indices
   
    return 0;
}
