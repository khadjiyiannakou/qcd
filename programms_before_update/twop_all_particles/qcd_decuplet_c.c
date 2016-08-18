/* KYRIAKOS HADJIYIANNAKOU 
 twop for decuplet */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>

#ifndef MAX_PARTICLES_2
#define MAX_PARTICLES_2 30
#endif

int decuplet_c(qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator *sprop, qcd_propagator *cprop, qcd_geometry *geo,
          qcd_complex_16 *block[MAX_PARTICLES_2], qcd_uint_4 t, qcd_uint_4 nz_counter[4],qcd_int_2 Pcg5cg5b_ind[4][16*16*16][6],qcd_complex_16 Pcg5cg5b_val[4][16*16*16]){
 
   qcd_uint_2 gamma1, gamma, alpha1, alpha, beta1, beta;
   qcd_uint_2 a, b, c, a1, b1, c1;
   qcd_uint_4 v3,v;
   
   
        
        for( int icg = 0 ; icg < 3 ; icg ++)
        for(int inz = 0 ; inz < nz_counter[icg] ; inz++){
               
                  gamma1 = Pcg5cg5b_ind[icg][inz][0];
                  gamma = Pcg5cg5b_ind[icg][inz][1];
                  alpha = Pcg5cg5b_ind[icg][inz][2];
                  beta = Pcg5cg5b_ind[icg][inz][3];
                  beta1 = Pcg5cg5b_ind[icg][inz][4];
                  alpha1 = Pcg5cg5b_ind[icg][inz][5];
                  
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

        
 
                        
                ////// sigma_star_plus ////////////////////////////
                        block[24][v3] = qcd_CSUB( block[24][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                Pcg5cg5b_val[icg][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][beta1][a][b1] )
                                                ,uprop->D[v][beta][alpha1][b][a1] ), cprop->D[v][gamma][gamma1][c][c1] ), 1./3.));
                        
                        block[24][v3] = qcd_CADD( block[24][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                Pcg5cg5b_val[icg][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][alpha1][a][a1] )
                                                ,uprop->D[v][beta][beta1][b][b1] ), cprop->D[v][gamma][gamma1][c][c1] ), 1./3.));
                        // *2
                        block[24][v3] = qcd_CSUB( block[24][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                Pcg5cg5b_val[icg][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][gamma1][a][c1] )
                                                ,uprop->D[v][beta][beta1][b][b1] ), cprop->D[v][gamma][alpha1][c][a1] ), 2./3.));
                        
                        block[24][v3] = qcd_CADD( block[24][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                Pcg5cg5b_val[icg][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][beta1][a][b1] )
                                                ,uprop->D[v][beta][gamma1][b][c1] ), cprop->D[v][gamma][alpha1][c][a1] ), 2./3.));
                        
                        block[24][v3] = qcd_CSUB( block[24][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                Pcg5cg5b_val[icg][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),cprop->D[v][alpha][gamma1][a][c1] )
                                                ,uprop->D[v][beta][beta1][b][b1] ), uprop->D[v][gamma][alpha1][c][a1] ), 2./3.));
                        
                        block[24][v3] = qcd_CADD( block[24][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                Pcg5cg5b_val[icg][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),cprop->D[v][alpha][gamma1][a][c1] )
                                                ,uprop->D[v][beta][alpha1][b][a1] ), uprop->D[v][gamma][beta1][c][b1] ), 2./3.));
                        // * 4 
                        
                        block[24][v3] = qcd_CSUB( block[24][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                Pcg5cg5b_val[icg][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),cprop->D[v][alpha][alpha1][a][a1] )
                                                ,uprop->D[v][beta][gamma1][b][c1] ), uprop->D[v][gamma][beta1][c][b1] ), 4./3.));
                        
                        block[24][v3] = qcd_CADD( block[24][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                Pcg5cg5b_val[icg][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),cprop->D[v][alpha][alpha1][a][a1] )
                                                ,uprop->D[v][beta][beta1][b][b1] ), uprop->D[v][gamma][gamma1][c][c1] ), 4./3.));
                //////////////////////////////////////////////////////
                        
                ///////////////////// sigma_star_zero ///////////////
                        block[25][v3] = qcd_CADD( block[25][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                Pcg5cg5b_val[icg][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),cprop->D[v][alpha][beta1][a][b1] )
                                                ,uprop->D[v][beta][gamma1][b][c1] ), dprop->D[v][gamma][alpha1][c][a1] ), 2./3.));
                        
                        block[25][v3] = qcd_CADD( block[25][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                Pcg5cg5b_val[icg][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][beta1][a][b1] )
                                                ,dprop->D[v][beta][gamma1][b][c1] ), cprop->D[v][gamma][alpha1][c][a1] ), 2./3.));
                        
                        block[25][v3] = qcd_CADD( block[25][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                Pcg5cg5b_val[icg][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),dprop->D[v][alpha][beta1][a][b1] )
                                                ,cprop->D[v][beta][gamma1][b][c1] ), uprop->D[v][gamma][alpha1][c][a1] ), 2./3.));
                        
                        block[25][v3] = qcd_CADD( block[25][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                Pcg5cg5b_val[icg][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),cprop->D[v][alpha][gamma1][a][c1] )
                                                ,uprop->D[v][beta][alpha1][b][a1] ), dprop->D[v][gamma][beta1][c][b1] ), 2./3.));
                        
                        block[25][v3] = qcd_CADD( block[25][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                Pcg5cg5b_val[icg][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][gamma1][a][c1] )
                                                ,dprop->D[v][beta][alpha1][b][a1] ), cprop->D[v][gamma][beta1][c][b1] ), 2./3.));
                        
                        block[25][v3] = qcd_CADD( block[25][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                Pcg5cg5b_val[icg][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),dprop->D[v][alpha][gamma1][a][c1] )
                                                ,cprop->D[v][beta][alpha1][b][a1] ), uprop->D[v][gamma][beta1][c][b1] ), 2./3.));
                        
                        block[25][v3] = qcd_CADD( block[25][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                Pcg5cg5b_val[icg][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),cprop->D[v][alpha][alpha1][a][a1] )
                                                ,uprop->D[v][beta][beta1][b][b1] ), dprop->D[v][gamma][gamma1][c][c1] ), 2./3.));
                        
                        block[25][v3] = qcd_CADD( block[25][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                Pcg5cg5b_val[icg][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][alpha1][a][a1] )
                                                ,dprop->D[v][beta][beta1][b][b1] ), cprop->D[v][gamma][gamma1][c][c1] ), 2./3.));
                        
                        block[25][v3] = qcd_CADD( block[25][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                Pcg5cg5b_val[icg][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),dprop->D[v][alpha][alpha1][a][a1] )
                                                ,cprop->D[v][beta][beta1][b][b1] ), uprop->D[v][gamma][gamma1][c][c1] ), 2./3.));
                        
                ////////////////////////////////////////////////////////
                        
                /////////////////////// sigma_star_minus /////////////
                        block[26][v3] = qcd_CSUB( block[26][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                Pcg5cg5b_val[icg][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),dprop->D[v][alpha][beta1][a][b1] )
                                                ,dprop->D[v][beta][alpha1][b][a1] ), cprop->D[v][gamma][gamma1][c][c1] ), 1./3.));
                        
                        block[26][v3] = qcd_CADD( block[26][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                Pcg5cg5b_val[icg][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),dprop->D[v][alpha][alpha1][a][a1] )
                                                ,dprop->D[v][beta][beta1][b][b1] ), cprop->D[v][gamma][gamma1][c][c1] ), 1./3.));
                        // *2
                        block[26][v3] = qcd_CSUB( block[26][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                Pcg5cg5b_val[icg][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),dprop->D[v][alpha][gamma1][a][c1] )
                                                ,dprop->D[v][beta][beta1][b][b1] ), cprop->D[v][gamma][alpha1][c][a1] ), 2./3.));
                        
                        block[26][v3] = qcd_CADD( block[26][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                Pcg5cg5b_val[icg][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),dprop->D[v][alpha][beta1][a][b1] )
                                                ,dprop->D[v][beta][gamma1][b][c1] ), cprop->D[v][gamma][alpha1][c][a1] ), 2./3.));
                        
                        block[26][v3] = qcd_CSUB( block[26][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                Pcg5cg5b_val[icg][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),cprop->D[v][alpha][gamma1][a][c1] )
                                                ,dprop->D[v][beta][beta1][b][b1] ), dprop->D[v][gamma][alpha1][c][a1] ), 2./3.));
                        
                        block[26][v3] = qcd_CADD( block[26][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                Pcg5cg5b_val[icg][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),cprop->D[v][alpha][gamma1][a][c1] )
                                                ,dprop->D[v][beta][alpha1][b][a1] ), dprop->D[v][gamma][beta1][c][b1] ), 2./3.));
                        // * 4 
                        
                        block[26][v3] = qcd_CSUB( block[26][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                Pcg5cg5b_val[icg][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),cprop->D[v][alpha][alpha1][a][a1] )
                                                ,dprop->D[v][beta][gamma1][b][c1] ), dprop->D[v][gamma][beta1][c][b1] ), 4./3.));
                        
                        block[26][v3] = qcd_CADD( block[26][v3] ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                Pcg5cg5b_val[icg][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),cprop->D[v][alpha][alpha1][a][a1] )
                                                ,dprop->D[v][beta][beta1][b][b1] ), dprop->D[v][gamma][gamma1][c][c1] ), 4./3.));
                        
                /////////////////////////////////////////////////////
                        
                //////////////////////// ksi_star_zero ////////////////////////
                        block[27][v3] = qcd_CSUB( block[27][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                Pcg5cg5b_val[icg][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),cprop->D[v][alpha][gamma1][a][c1] )
                                                ,uprop->D[v][beta][beta1][b][b1] ), cprop->D[v][gamma][alpha1][c][a1] ));
                        
                        block[27][v3] = qcd_CADD( block[27][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                Pcg5cg5b_val[icg][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),cprop->D[v][alpha][alpha1][a][a1] )
                                                ,uprop->D[v][beta][beta1][b][b1] ), cprop->D[v][gamma][gamma1][c][c1] ));
                //////////////////////////////////////////////////////////////
                        
                //////////////////////// ksi_star_minus ////////////////////////
                        block[28][v3] = qcd_CSUB( block[28][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                Pcg5cg5b_val[icg][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),cprop->D[v][alpha][gamma1][a][c1] )
                                                ,dprop->D[v][beta][beta1][b][b1] ), cprop->D[v][gamma][alpha1][c][a1] ));
                        
                        block[28][v3] = qcd_CADD( block[28][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                Pcg5cg5b_val[icg][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),cprop->D[v][alpha][alpha1][a][a1] )
                                                ,dprop->D[v][beta][beta1][b][b1] ), cprop->D[v][gamma][gamma1][c][c1] ));
                //////////////////////////////////////////////////////////////
                        
                ////////// omega ///////////////////////////////////////////////////
                        
                        block[29][v3] = qcd_CSUB( block[29][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                Pcg5cg5b_val[icg][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),cprop->D[v][alpha][gamma1][a][c1] )
                                                ,cprop->D[v][beta][beta1][b][b1] ), cprop->D[v][gamma][alpha1][c][a1] ));
                        
                        block[29][v3] = qcd_CADD( block[29][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                Pcg5cg5b_val[icg][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),cprop->D[v][alpha][beta1][a][b1] )
                                                ,cprop->D[v][beta][gamma1][b][c1] ), cprop->D[v][gamma][alpha1][c][a1] ));
                        
                        block[29][v3] = qcd_CADD( block[29][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                Pcg5cg5b_val[icg][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),cprop->D[v][alpha][gamma1][a][c1] )
                                                ,cprop->D[v][beta][alpha1][b][a1] ), cprop->D[v][gamma][beta1][c][b1] ));
                        
                        block[29][v3] = qcd_CSUB( block[29][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                Pcg5cg5b_val[icg][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),cprop->D[v][alpha][alpha1][a][a1] )
                                                ,cprop->D[v][beta][gamma1][b][c1] ), cprop->D[v][gamma][beta1][c][b1] ));
                        
                        block[29][v3] = qcd_CSUB( block[29][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                Pcg5cg5b_val[icg][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),cprop->D[v][alpha][beta1][a][b1] )
                                                ,cprop->D[v][beta][alpha1][b][a1] ), cprop->D[v][gamma][gamma1][c][c1] ));
                        
                        block[29][v3] = qcd_CADD( block[29][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
                                                Pcg5cg5b_val[icg][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),cprop->D[v][alpha][alpha1][a][a1] )
                                                ,cprop->D[v][beta][beta1][b][b1] ), cprop->D[v][gamma][gamma1][c][c1] ));
                  /////////////////////////////////////////////////////////////////////////////
                    } // close spatial 
                  } // close color 1
               } // close color 2
               
             }  // close spin indices
   
    return 0;
}
