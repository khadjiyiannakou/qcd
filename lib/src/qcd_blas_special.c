#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>

void qcd_copyVectorPropagator_timerange(qcd_vector *vec, qcd_propagator *prop, qcd_uint_2 nu, qcd_uint_2 c2, qcd_uint_4 cur_time, qcd_uint_4 after_tcur, qcd_uint_4 number_tsinks)
{
  qcd_uint_8 i;
  qcd_uint_2 mu,c1;
  qcd_uint_4 t;
  qcd_uint_4 z,y,x;

  for(int it = 0; it < number_tsinks ; it++){
    t = ( (it + cur_time + after_tcur) % vec->geo->L[0]);

    for(z=0; z < vec->geo->lL[3] ;z++)
      for(y=0; y < vec->geo->lL[2] ;y++)
	for(x=0; x < vec->geo->lL[1] ;x++){
	  i=qcd_LEXIC(t,x,y,z,vec->geo->lL);
	  for(mu=0; mu<4; mu++)
	    for(c1=0; c1<3; c1++)
	      vec->D[i][mu][c1] = prop->D[i][mu][nu][c1][c2];
	}
  }
}

void qcd_copyPropagatorVector_timerange(qcd_propagator *prop, qcd_vector *vec, qcd_uint_2 nu, qcd_uint_2 c2, qcd_uint_4 cur_time, qcd_uint_4 after_tcur, qcd_uint_4 number_tsinks)
{
  qcd_uint_8 i;
  qcd_uint_2 mu,c1;
  qcd_uint_4 t;
  qcd_uint_4 z,y,x;

  for(int it = 0; it < number_tsinks ; it++){
    t = ( (it + cur_time + after_tcur) % vec->geo->L[0]);

    for(z=0; z < vec->geo->lL[3] ;z++)
      for(y=0; y < vec->geo->lL[2] ;y++)
	for(x=0; x < vec->geo->lL[1] ;x++){
          i=qcd_LEXIC(t,x,y,z,vec->geo->lL);
	  
	  for(mu=0; mu<4; mu++)
	    for(c1=0; c1<3; c1++)
	      prop->D[i][mu][nu][c1][c2] = vec->D[i][mu][c1];
	}
  }
}

void qcd_tranformPropagatorPhysicalPlus(qcd_propagator *prop, qcd_geometry *geo, qcd_uint_4 cur_time, qcd_uint_4 after_tcur, qcd_uint_4 number_tsinks){

 qcd_complex_16 onePg5[4][4];
 qcd_complex_16 one[4][4];
 qcd_complex_16 imag, C;
 qcd_int_2 onePg5Sq_ind[16*16][4];
 qcd_complex_16 onePg5Sq_val[16*16];
 qcd_uint_4 counter = 0;
 qcd_uint_2 t;
 qcd_uint_8 v;
 qcd_complex_16 tmp[4][4];
 qcd_uint_2 alpha, beta, gamma, delta;
 
 imag = (qcd_complex_16) {0,1};
 
 for(int i=0; i<4 ; i++)
  for(int j=0 ; j<4 ; j++)
   one[i][j] = (qcd_complex_16) {0,0};
  
  one[0][0] = (qcd_complex_16) {1,0};
  one[1][1] = (qcd_complex_16) {1,0};
  one[2][2] = (qcd_complex_16) {1,0};
  one[3][3] = (qcd_complex_16) {1,0};
 
 for(int i= 0 ; i < 4 ; i++)
  for(int j=0; j < 4 ; j++)
   onePg5[i][j] = qcd_CADD(one[i][j],qcd_CMUL(imag,qcd_GAMMA[5][i][j]));
  
  for(int alpha = 0 ; alpha < 4 ; alpha++)
   for(int beta = 0 ; beta < 4 ; beta++)
    for(int gamma = 0 ; gamma < 4 ; gamma++)
     for(int delta = 0 ; delta < 4 ; delta++)
     {
      C = qcd_CMUL(onePg5[alpha][beta], onePg5[gamma][delta]);
          if(qcd_NORM(C)>1e-3){
            onePg5Sq_val[counter].re = 0.5 * C.re;
            onePg5Sq_val[counter].im = 0.5 * C.im;
            onePg5Sq_ind[counter][0] = alpha ;
            onePg5Sq_ind[counter][1] = beta ;
            onePg5Sq_ind[counter][2] = gamma ;
            onePg5Sq_ind[counter][3] = delta ;
            counter++;
          }
     }
     
   
   for(int it = 0; it < number_tsinks ; it++){
    t = ( (it + cur_time + after_tcur) % geo->L[0]);

    for(int z=0; z < geo->lL[3] ;z++)
      for(int y=0; y < geo->lL[2] ;y++)
        for(int x=0; x < geo->lL[1] ;x++){
          v=qcd_LEXIC(t,x,y,z,geo->lL);
          
          for(int c1 = 0 ; c1 < 3 ; c1++)
           for(int c2 = 0 ; c2 < 3 ; c2++){
            
            for(int i = 0 ; i < 4 ; i++)
             for(int j = 0 ; j< 4 ; j++)
              tmp[i][j] = (qcd_complex_16) {0,0};
            
            for(int spins =0 ; spins < counter ; spins++){
              alpha = onePg5Sq_ind[spins][0];
              beta = onePg5Sq_ind[spins][1];
              gamma = onePg5Sq_ind[spins][2];
              delta = onePg5Sq_ind[spins][3];
              
              tmp[alpha][delta] = qcd_CADD(tmp[alpha][delta],qcd_CMUL(onePg5Sq_val[spins],prop->D[v][beta][gamma][c1][c2]));
            } // close spins
            
            for(int alpha = 0 ; alpha < 4 ; alpha++)
             for(int delta = 0 ; delta < 4 ; delta++)
              prop->D[v][alpha][delta][c1][c2] = tmp[alpha][delta];
             
           } // close colors
          
        } // close spatial local volume
        
   } // close time
  
 
}

void qcd_tranformPropagatorPhysicalMinus(qcd_propagator *prop, qcd_geometry *geo, qcd_uint_4 cur_time, qcd_uint_4 after_tcur, qcd_uint_4 number_tsinks){

 qcd_complex_16 oneMg5[4][4];
 qcd_complex_16 one[4][4];
 qcd_complex_16 imag, C;
 qcd_int_2 oneMg5Sq_ind[16*16][4];
 qcd_complex_16 oneMg5Sq_val[16*16];
 qcd_uint_4 counter = 0;
 qcd_uint_2 t;
 qcd_uint_8 v;
 qcd_complex_16 tmp[4][4];
 qcd_uint_2 alpha, beta, gamma, delta;
 
 imag = (qcd_complex_16) {0,1};
 
 for(int i=0; i<4 ; i++)
  for(int j=0 ; j<4 ; j++)
   one[i][j] = (qcd_complex_16) {0,0};
  
  one[0][0] = (qcd_complex_16) {1,0};
  one[1][1] = (qcd_complex_16) {1,0};
  one[2][2] = (qcd_complex_16) {1,0};
  one[3][3] = (qcd_complex_16) {1,0};
 
 for(int i= 0 ; i < 4 ; i++)
  for(int j=0; j < 4 ; j++)
   oneMg5[i][j] = qcd_CSUB(one[i][j],qcd_CMUL(imag,qcd_GAMMA[5][i][j]));
  
  for(int alpha = 0 ; alpha < 4 ; alpha++)
   for(int beta = 0 ; beta < 4 ; beta++)
    for(int gamma = 0 ; gamma < 4 ; gamma++)
     for(int delta = 0 ; delta < 4 ; delta++)
     {
      C = qcd_CMUL(oneMg5[alpha][beta], oneMg5[gamma][delta]);
          if(qcd_NORM(C)>1e-3){
            oneMg5Sq_val[counter].re = 0.5 * C.re;
            oneMg5Sq_val[counter].im = 0.5 * C.im;
            oneMg5Sq_ind[counter][0] = alpha ;
            oneMg5Sq_ind[counter][1] = beta ;
            oneMg5Sq_ind[counter][2] = gamma ;
            oneMg5Sq_ind[counter][3] = delta ;
            counter++;
          }
     }
     
   
   for(int it = 0; it < number_tsinks ; it++){
    t = ( (it + cur_time + after_tcur) % geo->L[0]);

    for(int z=0; z < geo->lL[3] ;z++)
      for(int y=0; y < geo->lL[2] ;y++)
        for(int x=0; x < geo->lL[1] ;x++){
          v=qcd_LEXIC(t,x,y,z,geo->lL);
          
          for(int c1 = 0 ; c1 < 3 ; c1++)
           for(int c2 = 0 ; c2 < 3 ; c2++){
            
            for(int i = 0 ; i < 4 ; i++)
             for(int j = 0 ; j< 4 ; j++)
              tmp[i][j] = (qcd_complex_16) {0,0};
            
            for(int spins =0 ; spins < counter ; spins++){
              alpha = oneMg5Sq_ind[spins][0];
              beta = oneMg5Sq_ind[spins][1];
              gamma = oneMg5Sq_ind[spins][2];
              delta = oneMg5Sq_ind[spins][3];
              
              tmp[alpha][delta] = qcd_CADD(tmp[alpha][delta],qcd_CMUL(oneMg5Sq_val[spins],prop->D[v][beta][gamma][c1][c2]));
            } // close spins
            
            for(int alpha = 0 ; alpha < 4 ; alpha++)
             for(int delta = 0 ; delta < 4 ; delta++)
              prop->D[v][alpha][delta][c1][c2] = tmp[alpha][delta];
             
           } // close colors
          
        } // close spatial local volume
        
   } // close time
  
 
}
