#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>
#ifndef MAX_PARTICLES_2
#define MAX_PARTICLES_2 40
#endif

#ifndef qcd_CMUL
#define qcd_CMUL(x,y)  ( (qcd_complex_16) {x.re * y.re - x.im * y.im, x.re * y.im + x.im * y.re } )
#endif

#define qcd_CS(x,a)  ( (qcd_complex_16) {x.re*(a), x.im*(a)})

#define qcd_CMUL_P(x,y)  ( (qcd_complex_16) {x.re * y.re - x.im * y.im, x.re * y.im + x.im * y.re } )
#define qcd_CMUL_M(x,y)  ( (qcd_complex_16) {x.im * y.im - x.re * y.re, - x.re * y.im - x.im * y.re} )

#define qcd_CMUL3_P(x,y,z)  ( qcd_CMUL_P(qcd_CMUL(x,y),z) )
#define qcd_CMUL3_M(x,y,z)  ( qcd_CMUL_M(qcd_CMUL(x,y),z) )

#define qcd_CMUL3_P_CS(x,y,z,c) ( qcd_CS(qcd_CMUL_P(qcd_CMUL(x,y),z),c) )
#define qcd_CMUL3_M_CS(x,y,z,c) ( qcd_CS(qcd_CMUL_M(qcd_CMUL(x,y),z),c) )

void p1_p2_p3_fun_oc(qcd_propagator *p1,qcd_propagator *p2,qcd_propagator *p3,qcd_geometry *geo,
		 qcd_complex_16 *block[MAX_PARTICLES_2][16][10], qcd_uint_4 t, qcd_uint_4 nz_counter[10],
		 qcd_int_2 cg5cg5b_ind[10][16*16][4],qcd_complex_16 cg5cg5b_val[10][16*16],qcd_uint_2 par ){
  
  qcd_uint_2 alpha1, alpha, beta1, beta;
  qcd_uint_2 a, b, c, a1, b1, c1;
  qcd_uint_4 v3,v;
  qcd_complex_16 ccg;
  qcd_complex_16 value;
       
  for(int lx=0; lx < geo->lL[1]; lx++)
    for(int ly=0; ly < geo->lL[2]; ly++)
      for(int lz=0; lz < geo->lL[3]; lz++)
	{
	  v3 = qcd_LEXIC0(lx,ly,lz,geo->lL);
	  v =  qcd_LEXIC(t,lx,ly,lz,geo->lL);

	  for(int gamma = 0 ; gamma < 4 ; gamma++)
	    for(int gamma1 = 0 ; gamma1 < 4 ; gamma1++){
	      
	      value = (qcd_complex_16) {0.,0.} ;

	      for(int inz = 0 ; inz < nz_counter[9] ; inz++){
		alpha = cg5cg5b_ind[9][inz][0];
		beta = cg5cg5b_ind[9][inz][1];
		beta1 = cg5cg5b_ind[9][inz][2];
		alpha1 = cg5cg5b_ind[9][inz][3];
		 
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

			ccg = qcd_CS(cg5cg5b_val[9][inz],qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]);

		  value = qcd_CADD(value,qcd_CMUL(ccg,qcd_CMUL3_P(p1->D[v][alpha][alpha1][a][a1],p2->D[v][beta][beta1][b][b1],p3->D[v][gamma][gamma1][c][c1]))) ;
   
  		      } // cc2
		  }	//cc1	       		       
	      } // suppressed dirac octet
	      block[par][gamma*4+gamma1][9][v3] = value ;
	      //tranfer from temp to block

	    } // close gamma gamma1
	} //space loop
  
  
}
