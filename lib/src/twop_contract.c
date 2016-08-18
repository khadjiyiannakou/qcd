#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>
#include <gamma_J32.h>
#include "twop_contract.h"
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


/*
#define p1_p2_p1(p1,p2)  ( qcd_CADD(qcd_CMUL3_M(p1->D[v][alpha][gamma1][a][c1],  p2->D[v][beta][beta1][b][b1],  p1->D[v][gamma][alpha1][c][a1]), \
				   qcd_CMUL3_P(p1->D[v][alpha][alpha1][a][a1],  p2->D[v][beta][beta1][b][b1],  p1->D[v][gamma][gamma1][c][c1]) ) )  

#define lambdas(p1,p2,p3)   ( qcd_CADD(qcd_CADD(qcd_CADD(qcd_CADD(qcd_CADD(qcd_CADD(qcd_CADD(qcd_CADD(\
					   qcd_CMUL3_P_CS(p1->D[v][alpha][gamma1][a][c1],p2->D[v][gamma][alpha1][c][a1],p3->D[v][beta][beta1][b][b1],1./6.)\
					   ,qcd_CMUL3_P_CS(p1->D[v][gamma][alpha1][c][a1],p2->D[v][alpha][gamma1][a][c1],p3->D[v][beta][beta1][b][b1],1./6.) ),\
					   qcd_CMUL3_P_CS(p1->D[v][alpha][alpha1][a][a1],p2->D[v][gamma][gamma1][c][c1],p3->D[v][beta][beta1][b][b1],1./6.) )\
					   ,qcd_CMUL3_P_CS(p1->D[v][gamma][gamma1][c][c1],p2->D[v][alpha][alpha1][a][a1],p3->D[v][beta][beta1][b][b1],1./6.) ),\
					   qcd_CMUL3_M_CS(p1->D[v][gamma][alpha1][c][a1],p2->D[v][alpha][beta1][a][b1],p3->D[v][beta][gamma1][b][c1],2./6.) )\
					   ,qcd_CMUL3_M_CS(p1->D[v][alpha][alpha1][a][a1],p2->D[v][gamma][beta1][c][b1],p3->D[v][beta][gamma1][b][c1],2./6.) ),\
					   qcd_CMUL3_M_CS(p1->D[v][alpha][gamma1][a][c1],p2->D[v][beta][alpha1][b][a1],p3->D[v][gamma][beta1][c][b1],2./6.) )\
					   ,qcd_CMUL3_M_CS(p1->D[v][alpha][alpha1][a][a1],p2->D[v][beta][gamma1][b][c1],p3->D[v][gamma][beta1][c][b1],2./6.) ),\
					   qcd_CMUL3_P_CS(p1->D[v][alpha][alpha1][a][a1],p2->D[v][beta][beta1][b][b1],p3->D[v][gamma][gamma1][c][c1],4./6.) ) )

#define p132_p231(p1,p2,p3) ( qcd_CADD(qcd_CADD(qcd_CADD(\
					  qcd_CMUL3_M_CS(p1->D[v][alpha][gamma1][a][c1],p2->D[v][gamma][alpha1][c][a1],p3->D[v][beta][beta1][b][b1],0.5)\
					  ,qcd_CMUL3_M_CS(p1->D[v][gamma][alpha1][c][a1],p2->D[v][alpha][gamma1][a][c1],p3->D[v][beta][beta1][b][b1],0.5) ),\
					  qcd_CMUL3_P_CS(p1->D[v][alpha][alpha1][a][a1],p2->D[v][gamma][gamma1][c][c1],p3->D[v][beta][beta1][b][b1],0.5))\
					  ,qcd_CMUL3_P_CS(p1->D[v][gamma][gamma1][c][c1],p2->D[v][alpha][alpha1][a][a1],p3->D[v][beta][beta1][b][b1],0.5))\
					  )

#define p1_p1_p1(p1)        ( qcd_CADD(qcd_CADD(qcd_CADD(qcd_CADD(qcd_CADD(\
					   qcd_CMUL3_M(p1->D[v][alpha][gamma1][a][c1],p1->D[v][beta][beta1][b][b1],p1->D[v][gamma][alpha1][c][a1])\
					  ,qcd_CMUL3_P(p1->D[v][alpha][beta1][a][b1],p1->D[v][beta][gamma1][b][c1],p1->D[v][gamma][alpha1][c][a1]))\
					  ,qcd_CMUL3_P(p1->D[v][alpha][gamma1][a][c1],p1->D[v][beta][alpha1][b][a1],p1->D[v][gamma][beta1][c][b1]))\
					  ,qcd_CMUL3_M(p1->D[v][alpha][alpha1][a][a1],p1->D[v][beta][gamma1][b][c1],p1->D[v][gamma][beta1][c][b1]))\
					  ,qcd_CMUL3_M(p1->D[v][alpha][beta1][a][b1],p1->D[v][beta][alpha1][b][a1],p1->D[v][gamma][gamma1][c][c1]))\
					  ,qcd_CMUL3_P(p1->D[v][alpha][alpha1][a][a1],p1->D[v][beta][beta1][b][b1],p1->D[v][gamma][gamma1][c][c1]))\
					  )
					  
#define p1_p2_p3(p1,p2,p3)  ( qcd_CMUL3_P(p1->D[v][alpha][alpha1][a][a1],p2->D[v][beta][beta1][b][b1],p3->D[v][gamma][gamma1][c][c1]) )

#define delta(p1,p2)        ( qcd_CADD(qcd_CADD(qcd_CADD(qcd_CADD(qcd_CADD(qcd_CADD(qcd_CADD(\
					     qcd_CMUL3_M_CS(p1->D[v][alpha][beta1][a][b1],p1->D[v][beta][alpha1][b][a1],p2->D[v][gamma][gamma1][c][c1],1./3.) \
					    ,qcd_CMUL3_P_CS(p1->D[v][alpha][alpha1][a][a1],p1->D[v][beta][beta1][b][b1],p2->D[v][gamma][gamma1][c][c1],1./3.))\
					    ,qcd_CMUL3_P_CS(p1->D[v][alpha][beta1][a][b1],p1->D[v][gamma][alpha1][c][a1],p2->D[v][beta][gamma1][b][c1],2./3.))\
					    ,qcd_CMUL3_P_CS(p1->D[v][alpha][gamma1][a][c1],p1->D[v][beta][alpha1][b][a1],p2->D[v][gamma][beta1][c][b1],2./3.))\
					    ,qcd_CMUL3_M_CS(p1->D[v][alpha][alpha1][a][a1],p1->D[v][beta][gamma1][b][c1],p2->D[v][gamma][beta1][c][b1],2./3.))\
					    ,qcd_CMUL3_M_CS(p1->D[v][alpha][alpha1][a][a1],p1->D[v][gamma][beta1][c][b1],p2->D[v][beta][gamma1][b][c1],2./3.))\
					    ,qcd_CMUL3_M_CS(p1->D[v][alpha][gamma1][a][c1],p1->D[v][gamma][alpha1][c][a1],p2->D[v][beta][beta1][b][b1],4./3.))\
					    ,qcd_CMUL3_P_CS(p1->D[v][alpha][alpha1][a][a1],p1->D[v][gamma][gamma1][c][c1],p2->D[v][beta][beta1][b][b1],4./3.))\
					    )
					    
#define sigmas4(p1,p2)     ( qcd_CADD(qcd_CADD(qcd_CADD(qcd_CADD(qcd_CADD(qcd_CADD(qcd_CADD(\
					     qcd_CMUL3_M_CS(p1->D[v][alpha][beta1][a][b1],p1->D[v][beta][alpha1][b][a1],p2->D[v][gamma][gamma1][c][c1],1./3.) \
					    ,qcd_CMUL3_P_CS(p1->D[v][alpha][alpha1][a][a1],p1->D[v][beta][beta1][b][b1],p2->D[v][gamma][gamma1][c][c1],1./3.))\
					    ,qcd_CMUL3_M_CS(p1->D[v][alpha][gamma1][a][c1],p1->D[v][beta][beta1][b][b1],p2->D[v][gamma][alpha1][c][a1],2./3.))\
					    ,qcd_CMUL3_P_CS(p1->D[v][alpha][beta1][a][b1],p1->D[v][beta][gamma1][b][c1],p2->D[v][gamma][alpha1][c][a1],2./3.))\
					    ,qcd_CMUL3_M_CS(p1->D[v][beta][beta1][b][b1],p1->D[v][gamma][alpha1][c][a1],p2->D[v][alpha][gamma1][a][c1],2./3.))\
					    ,qcd_CMUL3_P_CS(p1->D[v][beta][alpha1][b][a1],p1->D[v][gamma][beta1][c][b1],p2->D[v][alpha][gamma1][a][c1],2./3.))\
					    ,qcd_CMUL3_M_CS(p1->D[v][beta][gamma1][b][c1],p1->D[v][gamma][beta1][c][b1],p2->D[v][alpha][alpha1][a][a1],4./3.))\
					    ,qcd_CMUL3_P_CS(p1->D[v][beta][beta1][b][b1],p1->D[v][gamma][gamma1][c][c1],p2->D[v][alpha][alpha1][a][a1],4./3.))\
					    )
					    
#define sigmas2(p1,p2,p3)  ( qcd_CADD(qcd_CADD(qcd_CADD(qcd_CADD(qcd_CADD(qcd_CADD(qcd_CADD(qcd_CADD(\
					      qcd_CMUL3_P_CS(p1->D[v][beta][gamma1][b][c1],p2->D[v][gamma][alpha1][c][a1],p3->D[v][alpha][beta1][a][b1],2./3.) \
					     ,qcd_CMUL3_P_CS(p1->D[v][alpha][beta1][a][b1],p2->D[v][beta][gamma1][b][c1],p3->D[v][gamma][alpha1][c][a1],2./3.))\
					     ,qcd_CMUL3_P_CS(p1->D[v][gamma][alpha1][c][a1],p2->D[v][alpha][beta1][a][b1],p3->D[v][beta][gamma1][b][c1],2./3.))\
					     ,qcd_CMUL3_P_CS(p1->D[v][beta][alpha1][b][a1],p2->D[v][gamma][beta1][c][b1],p3->D[v][alpha][gamma1][a][c1],2./3.))\
					     ,qcd_CMUL3_P_CS(p1->D[v][alpha][gamma1][a][c1],p2->D[v][beta][alpha1][b][a1],p3->D[v][gamma][beta1][c][b1],2./3.))\
					     ,qcd_CMUL3_P_CS(p1->D[v][gamma][beta1][c][b1],p2->D[v][alpha][gamma1][a][c1],p3->D[v][beta][alpha1][b][a1],2./3.))\
					     ,qcd_CMUL3_P_CS(p1->D[v][beta][beta1][b][b1],p2->D[v][gamma][gamma1][c][c1],p3->D[v][alpha][alpha1][a][a1],2./3.))\
					     ,qcd_CMUL3_P_CS(p1->D[v][alpha][alpha1][a][a1],p2->D[v][beta][beta1][b][b1],p3->D[v][gamma][gamma1][c][c1],2./3.))\
					     ,qcd_CMUL3_P_CS(p1->D[v][gamma][gamma1][c][c1],p2->D[v][alpha][alpha1][a][a1],p3->D[v][beta][beta1][b][b1],2./3.))\
					     )
*/
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////














//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void project32(qcd_complex_16 *block[MAX_PARTICLES_2][16][10],qcd_geometry *geo){
  
  qcd_complex_16 TermJ32[3][4][4], Term3[4][4], Term6[4][4];  
  qcd_uint_4 v3;	
 
  for(int lx=0; lx < geo->lL[1]; lx++)
    for(int ly=0; ly < geo->lL[2]; ly++)
      for(int lz=0; lz < geo->lL[3]; lz++)
	{	 
	  v3 = qcd_LEXIC0(lx,ly,lz,geo->lL);
	  for(int iparticle = 0 ; iparticle < MAX_PARTICLES_2 ; iparticle++)
	    {
	      if( ((iparticle >= 8) && (iparticle <= 17)) || ((iparticle >= 24) && (iparticle <= 29)) || ((iparticle >= 36) && (iparticle <= 39)) ){
		  
                for(int i = 0 ; i < 3 ; i++)
		  for(int mu = 0 ; mu < 4 ; mu++)
		    for(int nu = 0 ; nu < 4 ; nu++)
		      TermJ32[i][mu][nu] = (qcd_complex_16) {0.,0.} ;

		for(int mu = 0 ; mu < 4 ; mu++)
		  for(int nu = 0 ; nu < 4 ; nu++)
		    for(int lu = 0 ; lu < 4 ; lu++){
		      TermJ32[0][mu][nu] = qcd_CADD(TermJ32[0][mu][nu] , qcd_CMUL( qcd_GAMMA_J32[0][mu][lu] , qcd_CSUB(block[iparticle][lu*4+nu][1*3+0][v3] , block[iparticle][lu*4+nu][0*3+1][v3] ) ) );
		      TermJ32[1][mu][nu] = qcd_CADD(TermJ32[1][mu][nu] , qcd_CMUL( qcd_GAMMA_J32[1][mu][lu] , qcd_CSUB(block[iparticle][lu*4+nu][2*3+0][v3] , block[iparticle][lu*4+nu][0*3+2][v3] ) ) );
		      TermJ32[2][mu][nu] = qcd_CADD(TermJ32[2][mu][nu] , qcd_CMUL( qcd_GAMMA_J32[2][mu][lu] , qcd_CSUB(block[iparticle][lu*4+nu][2*3+1][v3] , block[iparticle][lu*4+nu][1*3+2][v3] ) ) );
		    }

		for(int mu = 0 ; mu < 4 ; mu++)
		  for(int nu = 0 ; nu < 4 ; nu++){
		    Term3[mu][nu] = qcd_CSCALE( qcd_CADD(qcd_CADD( block[iparticle][mu*4+nu][0*3+0][v3], block[iparticle][mu*4+nu][1*3+1][v3] ) , block[iparticle][mu*4+nu][2*3+2][v3] ) , 1./3. );
		    Term6[mu][nu] = qcd_CSCALE( qcd_CADD(qcd_CADD(TermJ32[0][mu][nu] , TermJ32[1][mu][nu] ) , TermJ32[2][mu][nu] ) , 1./6.) ;
		    block[iparticle][mu*4+nu][9][v3] = qcd_CSUB( Term3[mu][nu] , Term6[mu][nu] )  ;
		    
		  }
		  
	      } // close if
	    } // close particle
	} // close space
}


void contract2pf(qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator *sprop, qcd_propagator *cprop, qcd_geometry *geo,
		 qcd_complex_16 *block[MAX_PARTICLES_2][16][10], qcd_uint_4 t, qcd_uint_4 nz_counter[10],qcd_int_2 cg5cg5b_ind[10][16*16][4],qcd_complex_16 cg5cg5b_val[10][16*16])
{
  
 

			///////////////////////////////////////////////////////////////////			
			// Proton
			p1_p2_p1_fun_oc(uprop,dprop,geo,block,t,nz_counter,cg5cg5b_ind,cg5cg5b_val,0);
			
			//Neutron
			p1_p2_p1_fun_oc(dprop,uprop,geo,block,t,nz_counter,cg5cg5b_ind,cg5cg5b_val,1);	      
			
			//Lambda
			 lambdas_fun(uprop,dprop,sprop,geo,block,t,nz_counter,cg5cg5b_ind,cg5cg5b_val,2);
			 
			//SIGMA_PLUS
			 p1_p2_p1_fun_oc(uprop,sprop,geo,block,t,nz_counter,cg5cg5b_ind,cg5cg5b_val,3);
		      
			//SIGMA_ZERO
			 p132_p231_fun(uprop,dprop,sprop,geo,block,t,nz_counter,cg5cg5b_ind,cg5cg5b_val,4);
			
			//SIGMA_MINUS
			 p1_p2_p1_fun_oc(dprop,sprop,geo,block,t,nz_counter,cg5cg5b_ind,cg5cg5b_val,5);
				      
			//KSI_ZERO
			p1_p2_p1_fun_oc(sprop,uprop,geo,block,t,nz_counter,cg5cg5b_ind,cg5cg5b_val,6);
			  
			//KSI_MINUS
			p1_p2_p1_fun_oc(sprop,dprop,geo,block,t,nz_counter,cg5cg5b_ind,cg5cg5b_val,7);
		      
			//LAMBDA_PLUS_C
			lambdas_fun(uprop,dprop,cprop,geo,block,t,nz_counter,cg5cg5b_ind,cg5cg5b_val,18);
		      
			//SIGMA_PLUS_PLUS_C
			 p1_p2_p1_fun_oc(uprop,cprop,geo,block,t,nz_counter,cg5cg5b_ind,cg5cg5b_val,19);
		      
			//SIGMA_PLUS_C
			  p132_p231_fun(uprop,dprop,cprop,geo,block,t,nz_counter,cg5cg5b_ind,cg5cg5b_val,20);
		      
			//SIGMA_ZERO_C
			  p1_p2_p1_fun_oc(dprop,cprop,geo,block,t,nz_counter,cg5cg5b_ind,cg5cg5b_val,21);
		      
			//KSI_PLUS_PLUS_C_C
			  p1_p2_p1_fun_oc(cprop,uprop,geo,block,t,nz_counter,cg5cg5b_ind,cg5cg5b_val,22);
			  
			//KSI_PLUS_C_C
			  p1_p2_p1_fun_oc(cprop,dprop,geo,block,t,nz_counter,cg5cg5b_ind,cg5cg5b_val,23);
		      
			//OMEGA_PLUS_C_C
			   p1_p2_p1_fun_oc(cprop,sprop,geo,block,t,nz_counter,cg5cg5b_ind,cg5cg5b_val,35);
		      
			//KSI_PLUS_C
			   p1_p2_p3_fun_oc(uprop,sprop,cprop,geo,block,t,nz_counter,cg5cg5b_ind,cg5cg5b_val,34);
			
			//KSI_PRIME_PLUS_C
			   p132_p231_fun(uprop,sprop,cprop,geo,block,t,nz_counter,cg5cg5b_ind,cg5cg5b_val,33);
		      
			//KSI_ZERO_C
			   p1_p2_p3_fun_oc(dprop,sprop,cprop,geo,block,t,nz_counter,cg5cg5b_ind,cg5cg5b_val,32);
		      
			//KSI_PRIME_ZERO_C
			   p132_p231_fun(dprop,sprop,cprop,geo,block,t,nz_counter,cg5cg5b_ind,cg5cg5b_val,31);
		      
			//OMEGA_ZERO_C
			   p1_p2_p1_fun_oc(sprop,cprop,geo,block,t,nz_counter,cg5cg5b_ind,cg5cg5b_val,30);
			

			    
			    //DELTA_PLUS_PLUS
			    p1_p1_p1_fun(uprop,geo,block,t,nz_counter,cg5cg5b_ind,cg5cg5b_val,8);
			    
			    //DELTA_PLUS
			    delta_fun(uprop,dprop,geo,block,t,nz_counter,cg5cg5b_ind,cg5cg5b_val,9);
			    
			    //DELTA_ZERO
			    delta_fun(dprop,uprop,geo,block,t,nz_counter,cg5cg5b_ind,cg5cg5b_val,10);
			    
			    //DELTA_MINUS
			     p1_p1_p1_fun(dprop,geo,block,t,nz_counter,cg5cg5b_ind,cg5cg5b_val,11);

			    //SIGMA_STAR_PLUS
			     sigmas4_fun(uprop,sprop,geo,block,t,nz_counter,cg5cg5b_ind,cg5cg5b_val,12);
			    
			    //SIGMA_STAR_ZERO
			    sigmas2_fun(uprop,dprop,sprop,geo,block,t,nz_counter,cg5cg5b_ind,cg5cg5b_val,13);
			    
			    //SIGMA_STAR_MINUS
			    sigmas4_fun(dprop,sprop,geo,block,t,nz_counter,cg5cg5b_ind,cg5cg5b_val,14);

			    //KSI_STAR_ZERO
			     p1_p2_p1_fun_dec(sprop,uprop,geo,block,t,nz_counter,cg5cg5b_ind,cg5cg5b_val,15);
			    
			    //KSI_STAR_MINUS
			     p1_p2_p1_fun_dec(sprop,dprop,geo,block,t,nz_counter,cg5cg5b_ind,cg5cg5b_val,16);
			    
			    //OMEGA
			     p1_p1_p1_fun(sprop,geo,block,t,nz_counter,cg5cg5b_ind,cg5cg5b_val,17);
			    
			    //SIGMA_STAR_PLUS_PLUS_C
			      sigmas4_fun(uprop,cprop,geo,block,t,nz_counter,cg5cg5b_ind,cg5cg5b_val,24);
			    
			    //SIGMA_STAR_PLUS_C
			      sigmas2_fun(uprop,dprop,cprop,geo,block,t,nz_counter,cg5cg5b_ind,cg5cg5b_val,25);
			  
			    //SIGMA_STAR_ZERO_C
			      sigmas4_fun(dprop,cprop,geo,block,t,nz_counter,cg5cg5b_ind,cg5cg5b_val,26);
			    
			    //KSI_STAR_PLUS_PLUS_C_C
			      p1_p2_p1_fun_dec(cprop,uprop,geo,block,t,nz_counter,cg5cg5b_ind,cg5cg5b_val,27);
			      
			    //KSI_STAR_PLUS_C_C
			      p1_p2_p1_fun_dec(cprop,dprop,geo,block,t,nz_counter,cg5cg5b_ind,cg5cg5b_val,28);
				    
			    //OMEGA_PLUS_PLUS_C_C_C
			      p1_p1_p1_fun(cprop,geo,block,t,nz_counter,cg5cg5b_ind,cg5cg5b_val,29);
			      
			    //KSI_STAR_PLUS_C
			      p1_p2_p3_fun_dec(sprop,uprop,cprop,geo,block,t,nz_counter,cg5cg5b_ind,cg5cg5b_val,38);
			    
			    
			    //KSI_STAR_ZERO_C
			       p1_p2_p3_fun_dec(sprop,dprop,cprop,geo,block,t,nz_counter,cg5cg5b_ind,cg5cg5b_val,37);
			    
			    //OMEGA_STAR_ZERO_C
			       p1_p2_p1_fun_dec(sprop,cprop,geo,block,t,nz_counter,cg5cg5b_ind,cg5cg5b_val,36);
			    
			    //OMEGA_STAR_PLUS_C_C
			        p1_p2_p1_fun_dec(cprop,sprop, geo,block,t,nz_counter,cg5cg5b_ind,cg5cg5b_val,39);
			    
  project32(block,geo);
}
