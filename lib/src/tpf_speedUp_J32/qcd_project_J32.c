#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>
#include "gamma_J32.h"

void Project_J32(qcd_complex_16 *block_in[4][4][9],qcd_complex_16 *block_out[4][4], qcd_geometry *geo)
{
  
  qcd_uint_4 v3;
  qcd_complex_16 TermJ32[3][4][4];
  qcd_complex_16 Term3[4][4],Term6[4][4];
  
   for(int lx=0; lx < geo->lL[1]; lx++)
      for(int ly=0; ly < geo->lL[2]; ly++)
        for(int lz=0; lz < geo->lL[3]; lz++)
           {
            v3 = qcd_LEXIC0(lx,ly,lz,geo->lL);
                        
  
  	      for(int i = 0 ; i < 3 ; i++)
		for(int mu = 0 ; mu < 4 ; mu++)
		  for(int nu = 0 ; nu < 4 ; nu++){
		   TermJ32[i][mu][nu] = (qcd_complex_16) {0.,0.} ;
		  }

	      for(int mu = 0 ; mu < 4 ; mu++)
		for(int nu = 0 ; nu < 4 ; nu++)
		  for(int lu = 0 ; lu < 4 ; lu++){
		    
		    TermJ32[0][mu][nu] = qcd_CADD(TermJ32[0][mu][nu] , qcd_CMUL( qcd_GAMMA_J32[0][mu][lu] , qcd_CSUB(block_in[lu][nu][1*3+0][v3] , block_in[lu][nu][0*3+1][v3] ) ) );
		    TermJ32[1][mu][nu] = qcd_CADD(TermJ32[1][mu][nu] , qcd_CMUL( qcd_GAMMA_J32[1][mu][lu] , qcd_CSUB(block_in[lu][nu][2*3+0][v3] , block_in[lu][nu][0*3+2][v3] ) ) );
		    TermJ32[2][mu][nu] = qcd_CADD(TermJ32[2][mu][nu] , qcd_CMUL( qcd_GAMMA_J32[2][mu][lu] , qcd_CSUB(block_in[lu][nu][2*3+1][v3] , block_in[lu][nu][1*3+2][v3] ) ) );

		  }



	      for(int mu = 0 ; mu < 4 ; mu++)
		for(int nu = 0 ; nu < 4 ; nu++){
		  ///////////// upart //////////////
		  Term3[mu][nu] = qcd_CSCALE( qcd_CADD(qcd_CADD( block_in[mu][nu][0*3+0][v3], block_in[mu][nu][1*3+1][v3] ) , block_in[mu][nu][2*3+2][v3] ) , 1./3. );
		  Term6[mu][nu] = qcd_CSCALE( qcd_CADD(qcd_CADD(TermJ32[0][mu][nu] , TermJ32[1][mu][nu] ) , TermJ32[2][mu][nu] ) , 1./6.) ;
		  block_out[mu][nu][v3] = qcd_CSUB( Term3[mu][nu] , Term6[mu][nu] ) ;

		}	      
  
	}
  
}

