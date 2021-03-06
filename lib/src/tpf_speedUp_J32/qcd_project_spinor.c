#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>
#include "projectors_pb.h"

void Project_Spinor(qcd_complex_16 *block_in[4][4],qcd_complex_16 *block_out[5],qcd_int_2 list_projectors[5],qcd_geometry *geo){
  
  qcd_uint_4 v3;
  
  
 for(int lx=0; lx < geo->lL[1]; lx++)
      for(int ly=0; ly < geo->lL[2]; ly++)
        for(int lz=0; lz < geo->lL[3]; lz++)
           {
            v3 = qcd_LEXIC0(lx,ly,lz,geo->lL);
	    
	    for(int iproj = 0 ; iproj < 5 ; iproj++)
	      for(int alpha = 0 ; alpha < 4 ; alpha++)
		  for(int beta = 0 ; beta < 4 ; beta++){
		    block_out[iproj][v3] = qcd_CADD(block_out[iproj][v3],qcd_CMUL(PROJECTOR[list_projectors[iproj]][alpha][beta],block_in[beta][alpha][v3]));
		  }
	    
	   }
}
