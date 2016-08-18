#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>
#include <qcd_gamma.h>
#include <projectors_pb_static.h>

static void projector_32(qcd_complex_16 Proj[3][3][4][4]){
  qcd_complex_16 Proj_temp[3][3][4][4];

  for(int i = 0 ; i < 3 ; i++)
    for(int j = 0 ; j < 3 ; j++)
      for(int mu = 0 ; mu < 4 ; mu++)
	for(int nu = 0 ; nu < 4 ; nu++){
	  Proj_temp[i][j][mu][nu] = (qcd_complex_16) {0.,0.};
	  for(int lu = 0 ; lu < 4 ; lu++)
	    Proj_temp[i][j][mu][nu] = qcd_CADD(Proj_temp[i][j][mu][nu],qcd_CMUL(qcd_GAMMA[i+1][mu][lu],qcd_GAMMA[j+1][lu][nu]));
	}
  qcd_complex_16 delta;
  for(int i = 0 ; i < 3 ; i++)
    for(int j = 0 ; j < 3 ; j++)
      for(int mu = 0 ; mu < 4 ; mu++)
	for(int nu = 0 ; nu < 4 ; nu++){
	  if( (i==j) && (mu==nu) )
	    delta = (qcd_complex_16) {1.,0.};
	  else
	    delta = (qcd_complex_16) {0.,0.};
	  Proj[i][j][mu][nu] = qcd_CSUB(delta,qcd_CSCALE(Proj_temp[i][j][mu][nu],1./3.));
	}
}

#include "threep_seq_all_J32/qcd_proton.c"
#include "threep_seq_all_J32/qcd_neutron.c"
#include "threep_seq_all_J32/qcd_sigma_plus.c"
#include "threep_seq_all_J32/qcd_sigma_minus.c"
#include "threep_seq_all_J32/qcd_sigma_zero.c"
#include "threep_seq_all_J32/qcd_lambda.c"
#include "threep_seq_all_J32/qcd_octet_uds.c"

#include "threep_seq_all_J32/qcd_ksi_zero.c"
#include "threep_seq_all_J32/qcd_ksi_minus.c"
#include "threep_seq_all_J32/qcd_delta_plus_plus_J32.c"
#include "threep_seq_all_J32/qcd_delta_minus_J32.c"
#include "threep_seq_all_J32/qcd_delta_plus_J32.c"
#include "threep_seq_all_J32/qcd_delta_zero_J32.c"
#include "threep_seq_all_J32/qcd_sigma_star_plus_J32.c"
#include "threep_seq_all_J32/qcd_sigma_star_zero_J32.c"
#include "threep_seq_all_J32/qcd_sigma_star_minus_J32.c"
#include "threep_seq_all_J32/qcd_ksi_star_zero_J32.c"
#include "threep_seq_all_J32/qcd_ksi_star_minus_J32.c"
#include "threep_seq_all_J32/qcd_omega_J32.c"
#include "threep_seq_all_J32/qcd_decuplet_uds_J32.c"
