#ifndef H_QCD_ALL_PARTICLES_J12
#define H_QCD_ALL_PARTICLES_J12 1



int tpf_delta_plus_plus_J12(qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator *sprop, qcd_propagator *seq_uprop, qcd_propagator *seq_dprop,
		  qcd_propagator *seq_sprop, qcd_geometry *geo, qcd_int_2 cur_time, qcd_int_2 number_tsinks, qcd_complex_16 *U_part_deltaPP[4][4], 
		  qcd_complex_16 *D_part_deltaPP[4][4], qcd_complex_16 *S_part_deltaPP[4][4],qcd_int_4(*mom)[3], int num_momenta, qcd_uint_4 *x_src, qcd_int_2 after_tcur);


int tpf_delta_plus_J12(qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator *sprop, qcd_propagator *seq_uprop, qcd_propagator *seq_dprop,
		   qcd_propagator *seq_sprop, qcd_geometry *geo, qcd_int_2 cur_time, qcd_int_2 number_tsinks, qcd_complex_16 *U_part_deltaP[4][4], 
		   qcd_complex_16 *D_part_deltaP[4][4], qcd_complex_16 *S_part_deltaP[4][4],qcd_int_4(*mom)[3], int num_momenta, qcd_uint_4 *x_src, qcd_int_2 after_tcur);

int tpf_delta_zero_J12(qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator *sprop, qcd_propagator *seq_uprop, qcd_propagator *seq_dprop,
		   qcd_propagator *seq_sprop, qcd_geometry *geo, qcd_int_2 cur_time, qcd_int_2 number_tsinks, qcd_complex_16 *U_part_deltaZ[4][4], 
		   qcd_complex_16 *D_part_deltaZ[4][4], qcd_complex_16 *S_part_deltaZ[4][4],qcd_int_4(*mom)[3], int num_momenta, qcd_uint_4 *x_src, qcd_int_2 after_tcur);

int tpf_delta_minus_J12(qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator *sprop, qcd_propagator *seq_uprop, qcd_propagator *seq_dprop,
		    qcd_propagator *seq_sprop, qcd_geometry *geo, qcd_int_2 cur_time, qcd_int_2 number_tsinks, qcd_complex_16 *U_part_deltaM[4][4], 
		    qcd_complex_16 *D_part_deltaM[4][4], qcd_complex_16 *S_part_deltaM[4][4],qcd_int_4(*mom)[3], int num_momenta, qcd_uint_4 *x_src, qcd_int_2 after_tcur);

int tpf_sigma_star_plus_J12(qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator *sprop, qcd_propagator *seq_uprop, qcd_propagator *seq_dprop,
		qcd_propagator *seq_sprop, qcd_geometry *geo, qcd_int_2 cur_time, qcd_int_2 number_tsinks, qcd_complex_16 *U_part_sigmaSP[4][4], 
		qcd_complex_16 *D_part_sigmaSP[4][4], qcd_complex_16 *S_part_sigmaSP[4][4],qcd_int_4(*mom)[3], int num_momenta, qcd_uint_4 *x_src, qcd_int_2 after_tcur);

int tpf_sigma_star_zero_J12(qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator *sprop, qcd_propagator *seq_uprop, qcd_propagator *seq_dprop,
		   qcd_propagator *seq_sprop, qcd_geometry *geo, qcd_int_2 cur_time, qcd_int_2 number_tsinks, qcd_complex_16 *U_part_sigmaSZ[4][4], 
		   qcd_complex_16 *D_part_sigmaSZ[4][4], qcd_complex_16 *S_part_sigmaSZ[4][4],qcd_int_4(*mom)[3], int num_momenta, qcd_uint_4 *x_src, qcd_int_2 after_tcur);

int tpf_sigma_star_minus_J12(qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator *sprop, qcd_propagator *seq_uprop, qcd_propagator *seq_dprop,
		qcd_propagator *seq_sprop, qcd_geometry *geo, qcd_int_2 cur_time, qcd_int_2 number_tsinks, qcd_complex_16 *U_part_sigmaSM[4][4], 
		qcd_complex_16 *D_part_sigmaSM[4][4], qcd_complex_16 *S_part_sigmaSM[4][4],qcd_int_4(*mom)[3], int num_momenta, qcd_uint_4 *x_src, qcd_int_2 after_tcur);

int tpf_ksi_star_zero_J12(qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator *sprop, qcd_propagator *seq_uprop, qcd_propagator *seq_dprop,
		      qcd_propagator *seq_sprop, qcd_geometry *geo, qcd_int_2 cur_time, qcd_int_2 number_tsinks, qcd_complex_16 *U_part_ksiSZ[4][4], 
		      qcd_complex_16 *D_part_ksiSZ[4][4], qcd_complex_16 *S_part_ksiSZ[4][4],qcd_int_4(*mom)[3], int num_momenta, qcd_uint_4 *x_src, qcd_int_2 after_tcur);

int tpf_ksi_star_minus_J12(qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator *sprop, qcd_propagator *seq_uprop, qcd_propagator *seq_dprop,
		       qcd_propagator *seq_sprop, qcd_geometry *geo, qcd_int_2 cur_time, qcd_int_2 number_tsinks, qcd_complex_16 *U_part_ksiSM[4][4], 
		       qcd_complex_16 *D_part_ksiSM[4][4], qcd_complex_16 *S_part_ksiSM[4][4],qcd_int_4(*mom)[3], int num_momenta, qcd_uint_4 *x_src, qcd_int_2 after_tcur);

int tpf_omega_J12(qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator *sprop, qcd_propagator *seq_uprop, qcd_propagator *seq_dprop,
	      qcd_propagator *seq_sprop, qcd_geometry *geo, qcd_int_2 cur_time, qcd_int_2 number_tsinks, qcd_complex_16 *U_part_omega[4][4], 
	      qcd_complex_16 *D_part_omega[4][4], qcd_complex_16 *S_part_omega[4][4],qcd_int_4(*mom)[3], int num_momenta, qcd_uint_4 *x_src, qcd_int_2 after_tcur);

int tpf_decuplet_uds_J12(qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator *sprop, qcd_propagator *seq_uprop, qcd_propagator *seq_dprop,
		     qcd_propagator *seq_sprop, qcd_geometry *geo, qcd_int_2 cur_time, qcd_int_2 number_tsinks, qcd_complex_16 *U_part_uds[4][4],
		     qcd_complex_16 *D_part_uds[4][4], qcd_complex_16 *S_part_uds[4][4],qcd_int_4(*mom)[3], int num_momenta, qcd_uint_4 *x_src, qcd_int_2 after_tcur);

#endif
