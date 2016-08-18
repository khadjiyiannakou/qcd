#ifndef H_TEST
#define H_TEST
int tpf_sigma_zero_test(qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator *sprop, qcd_propagator *seq_uprop, qcd_propagator *seq_dprop,
			qcd_propagator *seq_sprop, qcd_geometry *geo, qcd_int_2 cur_time, qcd_int_2 number_tsinks, qcd_complex_16 *U_part_sigmaZ[4][4], 
			qcd_complex_16 *D_part_sigmaZ[4][4], qcd_complex_16 *S_part_sigmaZ[4][4],qcd_int_4(*mom)[3], int num_momenta, qcd_uint_4 *x_src, qcd_int_2 after_tcur);

#endif
