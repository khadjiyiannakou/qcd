#ifndef H_QCD_TPF_SPEEDUP
#define H_QCD_TPF_SPEEDUP 1

void Project_J32(qcd_complex_16 *block_in[4][4][9],qcd_complex_16 *block_out[4][4], qcd_geometry *geo);
void Project_Spinor(qcd_complex_16 *block_in[4][4],qcd_complex_16 *block_out[5],qcd_int_2 list_projectors[5],qcd_geometry *geo);

int tpf_sigma_star_zero_J32_sp(qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator *sprop, qcd_propagator *seq_uprop, qcd_propagator *seq_dprop,
			       qcd_propagator *seq_sprop, qcd_geometry *geo, qcd_int_2 cur_time, qcd_int_2 number_tsinks, qcd_complex_16 *U_part_sigmaSZ[5], 
			       qcd_complex_16 *D_part_sigmaSZ[5], qcd_complex_16 *S_part_sigmaSZ[5],qcd_int_4(*mom)[3], int num_momenta, qcd_uint_4 *x_src,
                               qcd_int_2 after_tcur,qcd_int_2 list_projectors[5]);

#endif
