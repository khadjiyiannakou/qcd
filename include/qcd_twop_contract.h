#ifndef MAX_PARTICLES_2
#define MAX_PARTICLES_2 40
#endif

void project32(qcd_complex_16 *block[MAX_PARTICLES_2][16][10],qcd_geometry *geo);

void contract2pf(qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator *sprop, qcd_propagator *cprop, qcd_geometry *geo,
		 qcd_complex_16 *block[MAX_PARTICLES_2][16][10], qcd_uint_4 t, qcd_uint_4 nz_counter[10],qcd_int_2 cg5cg5b_ind[10][16*16][4],qcd_complex_16 cg5cg5b_val[10][16*16]);
