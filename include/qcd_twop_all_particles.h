#ifndef MAX_PARTICLES_2
#define MAX_PARTICLES_2 40
#endif

int octet(qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator *sprop, qcd_propagator *cprop, qcd_geometry *geo,
          qcd_complex_16 *block[MAX_PARTICLES_2], qcd_uint_4 t, qcd_uint_4 nz_counter[4],qcd_int_2 Pcg5cg5b_ind[4][16*16*16][6],qcd_complex_16 Pcg5cg5b_val[4][16*16*16]);

int decuplet(qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator *sprop, qcd_propagator *cprop, qcd_geometry *geo,
            qcd_complex_16 *block[MAX_PARTICLES_2], qcd_uint_4 t, qcd_uint_4 nz_counter[4],qcd_int_2 Pcg5cg5b_ind[4][16*16*16][6],qcd_complex_16 Pcg5cg5b_val[4][16*16*16]);

int octet_c(qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator *sprop, qcd_propagator *cprop, qcd_geometry *geo,
	     qcd_complex_16 *block[MAX_PARTICLES_2], qcd_uint_4 t, qcd_uint_4 nz_counter[4],qcd_int_2 Pcg5cg5b_ind[4][16*16*16][6],qcd_complex_16 Pcg5cg5b_val[4][16*16*16]);

int decuplet_c(qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator *sprop, qcd_propagator *cprop, qcd_geometry *geo,
	     qcd_complex_16 *block[MAX_PARTICLES_2], qcd_uint_4 t, qcd_uint_4 nz_counter[4],qcd_int_2 Pcg5cg5b_ind[4][16*16*16][6],qcd_complex_16 Pcg5cg5b_val[4][16*16*16]);

int octet_final6(qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator *sprop, qcd_propagator *cprop, qcd_geometry *geo,
		 qcd_complex_16 *block[MAX_PARTICLES_2], qcd_uint_4 t, qcd_uint_4 nz_counter[4],qcd_int_2 Pcg5cg5b_ind[4][16*16*16][6],qcd_complex_16 Pcg5cg5b_val[4][16*16*16]);

int decuplet_final4(qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator *sprop, qcd_propagator *cprop, qcd_geometry *geo,
		    qcd_complex_16 *block[MAX_PARTICLES_2], qcd_uint_4 t, qcd_uint_4 nz_counter[4],qcd_int_2 Pcg5cg5b_ind[4][16*16*16][6],qcd_complex_16 Pcg5cg5b_val[4][16*16*16]);
