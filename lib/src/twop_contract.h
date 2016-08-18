void p1_p2_p1_fun_oc(qcd_propagator *p1,qcd_propagator *p2,qcd_geometry *geo,
		 qcd_complex_16 *block[MAX_PARTICLES_2][16][10], qcd_uint_4 t, qcd_uint_4 nz_counter[10],
		 qcd_int_2 cg5cg5b_ind[10][16*16][4],qcd_complex_16 cg5cg5b_val[10][16*16],qcd_uint_2 par);

void p1_p2_p1_fun_dec(qcd_propagator *p1,qcd_propagator *p2,qcd_geometry *geo,
		 qcd_complex_16 *block[MAX_PARTICLES_2][16][10], qcd_uint_4 t, qcd_uint_4 nz_counter[10],
		 qcd_int_2 cg5cg5b_ind[10][16*16][4],qcd_complex_16 cg5cg5b_val[10][16*16],qcd_uint_2 par);

void lambdas_fun(qcd_propagator *p1,qcd_propagator *p2,qcd_propagator *p3,qcd_geometry *geo,
		 qcd_complex_16 *block[MAX_PARTICLES_2][16][10], qcd_uint_4 t, qcd_uint_4 nz_counter[10],
		 qcd_int_2 cg5cg5b_ind[10][16*16][4],qcd_complex_16 cg5cg5b_val[10][16*16],qcd_uint_2 par);

void p132_p231_fun(qcd_propagator *p1,qcd_propagator *p2,qcd_propagator *p3,qcd_geometry *geo,
		 qcd_complex_16 *block[MAX_PARTICLES_2][16][10], qcd_uint_4 t, qcd_uint_4 nz_counter[10],
		 qcd_int_2 cg5cg5b_ind[10][16*16][4],qcd_complex_16 cg5cg5b_val[10][16*16],qcd_uint_2 par);

void p1_p1_p1_fun(qcd_propagator *p1,qcd_geometry *geo,
		 qcd_complex_16 *block[MAX_PARTICLES_2][16][10], qcd_uint_4 t,
		 qcd_uint_4 nz_counter[10],qcd_int_2 cg5cg5b_ind[10][16*16][4],qcd_complex_16 cg5cg5b_val[10][16*16],qcd_uint_2 par);

void p1_p2_p3_fun_oc(qcd_propagator *p1,qcd_propagator *p2,qcd_propagator *p3,qcd_geometry *geo,
		 qcd_complex_16 *block[MAX_PARTICLES_2][16][10], qcd_uint_4 t, qcd_uint_4 nz_counter[10],
		 qcd_int_2 cg5cg5b_ind[10][16*16][4],qcd_complex_16 cg5cg5b_val[10][16*16],qcd_uint_2 par );

void p1_p2_p3_fun_dec(qcd_propagator *p1,qcd_propagator *p2,qcd_propagator *p3,qcd_geometry *geo,
		 qcd_complex_16 *block[MAX_PARTICLES_2][16][10], qcd_uint_4 t, qcd_uint_4 nz_counter[10],
		 qcd_int_2 cg5cg5b_ind[10][16*16][4],qcd_complex_16 cg5cg5b_val[10][16*16],qcd_uint_2 par );

void delta_fun(qcd_propagator *p1,qcd_propagator *p2,qcd_geometry *geo,
		 qcd_complex_16 *block[MAX_PARTICLES_2][16][10], qcd_uint_4 t, qcd_uint_4 nz_counter[10],
		 qcd_int_2 cg5cg5b_ind[10][16*16][4],qcd_complex_16 cg5cg5b_val[10][16*16],qcd_uint_2 par);

void sigmas4_fun(qcd_propagator *p1,qcd_propagator *p2,qcd_geometry *geo,
		 qcd_complex_16 *block[MAX_PARTICLES_2][16][10], qcd_uint_4 t, qcd_uint_4 nz_counter[10],
		 qcd_int_2 cg5cg5b_ind[10][16*16][4],qcd_complex_16 cg5cg5b_val[10][16*16],qcd_uint_2 par);

void sigmas2_fun(qcd_propagator *p1,qcd_propagator *p2,qcd_propagator *p3,qcd_geometry *geo,
		 qcd_complex_16 *block[MAX_PARTICLES_2][16][10], qcd_uint_4 t, qcd_uint_4 nz_counter[10],
		 qcd_int_2 cg5cg5b_ind[10][16*16][4],qcd_complex_16 cg5cg5b_val[10][16*16],qcd_uint_2 par);