#ifndef H_QCD_BLAS_SPECIAL
#define H_QCD_BLAS_SPECIAL 1

void qcd_copyVectorPropagator_timerange(qcd_vector *vec, qcd_propagator *prop, qcd_uint_2 nu, qcd_uint_2 c2, qcd_uint_4 cur_time, qcd_uint_4 after_tcur, qcd_uint_4 number_tsinks);
void qcd_copyPropagatorVector_timerange(qcd_propagator *prop, qcd_vector *vec, qcd_uint_2 nu, qcd_uint_2 c2, qcd_uint_4 cur_time, qcd_uint_4 after_tcur, qcd_uint_4 number_tsinks);
void qcd_tranformPropagatorPhysicalPlus(qcd_propagator *prop, qcd_geometry *geo, qcd_uint_4 cur_time, qcd_uint_4 after_tcur, qcd_uint_4 number_tsinks);
void qcd_tranformPropagatorPhysicalMinus(qcd_propagator *prop, qcd_geometry *geo, qcd_uint_4 cur_time, qcd_uint_4 after_tcur, qcd_uint_4 number_tsinks);

#endif
