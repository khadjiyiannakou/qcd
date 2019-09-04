/* qcd_observables.h
 *
 * header file for qcd_observables.c
 *
 * Tomasz Korzec 2009
 *******************************************/
 
#ifndef H_QCD_OBSERVABLES
#define H_QCD_OBSERVABLES 1  

/* prototypes */
qcd_real_8 qcd_calculatePlaquette(qcd_gaugeField *u);
void qcd_calculateGluonLoops(qcd_gaugeField *u, qcd_complex_16* gLoops, qcd_int_4 (*mom)[3], int Nmom);
void qcd_calculateGluonLoops_T00g_Sg(qcd_gaugeField *u, qcd_complex_16* gLoops_T00g, qcd_complex_16* gLoops_Sg, qcd_int_4 (*mom)[3], int Nmom);
#endif
