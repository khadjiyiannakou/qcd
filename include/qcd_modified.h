#ifndef H_QCD_MODIFIED
#define H_QCD_MODIFIED 1

qcd_real_8 qcd_calculatePlaquette_modified(qcd_gaugeField *u, MPI_Comm local_comm);
int qcd_communicateGaugeM_modified(qcd_gaugeField *u, MPI_Comm local_comm);
int qcd_communicateGaugePM_modified(qcd_gaugeField *u, MPI_Comm local_comm);
int qcd_apeSmear3d_modified(qcd_gaugeField *apeu, qcd_gaugeField *u, qcd_real_8 alpha, MPI_Comm local_comm);
int qcd_communicatePropagatorP_modified(qcd_propagator *p, MPI_Comm local_comm);
int qcd_gaussIteration3d_modified(qcd_vector *v, qcd_gaugeField *u, qcd_real_8 alpha, qcd_uint_4 t, MPI_Comm local_comm);
int qcd_communicateVectorPMGaugeP_modified(qcd_vector *v, qcd_gaugeField *u, MPI_Comm local_comm);


#endif
