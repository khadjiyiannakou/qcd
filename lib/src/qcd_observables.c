/* qcd_observables.c
 * 
 * computes simple observables, e.g. plaquette
 *
 * Tomasz Korzec 2009
 ***************************************************/
 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>

 
qcd_real_8 qcd_calculatePlaquette(qcd_gaugeField *u)
{
   qcd_communicateGaugeM(u);
   qcd_waitall(u->geo);
   qcd_complex_16 plaq[3][3],tmp[3][3];
   qcd_real_8 meanplaq=0;
   qcd_real_8 result=0;
   qcd_uint_4 l;
   qcd_uint_2 mu,nu;
   
   for(l=0; l<u->geo->lV; l++)
   for(mu=0; mu<3; mu++)
   for(nu=mu+1; nu<4; nu++)
   {
      qcd_MUL3x3(plaq,u->D[l][mu],u->D[u->geo->plus[l][mu]][nu]);
      qcd_MULADJOINT3x3(tmp,plaq,u->D[u->geo->plus[l][nu]][mu]);
      qcd_MULADJOINT3x3(plaq,tmp,u->D[l][nu]);
      meanplaq += qcd_SU3TRACER(plaq);
   }
   MPI_Allreduce(&meanplaq, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   return(result/(u->geo->V * 3.0 * 6.0)); // volume * N_color * number of mu-nu combinations
}//end qcd_calculatePlaquette 

void qcd_calculateGluonLoops(qcd_gaugeField *u, qcd_complex_16* gLoops, qcd_int_4 (*mom)[3], int Nmom)
{
  /*
    Note that the size of the gLoops is  NT * Nmom where Nmom is the number of momenta and NT is the temporal extent
   */
   qcd_communicateGaugeM(u);
   qcd_waitall(u->geo);
   qcd_complex_16 plaq[3][3],tmp[3][3];
   qcd_complex_16 smu;
   qcd_int_8 l;
   qcd_uint_4 lt;
   qcd_uint_2 mu,nu;
   qcd_uint_2 lx,ly,lz;
   qcd_uint_2 x,y,z;
   int v3;
   qcd_complex_16 *block;
   qcd_complex_16 gl,gl2;
   qcd_complex_16 C2;
   qcd_real_8 argPhase;

   block = (qcd_complex_16*) malloc(u->geo->lV3*sizeof(qcd_complex_16));
   if(block == NULL){
     fprintf(stderr,"process %i: Error in qcd_initPropagator! Out of memory\n",u->geo->myid);
     exit(EXIT_FAILURE);
   }

   if(u->geo->L[0] != u->geo->lL[0]){
     fprintf(stderr,"Error! function which does the gluon loops accepts only one task in the temporal direction");
     exit(EXIT_FAILURE);
   }

   for(lt=0; lt < u->geo->lL[0]; lt++){

     for(v3=0; v3 < u->geo->lV3; v3++) //set blocks to zero for each time-slice
       block[v3]= (qcd_complex_16) {0,0};


     for(lx=0;lx<u->geo->lL[1];lx++)
       for(ly=0;ly<u->geo->lL[2];ly++)
	 for(lz=0;lz<u->geo->lL[3];lz++){
	   l=qcd_LEXIC(lt,lx,ly,lz,u->geo->lL);
	   v3 = qcd_LEXIC0(lx,ly,lz,u->geo->lL);

	   // term1
	   smu.re=0; smu.im=0;
	   for(mu=1; mu<4; mu++)
	     {
	       qcd_MUL3x3(plaq,u->D[l][mu],u->D[u->geo->plus[l][mu]][0]);
	       qcd_MULADJOINT3x3(tmp,plaq,u->D[u->geo->plus[l][0]][mu]);
	       qcd_MULADJOINT3x3(plaq,tmp,u->D[l][0]);
	       smu.re += qcd_SU3TRACER(plaq);
	     }
	   block[v3] = qcd_CADD(block[v3],smu);

	   //term2
	   smu.re=0; smu.im=0;
	   for(mu=1; mu<4; mu++)
	     for(nu=mu+1;nu<4;nu++)
	       {
		 qcd_MUL3x3(plaq,u->D[l][mu],u->D[u->geo->plus[l][mu]][nu]);
		 qcd_MULADJOINT3x3(tmp,plaq,u->D[u->geo->plus[l][nu]][mu]);
		 qcd_MULADJOINT3x3(plaq,tmp,u->D[l][nu]);
		 smu.re -= qcd_SU3TRACER(plaq);
	       }
	   block[v3] = qcd_CADD(block[v3],smu);	   
	 } // create the gluon loop for each lattice point

     // fourier transform
     for(int imom = 0 ; imom < Nmom ; imom++){
       gl = (qcd_complex_16) {0,0};
       for(lx=0;lx<u->geo->lL[1];lx++)
	 for(ly=0;ly<u->geo->lL[2];ly++)
	   for(lz=0;lz<u->geo->lL[3];lz++){
	     v3 = qcd_LEXIC0(lx,ly,lz,u->geo->lL);
	     x=lx+u->geo->Pos[1]*u->geo->lL[1];
	     y=ly+u->geo->Pos[2]*u->geo->lL[2];
	     z=lz+u->geo->Pos[3]*u->geo->lL[3];
	     argPhase = (((double) mom[imom][0]*x)/u->geo->L[1] + ((double) mom[imom][1]*y)/u->geo->L[2] + ((double) mom[imom][2]*z)/u->geo->L[3])*2*M_PI;
	     C2=(qcd_complex_16) {cos(argPhase), -sin(argPhase)}; // note the sign convention of the gluon loop FT is the same as the quark loops 
	     gl=qcd_CADD(gl, qcd_CMUL(block[v3],C2));
	   }
       MPI_Reduce(&(gl.re), &(gl2.re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
       gLoops[lt*Nmom+imom] = gl2;
     } // close momentum loops
   } // do all time slices
   free(block);
}


void qcd_calculateGluonLoops_T00g_Sg(qcd_gaugeField *u, qcd_complex_16* gLoops_T00g, qcd_complex_16* gLoops_Sg, qcd_int_4 (*mom)[3], int Nmom)
{
  /*
    Note that the size of the gLoops is  NT * Nmom where Nmom is the number of momenta and NT is the temporal extent
    In this function we compute both Sg = (E^2+B^2) and T00g=(-E^2+B2) 
    paper:0707.3225
   */
   qcd_communicateGaugeM(u);
   qcd_waitall(u->geo);
   qcd_complex_16 plaq[3][3],tmp[3][3];
   qcd_complex_16 smu;
   qcd_complex_16 nine0;
   qcd_int_8 l;
   qcd_uint_4 lt;
   qcd_uint_2 mu,nu;
   qcd_uint_2 lx,ly,lz;
   qcd_uint_2 x,y,z;
   int v3;
   qcd_complex_16 *block_T00g;
   qcd_complex_16 *block_Sg;
   qcd_complex_16 gl_T00g,gl2_T00g;
   qcd_complex_16 gl_Sg,gl2_Sg;
   qcd_complex_16 C2;
   qcd_real_8 argPhase;

   nine0 = (qcd_complex_16) {9,0};

   block_T00g = (qcd_complex_16*) malloc(u->geo->lV3*sizeof(qcd_complex_16));
   block_Sg = (qcd_complex_16*) malloc(u->geo->lV3*sizeof(qcd_complex_16));
   if(block_T00g == NULL || block_Sg == NULL){
     fprintf(stderr,"process %i: Error allocating memory! Out of memory\n",u->geo->myid);
     exit(EXIT_FAILURE);
   }

   if(u->geo->L[0] != u->geo->lL[0]){
     fprintf(stderr,"Error! function which does the gluon loops accepts only one task in the temporal direction");
     exit(EXIT_FAILURE);
   }

   for(lt=0; lt < u->geo->lL[0]; lt++){

     for(v3=0; v3 < u->geo->lV3; v3++){ //set block_T00gs to zero for each time-slice
       block_T00g[v3]= (qcd_complex_16) {0,0};
       block_Sg[v3] = (qcd_complex_16) {0,0};
     }
     


     for(lx=0;lx<u->geo->lL[1];lx++)
       for(ly=0;ly<u->geo->lL[2];ly++)
	 for(lz=0;lz<u->geo->lL[3];lz++){
	   l=qcd_LEXIC(lt,lx,ly,lz,u->geo->lL);
	   v3 = qcd_LEXIC0(lx,ly,lz,u->geo->lL);

	   // term1
	   smu.re=0; smu.im=0;
	   for(mu=1; mu<4; mu++)
	     {
	       qcd_MUL3x3(plaq,u->D[l][mu],u->D[u->geo->plus[l][mu]][0]);
	       qcd_MULADJOINT3x3(tmp,plaq,u->D[u->geo->plus[l][0]][mu]);
	       qcd_MULADJOINT3x3(plaq,tmp,u->D[l][0]);
	       smu.re += qcd_SU3TRACER(plaq);
	     }
	   block_T00g[v3] = qcd_CADD(block_T00g[v3],smu);
	   block_Sg[v3] = qcd_CADD(block_Sg[v3],qcd_CSUB(nine0,smu));

	   //term2
	   smu.re=0; smu.im=0;
	   for(mu=1; mu<4; mu++)
	     for(nu=mu+1;nu<4;nu++)
	       {
		 qcd_MUL3x3(plaq,u->D[l][mu],u->D[u->geo->plus[l][mu]][nu]);
		 qcd_MULADJOINT3x3(tmp,plaq,u->D[u->geo->plus[l][nu]][mu]);
		 qcd_MULADJOINT3x3(plaq,tmp,u->D[l][nu]);
		 smu.re += qcd_SU3TRACER(plaq);
	       }
	   block_T00g[v3] = qcd_CSUB(block_T00g[v3],smu);	   
	   block_Sg[v3] = qcd_CADD(block_Sg[v3],qcd_CSUB(nine0,smu));
	 } // create the gluon loop for each lattice point

     // fourier transform
     for(int imom = 0 ; imom < Nmom ; imom++){
       gl_T00g = (qcd_complex_16) {0,0};
       gl_Sg = (qcd_complex_16) {0,0};
       for(lx=0;lx<u->geo->lL[1];lx++)
	 for(ly=0;ly<u->geo->lL[2];ly++)
	   for(lz=0;lz<u->geo->lL[3];lz++){
	     v3 = qcd_LEXIC0(lx,ly,lz,u->geo->lL);
	     x=lx+u->geo->Pos[1]*u->geo->lL[1];
	     y=ly+u->geo->Pos[2]*u->geo->lL[2];
	     z=lz+u->geo->Pos[3]*u->geo->lL[3];
	     argPhase = (((double) mom[imom][0]*x)/u->geo->L[1] + ((double) mom[imom][1]*y)/u->geo->L[2] + ((double) mom[imom][2]*z)/u->geo->L[3])*2*M_PI;
	     C2=(qcd_complex_16) {cos(argPhase), -sin(argPhase)}; // note the sign convention of the gluon loop FT is the same as the quark loops 
	     gl_T00g=qcd_CADD(gl_T00g, qcd_CMUL(block_T00g[v3],C2));
	     gl_Sg=qcd_CADD(gl_Sg, qcd_CMUL(block_Sg[v3],C2));
	   }
       MPI_Reduce(&(gl_T00g.re), &(gl2_T00g.re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
       MPI_Reduce(&(gl_Sg.re), &(gl2_Sg.re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
       gLoops_T00g[lt*Nmom+imom] = gl2_T00g;
       gLoops_Sg[lt*Nmom+imom] = gl2_Sg;
     } // close momentum loops
   } // do all time slices
   free(block_T00g);
   free(block_Sg);
}
