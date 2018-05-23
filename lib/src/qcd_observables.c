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

void qcd_calculateGluonLoops(qcd_gaugeField *u, qcd_real_8* gLoops)
{
   qcd_communicateGaugeM(u);
   qcd_waitall(u->geo);
   qcd_complex_16 plaq[3][3],tmp[3][3];
   qcd_real_8 meanplaq;
   qcd_real_8 res1,res2;
   qcd_int_8 l;
   qcd_uint_4 lt;
   qcd_uint_2 mu,nu;
   qcd_uint_2 x,y,z;
   for(lt=0; lt < u->geo->lL[0]; lt++){
     res1=0;res2=0;

     //first term
     meanplaq=0;
     for(x=0;x<u->geo->lL[1];x++)
       for(y=0;y<u->geo->lL[2];y++)
	 for(z=0;z<u->geo->lL[3];z++){
	   l=qcd_LEXIC(lt,x,y,z,u->geo->lL);
	     for(mu=1; mu<4; mu++)
	       {
		 qcd_MUL3x3(plaq,u->D[l][mu],u->D[u->geo->plus[l][mu]][0]);
		 qcd_MULADJOINT3x3(tmp,plaq,u->D[u->geo->plus[l][0]][mu]);
		 qcd_MULADJOINT3x3(plaq,tmp,u->D[l][0]);
		 meanplaq += qcd_SU3TRACER(plaq);
	       }
	   MPI_Allreduce(&meanplaq, &res1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	 }

     //second term
     meanplaq=0;
     for(x=0;x<u->geo->lL[1];x++)
       for(y=0;y<u->geo->lL[2];y++)
	 for(z=0;z<u->geo->lL[3];z++){
	   l=qcd_LEXIC(lt,x,y,z,u->geo->lL);
	     for(mu=1; mu<4; mu++)
	       for(nu=mu+1;nu<4;nu++)
		 {
		   qcd_MUL3x3(plaq,u->D[l][mu],u->D[u->geo->plus[l][mu]][nu]);
		   qcd_MULADJOINT3x3(tmp,plaq,u->D[u->geo->plus[l][nu]][mu]);
		   qcd_MULADJOINT3x3(plaq,tmp,u->D[l][nu]);
		   meanplaq += qcd_SU3TRACER(plaq);
		 }
	   MPI_Allreduce(&meanplaq, &res2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	 }
     gLoops[lt]=res1-res2;
   }
}
