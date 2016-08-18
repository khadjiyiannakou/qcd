#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>

int qcd_communicateGaugeM_modified(qcd_gaugeField *u, MPI_Comm local_comm)
{
  qcd_uint_2 b,bb;
  qcd_uint_8 startpos;

  if(u->geo->numOfRequests != 0)
    {
      fprintf(stderr,"process %i: Error in qcd_communicateGaugeM! Previous communication not finished\n",u->geo->myid);
      return(1);
    }

  //start communication                                                                                                                             
  for(b=0; b<4; b++)
    {
      if(u->geo->lL[b] < u->geo->L[b])
	{
	  //start communication in b direction                                                                                                      
	  //send to - / recieve from +                                                                                                              
	 MPI_Isend(&(u->D[0][0][0][0]), 1, u->geo->stypeU[b], u->geo->Pminus[b], 3*b, local_comm, &(u->geo->requests[u->geo->numOfRequests++]));
	 MPI_Irecv(&(u->Bplus[b][0][0][0]), 1, u->geo->rtypeU[b], u->geo->Pplus[b], 3*b, local_comm, &(u->geo->requests[u->geo->numOfRequests++])); 
	}
    }//end communication-start loop                                                                                            
  return 0;
}//end qcd_communicateGaugeM  

qcd_real_8 qcd_calculatePlaquette_modified(qcd_gaugeField *u, MPI_Comm local_comm )
{
  qcd_communicateGaugeM_modified(u,local_comm);
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
  MPI_Allreduce(&meanplaq, &result, 1, MPI_DOUBLE, MPI_SUM, local_comm);
  return(result/(u->geo->V * 3.0 * 6.0)); // volume * N_color * number of mu-nu combinations                                                        
}//end qcd_calculatePlaquette                                                                                                      

int qcd_communicateGaugePM_modified(qcd_gaugeField *u, MPI_Comm local_comm)
{
  qcd_uint_2 b,bb;
  qcd_uint_8 startpos;

  if(u->geo->numOfRequests != 0)
    {
      fprintf(stderr,"process %i: Error in qcd_communicateGaugePM! Previous communication not finished\n",u->geo->myid);
      return(1);
    }

  //start communication                                                                                                                             
  for(b=0; b<4; b++)
    {
      if(u->geo->lL[b] < u->geo->L[b])
	{
	  //start communication in b direction                                                                                                        
	  //send to - / recieve from +                                                                                                              
	  MPI_Isend(&(u->D[0][0][0][0]), 1, u->geo->stypeU[b], u->geo->Pminus[b], 3*b, local_comm, &(u->geo->requests[u->geo->numOfRequests++]));
	  MPI_Irecv(&(u->Bplus[b][0][0][0]), 1, u->geo->rtypeU[b], u->geo->Pplus[b], 3*b, local_comm, &(u->geo->requests[u->geo->numOfRequests++]));
	  startpos=u->geo->lL[b]-1;
	  for(bb=0; bb<b; bb++)
            startpos *= u->geo->lL[bb];

	  //send to + / recieve from -                                                              
	  //send to + / recieve from -                                                                                
	  MPI_Isend(&(u->D[startpos][0][0][0]), 1, u->geo->stypeU[b], u->geo->Pplus[b], 3*b+2, local_comm, &(u->geo->requests[u->geo->numOfRequests++]));
	  MPI_Irecv(&(u->Bminus[b][0][0][0]), 1, u->geo->rtypeU[b], u->geo->Pminus[b], 3*b+2, local_comm, &(u->geo->requests[u->geo->numOfRequests++]));
	}
    }//end communication-start loop                                                                                                                 
  return 0;
}//end qcd_communicateGaugePM  

int qcd_communicatePropagatorP_modified(qcd_propagator *p, MPI_Comm local_comm)
{
   qcd_uint_2 b,bb;
   qcd_uint_8 startpos;
   
   if(p->geo->numOfRequests != 0)
   {
      fprintf(stderr,"process %i: Error in qcd_communicatePropagatorP! Previous communication not finished\n",p->geo->myid);
      return(1);
   }
   
   //start communication   
   for(b=0; b<4; b++)
   {
      if(p->geo->lL[b] < p->geo->L[b])
      {
      
         startpos=p->geo->lL[b]-1;
         for(bb=0; bb<b; bb++)
            startpos *= p->geo->lL[bb];
            
         //send to + / recieve from -
         MPI_Isend(&(p->D[startpos][0][0][0][0]), 1, p->geo->stypeP[b], p->geo->Pplus[b], 3*b+1, local_comm, &(p->geo->requests[p->geo->numOfRequests++]));
         MPI_Irecv(&(p->Bminus[b][0][0][0][0]), 1, p->geo->rtypeP[b], p->geo->Pminus[b], 3*b+1, local_comm, &(p->geo->requests[p->geo->numOfRequests++]));
      }
      
   }//end communication-start loop  
      
   return 0;
}//end qcd_communicatePropagatorP

int qcd_communicateVectorPMGaugeP_modified(qcd_vector *v, qcd_gaugeField *u,MPI_Comm local_comm)
{
   qcd_uint_2 b,bb;
   qcd_uint_8 startpos;
   
   if(v->geo->numOfRequests != 0)
   {
      fprintf(stderr,"process %i: Error in qcd_communicateVectorPMGaugeP! Previous communication not finished\n",v->geo->myid);
      return(1);
   }
   
   //start communication   
   for(b=0; b<4; b++)
   {
      if(v->geo->lL[b] < v->geo->L[b])
      {
         //start communication in b direction
         //send to - / recieve from +
         MPI_Isend(&(v->D[0][0][0]), 1, v->geo->stypeV[b], v->geo->Pminus[b], 3*b, local_comm, &(v->geo->requests[v->geo->numOfRequests++]));
         MPI_Irecv(&(v->Bplus[b][0][0]), 1, v->geo->rtypeV[b], v->geo->Pplus[b], 3*b, local_comm, &(v->geo->requests[v->geo->numOfRequests++]));
         
         startpos=v->geo->lL[b]-1;
         for(bb=0; bb<b; bb++)
            startpos *= v->geo->lL[bb];
            
         //send to + / recieve from -
         MPI_Isend(&(v->D[startpos][0][0]), 1, v->geo->stypeV[b], v->geo->Pplus[b], 3*b+1, local_comm, &(v->geo->requests[v->geo->numOfRequests++]));
         MPI_Irecv(&(v->Bminus[b][0][0]), 1, v->geo->rtypeV[b], v->geo->Pminus[b], 3*b+1, local_comm, &(v->geo->requests[v->geo->numOfRequests++]));
         MPI_Isend(&(u->D[startpos][0][0][0]), 1, v->geo->stypeU[b], v->geo->Pplus[b], 3*b+2, local_comm, &(v->geo->requests[v->geo->numOfRequests++]));
         MPI_Irecv(&(u->Bminus[b][0][0][0]), 1, v->geo->rtypeU[b], v->geo->Pminus[b], 3*b+2, local_comm, &(v->geo->requests[v->geo->numOfRequests++]));
      }
      
   }//end communication-start loop  
      
   return 0;
}//end qcd_communicateVectorPMGaugeP

/* perform 1 iteration of 3d APE-smearing 
 * with parameter alpha.
 *
 * u -> SU3-projection( u + alpha * sum spatial staples)
 */
int qcd_apeSmear3d_modified(qcd_gaugeField *apeu, qcd_gaugeField *u, qcd_real_8 alpha, MPI_Comm local_comm)
{
   qcd_propagator edge;
   qcd_complex_16 stapleForward[3][3];
   qcd_complex_16 stapleBackward[3][3];
   qcd_uint_2 mu,nu,i=0,c1,c2;
   qcd_uint_4 l;
   qcd_complex_16 tmp[3][3];

   qcd_initPropagator(&edge,u->geo); // store edges in a propagator-structure. 

   /* since staples need next-to-nearest neighbors like U(x+mu-nu), this is done in 2 steps
      a) communicate U & calculate edges U_mu(x)U_nu(x+mu)
      b) communicate edges and put them together to staples.
   */
   
   qcd_communicateGaugePM_modified(u,local_comm);
   qcd_zeroGaugeField(apeu);   
   qcd_waitall(u->geo);

   for(mu=1; mu<4; mu++)
   for(nu=1;nu<4; nu++)
   if(mu!=nu)   
   for(l=0;l<u->geo->lV;l++)
   {
      qcd_MUL3x3(edge.D[l][mu][nu], u->D[l][mu], u->D[u->geo->plus[l][mu]][nu]);
   }
   
   qcd_communicatePropagatorP_modified(&edge,local_comm);

   //the forward staple doesn't need the edges
   for(mu=1; mu<4; mu++)
   for(nu=1;nu<4; nu++)
   if(mu!=nu)
   for(l=0;l<u->geo->lV;l++)
   {
      qcd_MUL3x3(tmp, u->D[l][nu],u->D[u->geo->plus[l][nu]][mu]);
      qcd_MULADJOINT3x3(stapleForward, tmp, u->D[u->geo->plus[l][mu]][nu]);
      for(c1=0;c1<3;c1++)
      for(c2=0;c2<3;c2++)
         apeu->D[l][mu][c1][c2] = qcd_CADD(apeu->D[l][mu][c1][c2],stapleForward[c1][c2]);
   }   
   
   qcd_waitall(u->geo);
   
   for(mu=1; mu<4; mu++)
   for(nu=1;nu<4; nu++)
   if(mu!=nu)
   for(l=0;l<u->geo->lV;l++)
   {
      qcd_ADJOINTMUL3x3(stapleBackward, u->D[u->geo->minus[l][nu]][nu], edge.D[u->geo->minus[l][nu]][mu][nu]);
      for(c1=0;c1<3;c1++)
      for(c2=0;c2<3;c2++)
         apeu->D[l][mu][c1][c2] = qcd_CADD(apeu->D[l][mu][c1][c2],stapleBackward[c1][c2]);
   }
   
   
   qcd_scaleGaugeField(apeu,alpha);
   qcd_addGaugeField(apeu,u,apeu);
   
   qcd_projectSU33d(apeu);
   
   qcd_destroyPropagator(&edge);
   return(0);
}//end qcd_apeSmear3d

int qcd_gaussIteration3d_modified(qcd_vector *v, qcd_gaugeField *u, qcd_real_8 alpha, qcd_uint_4 t, MPI_Comm local_comm)
{
   qcd_uint_8 i,j;
   qcd_uint_2 c1,mu,nu,b; 
   qcd_uint_2 bl1[4] = {2,2,1,1};
   qcd_uint_2 bl2[4] = {3,3,3,2};
   qcd_uint_4 x,y,z,b0,b1,b2,b3,tt=0;
   qcd_complex_16 tmp[3];
   qcd_real_8 nrm = 1.0/(1.0 + alpha*6.0);
   qcd_vector v2;
   qcd_real_8 *uu;
   qcd_real_8 *psi;
   qcd_real_8 upsi[24];
   qcd_real_8 udaggerpsi[24];
   qcd_real_8 *total;

   if(!v->initialized)
   {
      fprintf(stderr,"Error in qcd_gaussIteration3d! Vector mus be properly initialized\n");
      return(1);
   }
   if(qcd_initVector(&v2, v->geo))
   {
     fprintf(stderr,"process %i: Error in qcd_gaussIteration3d! Could not initialize vector", v->geo->myid);
     return(1);
   }   

   //start communication (4d comm is in principle too much, but hey...
   qcd_communicateVectorPMGaugeP_modified(v,u,local_comm);

   //smear inner points:   
   qcd_zeroVector(&v2);

   if(v->geo->lL[1]>2 && v->geo->lL[2]>2 && v->geo->lL[3]>2 
      && t>=v->geo->Pos[0]*v->geo->lL[0] && t<(v->geo->Pos[0]+1)*v->geo->lL[0])
   { 
      tt=t-v->geo->Pos[0]*v->geo->lL[0];
      for(z=1;z<v->geo->lL[3]-1;z++)
      for(y=1;y<v->geo->lL[2]-1;y++)
      for(x=1;x<v->geo->lL[1]-1;x++)
      {
         i=qcd_LEXIC(tt,x,y,z,v->geo->lL);
         for(nu=1;nu<4;nu++)
         {
             uu  = (qcd_real_8*) &(u->D[i][nu][0][0].re);
             psi = (qcd_real_8*) &(v->D[v->geo->plus[i][nu]][0][0].re);
             qcd_APPLY_U(uu,upsi,psi);

             uu  = (qcd_real_8*) &(u->D[v->geo->minus[i][nu]][nu][0][0].re);
             psi = (qcd_real_8*) &(v->D[v->geo->minus[i][nu]][0][0].re); 
             qcd_APPLY_U_DAGGER(uu,udaggerpsi,psi);

             total = (qcd_real_8*) &(v2.D[i][0][0].re);
             qcd_SUM_UP_HOPP(total,upsi,udaggerpsi);
         }
      }//end inner-point loop
   }//end inner-points condition

   qcd_waitall(v->geo);

   //now boundary points
   if(t>=v->geo->Pos[0]*v->geo->lL[0] && t<(v->geo->Pos[0]+1)*v->geo->lL[0])
   for(j=0; j<v->geo->edge0Points; j++)
   {
      i=v->geo->edge0[j]*v->geo->lL[0]+tt; // works only with present lexic
      for(nu=1;nu<4;nu++)
      {
         uu  = (qcd_real_8*) &(u->D[i][nu][0][0].re);
         psi = (qcd_real_8*) &(v->D[v->geo->plus[i][nu]][0][0].re);
         qcd_APPLY_U(uu,upsi,psi);

         uu  = (qcd_real_8*) &(u->D[v->geo->minus[i][nu]][nu][0][0].re);
         psi = (qcd_real_8*) &(v->D[v->geo->minus[i][nu]][0][0].re); 
         qcd_APPLY_U_DAGGER(uu,udaggerpsi,psi);

         total = (qcd_real_8*) &(v2.D[i][0][0].re);
         qcd_SUM_UP_HOPP(total,upsi,udaggerpsi);
      }           
   }//end boundaries loop                  
   
   qcd_scaleVector3d(&v2,alpha,t);
   qcd_addVector(v,v,&v2);
   qcd_scaleVector3d(v,nrm,t); 

   qcd_destroyVector(&v2);
   return 0;
} 
