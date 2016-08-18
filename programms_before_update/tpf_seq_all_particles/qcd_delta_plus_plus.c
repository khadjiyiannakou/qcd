/* KYRIAKOS HADJIYIANNAKOU 
 delta++ 3pf using sequential throw the current */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>

int tpf_delta_plus_plus(qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator *sprop, qcd_propagator *seq_uprop, qcd_propagator *seq_dprop,
	       qcd_propagator *seq_sprop, qcd_geometry *geo, qcd_int_2 cur_time, qcd_int_2 number_tsinks, qcd_complex_16 *U_part_deltaPP[4][4], 
	       qcd_complex_16 *D_part_deltaPP[4][4], qcd_complex_16 *S_part_deltaPP[4][4],qcd_int_4(*mom)[3], int num_momenta, qcd_uint_4 *x_src, qcd_int_2 after_tcur){
     
     qcd_complex_16 *block_Upart[4][4], *block_Dpart[4][4], *block_Spart[4][4];
     int malloc_flag, malloc_flag_check;
     qcd_uint_4 nz_counter[3];
     qcd_complex_16 C, expon;
     qcd_int_2 cg5cg5b_ind[3][16*16][4];
     qcd_complex_16 cg5cg5b_val[3][16*16];
     qcd_uint_4 t,x,y,z;
     qcd_real_8 tmp;
     qcd_uint_2 alpha, beta, alpha1, beta1;
     qcd_uint_2 a, b, c, a1, b1, c1;
     qcd_uint_4 v3,v;                                 // variables for v3 = lx*ly*lz , v = lx*ly*lz*lt
     qcd_complex_16 tmp_corr[3][4][4] ;
     qcd_uint_4 itmom;
     
     for(int mu=0; mu < 4 ; mu++)
       for(int nu=0; nu < 4; nu++){
	 block_Upart[mu][nu] = (qcd_complex_16*) malloc(geo->lV3*sizeof(qcd_complex_16));  // allocate local volume memory for each part   
	 block_Dpart[mu][nu] = (qcd_complex_16*) malloc(geo->lV3*sizeof(qcd_complex_16));
	 block_Spart[mu][nu] = (qcd_complex_16*) malloc(geo->lV3*sizeof(qcd_complex_16));
       }
	 
	 
     malloc_flag = 0;
     
     for(int mu=0; mu < 4 ; mu++)
	for(int nu=0; nu < 4; nu++){
          if(block_Upart[mu][nu] == NULL)malloc_flag++;  // check if the allocation is ok or if we are out of memory
          if(block_Dpart[mu][nu] == NULL)malloc_flag++;
          if(block_Spart[mu][nu] == NULL)malloc_flag++;
	}

     MPI_Allreduce(&malloc_flag, &malloc_flag_check, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
     
    if(malloc_flag_check>0) // check for memory allocation flag
    {
      if(geo->myid == 0) printf("not enough memory\n");
      exit(EXIT_FAILURE);
    }
    
    /////////// tabulate non zero elements ///////////////////
    
   for(int idec=0;idec<3;idec++) nz_counter[idec] = 0;
   
    for(int idec = 0 ; idec < 3 ; idec++)
      for(int alpha=0 ; alpha < 4; alpha++)
	for(int beta=0; beta < 4 ; beta++)
	  for(int beta1=0; beta1 < 4; beta1++)
	    for(int alpha1=0; alpha1 < 4 ; alpha1++)
	    {	
	      C = qcd_CMUL(qcd_CGAMMA[idec+1][alpha][beta],qcd_BAR_CGAMMA[idec+1][beta1][alpha1]);
	      if(qcd_NORM(C)>1e-3)
	      {
		cg5cg5b_val[idec][nz_counter[idec]].re = C.re;
		cg5cg5b_val[idec][nz_counter[idec]].im = C.im;
		cg5cg5b_ind[idec][nz_counter[idec]][0] = alpha;
		cg5cg5b_ind[idec][nz_counter[idec]][1] = beta;
		cg5cg5b_ind[idec][nz_counter[idec]][2] = beta1;
		cg5cg5b_ind[idec][nz_counter[idec]][3] = alpha1;                                                            
		nz_counter[idec]++;
	      }
	    }
    //////////////////////////start contractions ////////////////////////
    
    for(int it = 0; it < number_tsinks ; it++){
	  t = ((it + cur_time + after_tcur)%geo->L[0]);
	  if(geo->myid == 0)printf("Time step %d\n",it);

	  for(int mu = 0; mu < 4 ; mu++)
	  for(int nu = 0;nu < 4 ; nu++)
	  for(int v3=0; v3 < geo->lV3; v3++){
            block_Upart[mu][nu][v3]= (qcd_complex_16) {0,0};
	    block_Dpart[mu][nu][v3]= (qcd_complex_16) {0,0};
            block_Spart[mu][nu][v3]= (qcd_complex_16) {0,0};
	  }
	  
	  for(int idec =0 ; idec < 3 ; idec++)
	  for(int gamma=0 ; gamma < 4 ; gamma++)
	  for(int gamma1=0; gamma1 < 4 ; gamma1++)
	  for(int inz =0 ; inz < nz_counter[idec] ; inz++){
	    
	    alpha = cg5cg5b_ind[idec][inz][0];
	    beta = cg5cg5b_ind[idec][inz][1];
	    beta1 = cg5cg5b_ind[idec][inz][2];
	    alpha1 =cg5cg5b_ind[idec][inz][3];
	    
	    for(int cc1=0;cc1<6;cc1++)
               {
                  a = qcd_EPS[cc1][0];
                  b = qcd_EPS[cc1][1];
                  c = qcd_EPS[cc1][2];
                  for(int cc2=0;cc2<6;cc2++)
                  {          
                     a1 = qcd_EPS[cc2][0];
                     b1 = qcd_EPS[cc2][1];
                     c1 = qcd_EPS[cc2][2];
		     
                     for(int lx=0; lx < geo->lL[1]; lx++)
                     for(int ly=0; ly < geo->lL[2]; ly++)
                     for(int lz=0; lz < geo->lL[3]; lz++)
                     {
                        v3 = qcd_LEXIC0(lx,ly,lz,geo->lL);
                        v =  qcd_LEXIC(t,lx,ly,lz,geo->lL);
		     
			////////// Upart ////////////////
			block_Upart[gamma][gamma1][v3] = qcd_CSUB( block_Upart[gamma][gamma1][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
			                                 cg5cg5b_val[idec][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][beta][beta1][b][b1] )
			                                 ,uprop->D[v][gamma][alpha1][c][a1] ), seq_uprop->D[v][alpha][gamma1][a][c1] ));
			
			block_Upart[gamma][gamma1][v3] = qcd_CADD( block_Upart[gamma][gamma1][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
			                                 cg5cg5b_val[idec][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][beta][gamma1][b][c1] )
			                                 ,uprop->D[v][gamma][alpha1][c][a1] ), seq_uprop->D[v][alpha][beta1][a][b1] ));
			
			block_Upart[gamma][gamma1][v3] = qcd_CADD( block_Upart[gamma][gamma1][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
			                                 cg5cg5b_val[idec][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][beta][alpha1][b][a1] )
			                                 ,uprop->D[v][gamma][beta1][c][b1] ), seq_uprop->D[v][alpha][gamma1][a][c1] ));
			
			block_Upart[gamma][gamma1][v3] = qcd_CSUB( block_Upart[gamma][gamma1][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
			                                 cg5cg5b_val[idec][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][beta][gamma1][b][c1] )
			                                 ,uprop->D[v][gamma][beta1][c][b1] ), seq_uprop->D[v][alpha][alpha1][a][a1] ));
			
			block_Upart[gamma][gamma1][v3] = qcd_CSUB( block_Upart[gamma][gamma1][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
			                                 cg5cg5b_val[idec][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][beta][alpha1][b][a1] )
			                                 ,uprop->D[v][gamma][gamma1][c][c1] ), seq_uprop->D[v][alpha][beta1][a][b1] ));
			
			block_Upart[gamma][gamma1][v3] = qcd_CADD( block_Upart[gamma][gamma1][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
			                                 cg5cg5b_val[idec][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][beta][beta1][b][b1] )
			                                 ,uprop->D[v][gamma][gamma1][c][c1] ), seq_uprop->D[v][alpha][alpha1][a][a1] ));
			
			block_Upart[gamma][gamma1][v3] = qcd_CADD( block_Upart[gamma][gamma1][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
			                                 cg5cg5b_val[idec][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][beta1][a][b1] )
			                                 ,uprop->D[v][gamma][alpha1][c][a1] ), seq_uprop->D[v][beta][gamma1][b][c1] ));
			
			block_Upart[gamma][gamma1][v3] = qcd_CSUB( block_Upart[gamma][gamma1][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
			                                 cg5cg5b_val[idec][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][gamma1][a][c1] )
			                                 ,uprop->D[v][gamma][alpha1][c][a1] ), seq_uprop->D[v][beta][beta1][b][b1] ));
			
			block_Upart[gamma][gamma1][v3] = qcd_CSUB( block_Upart[gamma][gamma1][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
			                                 cg5cg5b_val[idec][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][alpha1][a][a1] )
			                                 ,uprop->D[v][gamma][beta1][c][b1] ), seq_uprop->D[v][beta][gamma1][b][c1] ));
			
			block_Upart[gamma][gamma1][v3] = qcd_CADD( block_Upart[gamma][gamma1][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
			                                 cg5cg5b_val[idec][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][gamma1][a][c1] )
			                                 ,uprop->D[v][gamma][beta1][c][b1] ), seq_uprop->D[v][beta][alpha1][b][a1] ));
			
			block_Upart[gamma][gamma1][v3] = qcd_CADD( block_Upart[gamma][gamma1][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
			                                 cg5cg5b_val[idec][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][alpha1][a][a1] )
			                                 ,uprop->D[v][gamma][gamma1][c][c1] ), seq_uprop->D[v][beta][beta1][b][b1] ));
			
			block_Upart[gamma][gamma1][v3] = qcd_CSUB( block_Upart[gamma][gamma1][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
			                                 cg5cg5b_val[idec][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][beta1][a][b1] )
			                                 ,uprop->D[v][gamma][gamma1][c][c1] ), seq_uprop->D[v][beta][alpha1][b][a1] ));

			block_Upart[gamma][gamma1][v3] = qcd_CSUB( block_Upart[gamma][gamma1][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
			                                 cg5cg5b_val[idec][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][beta1][a][b1] )
			                                 ,uprop->D[v][beta][alpha1][b][a1] ), seq_uprop->D[v][gamma][gamma1][c][c1] ));
			
			block_Upart[gamma][gamma1][v3] = qcd_CADD( block_Upart[gamma][gamma1][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
			                                 cg5cg5b_val[idec][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][gamma1][a][c1] )
			                                 ,uprop->D[v][beta][alpha1][b][a1] ), seq_uprop->D[v][gamma][beta1][c][b1] ));
			
			block_Upart[gamma][gamma1][v3] = qcd_CADD( block_Upart[gamma][gamma1][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
			                                 cg5cg5b_val[idec][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][alpha1][a][a1] )
			                                 ,uprop->D[v][beta][beta1][b][b1] ), seq_uprop->D[v][gamma][gamma1][c][c1] ));
			
			block_Upart[gamma][gamma1][v3] = qcd_CSUB( block_Upart[gamma][gamma1][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
			                                 cg5cg5b_val[idec][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][gamma1][a][c1] )
			                                 ,uprop->D[v][beta][beta1][b][b1] ), seq_uprop->D[v][gamma][alpha1][c][a1] ));
			
			block_Upart[gamma][gamma1][v3] = qcd_CSUB( block_Upart[gamma][gamma1][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
			                                 cg5cg5b_val[idec][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][alpha1][a][a1] )
			                                 ,uprop->D[v][beta][gamma1][b][c1] ), seq_uprop->D[v][gamma][beta1][c][b1] ));
			
			block_Upart[gamma][gamma1][v3] = qcd_CADD( block_Upart[gamma][gamma1][v3] ,qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( 
			                                 cg5cg5b_val[idec][inz], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),uprop->D[v][alpha][beta1][a][b1] )
			                                 ,uprop->D[v][beta][gamma1][b][c1] ), seq_uprop->D[v][gamma][alpha1][c][a1] ));
			
		     } // space
		  } // color
	       }  // color
	  } // close non vanishing spin elements
	  
	  for(int imom = 0; imom < num_momenta; imom++){
	   
	    for(int part = 0 ; part < 3 ; part++)
	     for(int mu = 0 ; mu < 4 ; mu++)
	       for(int nu = 0 ; nu < 4 ; nu++)
		 tmp_corr[part][mu][nu] = (qcd_complex_16) {0,0} ;
	    
	    for(int lx=0; lx < geo->lL[1]; lx++)
            for(int ly=0; ly < geo->lL[2]; ly++)
            for(int lz=0; lz < geo->lL[3]; lz++)
            {
               v3 = qcd_LEXIC0(lx,ly,lz,geo->lL);
               x=lx+geo->Pos[1]*geo->lL[1] - x_src[1];
               y=ly+geo->Pos[2]*geo->lL[2] - x_src[2];
               z=lz+geo->Pos[3]*geo->lL[3] - x_src[3];
               tmp = (((double) mom[imom][0]*x)/geo->L[1] + ((double) mom[imom][1]*y)/geo->L[2] + ((double) mom[imom][2]*z)/geo->L[3])*2*M_PI;
               expon=(qcd_complex_16) {cos(tmp), -sin(tmp)}; //TABULATE FOR LARGE SPEEDUP!!!
//               corr=qcd_CADD(corr, qcd_CMUL(block[v3],C2));
	       for(int gamma =0 ; gamma < 4 ; gamma++)
		 for(int gamma1 =0 ;gamma1 < 4 ; gamma1++){
		   tmp_corr[0][gamma][gamma1] = qcd_CADD(tmp_corr[0][gamma][gamma1],qcd_CMUL(block_Upart[gamma][gamma1][v3],expon));
		     // no Dpart here
		    // no Spart here
		 } // close spin indices
            } // close space sum
            
	      itmom = it * num_momenta + imom ;
            for(int gamma =0 ; gamma < 4 ; gamma++)
	      for(int gamma1 =0 ;gamma1 < 4 ; gamma1++){
	      MPI_Allreduce(&(tmp_corr[0][gamma][gamma1].re), &(U_part_deltaPP[gamma][gamma1][itmom].re), 2, MPI_DOUBLE, MPI_SUM,  MPI_COMM_WORLD );
		// no Dpart here
		// no Spart here
	      } // close MPI_Reduce
	      
	  } // close momenta
	  
	  
    } //close time
     
 return 0; 
}
