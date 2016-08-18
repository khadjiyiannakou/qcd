#include <stdio.h>

typedef struct{
  double re;
  double im;
}qcd_complex_16;

#define qcd_CMUL(x,y)  ( (qcd_complex_16) {x.re * y.re - x.im * y.im, x.re * y.im + x.im * y.re } )
#define qcd_CADD(x,y)  ( (qcd_complex_16) {x.re+y.re, x.im+y.im})
#define qcd_CSUB(x,y)  ( (qcd_complex_16) {x.re-y.re, x.im-y.im})
#define qcd_CSCALE(x,a)  ( (qcd_complex_16) {x.re*(a), x.im*(a)})

const qcd_complex_16 qcd_GAMMA[8][4][4]=
  {
    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=-1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=-1.0,.im=0.0} },
      { {.re=-1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=-1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=-1.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=-1.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=-1.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=-1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=-1.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=1.0} },
      { {.re=0.0,.im=1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=-1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0} },
      { {.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=-1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=-1.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=-1.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=-1.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    }
  };




static void projector_32(qcd_complex_16 Proj[3][3][4][4]){
  qcd_complex_16 Proj_temp[3][3][4][4];

  for(int i = 0 ; i < 3 ; i++)
    for(int j = 0 ; j < 3 ; j++)
      for(int mu = 0 ; mu < 4 ; mu++)
	for(int nu = 0 ; nu < 4 ; nu++){
	  Proj_temp[i][j][mu][nu] = (qcd_complex_16) {0.,0.};
	  for(int lu = 0 ; lu < 4 ; lu++)
	    Proj_temp[i][j][mu][nu] = qcd_CADD(Proj_temp[i][j][mu][nu],qcd_CMUL(qcd_GAMMA[i+1][mu][lu],qcd_GAMMA[j+1][lu][nu]));
	}
  qcd_complex_16 delta;
  for(int i = 0 ; i < 3 ; i++)
    for(int j = 0 ; j < 3 ; j++)
      for(int mu = 0 ; mu < 4 ; mu++)
	for(int nu = 0 ; nu < 4 ; nu++){
	  if( (i==j) && (mu==nu) )
	    delta = (qcd_complex_16) {1.,0.};
	  else
	    delta = (qcd_complex_16) {0.,0.};
	  Proj[i][j][mu][nu] = qcd_CSUB(delta,qcd_CSCALE(Proj_temp[i][j][mu][nu],1./3.));
	}
}

int main(){
  qcd_complex_16 P[3][3][4][4];
  projector_32(P);
  for(int i = 0 ; i < 3 ; i++)
    for(int j = 0 ; j < 3 ; j++)
      for(int mu = 0 ; mu <4 ; mu++)
	for(int nu = 0 ; nu < 4 ; nu++)
	  printf("%d %d %d %d \t %+f %+f\n",i,j,mu,nu,P[i][j][mu][nu].re,P[i][j][mu][nu].im);
  return 0;
}
