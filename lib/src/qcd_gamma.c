/* qcd_gamma.h
 * collection of gamma matrices
 * ETMC conventions
 *
 * Tomasz Korzec 2009
 **************************************/
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>

const qcd_uint_2 qcd_EPS[6][3]=
  {
    {0,1,2},
    {2,0,1},
    {1,2,0},
    {2,1,0},
    {0,2,1},
    {1,0,2}
  };

const qcd_int_2 qcd_SGN_EPS[6]=
  {
    +1,+1,+1,-1,-1,-1
  };

const qcd_complex_16 qcd_ONE[4][4]=
     {
      { {.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0} }
    };

/*
 *
 *    qcd_GAMMA[0,1,2,3,4,5,6,7] -> \gamma_{0,1,2,3,4,5,+,-}
 *
 */
/*
 *
 *    gamma_4=-gamma_0
 *
 */
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

const qcd_complex_16 qcd_CGAMMA[6][4][4]=
  {
    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=-1.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=-1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=-1.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=1.0} },
      { {.re=0.0,.im=-1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0} },
      { {.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=1.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=1.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=-1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=-1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=-1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=-1.0,.im=0.0},{.re=0.0,.im=0.0} }
    }
  };

const qcd_complex_16 qcd_BAR_CGAMMA[6][4][4]=
  {
    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=-1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=-1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=1.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=-1.0} },
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
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=-1.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=-1.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=-1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=-1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=-1.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=-1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=-1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=-1.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0} }
    }
  };

const qcd_complex_16 qcd_G5GAMMA[8][4][4]=
  {
    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=-1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=-1.0,.im=0.0} },
      { {.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=-1.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=-1.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=-1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=-1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=-1.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=-1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=-1.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=1.0} },
      { {.re=0.0,.im=-1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0} },
      { {.re=-1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=-1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=-1.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=-1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=-1.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=-1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    }
  };

const qcd_complex_16 qcd_GAMMAG5[8][4][4]=
  {
    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0} },
      { {.re=-1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=-1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=1.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=1.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=-1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=-1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=1.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=-1.0} },
      { {.re=0.0,.im=1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=-1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=-1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=-1.0,.im=0.0} },
      { {.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=1.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=1.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    }
  };

const qcd_complex_16 qcd_G5GAMMAG5[8][4][4]=
  {
    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0} },
      { {.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=1.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=1.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=-1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=-1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=-1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=-1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=1.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=-1.0} },
      { {.re=0.0,.im=-1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=-1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=-1.0,.im=0.0} },
      { {.re=-1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=-1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=1.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=-1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=1.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=-1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    }
  };

const qcd_complex_16 qcd_BAR_G5GAMMA[8][4][4]=
  {
    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=-1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=-1.0,.im=0.0} },
      { {.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=1.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=1.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=-1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=-1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=1.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=-1.0} },
      { {.re=0.0,.im=1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=-1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0} },
      { {.re=-1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=-1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=1.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=1.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    }
  };

const qcd_complex_16 qcd_ONE_PLUS_GAMMA[6][4][4]=
  {
    {
      { {.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=-1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=-1.0,.im=0.0} },
      { {.re=-1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=-1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0} }
    },

    {
      { {.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=-1.0} },
      { {.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=-1.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=1.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0} }
    },

    {
      { {.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=-1.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=-1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0} }
    },

    {
      { {.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=-1.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=1.0} },
      { {.re=0.0,.im=1.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=-1.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0} }
    },

    {
      { {.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0} },
      { {.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0} }
    },

    {
      { {.re=2.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=2.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    }
  };

const qcd_complex_16 qcd_ONE_MINUS_GAMMA[6][4][4]=
  {
    {
      { {.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0} },
      { {.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0} }
    },

    {
      { {.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=1.0} },
      { {.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=1.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=-1.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=-1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0} }
    },

    {
      { {.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=-1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=-1.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0} }
    },

    {
      { {.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=1.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=-1.0} },
      { {.re=0.0,.im=-1.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=1.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0} }
    },

    {
      { {.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=-1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=-1.0,.im=0.0} },
      { {.re=-1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=-1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=2.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=2.0,.im=0.0} }
    }
  };

void qcd_gamma5Vector(qcd_vector *v)
{
   qcd_uint_8 l;
   qcd_uint_2 c;
   for(l=0; l<v->geo->lV; l++)
   for(c=0; c<3; c++)
   {
      v->D[l][2][c].re = -v->D[l][2][c].re;
      v->D[l][2][c].im = -v->D[l][2][c].im;
      v->D[l][3][c].re = -v->D[l][3][c].re;
      v->D[l][3][c].im = -v->D[l][3][c].im;
   }
}//end qcd_gamma5Vector 


void qcd_gamma5Propagator(qcd_propagator *p)
{
   qcd_uint_8 l;
   qcd_uint_2 c1,c2,mu;
   for(l=0; l<p->geo->lV; l++)
   for(mu=0; mu<4; mu++)
   for(c1=0; c1<3; c1++)
   for(c2=0; c2<3; c2++)
   {
      p->D[l][2][mu][c1][c2].re = -p->D[l][2][mu][c1][c2].re;
      p->D[l][2][mu][c1][c2].im = -p->D[l][2][mu][c1][c2].im;
      p->D[l][3][mu][c1][c2].re = -p->D[l][3][mu][c1][c2].re;
      p->D[l][3][mu][c1][c2].im = -p->D[l][3][mu][c1][c2].im;
   }
}//end qcd_gamma5Vector 

