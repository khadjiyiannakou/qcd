/**** Kyriakos Hadjiyiannakou ***/ 
 /*                physical basis                                   |  twisted basis (changes only) 
   // --------------------------------------------------------------|----------------------------------- 
   op_G_0 = 0,   //    gamma5                                       |  i identity tau_3 
   op_G_1 = 1,   //    gamma1                                       | 
   op_G_2 = 2,   //    gamma2                                       | 
   op_G_3 = 3,   //    gamma3                                       | 
   op_G_4 = 4,   // -i gamma0 gamma5                                | 
   op_G_5 = 5,   // -i gamma0 gamma1                                |    gamma5 gamma0 gamma1 tau_3 
   op_G_6 = 6,   // -i gamma0 gamma2                                |    gamma5 gamma0 gamma2 tau_3 
   op_G_7 = 7,   // -i gamma0 gamma3                                |    gamma5 gamma0 gamma3 tau_3 
   op_G_8 = 8,   //    identity                                     |  i gamma5 tau_3 
   op_G_9 = 9,   // -i gamma_5 gamma1                               | 
   op_G_10 = 10, // -i gamma_5 gamma2                               | 
   op_G_11 = 11, // -i gamma_5 gamma3                               | 
   op_G_12 = 12, //    gamma_0                                      | 
   op_G_13 = 13, // -i gamma_5 gamma_0 gamma_1                      |    gamma_0 gamma_1 tau_3 
   op_G_14 = 14, // -i gamma_5 gamma_0 gamma_2                      |    gamma_0 gamma_2 tau_3 
   op_G_15 = 15  // -i gamma_5 gamma_0 gamma_3                      |    gamma_0 gamma_3 tau_3 
   op_G_16 = 16,   // -i gamma_5 ( gamma1 + gamma2 + gamma3)                             | 
 */ 
qcd_complex_16 qcd_GAMMA_OPERATORS_UP[17][4][4]=
  {
    {
      { {.re=0.0,.im=1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=1.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=1.0} }
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
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=-1.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=-1.0} },
      { {.re=0.0,.im=1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=-1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=-1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=-1.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=-1.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=-1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=-1.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=-1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=-1.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=1.0} }
    },

    {
      { {.re=0.0,.im=1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=-1.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=-1.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=-1.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=-1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=-1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=-1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=1.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=-1.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=-1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=-1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0} },
      { {.re=-1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=-1.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=-1.0,.im=0.0} },
      { {.re=-1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=-1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=-1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=-1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=1.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=1.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=-1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=1.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=-1.0,.im=0.0},{.re=0.0,.im=0.0} }
    },

    {
      { {.re=0.0,.im=-1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=1.0},{.re=0.0,.im=0.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=-1.0} }
    },

    {
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=-1.0,.im=0.0},{.re=-1.0,.im=1.0} },
      { {.re=0.0,.im=0.0},{.re=0.0,.im=0.0},{.re=-1.0,.im=-1.0},{.re=1.0,.im=0.0} },
      { {.re=-1.0,.im=0.0},{.re=-1.0,.im=1.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} },
      { {.re=-1.0,.im=-1.0},{.re=1.0,.im=0.0},{.re=0.0,.im=0.0},{.re=0.0,.im=0.0} }
    },
  };

