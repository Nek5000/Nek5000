c Maximum number of real particles on a processor
#ifndef LPM_LPART
#define LPM_LPART 100000
#endif

c Number of secondary real properties for a particle
#define LPM_LRP2 4

c Number of integer properties for a particle
#ifndef LPM_LIP
#define LPM_LIP 11
#endif

c Maximum number of ghost particles
#define LPM_LPART_GP 26*LPM_LPART

#ifndef LPM_LRP_PRO
#define LPM_LRP_PRO 1
#endif
#if LPM_LRP_PRO == 0
#undef LPM_LRP_PRO
#define LPM_LRP_PRO 1
#endif

#ifndef LPM_LRP_GP
#define LPM_LRP_GP 4+LPM_LRP_PRO
#endif

#include "lpm_comm.f"
#include "lpm_solve.f"
#include "lpm_io.f"
