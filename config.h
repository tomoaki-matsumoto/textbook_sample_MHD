! Header file for configuration
#include "constant.h"

! flux and EOS
#define _FLUX_SCHEME_ _MHD_


! #define _DIRECTIONAL_SPLIT_ _UNSPLIT_
#define _DIRECTIONAL_SPLIT_ _SPLIT_


! #define _TIME_MARCHING_ _EULER_
! #define _TIME_MARCHING_ _PC2_
! #define _TIME_MARCHING_ _RK2_
#define _TIME_MARCHING_ _RK3_


! #define _RECONSTRUCTION_ _NONE_
! #define _RECONSTRUCTION_ _MUSCL2_
#define _RECONSTRUCTION_ _MUSCL3_


#define _MUSCL2_LIMITER_ _MINMOD_
! #define _MUSCL2_LIMITER_ _SUPERBEE_
