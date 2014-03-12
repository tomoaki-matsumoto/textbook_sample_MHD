! Header file for configuration

! flux and EOS
#define FLUX_SCHEME_SCALAR_ADVECTION
! #define FLUX_SCHEME_HD
! #define FLUX_SCHEME_MHD


#define DIRECTIONAL_UNSPLIT
! #define DIRECTIONAL_SPLIT


! #define TIME_MARCHING_EULER
! #define TIME_MARCHING_PC2
! #define TIME_MARCHING_RK2
#define TIME_MARCHING_RK3
! #define TIME_MARCHING_K3


! #define RECONSTRUCTION_NONE
! #define RECONSTRUCTION_MUSCL2
#define RECONSTRUCTION_MUSCL3


#define MUSCL2_LIMITER_MINMOD
! #define MUSCL2_LIMITER_SUPERBEE

