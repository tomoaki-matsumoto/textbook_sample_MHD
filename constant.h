! Definition of constants for configuration

! flux and EOS
#define _SCALAR_ADVECTION_  0
#define _HD_  1
#define _MHD_ 2

! Directional split
#define _UNSPLIT_ 0
#define _SPLIT_   1

! Time marching
#define _EULER_   0
#define _PC2_     1
#define _RK2_     2
#define _RK3_     3

! Reconstruction
#define _NONE_    0
#define _MUSCL2_  1
#define _MUSCL3_  2

! Slope limiter for MUSCL2
#define _MINMOD_   0
#define _SUPERBEE_ 1

