! Header file for configuration

! -----------------
! flux and EOS
! -----------------
! #define FLUX_SCHEME_SCALAR_ADVECTION
! #define FLUX_SCHEME_HD
#define FLUX_SCHEME_MHD

! ----------------------------
! Directinal split or unsplit
! ----------------------------
#define DIRECTIONAL_UNSPLIT
! #define DIRECTIONAL_SPLIT

! ------------------------
! Method of time marching
! ------------------------
! #define TIME_MARCHING_EULER
#define TIME_MARCHING_PC2
! #define TIME_MARCHING_RK2
! #define TIME_MARCHING_RK3
! #define TIME_MARCHING_K3

! --------------------------------
! Method of spatial reconstruction
! --------------------------------
! #define RECONSTRUCTION_NONE
#define RECONSTRUCTION_MUSCL2
! #define RECONSTRUCTION_MUSCL3
! #define RECONSTRUCTION_LIMO3

! -----------------------------
! Limiter for 2nd order MUSCL
! -----------------------------
! #define MUSCL2_WO_LIMITER
#define MUSCL2_LIMITER_MINMOD
! #define MUSCL2_LIMITER_SUPERBEE


! ------------------------------
! IO : binary (default) or text
! ------------------------------
! #define IO_TEXT
