program main
  real :: positive, negative, zero
  real :: nan, inf
  positive = 1.0
  negative = -1.0
  zero = 0.0
  inf = positive/zero
  nan = sqrt(negative)
  print *, inf, negative/zero, inf**2, nan, nan*0, inf*0

  if (nan /= nan) print *, nan , 'is NaN.'
  if (inf /= inf) print *, inf , 'is NaN.'
  if ( (nan*0) /= 0) print *, nan , 'is number'
  if ( (inf*0) /= 0) print *, inf , 'is number'


end program main
