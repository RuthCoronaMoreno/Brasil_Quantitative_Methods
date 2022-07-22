asym_gait = function(t,T) {
  # function [alpha1 alpha2 valpha1 valpha2]=asym_gait(t,T)
  #
  # Returns positions and velocities in angle space for
  # customized gait
  alpha1 = pi/4*cos(2*pi*t/T)
  alpha2 = pi/4*sin(2*pi*t/T)
  valpha1 = -2*pi/T*pi/4*sin(2*pi*t/T)
  valpha2 = 2*pi/T*pi/4*cos(2*pi*t/T)
  return(data.frame(alpha1,alpha2,valpha1,valpha2))
}
