#body_coords function

body_coords = function(alpha,L){
  # function pos = body_coords(alpha,L)
  #
  # Converts two body angles into the four positions of the
  # 3 link swimmer - in their own body frame
  dx=-cos(alpha[1])
  dy= sin(alpha[1])
  pos=matrix(NA,ncol=2,nrow=4)
  pos[1,] = 2*L*c(dx,dy)+c(-L,0)
  pos[2,] = c(-L,0)
  pos[3,] = c(L,0)
  dx=cos(alpha[2])
  dy=sin(alpha[2])
  pos[4,] = c(L,0) + 2*L*c(dx,dy)
  
  return(pos)
}