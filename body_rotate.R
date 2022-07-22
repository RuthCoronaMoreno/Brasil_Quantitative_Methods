body_rotate =function(pos,theta){
# function rpos = body_rotate(pos,theta)
  rmat = matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),nrow=2,byrow=T)
  rpos = rmat%*%t(pos)
  rpos = t(rpos)
return(data.frame(rpos))
}
