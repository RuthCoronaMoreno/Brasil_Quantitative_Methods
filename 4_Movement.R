library(ggplot2)

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

pos=body_coords(c(pi/4,-pi/4),2)
pos=data.frame(pos)
colnames(pos)=c("alpha1", "alpha2")
ggplot(pos,aes(alpha1,alpha2)) +
  geom_path(size=1.5) +
  xlab("") +
  ylab("") +
  theme_test() +
  theme(axis.text=element_text(size=20),
        aspect.ratio = 1)


#-----------------------------------------
#Challenge Problem: Visualizing Configuration Space
#-----------------------------------------
#Extend the code snippet to visualize a single body configuration as a configuration space from -pi/2<=alpha_1<=pi/2 and -pi/2<=alpha_2<=pi/2
alpha1_vec=seq(-pi/2,pi/2, by=pi/8)
alpha2_vec=seq(-pi/2,pi/2, by=pi/8)

p=ggplot() +xlim(-2,2) +ylim(-2,2)
for(i in 1:length(alpha1_vec)){
  for(j in 1:length(alpha2_vec)){
    pos=body_coords(c(alpha1_vec[i],alpha2_vec[j]),0.05)
    pos=data.frame(pos)
    colnames(pos)=c("alpha1", "alpha2")
    pos[,1] = pos[,1]+alpha1_vec[i]
    pos[,2] = pos[,2]+alpha2_vec[j]
    p = p + geom_path(data=pos,aes(alpha1,alpha2),size=1.5)
  }
}
p = p +
  xlab(expression(alpha[1])) +
  ylab(expression(alpha[2])) +
  ggtitle("Configuration space") +
  theme_test() +
  theme(axis.text=element_text(size=27),
        axis.title = element_text(size=30),
        title = element_text(size=30),
        aspect.ratio = 1)
p


#######################################################
# ORBITS IN CONFIGURATION SPACE
#######################################################

sym_gait = function(t,T) {
  # Returns positions and velocities in angle space for
  # customized gait
  alpha1 = pi/4*cos(2*pi*t/T)
  alpha2 = pi/4*cos(2*pi*t/T)
  valpha1 = -2*pi/T*pi/4*sin(2*pi*t/T)
  valpha2 = -2*pi/T*pi/4*sin(2*pi*t/T)
  return(data.frame(alpha1,alpha2,valpha1,valpha2))
}

#Overlay the gait
gait_dyn = sym_gait(seq(0,1,0.01),1)
p2=p + geom_path(data=gait_dyn,aes(alpha1,alpha2),size=1.5)
ggplot(gait_dyn,aes(alpha1,alpha2)) +
  geom_line()
p2


#-----------------------------------------
#Challenge Problem: Challenge Problem: Designing Gaits
#-----------------------------------------
#Develop an anti-symmetric gait that reproduces the motions above and (optionally) a square-wave shape of an
#anti-symmetric gait. If you have the time, write a code snippet that reproduces the conguration snapshots
#below the top panels in the prior figure.
asym_gait = function(t,T) {
  # Returns positions and velocities in angle space for
  # customized gait
  alpha1 = pi/4*cos(2*pi*t/T)
  alpha2 = pi/4*sin(2*pi*t/T)
  valpha1 = -2*pi/T*pi/4*sin(2*pi*t/T)
  valpha2 = 2*pi/T*pi/4*cos(2*pi*t/T)
  return(data.frame(alpha1,alpha2,valpha1,valpha2))
}

gait_dyn_2 = asym_gait(seq(0,1,0.01),1)
p3=p + geom_path(data=gait_dyn_2,aes(alpha1,alpha2),size=1.5)
ggplot(gait_dyn_2,aes(alpha1,alpha2)) +
  geom_line()
p3


########################################
#FROM BORELLI TO NEWTON AND BACK AGAIN
########################################
#-----------------------------------------
#Challenge Problem: Coasting Time and Distance
#-----------------------------------------
#Consider a non-motile organism of mass m with an initial velocity of v0 in an environment with drag k. How far
#will the organism coast before it reaches a complete stop? What is the effective time over which the organism
#remains in motion? In answering this problem, use a simulation and (if possible) mathematical analysis to
#verify your finding.

#Dynamic code
model_coast = function(t,y,pars){
  # function dydt = model_coast(t,y,pars)
  dydt=rep(0,2)
  dydt[1]=y[2]
  dydt[2]=-pars$k*y[2]/pars$m
  return(list(dydt))
}

#Simulation and visualization code
require(deSolve)
require(matlab)

# Simulation
mrange=logspace(-1,1,10)
v0=1
pars=list()
yf=c()
for(i in 1:length(mrange)){
  pars$m=mrange[i]
  pars$k=1
  y=ode(c(0,v0),c(0,100),model_coast,pars)
  yf[i]=y[nrow(y),2]
}
# Plot simulation against theory
m=logspace(-1,1,1000)
myf=yf
ggplot(grid(seq(0,10,1),seq(0,10,1))) +
  geom_line(aes(m,m*v0/pars$k),size=1.2) +
  geom_point(aes(mrange,myf),size=3) +
  scale_x_continuous(trans = "log",limits = c(.1,10),breaks=c(.1,1,10)) +
  scale_y_continuous(trans = "log",limits = c(.1,10),breaks=c(.1,1,10)) +
  xlab("Mass") +
  ylab('Coasting distance') +
  theme_test() +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=18),
        aspect.ratio = 1)

#--------------------------------
#Borelli: Connecting External Forces and Internal Body Movements
#--------------------------------
#Complete simulation code

model_4bead = function(t,y,pars){
  # function dydt = model_4bead(t,y,pars)
  # RFT dynamics for a 3 link swimmer
  # y -> bx, by, theta, vbx, vby, vtheta, alpha1, alpha2
  # A 8-coordinate dynamics, only 6 of which are considered here
  # the other two are prescribed by the initernal gait dynamics
  # Coords
  xb=y[1]
  yb=y[2]
  theta=y[3]
  vxb=y[4]
  vyb=y[5]
  vtheta=y[6]
  gait_dyn = pars$gait_function(t,pars$T)
  alpha1=gait_dyn$alpha1
  valpha1=gait_dyn$valpha1
  alpha2=gait_dyn$alpha2
  valpha2=gait_dyn$valpha2
  L=pars$L
  # Point positions
  # *----*----*----*
  # 1 2 3 4
  rotmat = function(theta){
    A = matrix(c(cos(theta), -sin(theta),
                 sin(theta), cos(theta)),ncol=2,byrow = T)
  }
  vec_to_angle = function(x,y){
    return(atan2(y,x))
  }
  angle_to_vec = function(theta){
    return(c(cos(theta),sin(theta)))
  }
  A = rotmat(theta+pi-alpha1)
  dpos1=A%*%c(2*L,0)
  A = rotmat(theta+alpha2)
  #dpos2 = A*c(2*L,0)
  dpos2=A%*%c(2*L,0)
  bpos= matrix(c(xb-L*cos(theta)+dpos1[1], yb-L*sin(theta)+dpos1[2],
                 xb-L*cos(theta), yb-L*sin(theta),
                 xb+L*cos(theta), yb+L*sin(theta),
                 xb+L*cos(theta)+dpos2[1], yb+L*sin(theta)+dpos2[2]),ncol=2,byrow = TRUE)
  # Link vectors for anisotropy
  link=matrix(NA,ncol=length(bpos[1,]),nrow=4)
  link[1,]=bpos[1,]-bpos[2,]
  link[2,]=bpos[2,]-bpos[3,]
  link[3,]=bpos[3,]-bpos[2,]
  link[4,]=bpos[4,]-bpos[3,]
  # Forces on each bead, due to drag
  bvel_vec = c(vxb,vyb)
  fparr=matrix(NA,nrow=4,ncol=ncol(link))
  fperp=matrix(NA,nrow=4,ncol=ncol(link))
  for(i in 1:4){
    vec_projection = dot(bvel_vec,link[i,])/(2*L)*link[i,]/(2*L)
    vec_rejection = bvel_vec-dot(bvel_vec,link[i,])/(2*L)^2*link[i,]
    fparr[i,]=-pars$kparr*vec_projection
    fperp[i,]=-pars$kperp*vec_rejection
  }
  fparr_bead=matrix(NA,ncol=2,nrow=2)
  fperp_bead=matrix(NA,ncol=2,nrow=2)
  # Forces on the movable beads
  bead_vangle1=vec_to_angle(link[1,1],link[1,2])-pi/2*sign(valpha1)
  bead_vvec1 = angle_to_vec(bead_vangle1)
  vec_projection = dot(bead_vvec1,link[1,])/(2*L)*link[1,]/(2*L)
  vec_rejection = bead_vvec1-dot(bead_vvec1,link[1,])/(2*L)^2*link[1,]
  fparr_bead[1,]=-pars$kparr*abs(valpha1)*2*L*vec_projection
  fperp_bead[1,]=-pars$kperp*abs(valpha1)*2*L*vec_rejection
  bead_vangle2=vec_to_angle(link[4,1],link[4,2])+pi/2*sign(valpha2)
  bead_vvec2 = angle_to_vec(bead_vangle2)
  vec_projection = dot(bead_vvec2,link[4,])/(2*L)*link[4,]/(2*L)
  vec_rejection = bead_vvec2-dot(bead_vvec2,link[4,])/(2*L)^2*link[4,]
  fparr_bead[2,]=-pars$kparr*abs(valpha2)*2*L*vec_projection
  fperp_bead[2,]=-pars$kperp*abs(valpha2)*2*L*vec_rejection
  # Total force
  FF = colSums(fparr)+colSums(fperp)+colSums(fparr_bead)+colSums(fperp_bead)
  
  # Equations of motion
  dydt=c()
  dydt[1] = vxb
  dydt[2] = vyb
  dydt[3] = vtheta
  dydt[4] = pars$damping*FF[1]/pars$M
  dydt[5] = pars$damping*FF[2]/pars$M
  dydt[6] = 0
  return(list(dydt))
}

#Master script
library(transformr)
require(pracma)
require(deSolve)
require(ggplot2)
require(gganimate)
source("asym_gait.R")
source("body_coords.R")
source("body_rotate.R")


# Parameters
pars=list()
pars$L=1
pars$kparr=0.2
pars$kperp=2
pars$max_angle = pi/4
pars$damping=1
pars$M=0.1
pars$gait_function = asym_gait
pars$T=40
# Simulate
y0=c(0,0,0,0,0,0.0) # 6 equations, all in world frame
t=seq(0,200,0.2)
y=ode(y0,t,model_4bead,pars)#,method="ode45")
# Recover the gaits
alpha1=c()
alpha2=c()
valpha1=c()
valpha2=c()
for (i in 1:length(y[,1])){
  gait_dyn=pars$gait_function(y[,1][i],pars$T)
  alpha1[i]=gait_dyn[[1]]
  alpha2[i]=gait_dyn[[2]]
  valpha1[i]=gait_dyn[[3]]
  valpha2[i]=gait_dyn[[4]]
}
#Visualize
construct_df = function(i){
  pos = body_coords(c(alpha1[[i]],alpha2[[i]]),pars$L)
  rpos = body_rotate(pos,y[i,4])
  vis_pos = data.frame("x"=rpos[,1] + y[i,2],"y"=rpos[,2] + y[i,3],"time"=i)
  return(vis_pos)
}
df=lapply(1:length(t),construct_df)
df=do.call(rbind,df)
p=ggplot(df,aes(x,y,group=time)) +
  geom_line() +
  transition_time(time) +
  theme(axis.text = element_text(size=10),
        axis.title = element_text(size=15),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linetype="solid",fill=NA),
        aspect.ratio = 1) +
  xlim(-10,10) +
  ylim(-10,10) +
  xlab("") +
  ylab("") +
  theme_test() +
  theme(axis.text=element_text(size=20),
        aspect.ratio = 1)
print(p)



