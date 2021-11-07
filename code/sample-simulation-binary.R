
#-----------------------------------------------------------------------------------------------------#
#   Simulation studies for binary response variable                                                   #
#-----------------------------------------------------------------------------------------------------#

### Note:
  # The following represents a sample simulation script, which returns .csv files summarizing:
  #     1) the estimates of fixed- and random-effect parameters and variance component of the random 
  #        effect parameter under the logistic quasi-linear mixed model; and
  #     2) the proposed conditional AIC for the logistic quasi-linear mixed model
  # for candidate values of the parameter tau. The minumum value of the conditional AIC is provided
  # among the settings with convergence.

library(MASS)
library(lme4)
source("code/functions.R")



#-----------------------------------------------------------------------------------------------------
#   Function to conduct a simulation study															  
#-----------------------------------------------------------------------------------------------------

### Parameters of function
  # ni          vector          sample sizes for clusters
  # tau.true    numeric         true value of the parameter tau
  # be.true     vector          true value of the fixed effect parameter(s)
  # b.dist      character       type of the random effect distribution ("norm"/"exp"/"lnorm"/"mixture")
  # b.para      vector          specified parameter(s) of the random effect distribution. The maximum
  #								length is 2. When b.dist is "norm", b.para means variances of the 
  #           	                normal distribution. When b.dist is "exp", b.para means rate of the 
  #								exponential distribution. When b.dist is "lnorm", b.para means variances
  #								of the log-normal distribution with mean 2. When b.dist is "mixture",
  #								the random effect distribution is specified by means (0, b.para) and 
  #								variances (1,1).
  # b.vary      boolean         should the random effect parameter be generated for each simulation?
  # var.x       vector          variances of the fixed effect variable(s)
  # var.z       numeric         variances of the random effect variable
  # Tau         vector          candidate values of the parameter tau for the grid search
  # n.sim       numeric         number of simulations
  # iter.step   numeric         maximum number of the estimation steps 2 and 3 in the manuscript
  # iter.cal    numeric         maximum number of the iteration in each estimation step

sim_bin = function(ni,tau.true,be.true,b.dist,b.para,b.vary,var.x,var.z,Tau,n.sim,iter.step,iter.cal){

family = "binom"
m <<- length(ni)
n <<- sum(ni)
id = NULL
for(i in 1:m)
	id = c(id,rep(i,ni[i]))

set.seed(1)
dimb <<- length(b.para)
b.true = matrix(nrow=m,ncol=0)
for(i in 1:dimb){
	if(b.dist=="norm")
		b.true = cbind(b.true,rnorm(m,0,sqrt(b.para[i])))
	if(b.dist=="exp")
		b.true = cbind(b.true,rexp(n=m,b.para[i]))
	if(b.dist=="lnorm")
		b.true = cbind(b.true,rlnorm(n=m,2,sqrt(b.para[i])))
	if(b.dist=="mixture"){
		A = rbinom(m,1,prob=1/2)*rnorm(m,b.para,1); A[A==0] = rnorm(sum(A==0),0,1)
		b.true = cbind(b.true,A)
	}
}
for(i in 1:ncol(b.true))
	b.true[,i] = (b.true[,i]-mean(b.true[,i]))
bb.true = bb(b.true)

dv = list()
if(dimb==1){
	dv[[1]] = matrix(1,1,1)
}else{
	dv[[1]] = as.matrix(rbind(c(1,0),c(0,0)))
	dv[[2]] = as.matrix(rbind(c(0,1),c(1,0)))
	dv[[3]] = as.matrix(rbind(c(0,0),c(0,1)))
}

mit <<- 1

varyt = ifelse(b.vary,"",paste("_vary=",b.vary,sep=""))
dir0 = paste(format(Sys.time(), '%y%m%d%H%M%S'),"_",family,"_tau.true=",tau.true,"_b.dist=",b.dist,"_n.sim=",n.sim,sep="")
dir = paste(getwd(),"/result/",dir0,sep="")
dir.create(dir)

fn = paste(dir,"/true.param.csv",sep="")
write.table(t(c("tau",paste("be",0:(length(be.true)-1),sep=""),paste("b",1:(m*dimb),sep=""),"b.dist","b.para","b.vary","m","ni","var.x","var.z")),fn,sep=",",row.name=F,col.name=F,append=F)
write.table(t(c(tau.true,be.true,c(b.true),b.dist,b.para,as.character(b.vary),m,unique(ni)[1],var.x,var.z)),fn,sep=",",row.name=F,col.name=F,append=T)

for(it in 1:n.sim){
set.seed(it)
if(it%%100==0)print(it)
	if(b.vary){
		b.true = matrix(nrow=m,ncol=0)
		for(i in 1:dimb){
			if(b.dist=="norm")
				b.true = cbind(b.true,rnorm(m,0,sqrt(b.para[i])))
			if(b.dist=="exp")
				b.true = cbind(b.true,rexp(n=m,b.para[i]))
			if(b.dist=="lnorm")
				b.true = cbind(b.true,rlnorm(n=m,2,sqrt(b.para[i])))
			if(b.dist=="mixture"){
				A = rbinom(m,1,prob=1/2)*rnorm(m,b.para,1); A[A==0] = rnorm(sum(A==0),0,1)
				b.true = cbind(b.true,A)
			}
		}
		for(i in 1:ncol(b.true))
			b.true[,i] = (b.true[,i]-mean(b.true[,i]))
		bb.true = bb(b.true)
	}
	x = z = matrix(1,nrow=n,ncol=1)
	if(length(be.true)>1)for(i in 2:length(be.true)){
		x = cbind(x,rnorm(n,1,var.x))
	}
	if(dimb>1)for(i in 2:dimb)
		z = cbind(z,rnorm(n,0,var.z))
	zij <<- function(i,j){
		A = matrix(0,nrow=m*dimb,ncol=1)
		A[((i-1)*dimb+1):(i*dimb),1] = z[rij(i,j),]
		return(A)
	}
	if(tau.true==0){
		ff = x%*%be.true
		rr = apply(z*bb.true,1,sum)
		eta.true = ff + rr
		pi.true = logit.inv( eta.true )
	}else{
		ff = exp(tau.true*x%*%be.true)
		rr = exp(tau.true*apply(z*bb.true,1,sum))
		eta.true = 2/tau.true*log( 1/2*ff+1/2*rr)
		pi.true = logit.inv( eta.true )
	}
	y = matrix(rbinom(length(pi.true),1,pi.true),ncol=1)
	A = try(glmer.binom(y,x,z,id),silent=T)
	if(class(A)[1]!="try-error")
		out.glmer = A
	
	cc = NULL
	for(tau in Tau){
		fn = paste(dir,"/tau=",tau,sep="")
		be=be.glm=out.glmer$be.glm; b=b.glm=out.glmer$b.glm; th=th.glm=out.glmer$th.glm
		if(dimb==2){
			co=out.glmer$co.glm
			th=th.glm=c(th[1],sqrt(th[1]*th[2])*co,th[2])
		}
		ie = 1; dif = Inf; sdif = 10^-2

		while((!is.nan(dif) & dif>sdif) & ie<iter.step){
			be.pre=be; b.pre=b; th.pre=th
			A=try(beb.est.binom(tau,be,b,th,y,x,z,dv),silent=T)
			if(class(A)[1]!="try-error"){
				be=A$be; b=A$b
			}
			B=try(d.est.binom(tau,be,b,th,y,x,z,dv),silent=T)
			if(class(B)[1]!="try-error")
				th=B$th
			ie=ie+1
			dif = sum(abs(c(be-be.pre,c(b-b.pre),th-th.pre)))
		}
		conv = F
		if(class(A)[1]!="try-error" & class(B)[1]!="try-error")
			conv = A$upd & B$upd & !is.nan(dif) & dif<sdif
		ci = caic.binom(tau,be,b,th,y,x,z,dv)
		write.table(t(c(be)),paste(fn,"_be.csv",sep=""),sep=",",row.name=F,col.name=F,append=T)
		write.table(t(c(b)),paste(fn,"_b.csv",sep=""),sep=",",row.name=F,col.name=F,append=T)
		write.table(t(c(th)),paste(fn,"_th.csv",sep=""),sep=",",row.name=F,col.name=F,append=T)
		write.table(conv,paste(fn,"_convergence.csv",sep=""),sep=",",row.name=F,col.name=F,append=T)
		write.table(ci$caic,paste(fn,"_cAIC.csv",sep=""),sep=",",row.name=F,col.name=F,append=T)
	} # tau
} # it

Caic = Conv = NULL
for(tau in Tau){
	ct = read.table(paste(dir,"/tau=",tau,"_caic.csv",sep=""),sep=",")[1:n.sim,]
	Caic = cbind(Caic,ct)
	co = read.table(paste(dir,"/tau=",tau,"_convergence.csv",sep=""),sep=",")[1:n.sim,]
	Conv = cbind(Conv,co)
}

caic.ql = NULL; for(i in 1:nrow(Caic))caic.ql = c(caic.ql,min(Caic[i,][Conv[i,]]))
caic.gl = Caic[,which(Tau==0)]

fn = paste(dir,"/caic.csv",sep="")
write.table(t(c("GLMM(tau=0)","QLMM(tau.estimated)")),fn,sep=",",row.name=F,col.name=F,append=F)
write.table(cbind(caic.gl,caic.ql),fn,sep=",",row.name=F,col.name=F,append=T)
} # sim_binom



#-----------------------------------------------------------------------------------------------------
#   (Sample) Simulation parameter values															  
#-----------------------------------------------------------------------------------------------------

# Cluster size
m = 10

# Sample sizes for clusters
ni = rep(20,m)

# True value of tau
tau.true = -1

# True value of fixed effect parameter
be.true = matrix(c(0,1),ncol=1)

# Random effect distribution (type of the distribution and specified parameter of the distribution)
# Length of the specified parameter equals the dimension of the random effect distribution (<3)
b.dist = "norm"; b.para = 2
#b.dist = "exp"; b.para = 0.2
#b.dist = "lnorm"; b.para = 1
#b.dist = "mixture"; b.para = 3

# Should the random effect parameter be generated for each simulation?
b.vary = T

# Variances of fixed and random effect variables
var.x = rep(5,length(be.true)-1)
var.z = NA

# Candidate values of tau for the grid search
A = c(0.1,0.2,0.4,0.7,1,3)
Tau = sort(c(-A,0,A))
#Tau = c(0,tau.true)

# Number of simulations
n.sim = 1000

# Maximum number of the estimation steps 2 and 3
iter.step = 20

# Maximum number of the iteration in each estimation step
iter.cal = 100



#-----------------------------------------------------------------------------------------------------
#   Simulation run																					  
#-----------------------------------------------------------------------------------------------------

sim_bin(ni,tau.true,be.true,b.dist,b.para,b.vary,var.x,var.z,Tau,n.sim,iter.step,iter.cal)






