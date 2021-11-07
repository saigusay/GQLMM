
#-----------------------------------------------------------------------------------------------------#
#   Simulation studies for binary response variable in which the random slope parameter is omitted    #
#-----------------------------------------------------------------------------------------------------#

### Note:
  # The following represents a sample simulation script, which returns .csv files summarizing:
  #     1) the estimates of fixed- and random-effect parameters and variance component of the random 
  #        effect parameter under the logistic quasi-linear mixed model without the random slope 
  #        and the (correctly-specified) logistic linear mixed model with random slope; and
  #     2) the proposed conditional AIC for the logistic quasi-linear mixed model and the 
  #        (correctly-specified) logistic linear mixed model
  # for candidate values of the parameter tau. The minumum value of the conditional AIC is provided
  # among the settings with convergence for the quasi-linear model.

library(MASS)
library(lme4)
source("code/functions.R")



#-----------------------------------------------------------------------------------------------------
#   Function to conduct a simulation study															  
#-----------------------------------------------------------------------------------------------------

### Parameters of function
  # ni          vector          sample sizes for clusters
  # b.para      vector          variance-covariance of the bivariate normal random effect distribution
  # b.vary      boolean         should the random effect parameter be generated for each simulation?
  # var.x       vector          variances of the fixed effect variable(s). Note: The first fixed effect 
  #								variable is also used as the random slope variable.
  # Tau         vector          candidate values of the parameter tau for the grid search
  # n.sim       numeric         number of simulations
  # iter.step   numeric         maximum number of the estimation steps 2 and 3 in the manuscript
  # iter.cal    numeric         maximum number of the iteration in each estimation step

sim_bin_misspe = function(ni,be.true,b.para,b.vary,var.x,Tau,n.sim,iter.step,iter.cal){

family = "binom"
m <<- length(ni)
n <<- sum(ni)
id = NULL
for(i in 1:m)
	id = c(id,rep(i,ni[i]))

set.seed(1)
dimb <<- 2
cov.true = matrix(rbind(c(b.para[1],b.para[2]),c(b.para[2],b.para[3])),2,2)
b.true = mvrnorm(m,rep(0,2),cov.true)

for(i in 1:ncol(b.true))
	b.true[,i] = (b.true[,i]-mean(b.true[,i]))
bb.true = bb(b.true)

dv0 = dvr = list()
dv0[[1]] = matrix(1,1,1)
dvr[[1]] = as.matrix(rbind(c(1,0),c(0,0))); dvr[[2]] = as.matrix(rbind(c(0,1),c(1,0))); dvr[[3]] = as.matrix(rbind(c(0,0),c(0,1)))

mit <<- 1/2

varyt = ifelse(b.vary,"",paste("_vary=",b.vary,sep=""))
dir0 = paste(format(Sys.time(), '%y%m%d%H%M%S'),"_",family,"-misspe_n.sim=",n.sim,sep="")
dir = paste(getwd(),"/result/",dir0,sep="")
dir.create(dir)

fn = paste(dir,"/true.param.csv",sep="")
write.table(t(c(paste("be",0:(length(be.true)-1),sep=""),paste("b",1:(m*dimb),sep=""),"b.vary","b.para")),fn,sep=",",row.name=F,col.name=F,append=F)
write.table(t(c(be.true,c(b.true),as.character(b.vary),b.para)),fn,sep=",",row.name=F,col.name=F,append=T)
write.table(t(c("var.x")),fn,sep=",",row.name=F,col.name=F,append=T)
write.table(t(c(var.x)),fn,sep=",",row.name=F,col.name=F,append=T)

write.table(t(c("tau",paste("be",0:(length(be.true)-1),sep=""),paste("b",1:(m*dimb),sep=""),"b.para","","","b.vary","m","ni","var.x")),fn,sep=",",row.name=F,col.name=F,append=F)
write.table(t(c(0,be.true,c(b.true),b.para,as.character(b.vary),m,unique(ni)[1],var.x)),fn,sep=",",row.name=F,col.name=F,append=T)

for(it in 1:n.sim){
set.seed(it)
	if(b.vary){
		b.true = mvrnorm(m,rep(0,2),cov.true)
		for(i in 1:ncol(b.true))
			b.true[,i] = (b.true[,i]-mean(b.true[,i]))
		bb.true = bb(b.true)
	}
	x = z = matrix(1,nrow=n,ncol=1)
	if(length(be.true)>1)for(i in 2:length(be.true)){
		x = cbind(x,rnorm(n,1,var.x))
	}
	z = cbind(z,x[,2])
	zij <<- function(i,j){
		A = matrix(0,nrow=m*dimb,ncol=1)
		A[((i-1)*dimb+1):(i*dimb),1] = z[rij(i,j),]
		return(A)
	}
	ff = x%*%be.true
	rr = apply(z*bb.true,1,sum)
	eta.true = ff + rr
	pi.true = logit.inv( eta.true )
	y = matrix(rbinom(length(pi.true),1,pi.true),ncol=1)

	x0 = x; z0 = as.matrix(z[,1],ncol=1)
	out.glmer = glmer.binom(y,x0,z0,id)
	xr = x; zr = z; dv = dvr
	outr.glmer = glmer.binom(y,xr,zr,id)
	
	cc = NULL; dif = Inf; sdif = 10^-2
	for(tau in Tau){
		fn = paste(dir,"/tau=",tau,sep="")
		be=out.glmer$be.glm; b=out.glmer$b.glm; th=out.glmer$th.glm
		x = x0; z = z0; dv = dv0; dimb <<- ncol(b.true)-1
		ie = 1
		while((!is.nan(dif) & dif>sdif) & ie<iter.step){
			be.pre=be; b.pre=b; th.pre=th
			A=beb.est.binom(tau,be,b,th,y,x,z,dv)
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

		fn = paste(dir,"/tau=",tau,"slope",sep="")
		be=outr.glmer$be.glm; b=matrix(outr.glmer$b.glm,nrow=m); th=outr.glmer$th.glm; co=outr.glmer$co.glm;
		th[th<10^-10] = 10^-10; th[10^10<th] = 10^10; co[is.nan(co)] = 0
		th=c(th[1],sqrt(th[1]*th[2])*co,th[2])
		x=xr; z = zr; dv = dvr; dimb <<- ncol(b.true)
		ie = 1
		while((!is.nan(dif) & dif>sdif) & ie<iter.step){
			be.pre=be; b.pre=b; th.pre=th
			A=beb.est.binom(tau,be,b,th,y,x,z,dv)
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
if(1==0){
Caic = Conv = Caicr = Convr = NULL
for(tau in Tau){
	ct = read.table(paste(dir,"/tau=",tau,"_caic.csv",sep=""),sep=",")[1:n.sim,]
	Caic = cbind(Caic,ct)
	co = read.table(paste(dir,"/tau=",tau,"_convergence.csv",sep=""),sep=",")[1:n.sim,]
	Conv = cbind(Conv,co)
	ct = read.table(paste(dir,"/tau=",tau,"slope_caic.csv",sep=""),sep=",")[1:n.sim,]
	Caicr = cbind(Caicr,ct)
	co = read.table(paste(dir,"/tau=",tau,"slope_convergence.csv",sep=""),sep=",")[1:n.sim,]
	Convr = cbind(Convr,co)
}

caic.gl = Caic[,which(Tau==0)]
caic.ql = NULL; for(i in 1:nrow(Caic))caic.ql = c(caic.ql,min(Caic[i,][Conv[i,]]))
caic.glr = Caicr[,which(Tau==0)]
caic.qlr = NULL; for(i in 1:nrow(Caicr))caic.qlr = c(caic.qlr,min(Caicr[i,][Convr[i,]]))

fn = paste(dir,"/caic.csv",sep="")
write.table(t(c("GLMM(tau=0)","GLMM(random-slope)","QLMM","QLMM(random-slope)")),fn,sep=",",row.name=F,col.name=F,append=F)
write.table(cbind(caic.gl,caic.glr,caic.ql,caic.qlr),fn,sep=",",row.name=F,col.name=F,append=T)
}
} # sim_binom



#-----------------------------------------------------------------------------------------------------
#   (Sample) Simulation parameter values															  
#-----------------------------------------------------------------------------------------------------

# Cluster size
m = 5

# Sample sizes for clusters
ni = rep(20,m)

# True value of fixed effect parameter
be.true = matrix(c(1,-1),ncol=1)

# True value of parameter of bivariate normal distribution of random effect (variance1, covariance, variance2)
b.para = c(2,sqrt(2*2)*0,2)

# Should the random effect parameter be generated for each simulation?
b.vary = T

# Variances of fixed and random effect variables
var.x = rep(3,length(be.true)-1)

# Candidate values of tau for the grid search
# Including -Inf indicates applying the correctly-specified model
A = c(0.1,0.2,0.4,0.7,1,3)
Tau = sort(c(-A,0,A))

# Number of simulation
n.sim = 1000

# Maximum number of the estimation steps 2 and 3
iter.step = 20

# Maximum number of the iteration in each estimation step
iter.cal = 200



#-----------------------------------------------------------------------------------------------------
#   Simulation run																					  
#-----------------------------------------------------------------------------------------------------

sim_bin_misspe(ni,be.true,b.para,b.vary,var.x,Tau,n.sim,iter.step,iter.cal)





