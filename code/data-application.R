
#-----------------------------------------------------------------------------------------------------#
#   Applied Example: Respiratory data in Davis (1991)										          #
#-----------------------------------------------------------------------------------------------------#

library(MASS)
library(lme4)
source("code/functions.R")



#-----------------------------------------------------------------------------------------------------
#   Function to select the variable and to estimate the regression parameters						  
#-----------------------------------------------------------------------------------------------------

### Parameters of function
  # dat         dataframe       data set
  # vy			character		response variable
  # ydist		character		distribution of response variable can be "binomial" or "poisson"
  # vc			character		cluster variable
  # vf          vector          fixed effect variable
  # vr          character       random effect variable
  # force		vector			variable to be retained
  # Tau         vector          candidate values of the parameter tau for the grid search
  # selvari     boolean			should variables be selected by backward stepwise algorithm based on 
  # 							the conditional AIC? When `selvari' is FALSE, the median odds ratios 
  #								for clusters are estimated.
  # vor			character		variable for the odds ratio
  # iter.step   numeric         maximum number of the estimation steps 2 and 3 in the manuscript
  # iter.cal    numeric         maximum number of the iteration in each estimation step
app_resp = function(dat,vy,ydist,vc,vf,vr,force,Tau,selvari,vor,iter.step,iter.cal){

dir = paste(getwd(),"/result/",format(Sys.time(), "%y%m%d%H%M%S"),"_tau=0_no-interaction",sep="")
dir = paste(getwd(),"/result/201009164348_Davis1991_no-interaction",sep="")
dir.create(dir)

if(ydist=="binomial"){
	beb.est <<- beb.est.binom; d.est <<- d.est.binom; caicf <<- caic.binom
}else if(ydist=="poisson"){
	beb.est <<- beb.est.pois; d.est <<- d.est.pois; caicf <<- caic.pois
}

A = unique(dat[,vc])
ni = NULL; for(i in A)ni = c(ni,sum(dat[,vc]==i))
ni <<- ni
m <<- length(ni)
n <<- sum(ni)
id = NULL
for(i in 1:m)
	id = c(id,rep(i,ni[i]))

A = setdiff(vf,colnames(dat))
if(length(A)>0){
	for(v in A){
		vv = strsplit(v,"\\.")[[1]]
		xv = dat[,vv[1]]*dat[,vv[2]]
		dat = cbind(dat,xv); colnames(dat)[ncol(dat)] = v
	}
}
y = as.matrix(dat[,vy],ncol=1); colnames(y) = "y"
x = as.matrix(cbind(1,dat[,vf])); colnames(x) = c(1,vf)
if(vr!=""){
	z = as.matrix(cbind(1,dat[,vr])); colnames(z) = c(1,paste("random.",vr,sep=""))
}else{
	z = matrix(1,nrow=n); colnames(z) = c(1)
}
x0 = x; z0 = z
vf0 = colnames(x0)[-1]; vr0 = colnames(z0)[-1]

mit <<- 1/2

set.seed(1)
for(tau in Tau){
	if(selvari)
		write.table(t(c(vf0,vr0,"cAIC","effective df","-logL","conv")),paste(dir,"/summary_tau=",tau,".csv",sep=""),sep=",",row.name=F,col.name=F,append=F)
	step = 1; cc = min.caic = Inf
	while(min(cc)==min.caic&(selvari|step==1)){
		if(step==1){
			vfs = list(); vfs[[1]] = vf0
			vrs = list(); vrs[[1]] = vr0
		}else{
			vf = vfs[[which(cc==min(cc))]]; vr = vrs[[which(cc==min(cc))]]; il = 1
			vfs = list(); vrs = list()
			if(length(vf)>0)for(jf in 1:length(vf)){
				if(!is.element(vf[jf],force)){
					vfs[[il]] = vf[-jf]; vrs[[il]] = vr; il=il+1
				}
			}
			if(length(vr)>0)for(jr in 1:length(vr)){
				vfs[[il]] = vf; vrs[[il]] = vr[-jr]; il=il+1
			}
		}
		cc = NULL
		if(length(vfs)>0)for(iv in 1:length(vfs)){
			vf = vfs[[iv]]; vr = vrs[[iv]]
			x = as.matrix(cbind(1,x0[,vf])); colnames(x) = c(1,vf)
			z = as.matrix(cbind(1,z0[,vr])); colnames(z) = c(1,paste(vr,sep=""))
			zij <<- function(i,j){ #m ni
				A = matrix(0,nrow=m*dimb,ncol=1)
				A[((i-1)*dimb+1):(i*dimb),1] = z[rij(i,j),]
				return(A)
			}
			dimb <<- ncol(z)
			model.glm = "y~"
			for(v in colnames(x)[-1])
				model.glm = paste(model.glm,"+",v,sep="")
			rr = "1"
			for(v in colnames(z)[-1])
				rr = paste(rr,"+",v,sep="")
			model.glm = paste(model.glm,"+(",rr,"|id)",sep="")
			dat = data.frame(as.matrix(cbind(id,y,x[,-1],z[,-1]))); colnames(dat) = c("id","y",colnames(x)[-1],colnames(z)[-1])
			A = glmer(as.formula(model.glm),data=dat,family=ydist)
			be.glm = matrix(summary(A)$coef[,1],ncol=1)
			b.glm = matrix(unlist(ranef(A)),ncol=dimb)
			sd.glm = matrix(attr(summary(A)$varcor$id,"stddev"),ncol=1)
			cr.glm = as.matrix(attr(summary(A)$varcor$id,"correlation"))
			d.glm = sd.glm%*%t(sd.glm)*cr.glm
			if(nrow(d.glm)==2)d.glm[1,2] = d.glm[2,1]
			th.glm = unique(c(d.glm))

			dv = list()
			if(dimb==1){
				dv[[1]] = matrix(1,1,1)
			}else{
				dv[[1]] = as.matrix(rbind(c(1,0),c(0,0)))
				dv[[2]] = as.matrix(rbind(c(0,1),c(1,0)))
				dv[[3]] = as.matrix(rbind(c(0,0),c(0,1)))
			}
			be=be.glm; b=b.glm; th=th.glm
			ie = 1; dif = Inf; sdif = 10^-4
			while((!is.nan(dif) & dif>sdif) & ie<iter.step){
				be.pre=be; b.pre=b; th.pre=th
				A=try(beb.est(tau,be,b,th,y,x,z,dv),silent=T)
				if(class(A)!="try-error"){
					be=A$be; b=A$b
				}
				B=try(d.est(tau,be,b,th,y,x,z,dv),silent=T)
				if(class(B)!="try-error")
					th=B$th
				ie=ie+1
				dif = abs(sum(c(be-be.pre,c(b-b.pre),th-th.pre)))
			}
			conv = F
			if(class(A)!="try-error" & class(B)!="try-error")
				conv = A$upd & B$upd & !is.nan(dif) & dif<sdif
			A = caicf(tau,be,b,th,y,x,z,dv)
			caic=A$caic; edf = A$ef.df; mL = A$L
			out = NULL
			for(v in vf0)out = c(out,ifelse(is.element(v,vf),1,0))
			for(v in vr0)out = c(out,ifelse(is.element(v,vr),1,0))
			out = c(out,caic,edf,mL,conv)
			if(selvari)
				write.table(t(out),paste(dir,"/summary_tau=",tau,".csv",sep=""),sep=",",row.name=F,col.name=F,append=T)
			cc = c(cc,ifelse(!is.na(caic),caic,Inf))
			min.caic = ifelse(caic<min.caic,caic,min.caic)
			if(!selvari&!is.null(vor)&ydist=="binomial"){
				lor = p = q = rep(NA,m)
				uid = unique(id)
				for(i in uid){
					xi1 = xi0 = apply(x[dat$id==i,],2,median)
					zi = apply(as.matrix(z[dat$id==i,],ncol=1),2,median)
					if(tau==0){
						p[i] = logit.inv(xi1%*%be+zi%*%b[i,])
					}else{
						p[i] = logit.inv(2/tau*log(1/2*exp(tau*xi1%*%be)+1/2*exp(tau*zi%*%b[i,])))
					}
					xi1[vor] = 1; xi0[vor] = 0
					A = grep(vor,colnames(x)); A = setdiff(A,which(colnames(x)==vor))
					if(length(A)>0){
						for(v in colnames(x)[A]){
							vv = strsplit(v,"\\.")[[1]]
							xi1[v] = xi1[vv[1]]*xi1[vv[2]]
							xi0[v] = xi0[vv[1]]*xi0[vv[2]]
						}
					}
					zi1=zi0=zi
					if(is.element(paste("random.",vor,sep=""),colnames(z))){
						zi1[paste("random.",vor,sep="")] = 1; zi0[paste("random.",vor,sep="")] = 0
					}
					if(tau==0){
						q[i] = 1/2
						lor[i] = (xi1%*%be+zi1%*%b[i,]) - (xi0%*%be+zi0%*%b[i,])
					}else{
						q[i] = exp(tau*xi0%*%be)/(exp(tau*xi0%*%be)+exp(tau*zi0%*%b[i,]))
						lor[i] = 2/tau*log(1/2*exp(tau*xi1%*%be)+1/2*exp(tau*zi1%*%b[i,])) - 2/tau*log(1/2*exp(tau*xi0%*%be)+1/2*exp(tau*zi0%*%b[i,]))
					}
				}
				write.table(t(c("intercept",vf)),paste(dir,"/vf_tau=",tau,".csv",sep=""),sep=",",row.name=F,col.name=F,append=F)
				write.table(t(c("random.intercept",vr)),paste(dir,"/vr_tau=",tau,".csv",sep=""),sep=",",row.name=F,col.name=F,append=F)
				write.table(t(c(be)),paste(dir,"/be_tau=",tau,".csv",sep=""),sep=",",row.name=F,col.name=F,append=F)
				write.table(t(c(b)),paste(dir,"/b_tau=",tau,".csv",sep=""),sep=",",row.name=F,col.name=F,append=F)
				write.table(t(c(th)),paste(dir,"/th_tau=",tau,".csv",sep=""),sep=",",row.name=F,col.name=F,append=F)
				write.table(q,paste(dir,"/q_tau=",tau,".csv",sep=""),sep=",",row.name=F,col.name=F,append=F)
				write.table(p,paste(dir,"/p_tau=",tau,".csv",sep=""),sep=",",row.name=F,col.name=F,append=F)
				write.table(lor,paste(dir,"/log-odds-ratio_tau=",tau,".csv",sep=""),sep=",",row.name=F,col.name=F,append=F)
			}
		} # vfs vrs
		step = step+1
	} #step
} #tau

} # app_resp



#-----------------------------------------------------------------------------------------------------
#   Parameter settings																				  
#-----------------------------------------------------------------------------------------------------

# Importing the data
# Note: the data structure are given as
#    center   subject   treat   sex   age   baseline   visit   status
#    1        1         P       M     46    0          1       0
#    1        1         P       M     46    0          2       0
#    1        1         P       M     46    0          3       0
#    1        1         P       M     46    0          4       0
#    1        2         P       M     28    0          1       0
#    ...
dat = read.csv("data/Davis1991_respiratory_stack.csv",sep=",")
dat[,"center"] = as.numeric(dat[,"center"]==2)
dat[,"treat"] = as.numeric(dat[,"treat"]=="A")
dat[,"sex"] = as.numeric(dat[,"sex"]=="M")

# Response and cluster variable
vy = "status"; ydist = "binomial"
vc = "subject"

# Fixed effect variable
# Note: interaction term must be given as "variable1.variable2"
vf = c("center","treat","age","baseline")
#vf = c("center","treat","age","baseline","treat.baseline")
vf = c("center","treat","age","baseline","visit")

# Random effect variable
vr = "treat"

# Variables to be retained
force = c("treat","baseline")

# Candidate values of tau for the grid search
A = c(0.1,0.2,0.4,0.7,1,3)
Tau = sort(c(A))
Tau = 0.7

# should variables be selected?
selvari = F

# Variable for odds ratio
vor = "treat"

# Maximum number of the estimation steps 2 and 3
iter.step = 20

# Maximum number of the iteration in each estimation step
iter.cal = 100



#-----------------------------------------------------------------------------------------------------
#   Data analysis run																					
#-----------------------------------------------------------------------------------------------------

app_resp(dat,vy,ydist,vc,vf,vr,force,Tau,selvari,vor,iter.step,iter.cal)





#-----------------------------------------------------------------------------------------------------
#   Function to obtain the bootstrap estimates of the regression parameters for binary response		  
#-----------------------------------------------------------------------------------------------------

### Parameters of function
  # dat         dataframe       data set
  # vy			character		response variable
  # vc			character		cluster variable
  # vf          vector          fixed effect variable
  # vr          character       random effect variable
  # dstr		dataframe		data set for stratified sampling
  # tau         numeric         value of tau parameter
  # vor			character		variable for the odds ratio
  # outdir		character		output directory including the estimates of regression parameters by
  #								the app_resp function
  # iter.boot   numeric         number of times of bootstrap replicates
  # iter.step   numeric         maximum number of the estimation steps 2 and 3 in the manuscript
  # iter.cal    numeric         maximum number of the iteration in each estimation step
app_boot = function(dat,vy,vc,vf,vr,dstr,tau,vor,outdir,iter.boot,iter.step,iter.cal){

#setwd(outdir)

A = unique(dat[,vc])
ni = NULL; for(i in A)ni = c(ni,sum(dat[,vc]==i))
ni <<- ni
m <<- length(ni)
n <<- sum(ni)
A = setdiff(vf,colnames(dat))
if(length(A)>0){
	for(v in A){
		vv = strsplit(v,"\\.")[[1]]
		xv = dat[,vv[1]]*dat[,vv[2]]
		dat = cbind(dat,xv); colnames(dat)[ncol(dat)] = v
	}
}

id = dstr[,1]; uid = unique(id)
if(!file.exists(paste(outdir,"/id_bootstrap.csv",sep=""))){
	rr = NULL; for(i in uid)rr = c(rr,which(id==i)[1])
	dstrc = as.matrix(dstr[rr,-1],nrow=nrow(dstr))
	cond = list()
	for(i in 1:ncol(dstrc))cond[[i]] = sort(unique(dstrc[,i]))
	cc = expand.grid(cond)
	ids = list()
	for(ic in 1:nrow(cc)){
		A = matrix(T,nrow=length(uid),ncol=1)
		for(i in 1:length(cond))A = A&dstrc[,i]==cc[ic,i]
		ids[[ic]] = which(A)
	}
	for(ib in 1:iter.boot){
set.seed(ib)
		idi = NULL
		for(i in 1:length(ids)){
			if(length(ids[[i]])==0){
			}else if(length(ids[[i]])==1){
				idi = c(idi,ids[[i]])
			}else{
				idi = c(idi,sample(ids[[i]],replace=T))
			}
		}
		write.table(t(idi),paste(outdir,"/id_bootstrap.csv",sep=""),sep=",",row.name=F,col.name=F,append=T)
	}
}
idd = read.table(paste(outdir,"/id_bootstrap.csv",sep=""),sep=",")

y = as.matrix(dat[,vy],ncol=1); colnames(y) = "y"
x = as.matrix(cbind(1,dat[,vf])); colnames(x) = c(1,vf)
if(vr!=""){
	z = as.matrix(cbind(1,dat[,vr])); colnames(z) = c(1,paste("random.",vr,sep=""))
}else{
	z = matrix(1,nrow=n); colnames(z) = c(1)
}
dimb <<- ncol(z)
dv = list()
if(dimb==1){
	dv[[1]] = matrix(1,1,1)
}else{
	dv[[1]] = as.matrix(rbind(c(1,0),c(0,0)))
	dv[[2]] = as.matrix(rbind(c(0,1),c(1,0)))
	dv[[3]] = as.matrix(rbind(c(0,0),c(0,1)))
}
x0 = x1 = x
x1[,vor] = 1; x0[,vor] = 0
ict = colnames(x)[setdiff(grep(vor,colnames(x)),which(colnames(x)==vor))]
if(length(ict)>0){
	for(v in ict){
		vv = strsplit(v,"\\.")[[1]]
		x1[,v] = x1[,vv[1]]*x1[,vv[2]]
		x0[,v] = x0[,vv[1]]*x0[,vv[2]]
	}
}
z0 = z1 = z
if(is.element(paste("random.",vor,sep=""),colnames(z))){
	z1[,paste("random.",vor,sep="")] = 1; z0[,paste("random.",vor,sep="")] = 0
}

#initiate values of regression parameters
be0 = matrix(unlist(read.csv(paste(outdir,"/be_tau=",tau,".csv",sep=""),header=F)),ncol=1)
b0 = matrix(unlist(read.csv(paste(outdir,"/b_tau=",tau,".csv",sep=""),header=F)),ncol=dimb)
th0 = unlist(read.csv(paste(outdir,"/th_tau=",tau,".csv",sep=""),header=F))

mit <<- 1/2
fn0 = paste("tau=",tau,"_bootstrap",sep="")

for(ib in 1:iter.boot){
set.seed(ib)
	idi = idd[ib,]
	di = matrix(nrow=0,ncol=ncol(dat))
	yi = NULL; xi = xim0 = xim1 = matrix(nrow=0,ncol=ncol(x)); zi = zim0 = zim1 = matrix(nrow=0,ncol=ncol(z))
	for(i in idi){
		ri = id==i
		di = rbind(di,dat[ri,])
		yi = matrix(c(yi,y[ri]),ncol=1)
		xi = rbind(xi,matrix(x[ri,],ncol=ncol(x)))
		zi = rbind(zi,matrix(z[ri,],ncol=ncol(z)))
		xim0 = rbind(xim0,apply(x0[ri,],2,median)); xim1 = rbind(xim1,apply(x1[ri,],2,median))
		zim0 = rbind(zim0,apply(z0[ri,],2,median)); zim1 = rbind(zim1,apply(z1[ri,],2,median))
	}
	zij <<- function(i,j){ #m ni
		A = matrix(0,nrow=m*dimb,ncol=1)
		A[((i-1)*dimb+1):(i*dimb),1] = zi[rij(i,j),]
		return(A)
	}
	di$id = id

	be=be0; b=b0; th=th0
	ie = 1; dif = Inf; sdif = 10^-2
	while((!is.nan(dif) & dif>sdif) & ie<iter.step){
		be.pre=be; b.pre=b; th.pre=th
		A=try(beb.est.binom(tau,be,b,th,yi,xi,zi,dv),silent=T)
		if(class(A)!="try-error"){
			be=A$be; b=A$b
		}
		B=try(d.est.binom(tau,be,b,th,yi,xi,zi,dv),silent=T)
		if(class(B)!="try-error")
			th=B$th
		ie=ie+1; dif = abs(sum(c(be-be.pre,c(b-b.pre),th-th.pre)))
	}
	conv = F
	if(class(A)!="try-error" & class(B)!="try-error")
		conv = A$upd & B$upd & !is.nan(dif) & dif<sdif
	if(tau==0){
		ori = (xim1%*%be+apply(zim1*b,1,sum)) - (xim0%*%be+apply(zim0*b,1,sum))
	}else{
		ori = 2/tau*log(1/2*exp(tau*xim1%*%be)+1/2*exp(tau*apply(zim1*b,1,sum))) - 2/tau*log(1/2*exp(tau*xim0%*%be)+1/2*exp(tau*apply(zim0*b,1,sum)))
	}
	orr = rep(NA,m); bi = rep(NA,m*dimb)
	for(i in idi){
		ii = dat[which(id==which(i==idi)[1])[1],vc]
		orr[i] = ori[ii]
		for(j in 1:dimb)
			bi[i+(j-1)*length(uid)] = b[ii,j]
	}
	app = ifelse(ib==1,F,T)
	write.table(t(c(be)),paste(outdir,"/be_",fn0,".csv",sep=""),sep=",",row.name=F,col.name=F,append=app)
	write.table(t(c(bi)),paste(outdir,"/b_",fn0,".csv",sep=""),sep=",",row.name=F,col.name=F,append=app)
	write.table(t(c(th)),paste(outdir,"/th_",fn0,".csv",sep=""),sep=",",row.name=F,col.name=F,append=app)
	write.table(conv,paste(outdir,"/convergence_",fn0,".csv",sep=""),sep=",",row.name=F,col.name=F,append=app)
	write.table(t(orr),paste(outdir,"/log-odds-ratio_",fn0,".csv",sep=""),sep=",",row.name=F,col.name=F,append=app)
}

}



#-----------------------------------------------------------------------------------------------------
#   Parameter settings																				  
#-----------------------------------------------------------------------------------------------------

# Importing the data
dat = read.csv("data/Davis1991_respiratory_stack.csv",sep=",")
dat[,"center"] = as.numeric(dat[,"center"]==2)
dat[,"treat"] = as.numeric(dat[,"treat"]=="A")
dat[,"sex"] = as.numeric(dat[,"sex"]=="M")

# Response and cluster variable
vy = "status"; ydist = "binomial"
vc = "subject"

# Fixed effect variable
# Note: interaction term must be given as "variable1.variable2"
vf = c("center","treat","age","baseline","visit")
vf = c("center","treat","age","baseline","treat.baseline")

# Random effect variable
vr = "treat"

# Stratification variable
id = dat[,vc]
uid = unique(id)
rr = y2 = NULL
for(i in uid){
	rr = c(rr,which(dat[,vc]==i)[1])
	y2 = c(y2,sum(dat[dat[,vc]==i,vy]))
}
y2 = as.numeric(y2>=3)%x%rep(1,4)
age2 = as.numeric(dat[,"age"]>=median(dat[,"age"]))
dstr = data.frame(cbind(id,y2,dat[,c("baseline","treat")],age2))

# Value of tau parameter
tau = 0.7

# Variable for odds ratio
vor = "treat"

# Output directory including the estimates of regression parameters by the app_resp function
outdir = "result/201009152003_Davis1991_tau=0_interaction"

# Number of times of bootstrap replicates
iter.boot = 1

# Maximum number of the estimation steps 2 and 3
iter.step = 20

# Maximum number of the iteration in each estimation step
iter.cal = 100



#-----------------------------------------------------------------------------------------------------
#   Data analysis run																					
#-----------------------------------------------------------------------------------------------------

app_boot(dat,vy,vc,vf,vr,dstr,tau,vor,outdir,iter.boot,iter.step,iter.cal)


# calculate the median of the log-odds-ratios for treatment and the 95% confidence intervals
dir = "result/201009164348_Davis1991_no-interaction"
tau = 0.7
l0 = read.table(paste(dir,"/log-odds-ratio_tau=",tau,".csv",sep=""),sep=",",header=F)
cc = read.table(paste(dir,"/convergence_tau=",tau,"_bootstrap.csv",sep=""),sep=",",header=F)
ll = read.table(paste(dir,"/log-odds-ratio_tau=",tau,"_bootstrap.csv",sep=""),sep=",",header=F)[unlist(cc),]
id0 = NULL; for(i in unique(dat[,"subject"]))id0 = c(id0,which(dat[,"subject"]==i)[1])
base0 = dat[id0,"baseline"]==0; base1 = dat[id0,"baseline"]==1
pfun = function(vec)median(vec[!is.na(vec)])
median(unlist(l0)[base0]); quantile(apply(ll[,base0],1,pfun),c(0.05/2,1-0.05/2))
median(unlist(l0)[base1]); quantile(apply(ll[,base1],1,pfun),c(0.05/2,1-0.05/2))



