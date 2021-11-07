
#-----------------------------------------------------------------------------------------------------#
#   Production of the figures and tables   											                  #
#-----------------------------------------------------------------------------------------------------#

library(latex2exp)

#-----------------------------------------------------------------------------------------------------#
#   Function to draw the boxplot of the estimates of fixed and random effect parameters				  #
#-----------------------------------------------------------------------------------------------------#

### Parameters of function
  # dir		character	subfolder including the simulation results to draw the boxplot
  # tau		numeric		value of tau of the model to draw the boxplot

boxplot_estimate = function(dir,tau){
	P = read.table(paste(dir,"/true.param.csv",sep=""),sep=",",nrows=1,header=T)
	tau.true = P[1,"tau"]
	be.true = P[grep("be",colnames(P))]
	b.true = P[(max(grep("be",colnames(P)))+1):(max(grep("b",colnames(P)))-3)]
	b.vary = P[1,"b.vary"]

	minf = function(v){
		if(sum(!is.na(v))>0){
			return(which(v==min(v,na.rm=T))[1])
		}else{
			return(NA)
		}
	}
	if(b.vary){
		A = grep("tau=",list.files(dir)); di = list.files(dir)[A]
		Taui = NULL
		for(j in 1:length(di)){
			st = strsplit(strsplit(di[j],"tau=")[[1]][2],"_")[[1]][1]
			if(length(grep("slope",st))==1)st = strsplit(st,"slope")[[1]][1]
			Taui = c(Taui,as.numeric(st))
		}
		Taui = sort(unique(Taui))
	}else{
		Taui = tau
	}
	if(tau==0)Taui = 0
	n.sim = nrow(read.table(paste(dir,"/tau=",0,"_cAIC.csv",sep=""),sep=","))
	ll=10^-10; ul = 10^10
	Caic = matrix(nrow=n.sim,ncol=0)
	for(taui in Taui){
		caic = read.table(paste(dir,"/tau=",taui,"_cAIC.csv",sep=""),sep=",")
		co = read.table(paste(dir,"/tau=",taui,"_convergence.csv",sep=""),sep=",")
		th = read.table(paste(dir,"/tau=",taui,"_th.csv",sep=""),sep=",")
		dimth = ncol(th)
		if(dimth==1){
			caic[!(co&ll<th&th<ul)] = NA
		}else{
			caic[!(co&ll<apply(th[,-2],1,min)&apply(th[,-2],1,max)<ul)] = NA
		}
		Caic = cbind(Caic,caic)
	}
	tau.best = Taui[apply(Caic,1,minf)]
	
	belist = blist = list()
	for(taui in Taui){
		fn = paste(dir,"/tau=",taui,sep="")
		belist[[which(taui==Taui)]] = read.table(paste(fn,"_be.csv",sep=""),sep=",")
		blist[[which(taui==Taui)]] = read.table(paste(fn,"_b.csv",sep=""),sep=",")
	}
	BE = matrix(nrow=0,ncol=ncol(belist[[1]]))
	BB = matrix(nrow=0,ncol=ncol(blist[[1]]))
	for(it in 1:n.sim){
		taui = tau.best[it]
		if(!is.na(taui)){
			BE = rbind(BE,belist[[which(taui==Taui)]][it,])
			BB = rbind(BB,blist[[which(taui==Taui)]][it,])
		}
	}
	if(!b.vary){
		EE = cbind(BE,BB)
		ylim = c(-5,5); cex.axis = 3
		#par(mar=c(6,4,4,1))
		boxplot(EE,ylim=ylim,xaxt="n",las=2,cex.axis=cex.axis,pch=20)
		axis(1,at=1:ncol(EE),tick=F,line=1,labels=c(TeX(paste("$\\beta_",0:(ncol(BE)-1),"$",sep="")),TeX(paste("$\\mathit{b}_{",1:ncol(BB),"}$",sep=""))),cex.axis=cex.axis)
		points(1:length(c(be.true,b.true)),c(be.true,b.true),col="red",pch=20,cex=3)
	}else{
		EE = BE
		ylim = c(-5,5); cex.axis = 3
		#par(mar=c(6,4,4,1))
		boxplot(EE,ylim=ylim,xaxt="n",las=2,cex.axis=cex.axis,pch=20)
		axis(1,at=1:ncol(EE),tick=F,line=1,labels=c(TeX(paste("$\\beta_",0:(ncol(BE)-1),"$",sep=""))),cex.axis=cex.axis)
		points(1:length(c(be.true)),c(be.true),col="red",pch=20,cex=3)
	}
}

#-----------------------------------------------------------------------------------------------------#
#   Draw boxplot																					  #
#-----------------------------------------------------------------------------------------------------#

boxplot_estimate("result/210912230458_binom_tau.true=-1_b.dist=norm_n.sim=1000",tau=0)
boxplot_estimate("result/210912230458_binom_tau.true=-1_b.dist=norm_n.sim=1000",tau=-1)



#-----------------------------------------------------------------------------------------------------#
#   Function to make the table of median of the conditional AIC										  #
#-----------------------------------------------------------------------------------------------------#

### Parameters of function
  # dd		vector 		subfolders including the simulation results to make the table
  # Tau	    vector      candidate values of the parameter tau

summ_sim = function(dd,Tau){
	is.misspe = length(grep("misspe",dd))>0
	tabcaic = matrix(nrow=0,ncol=3+ifelse(is.misspe,2,0))
	tabconv = tabtau = tabconvr = tabtaur = matrix(nrow=0,ncol=length(Tau))
	colnames(tabconv) = colnames(tabtau) = colnames(tabconvr) = colnames(tabtaur) = Tau
	minf = function(v){
		if(sum(!is.na(v))>0){
			return(which(v==min(v,na.rm=T))[1])
		}else{
			return(NA)
		}
	}
	for(idir in 1:length(dd)){
		A = grep("tau=",list.files(dd[idir])); di = list.files(dd[idir])[A]
		Taui = NULL
		for(j in 1:length(di)){
			st = strsplit(strsplit(di[j],"tau=")[[1]][2],"_")[[1]][1]
			if(length(grep("slope",st))==1)st = strsplit(st,"slope")[[1]][1]
			Taui = c(Taui,as.numeric(st))
		}
		Taui = sort(unique(Taui))

		n.sim = nrow(read.table(paste(dd[1],"/tau=",Taui[1],"_cAIC.csv",sep=""),sep=","))
		ll=10^-10; ul = 10^10
		Caic = matrix(nrow=n.sim,ncol=0)
		for(tau in Taui){
			caic = read.table(paste(dd[idir],"/tau=",tau,"_cAIC.csv",sep=""),sep=",")
			co = read.table(paste(dd[idir],"/tau=",tau,"_convergence.csv",sep=""),sep=",")
			th = read.table(paste(dd[idir],"/tau=",tau,"_th.csv",sep=""),sep=",")
			dimth = ncol(th)
			if(dimth==1){
				caic[!(co&ll<th&th<ul)] = NA
				#caic[!(co)] = NA
			}else{
				caic[!(co&ll<apply(th[,-2],1,min)&apply(th[,-2],1,max)<ul)] = NA
				#caic[!(co)] = NA
			}
if(tau==Taui[1])nss = nrow(caic)
caic2 = rep(NA,n.sim); caic2[1:nrow(caic)] = unlist(caic); caic = caic2
			Caic = cbind(Caic,caic)
		}
		tau.best = Taui[apply(Caic,1,minf)]

		#Proportions of convergence
		cvg = tauest = rep("-",length(Tau))
		#for(ti in Taui)cvg[ti==Tau] = round(sum(!is.na(Caic[,Taui==ti]))/nrow(Caic),2)
for(ti in Taui)cvg[ti==Tau] = round(sum(!is.na(Caic[,Taui==ti]))/nss,2)
		tabconv = rbind(tabconv,cvg)

		#Proportions of estimated parameter of tau
		for(ti in Taui)tauest[ti==Tau] = round(sum(tau.best==ti,na.rm=T)/sum(!is.na(tau.best)),2)
		tabtau = rbind(tabtau,tauest)

		if(is.misspe){
			Caicr = matrix(nrow=n.sim,ncol=0)
			for(tau in Tau){
				caic = read.table(paste(dd[idir],"/tau=",tau,"slope_cAIC.csv",sep=""),sep=",")
				co = read.table(paste(dd[idir],"/tau=",tau,"slope_convergence.csv",sep=""),sep=",")
				th = read.table(paste(dd[idir],"/tau=",tau,"slope_th.csv",sep=""),sep=",")
				dimth = ncol(th)
				if(dimth==1){
					caic[!(co&ll<th&th<ul)] = NA
				}else{
					caic[!(co&ll<apply(th[,-2],1,min)&apply(th[,-2],1,max)<ul)] = NA
				}
				Caicr = cbind(Caicr,caic)
			}
			tau.bestr = Tau[apply(Caicr,1,minf)]

			#Proportions of convergence
			cvg = tauest = rep("-",length(Tau))
			for(ti in Taui)cvg[ti==Tau] = round(sum(!is.na(Caicr[,Taui==ti]))/nrow(Caicr),2)
			tabconvr = rbind(tabconvr,cvg)

			#Proportions of estimated parameter of tau
			for(ti in Taui)tauest[ti==Tau] = round(sum(tau.bestr==ti,na.rm=T)/sum(!is.na(tau.bestr)),2)
			tabtaur = rbind(tabtaur,tauest)
		}
		ci0 = Caic[,which(Taui==0)]; ciq = NULL
		for(j in 1:nrow(Caic))
			ciq = c(ciq,ifelse(is.na(tau.best[j]),NA,Caic[j,which(Taui==tau.best[j])]))
		mat.caic = cbind(ci0,ciq)
		if(is.misspe){
			ci0r = Caicr[,which(Taui==0)]; ciqr = NULL
			for(j in 1:nrow(Caicr))
				ciqr = c(ciqr,ifelse(is.na(tau.bestr[j]),NA,Caicr[j,which(Taui==tau.bestr[j])]))
			mat.caic = cbind(ci0,ciq,ci0r,ciqr)
		}
		median.na = function(x)median(x,na.rm=T)
		tabcaic = rbind(tabcaic,c(dd[idir],round(apply(mat.caic,2,median.na),1)))
	}
	if(is.misspe){
		colnames(tabcaic) = c("dir","GLMM","GQLMM","GLMMslope","GQLMMslope")
		return(list(tabcaic=tabcaic,tabconv=tabconv,tabtau=tabtau,tabconvr=tabconvr,tabtaur=tabtaur))
	}else{
		colnames(tabcaic) = c("dir","GLMM","GQLMM")
		return(list(tabcaic=tabcaic,tabconv=tabconv,tabtau=tabtau))
	}
}

#-----------------------------------------------------------------------------------------------------#
#   Make table																						  #
#-----------------------------------------------------------------------------------------------------#

library(Hmisc)

f = function(tab,dig=2){
	for(i in 1:nrow(tab))for(j in 1:ncol(tab)){
		if(tab[i,j]=="-"){
		}else{
			A = strsplit(tab[i,j],"[.]")[[1]]
			if(length(A)==1){
				tab[i,j] = paste(tab[i,j],".",sep="")
				for(k in 1:dig)tab[i,j] = paste(tab[i,j],"0",sep="")
			}
			if(length(A)==2){
				dij = A[2]
				if(nchar(dij)<dig){
					for(k in 1:(dig-nchar(dij)))
						tab[i,j] = paste(tab[i,j],"0",sep="")
				}
			}
		}
	}
	return(tab)
}
g = function(outsumm,set){
	caic = cbind(set,"",f(outsumm$tabcaic[,-1],1)); rownames(caic)=NULL
	latex(caic,file="tabcaic.tex")
	conv = cbind(set,"",f(outsumm$tabconv)); rownames(conv)=NULL
	latex(conv,file="tabconv.tex")
	if(!is.null(outsumm$tabconvr)){
		conv = cbind(set,"",f(outsumm$tabconvr)); rownames(conv)=NULL
		latex(conv,file="tabconvr.tex")
	}
	tauest = cbind(set,"",f(outsumm$tabtau)); rownames(tauest)=NULL
	latex(tauest,file="tabtau.tex")
	if(!is.null(outsumm$tabtaur)){
		tauest = cbind(set,"",f(outsumm$tabtaur)); rownames(tauest)=NULL
		latex(tauest,file="tabtaur.tex")
	}
}

# Candidate values of tau
A = c(0.1,0.2,0.4,0.7,1,3)
Tau = sort(c(-A,0,A))

summ_sim(dd="result/210912230458_binom_tau.true=-1_b.dist=norm_n.sim=1000",Tau=Tau)



