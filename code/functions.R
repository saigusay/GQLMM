

bv = function(mat){
	A = matrix(nrow=0,ncol=1)
	for(i in 1:nrow(mat))A = rbind(A,matrix(mat[i,],ncol=1))
	return(A)
}

diagm = function(mat,rep){
	dm = dim(mat)[1]
	A = matrix(0,nrow=dm*rep,ncol=dm*rep)
	for(i in 1:rep)
	A[(i-1)*dm+1:dm,(i-1)*dm+1:dm] = mat
	return(A)
}

sqmt = function(mat){
	A = svd(mat)
	return(A$u %*% diag(sqrt(A$d)) %*% t(A$v) )
}

logit = function(p)log(p/(1-p))

logit.inv = function(x)1/(1+exp(-x))

prd.bycol = function(mat1,mat2){
	if(nrow(mat1)!=nrow(mat2) | ncol(mat1)!=1)return(0)
	AA = matrix(rep(mat1,ncol(mat2)),ncol=ncol(mat2))
	return(AA*mat2)
}

sum.g = function(mat){
	A = matrix(nrow=m,ncol=ncol(mat))
	for(jm in 1:ncol(mat))for(im in 1:m)
	A[im,jm] = sum(mat[(sum(ni[0:(im-1)])+1):sum(ni[1:im]),jm])
	return(A)
}

bm = function(vec){ #m, dimb
	A = matrix(nrow=m,ncol=dimb)
	for(i in 1:m)
	A[i,] = vec[((i-1)*dimb+1):(i*dimb)]
	return(A)
}

bb = function(b){ #ni
	A = matrix(nrow=0,ncol=ncol(b))
	for(i in 1:m)
	A = rbind(A,matrix(b[i,],nrow=1)%x%rep(1,ni[i]))
	return(A)
}

rij = function(i,j) # ni
	ifelse(i==1,j,sum(ni[1:(i-1)])+j)

ka.binom = function(y,pi,b,D){
	bv = bv(b)
	return(- sum(y*log(pi)+(1-y)*log(1-pi)) + 1/2*t(bv)%*%ginv(D)%*%bv)
}

glmer.binom = function(y,x,z,id){
	dat = data.frame(cbind(id,y))
	colnames(dat) = c("id","y")
	model = "y~"
	if(ncol(x)>1){
	sc = ncol(dat)
	dat = data.frame(cbind(dat,x[,2:ncol(x)]))
	colnames(dat)[(sc+1):ncol(dat)] = paste("x",1:(ncol(x)-1),sep="")
	for(i in 2:ncol(x))
	model = paste(model,"+x",i-1,sep="")
	}
	model = paste(model,"+(1",sep="")
	if(ncol(z)>1){
	sc = ncol(dat)
	dat = data.frame(cbind(dat,z[,2:ncol(z)]))
	colnames(dat)[(sc+1):ncol(dat)] = paste("z",1:(ncol(z)-1),sep="")
	for(i in 2:ncol(z))
	model = paste(model,"+z",i-1,sep="")
	}
	model = paste(model,"|id)",sep="")

	A = glmer(as.formula(model),data=dat,family="binomial")
	be.glm = matrix(summary(A)$coef[,1],ncol=1)
	b.glm = matrix(unlist(ranef(A)),ncol=ncol(z))
	th.glm = matrix(attr(summary(A)$varcor$id,"stddev")^2,ncol=1)
	co.glm = attr(summary(A)$varcor$id,"correlation")
	if(dim(co.glm)[1]==2)co.glm = co.glm[1,2]
	return(list(be.glm=be.glm,b.glm=b.glm,th.glm=th.glm,co.glm=co.glm))
}

beb.est.binom = function(tau,be,b,th,y,x,z,dv){
	dimb = ncol(b)
	di = matrix(0,nrow=dimb,ncol=dimb)
	for(i in 1:length(th))
		di = di + dv[[i]]*th[i]
	D = diagm(di,m)
	be.pre=be; b.pre=b
	kk = NULL; dif=Inf; it2=0; ka0=NaN
	while((!is.nan(dif) & dif>10^-5) & it2<iter.cal){
		it2 = it2+1; upd = F
		if(tau==0){
			rr0 = apply(z*bb(b),1,sum)
			pi = 1/(1+ exp(-(x%*%be+rr0)) )
			q = matrix(1/2,nrow=nrow(pi),ncol=1)
		}else{
			ff = exp(tau*x%*%be)
			rr = exp(tau*apply(z*bb(b),1,sum))
			pi = 1/(1+ (1/2*ff+1/2*rr)^(-2/tau) )
			q = ff/(ff+rr)
		}
		g1 = - 2*t(x)%*%(q*(y-pi))
		h11_1 = 2*t(x)%*%(prd.bycol(q*(2*q*pi*(1-pi)-tau*(1-q)*(y-pi)),x))
		h11_2 = 4*t(x)%*%(prd.bycol(q^2*pi*(1-pi),x))
		h11_3 = diag(diag(h11_2))
		be1 = be2 = be3 = matrix(NaN,nrow=length(be))
		A=try(ginv(h11_1),silent=T); if(class(A)[1]!="try-error")be1 = be - A%*%g1 *mit
		A=try(ginv(h11_2),silent=T); if(class(A)[1]!="try-error")be2 = be - A%*%g1 *mit
		A=try(ginv(h11_3),silent=T); if(class(A)[1]!="try-error")be3 = be - A%*%g1 *mit

		if(tau==0){
			pi1 = 1/(1+ exp(-(x%*%be1+rr0)) )
			pi2 = 1/(1+ exp(-(x%*%be2+rr0)) )
			pi3 = 1/(1+ exp(-(x%*%be3+rr0)) )
		}else{
			pi1 = 1/(1+(1/2*exp(tau*x%*%be1)+1/2*rr)^(-2/tau))
			pi2 = 1/(1+(1/2*exp(tau*x%*%be2)+1/2*rr)^(-2/tau))
			pi3 = 1/(1+(1/2*exp(tau*x%*%be3)+1/2*rr)^(-2/tau))
		}
		ka1 = ka.binom(y,pi1,b,D)
		ka2 = ka.binom(y,pi2,b,D)
		ka3 = ka.binom(y,pi3,b,D)
		if(!is.nan(ka1)){
			be = be1; ka0 = ka1; upd = T
		}
		if(!is.nan(ka1)&!is.nan(ka2)&!is.nan(ka3)&(ka2<=ka1&ka2<=ka3)){
			be = be2; ka0 = ka2; upd = T
		}else if(!is.nan(ka1)&!is.nan(ka2)&!is.nan(ka3)&(ka3<=ka1&ka3<=ka2)){
			be = be3; ka0 = ka3; upd = T
		}
		kk = c(kk,ka0)

		if(tau==0){
			ff0 = x%*%be
			pi = 1/(1+ exp(-(ff0+rr0)) )
			q = matrix(1/2,nrow=nrow(pi),ncol=1)
		}else{
			ff = exp(tau*x%*%be)
			pi = 1/(1+(1/2*ff+1/2*rr)^(-2/tau))
			q = ff/(ff+rr)
		}
		A = matrix(nrow=nrow(b),ncol=0)
		for(i in 1:ncol(b))
			A = cbind(A,matrix(sum.g((1-q)*(y-pi)*z[,i]),ncol=1))
		g2 = - 2*bv(A) + ginv(D)%*%bv(b)
		h22_1 = h22_2 = ginv(D) #matrix(0,nrow=length(b),ncol=length(b))
		sa = 2*(1-q)^2*pi*(1-pi)+tau*q*(1-q)*(y-pi); sb = (1-q)^2*pi*(1-pi)
		for(i in 1:m)for(j in 1:ni[i]){
			h22_1 = h22_1 + 2*sa[rij(i,j)]*zij(i,j)%*%t(zij(i,j))
			h22_2 = h22_2 + 4*sb[rij(i,j)]*zij(i,j)%*%t(zij(i,j))
		}
		h22_3 = diag(diag(h22_2))
		b1 = b2 = b3 = matrix(NaN,nrow=nrow(b),ncol=ncol(b))
		A=try(ginv(h22_1),silent=T); if(class(A)[1]!="try-error")b1 = b - bm(A%*%g2) *mit
		A=try(ginv(h22_2),silent=T); if(class(A)[1]!="try-error")b2 = b - bm(A%*%g2) *mit
		A=try(ginv(h22_3),silent=T); if(class(A)[1]!="try-error")b3 = b - bm(A%*%g2) *mit

		if(tau==0){
			pi1 = 1/(1+ exp(-(ff0+apply(z*bb(b1),1,sum))) )
			pi2 = 1/(1+ exp(-(ff0+apply(z*bb(b2),1,sum))) )
			pi3 = 1/(1+ exp(-(ff0+apply(z*bb(b3),1,sum))) )
		}else{
			pi1 = 1/(1+(1/2*ff+1/2*exp(tau*apply(z*bb(b1),1,sum)))^(-2/tau))
			pi2 = 1/(1+(1/2*ff+1/2*exp(tau*apply(z*bb(b2),1,sum)))^(-2/tau))
			pi3 = 1/(1+(1/2*ff+1/2*exp(tau*apply(z*bb(b3),1,sum)))^(-2/tau))
		}
		ka1 = ka.binom(y,pi1,b1,D)
		ka2 = ka.binom(y,pi2,b2,D)
		ka3 = ka.binom(y,pi3,b3,D)
		if(!is.nan(ka1)){
			b = b1; ka0 = ka1; upd = T
		}
		if(!is.nan(ka1)&!is.nan(ka2)&!is.nan(ka3)&(ka2<=ka1&ka2<=ka3)){
			b = b2; ka0 = ka2; upd = T
		}else if(!is.nan(ka1)&!is.nan(ka2)&!is.nan(ka3)&(ka3<=ka1&ka3<=ka2)){
			b = b3; ka0 = ka3; upd = T
		}
		kk = c(kk,ka0)

		#be.i = rbind(be.i,c(be)); b.i = rbind(b.i,c(b))
		dif = sum(abs(c(be-be.pre,c(b-b.pre))))
		be.pre=be; b.pre=b
	} #while
	return(list(be=be,b=b,upd=upd))
}

d.est.binom = function(tau,be,b,th,y,x,z,dv){ #m, dimb, ni
	th.pre=th; dif=Inf
	it2=0
	TH = NULL
	while((!is.nan(dif) & dif>10^-5) & it2<iter.cal){
		it2 = it2+1; upd = T
		#if(it2==iter.cal)print("count.stop: d.est")
		if(tau==0){
			pi = 1/(1+ exp(-(x%*%be + apply(z*bb(b),1,sum))) )
			q = matrix(1/2,nrow=nrow(pi),ncol=1)
		}else{
			ff = exp(tau*x%*%be)
			rr = exp(tau*apply(z*bb(b),1,sum))
			pi = 1/(1+(1/2*ff+1/2*rr)^(-2/tau))
			q = ff/(ff+rr)
		}
		U = 2*(1-q)*apply(z*bb(b),1,sum) + (y-pi)/(pi*(1-pi)) #Y-2QXbe
		W = diag(c(pi*(1-pi)))
		Z = matrix(rbind(matrix(z[1:ni[1],],ncol=dimb),matrix(0,nrow=sum(ni)-ni[1],ncol=dimb)),ncol=dimb)
		for(i in 2:m)
			Z = cbind(Z,matrix(rbind(matrix(0,nrow=sum(ni[1:(i-1)]),ncol=dimb),matrix(z[(sum(ni[1:(i-1)])+1):sum(ni[1:i]),],ncol=dimb),matrix(0,nrow=n-sum(ni[1:i]),ncol=dimb)),ncol=dimb))

		di = matrix(0,nrow=dimb,ncol=dimb)
		for(i in 1:length(th))
			di = di + dv[[i]]*th[i]
		D = diagm(di,m)
		V = ginv(W) + 4*diag(c(1-q))%*%Z%*%D%*%t(diag(c(1-q))%*%Z); Vinv = ginv(V)
		P = Vinv-Vinv%*%diag(c(q))%*%x%*%ginv(t(diag(c(q))%*%x)%*%Vinv%*%diag(c(q))%*%x)%*%t(diag(c(q))%*%x)%*%Vinv
		dimt = length(th)
		g = matrix(NA,nrow=dimt,ncol=1); h = matrix(NA,dimt,dimt)
		dV =list()
		for(i in 1:length(th))
			dV[[i]] = 4*diag(c(1-q))%*%Z%*%diagm(dv[[i]],m)%*%t(diag(c(1-q))%*%Z)
		for(i in 1:dimt)
			g[i] = 1/2*( t(U)%*%Vinv%*%dV[[i]]%*%Vinv%*%U - sum(diag(P%*%dV[[i]])) )
		for(i in 1:dimt)for(j in 1:dimt)
			h[i,j] = -1/2*sum(diag(P%*%dV[[i]]%*%P%*%dV[[j]]))
		if(dimt==1){
			th = th - h^-1*g * mit
		}else{
			th = th - ginv(h)%*%g * mit
		}
		th[th>10^10] = 10^10
		if(dimb==1)th[th< 0] = 10^-10
		if(dimb==2)th[(th< 0)&c(T,F,T)] = 10^-10
		if(is.element(NaN,th)){
			th[is.nan(th)] = 10^10; upd = F
		}
		dif = sum(abs(th - th.pre))
		th.pre = th
	} #while
	D = diag(rep(th,m))
	return(list(th=th,upd=upd))
}

caic.binom = function(tau,be,b,th,y,x,z,dv){
	if(tau==0){
		pi = 1/(1+ exp(-(x%*%be+apply(z*bb(b),1,sum))) )
		q = matrix(1/2,nrow=nrow(pi),ncol=1)
	}else{
		ff = exp(tau*x%*%be)
		rr = exp(tau*apply(z*bb(b),1,sum))
		pi = 1/(1+(1/2*ff+1/2*rr)^(-2/tau))
		q = ff/(ff+rr)
	}
	Q = diag(c(q))
	iqz = 2*z*((1-q)%x%matrix(1,nrow=1,ncol=dimb))
	QZ = matrix(rbind(matrix(iqz[1:ni[1],],ncol=dimb),matrix(0,nrow=sum(ni)-ni[1],ncol=dimb)),ncol=dimb)
	for(i in 2:m)
		QZ = cbind(QZ,matrix(rbind(matrix(0,nrow=sum(ni[1:(i-1)]),ncol=dimb),matrix(iqz[(sum(ni[1:(i-1)])+1):sum(ni[1:i]),],ncol=dimb),matrix(0,nrow=n-sum(ni[1:i]),ncol=dimb)),ncol=dimb))
	M = as.matrix(cbind(2*Q%*%x,QZ))
	W = diag(c(pi*(1-pi)))
	di = matrix(0,nrow=dimb,ncol=dimb)
	for(i in 1:length(th))
		di = di + dv[[i]]*th[i]
	D = diagm(di,m)
	p = length(be)
	A = t(M)%*%W%*%M
	A[(p+1):(p+m*dimb),(p+1):(p+m*dimb)] = A[(p+1):(p+m*dimb),(p+1):(p+m*dimb)] + ginv(D)
	Ainv = try(ginv(A),silent=T)
	ef.df =  ifelse(class(Ainv)[1]!="try-error", sum(diag(t(M)%*%W%*%M%*%Ainv)), NA)
	caic = -2*(sum(y*log(pi)+(1-y)*log(1-pi))) + 2*ef.df
	return(list(caic=caic,ef.df=ef.df,L=sum(y*log(pi)+(1-y)*log(1-pi))))
}


ka.pois = function(y,lam,b,D){
	bv = bv(b)
	return(- sum(y*log(lam)-lam) + 1/2*t(bv)%*%ginv(D)%*%bv)
}

glmer.pois = function(y,x,z,id){
	dat = data.frame(cbind(id,y))
	colnames(dat) = c("id","y")
	model = "y~"
	if(ncol(x)>1){
	sc = ncol(dat)
	dat = data.frame(cbind(dat,x[,2:ncol(x)]))
	colnames(dat)[(sc+1):ncol(dat)] = paste("x",1:(ncol(x)-1),sep="")
	for(i in 2:ncol(x))
	model = paste(model,"+x",i-1,sep="")
	}
	model = paste(model,"+(1",sep="")
	if(ncol(z)>1){
	sc = ncol(dat)
	dat = data.frame(cbind(dat,z[,2:ncol(z)]))
	colnames(dat)[(sc+1):ncol(dat)] = paste("z",1:(ncol(z)-1),sep="")
	for(i in 2:ncol(z))
	model = paste(model,"+z",i-1,sep="")
	}
	model = paste(model,"|id)",sep="")

	A = glmer(as.formula(model),data=dat,family="poisson")
	be.glm = matrix(summary(A)$coef[,1],ncol=1)
	b.glm = matrix(unlist(ranef(A)),ncol=ncol(z))
	th.glm = matrix(attr(summary(A)$varcor$id,"stddev")^2,ncol=1)
	co.glm = attr(summary(A)$varcor$id,"correlation")
	if(dim(co.glm)[1]==2)co.glm = co.glm[1,2]
	return(list(be.glm=be.glm,b.glm=b.glm,th.glm=th.glm,co.glm=co.glm))
}

beb.est.pois = function(tau,be,b,th,y,x,z,dv){
di = matrix(0,nrow=dimb,ncol=dimb)
for(i in 1:length(th))
	di = di + dv[[i]]*th[i]
D = diagm(di,m)
be.pre=be; b.pre=b
dif=Inf; it2=0; kk = NULL; ka0=NaN
while((!is.nan(dif) & dif>10^-5) & it2<iter.cal){
	it2 = it2+1; upd = F

	if(tau==0){
		rr0 = apply(z*bb(b),1,sum)
		lam = exp(x%*%be + rr0)
		q = matrix(1/2,nrow=nrow(lam),ncol=1)
	}else{
		ff = exp(tau*x%*%be)
		rr = exp(tau*apply(z*bb(b),1,sum))
		lam = (1/2*ff+1/2*rr)^(2/tau)
		q = ff/(ff+rr)
	}
	g1 = - 2*t(x)%*%(q*(y-lam))
	h11_1 = 2*t(x)%*%(prd.bycol(lam*q*(2*q+tau-tau*q),x))
	h11_2 = 4*t(x)%*%(prd.bycol(lam*q^2,x))
	h11_3 = diag(diag(h11_2))
	be1 = be2 = be3 = matrix(NaN,nrow=length(be))
	A=try(ginv(h11_1),silent=T); if(class(A)[1]!="try-error")be1 = be - A%*%g1 *mit
	A=try(ginv(h11_2),silent=T); if(class(A)[1]!="try-error")be2 = be - A%*%g1 *mit
	A=try(ginv(h11_3),silent=T); if(class(A)[1]!="try-error")be3 = be - A%*%g1 *mit

	if(tau==0){
		lam1 = exp(x%*%be1 + rr0)
		lam2 = exp(x%*%be2 + rr0)
		lam3 = exp(x%*%be3 + rr0)
	}else{
		lam1 = (1/2*exp(tau*x%*%be1)+1/2*rr)^(2/tau)
		lam2 = (1/2*exp(tau*x%*%be2)+1/2*rr)^(2/tau)
		lam3 = (1/2*exp(tau*x%*%be3)+1/2*rr)^(2/tau)
	}
	ka1 = ka.pois(y,lam1,b,D)
	ka2 = ka.pois(y,lam2,b,D)
	ka3 = ka.pois(y,lam3,b,D)
	if(!is.nan(ka1)){
	be = be1; ka0 = ka1; upd = T
	}
	if(!is.nan(ka1)&!is.nan(ka2)&!is.nan(ka3)&(ka2<=ka1&ka2<=ka3)){
		be = be2; ka0 = ka2; upd = T
	}else if(!is.nan(ka1)&!is.nan(ka2)&!is.nan(ka3)&(ka3<=ka1&ka3<=ka2)){
		be = be3; ka0 = ka3; upd = T
	}

	if(tau==0){
		ff0 = x%*%be
		lam = exp(ff0 + rr0)
		q = matrix(1/2,nrow=nrow(lam),ncol=1)
	}else{
		ff = exp(tau*x%*%be)
		lam = (1/2*ff+1/2*rr)^(2/tau)
		q = ff/(ff+rr)
	}
	A = matrix(nrow=nrow(b),ncol=0)
	for(i in 1:ncol(b))
		A = cbind(A,matrix(sum.g((1-q)*(y-lam)*z[,i]),ncol=1))
	g2 = - 2*bv(A) + ginv(D)%*%bv(b)
	h22_1 = h22_2 = ginv(D) #matrix(0,nrow=length(b),ncol=length(b))
	sa = lam*(1-q)*(2-2*q+tau*q); sb = lam*(1-q)^2
	for(i in 1:m)for(j in 1:ni[i]){
		h22_1 = h22_1 + 2*sa[rij(i,j)]*zij(i,j)%*%t(zij(i,j))
		h22_2 = h22_2 + 4*sb[rij(i,j)]*zij(i,j)%*%t(zij(i,j))
	}
	h22_3 = diag(diag(h22_2))
	b1 = b2 = b3 = matrix(NaN,nrow=nrow(b),ncol=ncol(b))
	A=try(ginv(h22_1),silent=T); if(class(A)[1]!="try-error")b1 = b - bm(A%*%g2) *mit
	A=try(ginv(h22_2),silent=T); if(class(A)[1]!="try-error")b2 = b - bm(A%*%g2) *mit
	A=try(ginv(h22_3),silent=T); if(class(A)[1]!="try-error")b3 = b - bm(A%*%g2) *mit

	if(tau==0){
		lam1 = exp(ff0+apply(z*bb(b1),1,sum))
		lam2 = exp(ff0+apply(z*bb(b2),1,sum))
		lam3 = exp(ff0+apply(z*bb(b3),1,sum))
	}else{
		lam1 = (1/2*ff+1/2*exp(tau*apply(z*bb(b1),1,sum)))^(2/tau)
		lam2 = (1/2*ff+1/2*exp(tau*apply(z*bb(b2),1,sum)))^(2/tau)
		lam3 = (1/2*ff+1/2*exp(tau*apply(z*bb(b3),1,sum)))^(2/tau)
	}
	ka1 = ka.pois(y,lam1,b1,D)
	ka2 = ka.pois(y,lam2,b2,D)
	ka3 = ka.pois(y,lam3,b3,D)
	if(!is.nan(ka1)){
		b = b1; ka0 = ka1; upd = T
	}
	if(!is.nan(ka1)&!is.nan(ka2)&!is.nan(ka3)&(ka2<=ka1&ka2<=ka3)){
		b = b2; ka0 = ka2; upd = T
	}else if(!is.nan(ka1)&!is.nan(ka2)&!is.nan(ka3)&(ka3<=ka1&ka3<=ka2)){
		b = b3; ka0 = ka3; upd = T
	}

	dif = sum(abs(c(be-be.pre,c(b-b.pre))))
	be.pre=be; b.pre=b
} #while
return(list(be=be,b=b,upd=upd))
}

d.est.pois = function(tau,be,b,th,y,x,z,dv){
th.pre=th; dif=Inf
it2=0
TH = NULL
while((!is.nan(dif) & dif>10^-5) & it2<iter.cal){
	it2 = it2+1; upd = T
	if(tau==0){
		lam = exp(x%*%be + apply(z*bb(b),1,sum))
		q = matrix(1/2,nrow=nrow(lam),ncol=1)
	}else{
		ff = exp(tau*x%*%be)
		rr = exp(tau*apply(z*bb(b),1,sum))
		lam = (1/2*ff+1/2*rr)^(2/tau)
		q = ff/(ff+rr)
	}
	U = 2*(1-q)*apply(z*bb(b),1,sum) + (y-lam)/lam
	W = diag(c(lam))
	Z = matrix(rbind(matrix(z[1:ni[1],],ncol=dimb),matrix(0,nrow=sum(ni)-ni[1],ncol=dimb)),ncol=dimb)
	for(i in 2:m)
		Z = cbind(Z,matrix(rbind(matrix(0,nrow=sum(ni[1:(i-1)]),ncol=dimb),matrix(z[(sum(ni[1:(i-1)])+1):sum(ni[1:i]),],ncol=dimb),matrix(0,nrow=n-sum(ni[1:i]),ncol=dimb)),ncol=dimb))

	di = matrix(0,nrow=dimb,ncol=dimb)
	for(i in 1:length(th))
		di = di + dv[[i]]*th[i]
	D = diagm(di,m)
	V = ginv(W) + 4*diag(c(1-q))%*%Z%*%D%*%t(diag(c(1-q))%*%Z); Vinv = ginv(V)
	P = Vinv-Vinv%*%diag(c(q))%*%x%*%ginv(t(diag(c(q))%*%x)%*%Vinv%*%diag(c(q))%*%x)%*%t(diag(c(q))%*%x)%*%Vinv
	dimt = length(th)
	g = matrix(NA,nrow=dimt,ncol=1); h = matrix(NA,dimt,dimt)
	dV =list()
	for(i in 1:length(th))
		dV[[i]] = 4*diag(c(1-q))%*%Z%*%diagm(dv[[i]],m)%*%t(diag(c(1-q))%*%Z)
	for(i in 1:dimt)
		g[i] = 1/2*( t(U)%*%Vinv%*%dV[[i]]%*%Vinv%*%U - sum(diag(P%*%dV[[i]])) )
	for(i in 1:dimt)for(j in 1:dimt)
		h[i,j] = -1/2*sum(diag(P%*%dV[[i]]%*%P%*%dV[[j]]))
	if(dimt==1){
		th = th - h^-1*g * mit
	}else{
		th = th - ginv(h)%*%g * mit
	}
	th[th>10^10] = 10^10
	if(dimb==1)th[th< 0] = 10^-10
	if(dimb==2)th[(th< 0)&c(T,F,T)] = 10^-10
	if(is.element(NaN,th)){
		th[is.nan(th)] = 10^10; upd = F
	}
	dif = abs(sum(th - th.pre))
	th.pre = th
} #while
D = diag(rep(th,m))
return(list(th=th,upd=upd))
}

caic.pois = function(tau,be,b,th,y,x,z,dv){
	if(tau==0){
		lam = exp(x%*%be+apply(z*bb(b),1,sum))
		q = matrix(1/2,nrow=nrow(lam),ncol=1)
	}else{
		ff = exp(tau*x%*%be)
		rr = exp(tau*apply(z*bb(b),1,sum))
		lam = (1/2*ff+1/2*rr)^(2/tau)
		q = ff/(ff+rr)
	}
	Q = diag(c(q))
	iqz = 2*z*((1-q)%x%matrix(1,nrow=1,ncol=dimb))
	QZ = matrix(rbind(matrix(iqz[1:ni[1],],ncol=dimb),matrix(0,nrow=sum(ni)-ni[1],ncol=dimb)),ncol=dimb)
	for(i in 2:m)
		QZ = cbind(QZ,matrix(rbind(matrix(0,nrow=sum(ni[1:(i-1)]),ncol=dimb),matrix(iqz[(sum(ni[1:(i-1)])+1):sum(ni[1:i]),],ncol=dimb),matrix(0,nrow=n-sum(ni[1:i]),ncol=dimb)),ncol=dimb))
	M = as.matrix(cbind(2*Q%*%x,QZ))
	W = diag(c(lam))
	di = matrix(0,nrow=dimb,ncol=dimb)
	for(i in 1:length(th))
		di = di + dv[[i]]*th[i]
	D = diagm(di,m)
	p = length(be)
	A = t(M)%*%W%*%M
	A[(p+1):(p+m*dimb),(p+1):(p+m*dimb)] = A[(p+1):(p+m*dimb),(p+1):(p+m*dimb)] + ginv(D)
	Ainv = try(ginv(A),silent=T)
	ef.df =  ifelse(class(Ainv)[1]!="try-error", sum(diag(t(M)%*%W%*%M%*%Ainv)), NA)
	caic = -2*(sum(y*log(lam)-lam)) + 2*ef.df
	return(list(caic=caic,ef.df=ef.df,L=sum(y*log(lam)-lam)))
}

