#Optimierung mit korrekten Varianzen und korrekten Determinanten
# General Data
#load("C:\\Users\\LocalAdmin\\Desktop\\Diss\\2015-08-15 Optim m korrekten Varianzen\\M2.RData")
checklist<-0 #controlling measure from inside of loops/functions
R=8e6 # Budget
d=8# dimension
# bet<-0.9 OLD 
bet<-0.85 #Smoothing parameter, percentage of the OLD PARAMETER exploration
#explstart<-1.0
#explend<-1.3
#bet2v<-c(seq(explstart,explend,length.out=(40-exex)))
#bet2v<-c(rep(1.14,10),seq(explstart,explend,length.out=(40-exex-10)))
bet2<-1.13 #Smoothing for exploitation, >1?Happy Engineering
draw<-0
dline<-1
dyn<-0 #=1 if regression changes dynamically
excludezero<-0 #if 0, the origin is added as a regression point
exex<-16 #Border between exploration and expliotation
punish<-10 # punishment exponent
forcecoef<-1 #forcing negative quadratic parabola
elaind<-2
exadd<- 0
DVBIND<-0 # all changed dvec entries are to be bound here
DVBIND2<-0 # all changed dvec entries are to be bound here
flag<-rep(0,9) #set + 1 if "global" optimum was chosen from i-th algorithm
dtest<-1 # deterministic function test if = 1, random functions if not
norming<-1 # alternative normings: if=1 reutrns are normed by functionindex, if=2 by seedsample, if=3 by optimization instance
require("quadprog")
require("GenSA")
require("car")
require("multcomp")
source("Variance List New.txt")
# Turnover functions 
fmult<-function(x,i) {a[i]*x^b[i]}
fex<-function(x,i) {Mex[i]*(1-exp(x*h[i]))*1000000}
fad<-function(x,i) {Mad[i]*x^fc[i]/(g[i]+x^fc[i])*1000000}
fadS<-function(x,i) {Mad[i]*x^fcS[i]/(gS[i]+x^fcS[i])*1000000}
epsadS<-function(x) {fcS*gS/(gS+x^fcS)}
opfAD<- function(x) {
	v<-0
	s<-0
	for(i in 1:(d-1)){
		v<- v - functions[[fin[i]]](x[i],i)
		s<- s+ x[i]
		}
	v<- v- functions[[fin[i]]](R-s,d)
	return(v)
	}
# program functions
dv<-function(colm) {c(floor((colm-1)/6)%%2,floor((((colm-1)%%6))/3),(-1*(((colm-1)%%3)-1)),floor((colm-1)/12))}
revenew<-function(v){
pos1<-which(v %in% sort(v)[1])
return(pos1)
}
quickopt<- function(fv){
#preparations
{
RESMAT<-matrix(nrow=9,ncol=d+1,data=0)
dimnames(RESMAT)<-list(c("ALB","ROT1","ROT2","ROT3","ROT1P","ROT2P","ROT3P","GenSA","QUAD"))
ZSTqo=matrix(nrow=20,ncol=8,data=0)
ZSTxqo=ZSTqo
DFqo<<-as.list(1:d)
u2<-rep(0, times=d)
u1<-u2
x1<-rep(R/d, times=d)
#a<-ust/(x1^b)
#mult<-function(x,i) {a[i]*x^b[i]}
#mult<-function(x,i) {nmult(multi,c(x,i))}
rv<-sample(1:d)
PV<-1/(rv+1)
#if(d%%2==0){PV<-rv*d/sum(rv)}else{PV<-rv*(d+1)/sum(rv)}
#rv<- c(3,8,4,5,1,7,2,6) 
#PV<-(rv*0.1 - 0.45)+c(1,1,1,1,1,1,1,1)
x2<-x1*PV
#x2<-x1+c(1e-5,1e-5,1e-5,1e-5,-1e-5,-1e-5,-1e-5,-1e-5)
x2<-R*x2/sum(x2)
nem<-x2

f<-function(x,i) {functions[[fv[i]]](x,i)}
for(i in 1:d){
u2[i]<-nmult(f,c(x2[i],i),sig=0)
u1[i]<-nmult(f,c(x1[i],i),sig=0)
}
er<-0
STVEC<-c(x1,x2,u1,u2)
}
#### Albers
{
for(j in 1:150){

	while(er!=0){
	er<-scheck(x1,x2)
	#print(er)
	if(er!=0){neww<-snew(x1,x2,f,R,er,osig=rep(0,d))
	x1<-neww[[1]];x2<-neww[[2]];u1<-neww[[3]];u2<-neww[[4]]
		}
		}
		#print(c(j,x1,x2))
	#eps<-log(u2/u1)/log(x2/x1)
	eps<-((u2-u1)/u1)/((x2-x1)/x1)
	#print(c(j,eps))
	if(is.nan(prod(eps))){
	break
	#DA<-GenSA(rep(0,d-1),rep(0,d-1),rep(R,d-1),fn=opfAD,control=list(maxit=500))

	#DA<<-DA
	#nn<-c(DA[[2]],R-sum(DA[[2]]))
	#uu1<- -DA[[1]]

	#RESLIST<<-list(n,nn,sum(u1),uu1)
	#if(sum(u1)<uu1){n<-nn;u1<-uu1}
	
	
	
	#return(list(vecreverse(pri),pri,sum(u1),n,nem))}
	}
	#if(prod(fv==4)){eps<-epsadS(x1)}
		if(j==1){print(eps)}
		eps<-epscheck1(eps)
		
	n<- u1*eps/(t(u1)%*%eps)*R
	if(j==1){pri<-n}
	x2<-x1
	u2<-u1
	x1<-n
	for(i in 1:d){
	u1[i]<-nmult(f,c(n[i],i),sig=0)

		}
		
		#if(j<6){print(c(j,x1,x2,u1,u2,n,eps))}
		#nor<-sum(u2);subop<-abs(opt-nor)/opt
				}#loop50
RESMAT[1,]<-c(n,sum(u1))
}
#### ROT 1 2 3 NOPERT
{
x1<-STVEC[1:d]
x2<-STVEC[(d+1):(2*d)]
u1<-STVEC[(2*d+1):(3*d)]
u2<-STVEC[(3*d+1):(4*d)]

up1<-u1; urel1<-u1; bp1<-x1; brel1<-x1; upret<-u1; urelret<-u1;usret<-u1
SAET<-x1; Maxvec<-1:d

ZSTQ3<-matrix(nrow=150,ncol=d,data=NA)

for(j in 1:150){


v<- up1/sum(up1)
vrel<-(urel1/brel1)/sum(urel1/brel1)
bp1<- R*v
brel1<- R*vrel
for(i in 1:d){
	upret[i]<-nmult(f,c(bp1[i],i),sig=0)#sig=rep(R,d))
	urelret[i]<-nmult(f,c(brel1[i],i),sig=0)
	usret[i]<-nmult(f,c(SAET[i],i),sig=0)
	ZSTQ3[j,i]<-usret[i]
	Maxvec[i]<-max(ZSTQ3[(1:j),i])
		}
		up1<-upret
		urel1<-urelret
		SAET<-R*Maxvec/sum(Maxvec)
		

	
}
RESMAT[2,]<-c(bp1,sum(upret))
RESMAT[3,]<-c(brel1,sum(urelret))
RESMAT[4,]<-c(SAET,sum(usret))	
}		
#### ROT 1 2 3 PERT
{
x1<-STVEC[1:d]
x2<-STVEC[(d+1):(2*d)]
u1<-STVEC[(2*d+1):(3*d)]
u2<-STVEC[(3*d+1):(4*d)]

up1<-u1; urel1<-u1; bp1<-x1; brel1<-x1; upret<-u1; urelret<-u1;usret<-u1
SAET<-x1; Maxvec<-1:d

ZSTQ3<-matrix(nrow=40,ncol=d,data=NA)

for(j in 1:40){


v<- up1/sum(up1)
vrel<-(urel1/brel1)/sum(urel1/brel1)
bp1<- R*v
brel1<- R*vrel
for(i in 1:d){
	upret[i]<-nmult(f,c(bp1[i],i),sig=sigm[i])
	urelret[i]<-nmult(f,c(brel1[i],i),sig=sigm[i])
	usret[i]<-nmult(f,c(SAET[i],i),sig=sigm[i])

	ZSTQ3[j,i]<-usret[i]

	Maxvec[i]<-max(ZSTQ3[(1:j),i])
		}
		up1<-upret
		urel1<-urelret
		SAET<-R*Maxvec/sum(Maxvec)
		
		


		}
	

#last update without perturbations
for(i in 1:d){
	upret[i]<-nmult(f,c(bp1[i],i),sig=0)#sig=rep(R,d))
	urelret[i]<-nmult(f,c(brel1[i],i),sig=0)
	usret[i]<-nmult(f,c(SAET[i],i),sig=0)
}




RESMAT[5,]<-c(bp1,sum(upret))
RESMAT[6,]<-c(brel1,sum(urelret))
RESMAT[7,]<-c(SAET,sum(usret))	
ZSTQ3<<-ZSTQ3
}		
#### SIMANN
{
#if(1-prod(1-(fv==4))){
	DA<-GenSA(rep(0,d-1),rep(0,d-1),rep(R,d-1),fn=opfAD,control=list(maxit=500))

	DA<<-DA
	nn<-c(DA[[2]],R-sum(DA[[2]]))
	uu1<- -DA[[1]]
#				}
	RESLIST<<-list(n,nn,sum(u1),uu1)
	#if(sum(u1)<uu1){n<-nn;u1<-uu1}
#	print(n)
	RESMAT[8,]<-c(n,sum(uu1))
}
###QUAD
{
osig=rep(1e9,8)
for(j in 1:20){
er<-scheck(x1,x2)# 1 if negative entry, 2 if equal entries, 3 if entries close to 0
	if(er!=0){print(c(j,"PERTOPT1"));neww<-snew(x1,x2,f,R,er,osig)
	x1<-neww[[1]];x2<-neww[[2]];u1<-neww[[3]];u2<-neww[[4]]
		}
		
		if(elaind==1){eps<-log(u2/u1)/log(x2/x1)}
		if(elaind==2){eps<-((u2-u1)/u1)/((x2-x1)/x1)}
		if(elaind==3){eps<-((u2-u1)/u2)/((x2-x1)/x2)}
		if(elaind==4){eps<-((u2-u1)/((u1+u2)/2))/((x2-x1)/((x1+x2)/2))}
		if((elaind-1)*(elaind-2)*(elaind-3)*(elaind-4)!=0){return("ERROR")}
	
	
	
		eps<-epscheck(eps)
		#if(j!=1){epsti<-eps*(1-betv[j]) +betv[j]*epsol}else{epsti<-rep(0.5,times=d)}
		if(j!=1){epsti<-eps*(1-bet) +bet*epsol}else{epsti<-rep(0.5,times=d)}
	n<- u1*epsti/(t(u1)%*%epsti)*R
	epsol<-epsti
	x2<-x1
	u2<-u1
	er2<-scheck(n,x2)
		if(er2!=0){print(c(j,"PERTOPT2"));new<-snew(n,x2,f,R,er2,osig)
	n<-new[[1]];x2<-new[[2]];u1<-new[[3]];u2<-new[[4]]
		}

	x1<-n
	for(i in 1:d){
	u1[i]<-nmult(functions[[fin[i]]],c(n[i],i),sig=osig[i])
				ZSTqo[j,i]<-u1[i]
				ZSTxqo[j,i]<-x1[i]
		}







}

TABq<<-as.list(1:d) #quadratic coefficients
TAB2q<<-as.list(1:d) #linear coefficients
CMq<<-matrix(nrow=d,ncol=3,data=NA) #ncol should be 3 if regression is with constant 2 w/out
CM2q<<-matrix(nrow=d,ncol=2,data=NA)
qfuq<<-as.list(1:d)
xoldq<-rep(NA,d)
posi1q<-NA
unq<-rep(0,d)
for(i in 1:d){

	DFqo[[i]]<<-as.data.frame(matrix(nrow=20,ncol=2,data=c(ZSTqo[1:20,i],ZSTxqo[1:20,i])))
	xoldq[i]<-ZSTxqo[20,i]
	#print(DF[[1]])
	DFqo[[i]]<<-cbind(DFqo[[i]],DFqo[[i]][,2]^2)
	names(DFqo[[i]])<<-c("y","a","b")
	}
	
	for(i in 1:d){
		TABq[[i]]<<-lm(y~a+b,data=DFqo[[i]])
		TAB2q[[i]]<<-lm(y~a,data=DFqo[[i]])
		CMq[i,]<<-TABq[[i]]$coefficients
		CM2q[i,]<<-TAB2q[[i]]$coefficients
	}


	DIVq<<- (CMq[,3]-abs(CMq[,3]))/2




	dvecq<<- CMq[,2]#with const
	for(i in 1:d){if(DIVq[i]==0){dvecq[i]<<-TAB2q[[i]]$coefficients[2];DIVq[i]<<- -1e-15
			#DVBIND<<-cbind(DVBIND,dvec[i])
			if(dvecq[i]<0){#dvec[i]<<-Meean2 
			ABq<-DFqo[[i]][c(1:j),]
			AB2q<-ABq[(ABq$y>15)|(ABq$y+ABq$a<50),]
			AB3q<-lm(y~a -1,data=as.data.frame(ABq))
			dvecq[i]<<-AB3q$coefficients[1]
			#DVBIND2<<-cbind(DVBIND2,dvec[i])
			}
		}
	}
	
	Dmatq<<- -2*diag(DIVq)
	Amatq<<-matrix(nrow=d,ncol=d+1,data=c(rep(-1,d),diag(rep(1,d))))
	bvecq<<- c(-R,rep(0,d))
	SOLq<<-solve.QP(Dmatq, dvecq, Amatq, bvecq, meq=1, factorized=FALSE)
	
	xnq<-SOLq[[1]]
		for(i in 1:d){

	unq[i]<-nmult(functions[[fin[i]]],c(xnq[i],i),sig=0)
	}
	
	RESMAT[9,]<-c(round(xnq),sum(unq))
	
}

# Returning
{
RESMAT<<-RESMAT
bestm<-revenew(vecreverse(RESMAT[,(d+1)]))
best<-bestm[1]
flag[best]<<-flag[best]+1
n<-RESMAT[best,(1:d)]
u1<-RESMAT[best,(d+1)]
print(best)
	return(list(vecreverse(pri),pri,sum(u1),n,nem))
}
}
 
{#controls and updates
scheck <-function(x1,x2) {
if(abs(sum(c(x1,x2))-sum(abs(c(x1,x2))))>1e-20){return(1) # at least one entry is negative
	}else{
		if(min(abs(x1-x2))<1e-5){return(2) # x1 and x2 entries too close to each other
			}else{
					if(min(c(x1,x2))<1e-5){return(3)} # at least one entry too close to zero
							return(0)}
	}
}
vecreverse<-function(v) { #reverses ORDERING of a vector, =! rev
no<-rev(sort(v))
no[rank(v)]
}
snew<- function(x1,x2,f,R,er,osig){
print("error code")
if(er==2){
	print(2)
	dis<-abs(x1-x2)
	u2<-rep(0, times=d)
	u1<-u2
	print(dis)
		for(i in 1:d){
		if(dis[i]<1e-10){b1a200<-runif(1)*(R/1000);x1[i]<-x1[i]+b1a200;x2[i]<-x2[i]-b1a200}
		}
	x1<-R*x1/sum(x1)
	x2<-R*x2/sum(x2)
	print(c(x1,x2))
	repeat{
		for(i in 1:d){
		u2[i]<-nmult(f,c(x2[i],i),sig=osig[i])
		u1[i]<-nmult(f,c(x1[i],i),sig=osig[i])
		}
		if(min(abs(1-u1/u2))>1/1e12){break}
		}
return(list(x1,x2,u1,u2))
}
if(er==3){
print(3)
u2<-rep(0, times=d)
u1<-u2
		for(i in 1:d){
			if(x1[i]<1e-10){x1[i]<-R/100}
			if(x2[i]<1e-10){x2[i]<-R/200}
			}
x1<-R*x1/sum(x1)
x2<-R*x2/sum(x2)
		for(i in 1:d){
		u2[i]<-nmult(f,c(x2[i],i),sig=osig[i])
		u1[i]<-nmult(f,c(x1[i],i),sig=osig[i])
		
		}
	return(list(x1,x2,u1,u2))
}
if(er==1){
print(1)
u2<-rep(0, times=d)
u1<-u2
		for(i in 1:d){
			if(x1[i]<1e-10){print(x1);x1[i]<-R/100}
			if(x2[i]<1e-10){x2[i]<-R/200}
			}
x1<-R*x1/sum(x1)
x2<-R*x2/sum(x2)


		for(i in 1:d){
		u2[i]<-nmult(f,c(x2[i],i),sig=osig[i])
		u1[i]<-nmult(f,c(x1[i],i),sig=osig[i])
		}
		#if(min(abs(1-u1/u2))>1/1e12){break} #Nonsense
		
		
return(list(x1,x2,u1,u2))
}
print("ERROR invalid error code");return(NaN)
}
functions<-c(fmult,fex,fad,fadS) #LIST, not a vector
nmult <- function(f,x,mag=1,r=0,n=rnorm,mu=00,sig=1000,p=0) { #normally distributed perturbation
	if(r==1){set.seed(1789)}
	d=length(x[1])
	x[1]<-(abs(x[1])+x[1])/2 #x<-max(x,0)
	A<- f(x[1],x[2])
	val<-A+mag*(n(length(A),mu,sqrt(sig)))
	return(max(c(val,10))) 
	}
epscheck<-function(x) {
l<-0;h<-0
	for(i in 1:length(x)){
		if(x[i]<0.01){x[i]<-0.01;l<-l+1}
		if(x[i]>0.5){x[i]<-0.5;h<-h+1}
	}
			#print(c("low",l,"high",h))
	return(x)
}
epscheck1<-function(x) {
l<-0;h<-0
	for(i in 1:length(x)){
		if(x[i]<0.05){x[i]<-0.05;l<-l+1}
		if(x[i]>0.5){x[i]<-0.5;h<-h+1}
	}
			#print(c("low",l,"high",h))
	return(x)
}
}
pertopt<- function(fin,sigm,tup,maxit=40){
#preparations
{
exex<-15
ANOVyVec<-rep(NA,8)
names(ANOVyVec)<-c("A1subop","A2subop","A3subop","A4subop","A1DB","A2DB","A3DB","A4DB")
ZST<-matrix(nrow=maxit,ncol=d*5,data=NA) ##captures the u-Vectors in every period of every algorithm
ZSTx<<-matrix(nrow=maxit,ncol=d,data=NA) ##captures the x-Vectors in every period of our alg
zv<-rep(-Inf,times=14) # the first 4 won't be changed or passed on
osig<-sigm #Variances from regression
if(length(osig!=d)){osig<-rep(mean(osig),d)}
x1<-tup[[1]] #start values from quickopt
x2<-tup[[2]]
u2<-rep(0, times=d)
u1<-u2
opt<-tup[[3]] #unperturbed optimum from quickopt
Xoptv<-tup[[4]] #Uoptv generated in next loop
#print(Xoptv)
Uoptv<-0
ureturn<-rep(0,times=d)
f<-function(x,i) {functions[[fin[i]]](x,i)}
}
#first update

for(i in 1:d){ 
u2[i]<-nmult(f,c(x2[i],i),sig=osig[i])
u1[i]<-nmult(f,c(x1[i],i),sig=osig[i])
ureturn[i]<-nmult(f,c(x1[i],i),sig=0)
Uoptv[i]<-nmult(f,c(Xoptv[i],i),sig=0)
}

zv[5]<-(1-sum(ureturn)/opt)*100

# bad rules including Albers original and saturation rule

valmat<-matrix(ncol=5,nrow=maxit,data=NA)
up1<-u1; urel1<-u1; bp1<-x1; brel1<-x1; upret<-u1; urelret<-u1;ua1<-u1;ua2<-u2;xa1<-x1;xa2<-x2;uaret<-u1;usret<-u1;us<-u1
SAET<-x1; Maxvec<-1:d	
for(j in 1:maxit){
epsa<-log(ua2/ua1)/log(xa2/xa1)
na<- ua1*epsa/(t(ua1)%*%epsa)*R
	xa2<-xa1
	ua2<-ua1
	xa1<-na

		



v<- up1/sum(up1)
vrel<-(urel1/brel1)/sum(urel1/brel1)
bp1<- R*v
brel1<- R*vrel
for(i in 1:d){
	up1[i]<-nmult(f,c(bp1[i],i),sig=osig[i])
	urel1[i]<-nmult(f,c(brel1[i],i),sig=osig[i])
	us[i]<-nmult(f,c(SAET[i],i),sig=osig[i])
	ua1[i]<-nmult(functions[[fin[i]]],c(na[i],i),sig=osig[i])
	ZST[j,i+d]<-up1[i]
	ZST[j,i+2*d]<-urel1[i]
	ZST[j,i+3*d]<-ua1[i]
	ZST[j,i+4*d]<-us[i]
	uaret[i]<-nmult(f,c(na[i],i),sig=0)
	upret[i]<-nmult(f,c(bp1[i],i),sig=0)
	urelret[i]<-nmult(f,c(brel1[i],i),sig=0)
	usret[i]<-nmult(f,c(SAET[i],i),sig=0)
	Maxvec[i]<-max(ZST[(1:j),i+4*d])
		}
		SAET<-R*Maxvec/sum(Maxvec)
		#norp<-sum(upret);norrel<-sum(urelret);subopp<-(abs(opt-norp)/opt);suboprel<-(abs(opt-norrel)/opt) #SUBOPTIMALITIES
		#nora<-sum(uaret);subopa<-abs(opt-nora)/opt											#SUBOPTIMALITIES
		norp<-sum(upret);norrel<-sum(urelret);subopp<-norp/opt;suboprel<-norrel/opt #OPTIMALITIES
		nora<-sum(uaret);subopa<-nora/opt											#OPTIMALITIES
		nors<-sum(usret);subops<-nors/opt
		#norp<-sum(upret);norrel<-sum(urelret);subopp<-(opt-norp);suboprel<-(opt-norrel) # unnormed SUBOPTIMALITIES
		#nora<-sum(uaret);subopa<-opt-nora		#unnormed SUBOPTIMALITIES
		
		valmat[j,2:5]<-round(c(100*subopp,100*suboprel,100*subopa,100*subops),4)
	if(j==1){zv[c(7,8)]<-round(c(100*subopp,100*suboprel),4)}	
	if(j==25){zv[c(10,11)]<-round(c(100*subopp,100*suboprel),4)}	
	#if(sw[2]==1&subopp<0.001){zv[13]<-j; sw[2]<-0}
	#if(sw[3]==1&suboprel<0.001){zv[14]<-j; sw[3]<-0}
	#if(j==maxit){
	#		if(sw[2]==0&subopp>0.001){zv[13]<- -Inf; sw[2]<-1}
	#		if(sw[3]==0&suboprel>0.001){zv[14]<- -Inf; sw[3]<-1}
	#	}

}
#good rules
#for(j in 1:exex){
j=1
while(j<exex+1){  
	er<-scheck(x1,x2)# 1 if negative entry, 2 if equal entries, 3 if entries close to 0
	if(er!=0){print(c(j,"PERTOPT1"));neww<-snew(x1,x2,f,R,er,osig)
	x1<-neww[[1]];x2<-neww[[2]];u1<-neww[[3]];u2<-neww[[4]]
		}
		
		if(elaind==1){eps<-log(u2/u1)/log(x2/x1)}
		if(elaind==2){eps<-((u2-u1)/u1)/((x2-x1)/x1)}
		if(elaind==3){eps<-((u2-u1)/u2)/((x2-x1)/x2)}
		if(elaind==4){eps<-((u2-u1)/((u1+u2)/2))/((x2-x1)/((x1+x2)/2))}
		if((elaind-1)*(elaind-2)*(elaind-3)*(elaind-4)!=0){return("ERROR")}
	
	
	
		eps<-epscheck(eps)
		#if(j!=1){epsti<-eps*(1-betv[j]) +betv[j]*epsol}else{epsti<-rep(0.5,times=d)}
		if(j!=1){epsti<-eps*(1-bet) +bet*epsol}else{epsti<-rep(0.5,times=d)}
	n<- u1*epsti/(t(u1)%*%epsti)*R
	epsol<-epsti
	x2<-x1
	u2<-u1
	er2<-scheck(n,x2)
		if(er2!=0){print(c(j,"PERTOPT2"));new<-snew(n,x2,f,R,er2,osig)
	n<-new[[1]];x2<-new[[2]];u1<-new[[3]];u2<-new[[4]]
		}

	x1<-n
	for(i in 1:d){
	u1[i]<-nmult(functions[[fin[i]]],c(n[i],i),sig=osig[i])
	ureturn[i]<-nmult(functions[[fin[i]]],c(n[i],i),sig=0)
	ZST[j,i]<-u1[i]
	ZSTx[j,i]<-x1[i]
		}
		
		#nor<-sum(ureturn);subop<-abs(opt-nor)/opt #Suboptimality
		nor<-sum(ureturn);subop<-nor/opt #optimality
		#nor<-sum(ureturn);subop<- opt-nor # unnormed Suboptimality
		valmat[j,1]<-round(100*subop,4)
		print(valmat[j,]) #ausgabe periodenweise
	if(j==1){zv[6]<-round(100*subop,4)}
	if(j==25){zv[9]<-round(100*subop,4)}
	#if(sw[1]==1&subop<0.001){zv[12]<-j; sw[1]<-0}
		#if(j==maxit){
	#	if(sw[1]==0&subop>0.001){zv[12]<- -Inf; sw[1]<-1}
		#}
		if(j==9){exex<-min(max((round(log(mean(diag(var(ZST[1:9,1:8]))),10)*8-81)),10),25)+exadd;print(exex);checklist<<-ZST;chx<<-ZSTx}
	j=j+1
	}

DF<<-as.list(1:d)
TAB<<-as.list(1:d) #quadratic coefficients
TAB2<<-as.list(1:d) #linear coefficients
CM<<-matrix(nrow=d,ncol=3,data=NA) #ncol should be 3 if regression is with constant 2 w/out
CM2<<-matrix(nrow=d,ncol=2,data=NA)
qfu<<-as.list(1:d)
xold<-rep(NA,d)
posi1<-NA
un<-rep(0,d)
#nor<-0
#indic<-0#rewarding improvement if 1, not in exex
for(i in 1:d){

	DF[[i]]<<-as.data.frame(matrix(nrow=exex,ncol=2,data=c(ZST[1:exex,i],ZSTx[1:exex,i])))
	xold[i]<-ZSTx[exex,i]
	#print(DF[[1]])
	DF[[i]]<<-cbind(DF[[i]],DF[[i]][,2]^2)
	names(DF[[i]])<<-c("y","a","b")
	}
	
	#plot(y=ZST[1:exex,1],x=ZSTx[,1])
for(j in (exex+1):maxit){

	#print(c(j,x1,sum(u1)))
	for(i in 1:d){
	
	#if(forcecoef==0){
		if(excludezero==0){TAB[[i]]<<-lm(y~a+b,data=as.data.frame(rbind(DF[[i]],c(0,0,0))))}else{#const
			
		#	if(excludezero<0){ZMA<-matrix(nrow=ceiling((j-exex)/(-excludezero)),ncol=3,data=0);colnames(ZMA)=c("y","a","b");
		#	TAB[[i]]<<-lm(y~a+b,data=as.data.frame(rbind(DF[[i]],ZMA)))
		#	}
		
		TAB[[i]]<<-lm(y~a+b,data=DF[[i]])}#const
		TAB2[[i]]<<-lm(y~a,data=DF[[i]])
		
		#}else{
		#AMAT<-matrix(ncol=3,nrow=j,data=c(rep(1,j),DF[[i]][,2],DF[[i]][,3]))
		#BVEC<-DF[[i]][,1]
		#BUNDER <-c(-1e10,-1e10,-1e10)
		#BOVER<-c(1e10,1e10,0)
		#AM0<-rbind(AMAT,c(1,0,0))
		#B0<-c(BVEC,0)
		#print(BVEC)
		#		if(excludezero==0){TAB[[i]]<<-bvls(AM0,B0,BUNDER,BOVER)}else{
		#TAB[[i]]<<-bvls(AMAT,BVEC,BUNDER,BOVER)}
		#if(j==exex+1){checklist<<-TAB}
		#names(TAB[[i]])[1]<<-"coefficients"

		#}
		
	CM[i,]<<-TAB[[i]]$coefficients
	CM2[i,]<<-TAB2[[i]]$coefficients
	qfu<<-function(x,i) {CM[i,3]*x^2+CM[i,2]*x+CM[i,1]} #with const
	#qfu<<-function(x,i) {CM[i,2]*x^2+CM[i,1]*x}#w/o const
	}
	if(draw==1){
		par(mfrow=c(ceiling(sqrt(d)),ceiling(sqrt(d))))
		for(i in 1:d){
		tfx<-function(x) {functions[[fin[i]]](x,i)}
		tfu<-function(x) {qfu(x,i)}
		#plot(rbind(DF[[i]][,c(2,1)],c(0,0),c(1e7,1e7)))
		#plot(rbind(DF[[i]][,c(2,1)],c(0,0),c(1e7,functions[[fin[i]]](1e7,i))))
		plot(rbind(DF[[i]][,c(2,1)],c(0,0),c(R,functions[[fin[i]]](R,i))))
		points(Xoptv[i],Uoptv[i],col="green",pch=16)
		if(j>exex+1){points(xperter[i],uperter[i],col="red",pch=16)}
		if(j==exex+1){CM1<-CM;DF1<-DF
		PERTCUR<-function(x) {qfu1(x,i)};qfu1<<-function(x,i) {CM1[i,3]*x^2+CM1[i,2]*x+CM1[i,1]}}# with const
		#PERTCUR<-function(x) {qfu1(x,i)};qfu1<<-function(x,i) {CM1[i,2]*x^2+CM1[i,1]*x}} #w/o const
		#plot(rbind(DF[[5]][,c(2,1)],c(0,0),c(1e7,fex(1e7,8))))
		#if(j==exex+1|j==exex+10){
		#plot(DF1[[i]][,c(2,1)],col="red",add=T)
		points(DF1[[i]][,2],DF1[[i]][,1],col="red",pch=1)
		curve(tfu,0,R*1.2,add=TRUE)
		curve(PERTCUR,0,R*1.2,add=TRUE,col="red")
		curve(tfx,0,R*1.2,add=TRUE)
		if(dline==1){abline(TAB2[[i]][1],TAB2[[i]][2],col="blue")}
		}
		if(d==8){
			tfx<-function(x) {functions[[fin[6]]](x,6)}
			curve(tfx,0,R/8,add=FALSE)
		for(din in 1:8){
			tfx<-function(x) {functions[[fin[din]]](x,din)}
				curve(tfx,0,R/8,add=TRUE)
				
			}
		}
		
	}
	#}
	#curve(y=qfu(seq(0,1e8,1e6)),x=seq(0,1e8,1e6),add=T)
		#print(summary(TAB[[1]]))
	{#with const
	
	DIAVEC<<- (CM[,3]-abs(CM[,3]))/2
	DIAVECn<<- (CM[,2]+abs(CM[,2]))/2
	DIAVEC2<<- (CM2[,2]+abs(CM2[,2]))/2
	Meean<<-sum(DIAVEC)
	Meeann<<-mean(CM[,2])
	Meean2<<-mean(DIAVEC2)
	}
	{#without const
	
	#DIAVEC<<- (CM[,2]-abs(CM[,2]))/2
	#DIAVECn<<- (CM[,1]+abs(CM[,1]))/2
	#DIAVEC2<<- (CM2[,1]+abs(CM2[,1]))/2
	#Meean<<-sum(DIAVEC)
	#Meeann<<-mean(CM[,1])
	#Meean2<<-mean(DIAVEC2)
	}
	
	DIV<<-DIAVEC
	#for(i in 1:d){if(DIV[i]==0){DIV[i]<<-Meean}}

	dvec<<- CM[,2]#with const
	#dvec<<- CM[,1] #w/o
	for(i in 1:d){if(DIV[i]==0){dvec[i]<<-TAB2[[i]]$coefficients[2];DIV[i]<<- -1e-15
			#DVBIND<<-cbind(DVBIND,dvec[i])
			if(dvec[i]<0){#dvec[i]<<-Meean2 
			AB<-DF[[i]][c(1:j),]
			AB2<-AB[(AB$y>15)|(AB$y+AB$a<50),]
			AB3<-lm(y~a -1,data=as.data.frame(AB))
			dvec[i]<<-AB3$coefficients[1]
			#DVBIND2<<-cbind(DVBIND2,dvec[i])
			}
		}
	}
	
	Dmat<<- -2*diag(DIV)
	Amat<<-matrix(nrow=d,ncol=d+1,data=c(rep(-1,d),diag(rep(1,d))))
	bvec<<- c(-R,rep(0,d))
	SOL<<-solve.QP(Dmat, dvec, Amat, bvec, meq=1, factorized=FALSE)
	
	xn<-SOL[[1]]
	xn<<-xn
	xx<-bet2*xn +(1-bet2)*xold
	xold<-xx #newly added, wasn't used before, bad...
	for(i in 1:d){
	#un[i]<-nmult(functions[[fin[i]]],c(xn[i],i),sig=osig[i])
	un[i]<-nmult(functions[[fin[i]]],c(xx[i],i),sig=osig[i])
	ureturn[i]<-nmult(functions[[fin[i]]],c(xx[i],i),sig=0)
	}
	if(j==exex+1){xperter<-xn;uperter<-un}
	#if(nor<sum(un)){indic<-1}else{indic<-0}#checking for improvement
			#nor<-sum(un);subop<-abs(opt-nor)/opt
		#	nur<-sum(ureturn);subop<-(opt-nur)/opt #Suboptimality
			nur<-sum(ureturn);subop<-nur/opt #optimality
			#	nur<-sum(ureturn);subop<-opt-nur # unnormed Suboptimality
		valmat[j,1]<-round(100*subop,4)
		print(valmat[j,]) #ausgabe periodenweise
		if(j==25){zv[9]<-round(100*subop,4)}
		DIFMAT<<-matrix(nrow=exex,ncol=d,data=NA) #Finde x-vektor mit geringestem Abstand zu neuem Vektor
		MEANMAT<<-matrix(nrow=exex,ncol=d,data=NA)#Finde x-vektor mit geringstem Abstand zum Mittel
		for(i in 1:d){for(kindex in 1:exex){
		DIFMAT[kindex,i]<<-(DF[[i]][kindex,2]-xn[i])^punish
		MEANMAT
		}}
		posi1<-revenew(rowMeans(DIFMAT))
	#if(indic==0){
		for(i in 1:d){
			#DF[[i]]<<-rbind(DF[[i]][2:exex,],c(un[i],xn[i],xn[i]^2))
			if(dyn==1){
			DF[[i]][posi1,]<<-c(un[i],xn[i],xn[i]^2)}else{
			DF[[i]]<<-rbind(DF[[i]],c(un[i],xn[i],xn[i]^2))}
				ZST[j,i]<-un[i]	
				ZSTx[j,i]<-x1[i]
				#print(c(j,ZST[j,i],un[i]))
		}
		#{
		#		print(c(mean(rowSums(ZST[1:j,1:d])),mean(rowSums(ZST[1:j,(d+1):(2*d)])),mean(rowSums(ZST[1:j,(2*d+1):(3*d)])),mean(rowSums(ZST[1:j,(4*d+1):(5*d)]))))
		#}
	#}
}
	#print(ZST)	
	#return(zv[5:14])
	#print(valmat)
	#ANOVyVec[1:3]<-colMeans(valmat[((exex+1):maxit),])[1:3]
	#ANOVyVec[4]<-mean(rowSums(ZST[((exex+1):maxit),1:d]))
	#ANOVyVec[5]<-mean(rowSums(ZST[((exex+1):maxit),(d+1):(2*d)]))
	#ANOVyVec[6]<-mean(rowSums(ZST[((exex+1):maxit),(2*d+1):(3*d)]))
	
	#ANOVyVec[8]<-mean(rowSums(ZST[,(3*d+1):(4*d)]))
	
	ANOVyVec[1:4]<-colMeans(valmat)[c(1,2,3,5)]
	ANOVyVec[5]<-mean(rowSums(ZST[,1:d]))
	ANOVyVec[6]<-mean(rowSums(ZST[,(d+1):(2*d)]))
	ANOVyVec[7]<-mean(rowSums(ZST[,(2*d+1):(3*d)]))
	ANOVyVec[8]<-mean(rowSums(ZST[,(4*d+1):(5*d)]))
	checklist<<-ZST
	CLST<<-ZSTx
	return(ANOVyVec)
}

freqtest<-function(fin,colm,R2,oldsig,iter=20){
count<-0
for(i in 1:iter){
res<-omkv(fin,colm,R2,oldsig)
if(res[[9]]>10){count=count+1}
} 
return(count/iter)
}
coefnew<-function(fin,dv2,colm){# calculates all 9 new coefficient sets, returns the necessary ones
#dv2 has two entries, both are coded in binary.first: elasticity close 0 or far 1, second: M close 0 or far 1
	#setting DB distribution
	ust<-rep(2500000,d) #if(detvec[2]!=1){ust<-rep(2500000, times=d)}else{ust<-rep(c(1,1,4,4), times=ceiling(d/4))[1:d]*1000000}
	#setting b/fc distribution
	if(dv2[1]!=1){if(d==8){b<-c(0.26,0.27,0.28,0.29,0.31,0.32,0.33,0.34)}else{b<-seq(from=0.26, to=0.34, by=0.08/(d-1))}
						}else{b<-seq(from=0.11, to=0.18,length=d);b[(floor(d/2)+1):d]<-seq(from=0.47, to=0.5, length=floor((d+1)/2))}
	#fc<<-b
	b<<-b
	#fcS<<-fc+1
	#setting Maximum Values 
	if(dv2[2]!=1){if(d==8){Mex<<-c(6.1,6.2,6.3,6.4,6.6,6.7,6.8,6.9)}else{Mex<<-seq(from=6.1, to=6.9, by=0.8/(d-1))}
			}else{Mex<<-rep(c(4.5,10),times=ceiling(d/2))[1:d]}
	Mad<<-Mex
	#calculating fmult coefficient a
	###############aid<-rep(R/d, times=d)
	a<<-Mex*R^(-b)*1000000
	#calculating fex coefficient h
	h<<-rep(-0.1,d) #starting values
	for( i in 1:100){h<<- log((b*(1-exp(h*R/d)))/(-h*R/d))/(R/d)}

	#calculating fad coefficient g
	fc<<-rep(0.5,d) #starting values
	fcS<<-rep(1.5,d) #starting values
	g<<-rep(1e5,d) #starting values
	gS<<-rep(1e5,d) #starting values
	#Convmat<-matrix(ncol=d,byrow=TRUE,nrow=4,data=c(rep(12,8),rep(35,4),rep(5,4),rep(c(10,20),4),rep(30,4),rep(5,4),c(15,55,15,55,5,55,5,55)))
		CR1<-2*floor(Mad)
		CR2<- -10*ceiling(10*b)+55
		CR3<- 2*ceiling(Mad)
		CR4<- ((ceiling(b*10))*(-4/3)+23/3)*ceiling(Mad)
		Convmat2<-matrix(ncol=d,byrow=TRUE,nrow=4,data=c(CR1,CR2,CR3,CR4))
		
		for(i in 1:2000){fc<<-log((ust*g)/(Mad*1e6))/log(1e6);g<<-b*(g+1e6^fc)/fc
					fcS<<-log((ust*Convmat2[colm,]*gS)/(Mad*1e6))/log(1e6);gS<<-b*(gS+1e6^fcS)/fcS}
	#calculating fadS coefficient gS
#	gS<<-(Mad*1000000/ust -1 )*1000000^fcS
		#calculating starting solution
	v<-quickopt(fin)
	#if(detvec[3]==0){x1<-rep(R/d, times=d)}else{if(detvec[3]==1){x1<-v[[2]]}else{x1<-v[[1]]}}
	x1<-rep(R/d,d)
		return(list(x1,v[[5]],v[[3]],v[[4]]))
}

{#variance and R^2 calculations
newgen<-function(f1,i,n,sig=0){
x1<-sample(1:R,n,replace=T)
x<-matrix(ncol=1, nrow=n,data=x1)
y<-f1(x1,i) +rnorm(n,sd=sqrt(sig))
return(data.frame(x=x,y=y))
}
rad<- function(f=fad,i,n=2000,s){
r<<-1
NT<-T
while(NT==T){
DF<<-newgen(f,i,n,s)
nfunc<-function(dff) {nls(y~1/((a*x^b)+1), data=dff,start = list(b = -fc[i],a=g[i]),trace=F,alg="plinear", control=list(tol=1e-1,maxiter = 500,warnOnly=T))}
err<<-tryCatch(TAB<-nfunc(DF),finally=if(is.na(TAB[2])==T){print(4);rad(f,i,n,s)})
NT<-is.na(TAB[2])
}
bhat<-summary(TAB)$coefficients[1]
ahat<-summary(TAB)$coefficients[2]
Mhat<-summary(TAB)$coefficients[3]
yhat<-Mhat/(ahat*(DF$x)^bhat+1)
den<-0; num<-0
bar1<-mean(DF$y)
for(j in 1:n){
num<-num+(DF$y[j]-yhat[j])^2
den<-den+(DF$y[j]-bar1)^2
}

return(1-num/den)
} 
radS<- function(f=fadS,i,n=2000,s){
rS<<-1
NT<-T
while(NT==T){
DF<<-newgen(f,i,n,s)
nfunc<-function(dff) {nls(y~1/((a*x^b)+1), data=dff,start = list(b = -fcS[i],a=gS[i]),trace=F,alg="plinear", control=list(tol=1e-1,maxiter = 1000,warnOnly=T))}
err<<-tryCatch(TAB<<-nfunc(DF),finally=if(is.na(TAB[2])==T){print(4);radS(f,i,n,s)})
NT<-is.na(TAB[2])
}
bhat<-summary(TAB)$coefficients[1]
ahat<-summary(TAB)$coefficients[2]
Mhat<-summary(TAB)$coefficients[3]
yhat<-Mhat/(ahat*(DF$x)^bhat+1)
den<-0; num<-0
bar1<-mean(DF$y)
for(j in 1:n){
num<-num+(DF$y[j]-yhat[j])^2
den<-den+(DF$y[j]-bar1)^2
}

return(1-num/den)
} 
radown<-function(f=fad,i,n=2000,s){
r3<<-1
DATF<<-newgen(f,i,n,s)
yhatt<-f(DATF[,1],i)
ybarr<-mean(yhatt)
NUME<-sum((yhatt-ybarr)^2)
DENOM<-sum((DATF[,2]-ybarr)^2)
return(NUME/DENOM)
}
radSown<-function(f=fadS,i,n=2000,s){
r3<<-1
DATFS<<-newgen(f,i,n,s)
yhatt<-f(DATFS[,1],i)
ybarr<-mean(yhatt)
NUME<-sum((yhatt-ybarr)^2)
DENOM<-sum((DATFS[,2]-ybarr)^2)
return(NUME/DENOM)
}
rmult1<- function(f=fmult,i,n=2000,s){
DF<-newgen(f,i,n,s)
TAB<-nls(y~x^b,data=DF,start = list(b = 0.1),alg = "plinear",trace=F, control=list(tol=1,maxiter = 500,warnOnly=T))
bhat<-summary(TAB)$coefficients[1]
ahat<-summary(TAB)$coefficients[2]
yhat<-ahat*(DF$x)^bhat
den<-0; num<-0
bar1<-mean(DF$y)
for(j in 1:n){
num<-num+(DF$y[j]-yhat[j])^2
den<-den+(DF$y[j]-bar1)^2
}
return(1-num/den)
}
rexp1<- function(f=fex,i=1,n=2000,s=0){
DF<-newgen(f,i,n,s)
TAB<-nls(y~(1-exp(-x*b)),data=DF,start = list(b = -1e-5),trace=F,alg="plinear", control=list(tol=1,maxiter = 500,warnOnly=T))
bhat<-summary(TAB)$coefficients[1]
Mhat<-summary(TAB)$coefficients[2]
yhat<-Mhat*(1-exp((DF$x)*-bhat))
den<-0; num<-0
bar1<-mean(DF$y)
for(j in 1:n){
num<-num+(DF$y[j]-yhat[j])^2
den<-den+(DF$y[j]-bar1)^2
}
return(1-num/den)
}
rfunction<-c(rmult1,rexp1,radown,radSown) #LIST
var8test<-function(fin,R2,llimit=20,splits=10){
if(sum(fin)==fin[1]){fin<-rep(fin,times=d)}
if(R2>1||R2<0){return("R2-ERROR")}
if(splits<5){return("splERROR")}
varvec<-rep(0,times=d)
rvec<-rep(0,times=d)
for(i in 1:d){
	fi<-fin[i] 
	if(fi>2){n<-rep(5000,times=21)}else{n<-c(100,500,1000,2000,rep(5000,times=17))}
	e<-1:splits;ri<-1:splits;rdarrange<-1:splits;l<-0
	si<-10^(10-splits/2 +e)
	r1<- -1
		while((abs(R2-r1)>0.0001|l<5)&&l<llimit){
		l<-l+1
		#n<-l*n
	for(kk in 1:splits){
	ri[kk]<-rfunction[[fi]](functions[[fi]],i=i,n=n[l+1],s=si[kk])
	rdarrange[kk]<-abs(ri[kk]-R2)
	}
	pos1<-which(rdarrange %in% sort(rdarrange)[1])
	plumin<-ri[pos1]-R2
	if(plumin>0){
	if(pos1==splits){pos2=pos1-1}else{pos2<-pos1+1}
	}else{
	if(pos1==1){pos2=pos1+1}else{pos2<-pos1-1}
	}
	#print(c(plumin,pos1,pos2))
	#pos2<-which(rdarrange %in% sort(rdarrange)[2])
	dissss<-((log(si[pos1])+log(si[pos2]))/2)# does not appear?
	stp<-(e-1)/(splits-1)
	#si<-exp(log(si[pos1])*(1-(stp))+log(si[pos2])*stp)
	ma<-max(si[pos1],si[pos2])
	mi<-min(si[pos1],si[pos2])
	si<-ma*stp + (1-stp)*mi
	r1<-ri[pos1]
	print(c(l,r1,si[pos1]))
	}
	varvec[i]<-si[pos1]
	rvec[i]<-r1
	}
return(varvec)
}
#var8test(rep(4,times=8),0.9)
ctest<-function(e,d,co,M){ # input: binary-coded detvec, output: column
testvector<-c(e,d,co,M)
for(colm in 1:48){
detvec<-c(floor((colm-1)/6)%%2,floor((((colm-1)%%6))/3),(-1*(((colm-1)%%3)-1)),floor((colm-1)/12))
test<-0
for(tt in 1:4){
if(testvector[tt]==detvec[tt]){test<-test+1}
}
if(test==4){return(colm)}
}
return("ungültige Eingabe")
}
}

omkv<-function(fin,colm,R2,oldsig){
if(sum(fin)==fin[1]){fin<-rep(fin,times=d)}#usually fin is a vector
d<<-length(fin)
fin<<-fin
#zv<-rep(-Inf,times=14)
#sw<<-c(1,1,1)  #Switches, turned off, if suboptimality is less than.1%
dv2<-c(abs(colm%%2-1),ceiling(colm/2)-1)
#detvec[1]=0 if b/fc are close, =1 of they're far apart (0 if modexp)
#detvec[2]=0 if DB=2500000, 1 if it's 1Mio and 4Mio
#detvec[3]=-1,0,1 depending on the correlation of the starting solution to the optimum
#detvec[4]=0 if Maxvalues are close, 1 if they aren't. 
#zv[c(1,2,3,4)]<-detvec

RRR<-(-R2*10+11)/2
if(R2==1){sigm<<-rep(0,times=d)}else{
if(oldsig==0){sigm<<-var8test(fin=fin,R2=R2)}else{if(oldsig==2){
for(i in 1:d){
sigm[i]<<-VLN[[fin[i]]][[RRR]][i,colm]
}
}}}
tup<-coefnew(fin,dv2,colm)
ANOVV<-pertopt(fin,sigm,tup)
#vecccc<-pertopt(fin,sigm,tup)
#zv[5:14]<-vecccc
#MAAA<-as.matrix(zv)
#rownames(MAAA)<-c("Elastizitäten","Deckungsbeiträge","Korrelation","M bzw. Sättigung","Startfehler","1.Per DB*EL---->","1.Per DB","1.Per DB/Budg","25.Per DB*EL--->","25.Per DB","25.Per DB/Budg","Konverg. DB*EL","Konverg. DB", "Konverg. DB/Budg")
#colnames(MAAA)<-colm
#return(MAAA)
return(ANOVV)
}
#omkv(2,1,1,1)
{#generating tables
Tabgen<-function(fin){ #old, produces Table as in Albers' Paper
v<<-matrix(ncol=24, nrow=14, data=0)
for(i in c(1,2,3,4,5,6,13,14,15,16,17,18)){
 v[,i]<<-omkv(fin,(i),1,0)

}
return(v)
}
Varmatrix<- function(){ #Naive function trying to produce column 6 of every function's Variance matrix
M=matrix(nrow=6,ncol=4,dimnames=list(c("0.9 var","0.9 sd","0.7 var","0.7 sd","0.5 var","0.5 sd"),c("mult","exp","ADB K","ADB S")),data=0)
omkv(c(1,2,3,4),6,0.9,0)
M[1,]<-sigm
stab<-sqrt(sigm)
M[2,]<-stab
omkv(c(1,2,3,4),6,0.7,0)
M[3,]<-sigm
stab<-sqrt(sigm)
M[4,]<-stab
omkv(c(1,2,3,4),6,0.5,0)
M[5,]<-sigm
stab<-sqrt(sigm)
M[6,]<-stab
return(M)
}
#set.seed(1)
#sigm<-rep(1e11,times=d)


#VLNi<-list(matrix(nrow=8,ncol=4,data=NA),matrix(nrow=8,ncol=4,data=NA),matrix(nrow=8,ncol=4,data=NA))
#VLN<-list(VLNi,VLNi,VLNi,VLNi)

varlistnew<-function(){ #produces variance list 
for(fin in 1:4){ #1:4
	for(rr in c(0.9,0.7,0.5)){ #c(0.9,0.7,0.5)
			RRR<-(-rr*10+11)/2
				for(colm in 1:4){ #1:nc
			omkv(fin,colm,rr,0)
			VLN[[fin]][[RRR]][,colm]<<-sigm
			}
			dump("VLN",file="Variance List New.txt")
			}
		}
	}	
#varlistnew()

#R2LN<-list(VLNi,VLNi,VLNi,VLNi)
 #RL1<-list(matrix(nrow=7,ncol=12,data=0),matrix(nrow=7,ncol=12,data=0),matrix(nrow=7,ncol=12,data=0))
 #RL2<-list(matrix(nrow=7,ncol=12,data=0),matrix(nrow=7,ncol=12,data=0),matrix(nrow=7,ncol=12,data=0))
 #RL3<-list(matrix(nrow=7,ncol=24,data=0),matrix(nrow=7,ncol=24,data=0),matrix(nrow=7,ncol=24,data=0))
 #RL4<-list(matrix(nrow=7,ncol=24,data=0),matrix(nrow=7,ncol=24,data=0),matrix(nrow=7,ncol=24,data=0))
 #NRL<-list(matrix(nrow=7,ncol=12,data=0),matrix(nrow=7,ncol=12,data=0),matrix(nrow=7,ncol=24,data=0),matrix(nrow=7,ncol=24,data=0))
 #RL<-list(RL1,RL2,RL3,RL4)
# R2L<-list(RL1,RL2,RL3,RL4)
multreg<- function(fin,colm,len,sig){ #produces a column necessary for reglist, as a mean of optimizations
DATMAT<-matrix(nrow=7,ncol=len,data=0)
sigm<<-sig
for(i in 1:len){
set.seed(i)
ttt<-omkv(fin,colm,0.5,1)
DATMAT[,i]<-ttt[c(5:11)]
}
return(t(t(rowMeans(DATMAT))))
}
#multreg(2,18,50,VL[[2]][[3]][,4])
multR2new<-function(fin,colm,sig){#produces a column necessary for R2list
dv2<-c(abs(colm%%2-1),ceiling(colm/2)-1)
	coefnew(rep(fin,d),dv2)
	COL<-1:d
	for(i in 1:d){
	COL[i]<-rfunction[[fin]](i=i,n=2000,s=sig[i])
}
	return(COL)
}

noreglist<-function(){
sigm<<-rep(0,d) #Produces means of 200 omkv's with all Types
for(fin in 1:4){
	if(fin<3){nc<-12}else{nc<-24}
		rr<-1
		for(colm in 1:nc){
		if(fin!=2){
		cind<-ceiling(colm/3)
		v<-omkv(fin=fin,colm=colm,1,1)[c(5:11)]
		}else{
		if(colm<7){cind<-ceiling(colm/3);v<-omkv(fin=fin,colm=colm,1,1)[c(5:11)]
		}else{cind<-ceiling(colm/3);omkv(fin=fin,colm=colm,1,1)[c(5:11)]}
		}
		NRL[[fin]][,colm]<<-v
	}
}
}
#noreglist()
R2listnew<-function(){#Produces a list of all R^2-Values
source("Variance List New.txt")
for(fin in 1:4){
	for(rr in c(0.9,0.7,0.5)){
	RRR<-(-rr*10+11)/2 #RRR<-c(1,2,3)
		for(colm in 1:4){
		v<-multR2new(fin=fin,colm=colm,VLN[[fin]][[RRR]][,colm])
		R2LN[[fin]][[RRR]][,colm]<<-v
	}
}
}

}
#R2listnew()

funcc<-function(x,t1,t2){x^t1/(t2+x^t1)}
#funcm<-function(x,t1,t2){t1*x^t2}
manregad<-function(f,x0,n=3,itmax=50,eps=1e-10,sig=0){ #Manual determination of correct entries for adbudg
dax1<-sample(1:100,n,replace=T)
day<-f(dax1,0.5,50) +rnorm(n,sd=sqrt(sig))

fc1<-function(x,t1,t2){
(t2*log(x)*(x^(t1)))/((t2+x^t1)^2)
}
fc2<-function(x,t1,t2){
(-x^(t1))/(((x^t1)+t2)^2)
}
#fc1<-function(x,t1,t2){x^t2}
#fc2<-function(x,t1,t2){t1*(x^t2)*log(x)}
V<-matrix(ncol=2,byrow=F,nrow=n,data=c(fc1(dax1,x0[1],x0[2]),fc2(dax1,x0[1],x0[2])))
eta<-f(dax1,x0[1],x0[2])
QR<-qr(V)
Q<-qr.Q(QR);RR<-qr.R(QR)
Ri<-solve(RR)
z<-day-f(dax1,x0[1],x0[2])	
w<-t(Q)%*%z
d0<-Ri%*%w
ind<-0
fehler<-10
#while(ind<maxit|fehler>eps){
#}
#return(c(day,f(dax1,x0[1],x0[2]),f(dax1,x0[1]+d0[1],x0[2]+d0[2])))
print(x0)
print(x0+d0)
}
#manregad(funcc,c(0.6,60),sig=1,n=3)
# manregad(funcc,c(0.2,50),sig=1,n=3)
coeffeps<-function(fin,detvec){# calculates all 18 coefficient sets, returns the necessary ones

	#setting DB distribution
	if(detvec[2]!=1){ust<-rep(2500000, times=d)}else{ust<-rep(c(1,1,4,4), times=ceiling(d/4))[1:d]*1000000}
	#setting b/fc distribution
	if(detvec[1]!=1){if(d==8){b<-c(0.26,0.27,0.28,0.29,0.31,0.32,0.33,0.34)}else{b<-seq(from=0.26, to=0.34, by=0.08/(d-1))}
						}else{b<-seq(from=0.11, to=0.18,length=d);b[(floor(d/2)+1):d]<-seq(from=0.47, to=0.5, length=floor((d+1)/2))}
	fc<<-b
	b<<-b
	fcS<<-fc+1
	#setting Maximum Values 
	if(detvec[4]!=1){if(d==8){Mex<<-c(6.1,6.2,6.3,6.4,6.6,6.7,6.8,6.9)}else{Mex<<-seq(from=6.1, to=6.9, by=0.8/(d-1))}
			}else{Mex<<-rep(c(4.5,10),times=ceiling(d/2))[1:d]}
	Mad<<-Mex
	#calculating fmult coefficient a
	aid<-rep(R/d, times=d)
	a<<-ust/(aid^b)
	#calculating fex coefficient h
	h<<- log(1-ust/(Mex*1000000))/1000000
	#calculating fad coefficient g
	g<<-(Mad*1000000/ust -1 )*1000000^fc	
	#calculating fadS coefficient gS
	gS<<-(Mad*1000000/ust -1 )*1000000^fcS
		#calculating starting solution
	v<-quickopt(fin)
	if(detvec[3]==0){x1<-rep(R/d, times=d)}else{if(detvec[3]==1){x1<-v[[2]]}else{x1<-v[[1]]}}
		#return(list(x1,v[[5]],v[[3]],v[[4]]))
}
}

####NEWDRAW
nd<-function() {
pdf("Newcurve.pdf")
fiin=7
omkv(2,4,1,1)
drf=function(x) {fex(x,fiin)}
op<<-CLST[40,]
ov<<-checklist[40,1:8]
curve(drf,0,3e6,type="b",pch=20,axes=FALSE,xlab="x/10000",ylim=c(0,drf(3000000)),xlim=c(0,3000000),ylab="y/10000")
axis(1,at=seq(0,3000000,500000),labels=(c("0","50","100","150","200","250","300")))
axis(2,at=seq(0,4500000,500000),labels=(c("0","50","100","150","200","250","300","350","400","450")))

#points(x=op[fiin],y=ov[fiin],col="blue")
mi<-min(CLST[,fiin])
ma<-max(CLST[,fiin])
#points(x=mi,y=drf(mi),col="red")
#points(x=ma,y=drf(ma),col="red")
VEC<-TAB[[7]]$coefficients
curve(VEC[3]*x^2+VEC[2]*x+VEC[1],type="b",pch=5,add=TRUE)
legend(x=2000000,y=3000000,legend=c("true function","fitted parabola"),pch=c(20,5))
segments(mi,drf(mi),mi,drf(ma))
segments(ma,drf(mi),ma,drf(ma))
segments(mi,drf(mi),ma,drf(mi))
segments(mi,drf(ma),ma,drf(ma))
dev.off()
}

########
#Factors for ANOVA
set.seed(1)

algo<-c("Quadallo","Faust1","Faust2")
form<-c("Multip","ModExp","ADBUDG-K","ADBUDG-S")
pert<-c("0.9","0.7","0.5")
elas<-c("close","far")
satu<-c("equal","different")

#DATAFRAME<-data.frame(A=factor(algo),F=factor(form),P=factor(pert),E<-factor(elas),M<-factor(marg))
#for(i1 in algo){ #algorithms are evaluated simultaneously

if(dtest==1){
####### deterministic functions

DATA<-matrix(ncol=7,nrow=1,data=NA)
set.seed(1)
for(finn in 1:4){
	for(rr in c(0.9,0.7,0.5)){
		RRR<-(-rr*10+11)/2
	#	exex<-RRR*4+6 #Optimal exex-Borders. Need to find an estimate based on current x-Values
		for(dve in 0:1){
			for(dvs in 0:1){
				colm<-dve+2*dvs+1
				sigm<<-VLN[[finn]][[RRR]][,colm]
				for(seedvar in 1:5){
				print(c(finn,colm,RRR))
					yvec<-omkv(finn,colm,0.5,1)
					v1<-c(yvec[1],yvec[5],"Quadallo",finn,RRR,dve,dvs)
					v2<-c(yvec[c(2,6)],"Faust1",finn,RRR,dve,dvs)
					v3<-c(yvec[c(3,7)],"Faust2",finn,RRR,dve,dvs)
					v4<-c(yvec[c(4,8)],"Faust3",finn,RRR,dve,dvs)
					DATA<<-rbind(DATA,v1)
					DATA<<-rbind(DATA,v2)
					DATA<<-rbind(DATA,v3)
					DATA<<-rbind(DATA,v4)
				}
			}
		}
	}
}

dlen<-length(DATA)/7
ANOMAT<-DATA[2:dlen,]
yval1<-as.numeric(ANOMAT[,1])
yval2<-as.numeric(ANOMAT[,2])


if(norming==1){#leads to smaller R^2
	part=4
	inv=length(yval2)/part
	for(jump in 1:part){ 
		repl<-yval2[(1+inv*(jump-1)):(inv*jump)]
		repl<-repl/max(repl)
		yval2[(1+inv*(jump-1)):(inv*jump)]<-repl
	}
}

if(norming==2){#leads to extremely bad R^2
	inv=seedvar
	part=length(yval2)/inv
	for(jump in 1:part){ 
		repl<-yval2[(1+inv*(jump-1)):(inv*jump)]
		repl<-repl/max(repl)
		yval2[(1+inv*(jump-1)):(inv*jump)]<-repl
	}
}



#yval1<-yval1/max(yval1)  #No normalization of the optimalities
yval1<-yval1/100 # back into [0,1], might not include 1
yval2<-yval2/max(yval2) # will include 1
algo<-ANOMAT[,3]
form<-as.numeric(ANOMAT[,4])
pert<-as.numeric(ANOMAT[,5])
elas<-as.numeric(ANOMAT[,6])
satu<-as.numeric(ANOMAT[,7])
ANODF1<-data.frame(optimality=yval1,Algorithm=relevel(factor(algo),ref="Quadallo"),Form=factor(form),Perturbations=factor(pert), Elasticities=factor(elas), Saturation=factor(satu))
ANODF2<-data.frame(Returns=yval2,Algorithm=relevel(factor(algo),ref="Quadallo"),Form=factor(form),Perturbations=factor(pert), Elasticities=factor(elas), Saturation=factor(satu))
RES<-lm(Returns~Algorithm*Form*Perturbations*Elasticities*Saturation, data=ANODF2)
RES2<-lm(optimality~Algorithm*Form*Perturbations*Elasticities*Saturation, data=ANODF1)
PURE1<-lm(Returns~Algorithm+Form+Perturbations+Elasticities+Saturation, data=ANODF2)
PURE2<-lm(optimality~Algorithm+Form+Perturbations+Elasticities+Saturation, data=ANODF1)
PAIR1<-lm(Returns~Algorithm+Form+Perturbations+Elasticities+Saturation+Algorithm:Form+Algorithm:Perturbations+Algorithm:Elasticities+Algorithm:Saturation+Form:Perturbations+Form:Elasticities+Form:Saturation+Perturbations:Elasticities+Perturbations:Saturation+Elasticities:Saturation, data=ANODF2)
PAIR2<-lm(optimality~Algorithm+Form+Perturbations+Elasticities+Saturation+Algorithm:Form+Algorithm:Perturbations+Algorithm:Elasticities+Algorithm:Saturation+Form:Perturbations+Form:Elasticities+Form:Saturation+Perturbations:Elasticities+Perturbations:Saturation+Elasticities:Saturation, data=ANODF1)
#PURE1<-lm(Returns~Algorithm+Form+Elasticities+Saturation, data=ANODF2)
summary(PURE1)
RP1<-lm(Returns~Algorithm:Form+Algorithm:Perturbations+Algorithm:Elasticities+Algorithm:Saturation+Form:Perturbations+Form:Elasticities+Form:Saturation+Perturbations:Elasticities+Perturbations:Saturation+Elasticities:Saturation, data=ANODF2)
testreg<-lm(Returns~Algorithm:Form+Algorithm:Perturbations+Algorithm:Elasticities+Algorithm:Saturation+Form:Perturbations+Form:Elasticities+Form:Saturation+Perturbations:Elasticities+Perturbations:Saturation+Elasticities:Saturation -1, data=ANODF2)
PINT1<-lm(Returns~Algorithm+Form+Perturbations+Elasticities+Saturation+Algorithm:Form+Algorithm:Perturbations+Algorithm:Elasticities+Algorithm:Saturation, data=ANODF2)
PINT2<-lm(optimality~Algorithm+Form+Perturbations+Elasticities+Saturation+Algorithm:Form+Algorithm:Perturbations+Algorithm:Elasticities+Algorithm:Saturation, data=ANODF1)


ANODF2<-data.frame(Returns=yval2,Algorithm=relevel(factor(algo),ref="Faust1"),Form=factor(form),Perturbations=factor(pert), Elasticities=factor(elas), Saturation=factor(satu))
#####Contrasts/F-Tests
ZM<-rep(0,32)
K1<-t(as.matrix(ZM))
K2<-t(as.matrix(ZM))
K3<-t(as.matrix(ZM))
KADB<-t(as.matrix(ZM))

##Unser Alg: 0+0+5+0+6+0+7+0+8+0+9+0+10+0+11
##Alg1:0+1+0+1+12+0+1+15+0+1+18+0+1+21+0+1+24+0+1+27+0+1+30
## Differenz: 5+6+7+8+9+10+11-12-15-18-21-24-27-30-8*1
K1[c(2,5,6,7,8,9,10,11,12,15,18,21,24,27,30)]<-c(-8,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1)
K2[c(3,5,6,7,8,9,10,11,13,16,19,22,25,28,31)]<-c(-8,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1)
K3[c(4,5,6,7,8,9,10,11,14,17,20,23,26,29,32)]<-c(-8,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1)
KADB[c(5,7,12,13,14,18,19,20)]<-c(1,-1,1,1,1,-1,-1,-1)
betpint1<-PINT1$coefficients
betpint2<-PINT1$coefficients
K1%*%betpint1
K2%*%betpint1
K3%*%betpint1
K1%*%betpint2
K2%*%betpint2
K3%*%betpint2

linearHypothesis(PINT1,hypothesis.matrix=K1)




##### Perttest, nur perturbations im Modell.
PERT1<-lm(Returns~Algorithm*Perturbations,data=ANODF2 )
PERT2<-lm(optimality~Algorithm*Perturbations,data=ANODF1 )
summary(PERT1)
K<-matrix(nrow=1,data=c(0,-4,0,0,0,0,2,0,0,2,0,0))
C<-matrix(nrow=1,data=c(0,-3,0,0,0,0,-1,0,0,-1,0,0))
t <- glht(PERT1, linfct = K)
summary(t)
betav<-PERT1$coefficients
REEEEE<-PERT1[[2]]
estvar<-sum(REEEEE^2)
DESMAT<-matrix(nrow=960,ncol=12,data=NA)
for(i in 1:960){
VV<-PERT1[[13]][i,]
DESMAT[i,1]<-1
DESMAT[i,2]<-if(VV[2]=="Faust1"){1}else{0}
DESMAT[i,3]<-if(VV[2]=="Faust2"){1}else{0}
DESMAT[i,4]<-if(VV[2]=="Faust3"){1}else{0}
DESMAT[i,5]<-if(VV[3]==2){1}else{0}
DESMAT[i,6]<-if(VV[3]==3){1}else{0}
DESMAT[i,7]<-DESMAT[i,2]*DESMAT[i,5]
DESMAT[i,8]<-DESMAT[i,3]*DESMAT[i,5]
DESMAT[i,9]<-DESMAT[i,4]*DESMAT[i,5]
DESMAT[i,10]<-DESMAT[i,2]*DESMAT[i,6]
DESMAT[i,11]<-DESMAT[i,3]*DESMAT[i,6]
DESMAT[i,12]<-DESMAT[i,4]*DESMAT[i,6]
}
FST<-t(C%*%betav)%*%solve(C%*%solve(t(DESMAT)%*%DESMAT)%*%t(C))%*%(C%*%betav)
pf(FST,1,960-11-1)
DATAAA<-linearHypothesis(PERT1,hypothesis.matrix=C)

betgigvec<-RES2$coefficients
KMGV<-c(rep(0,26),rep(1,6),rep(0,3),rep(1,5),rep(0,3),rep(1,149))
KMGVP<-rep(0,sum(KMGV));kt=1
for(i in 1:192){if(KMGV[i]==1){KMGVP[kt]<-i;kt=kt+1}}

KMG<-matrix(ncol=192,nrow=160,data=0)
for(i in 1:160){KMG[i,KMGVP[i]]<-1}
linearHypothesis(RES2,hypothesis.matrix=KMG)
##VGL Pair&ALGPAIR
KP<-matrix(nrow=17,ncol=49,data=0)
for(i in 1:17){KP[i,i+32]<-1}
linearHypothesis(PAIR2,hypothesis.matrix=KP)
###VGL Algo und Faust1 im RES Modell
#Diff:(falsch)

erv<-c(-1,rep(1,7),rep(-1,5),rep(1,6),-1,rep(1,5),-1,rep(1,6),rep(-1,11),rep(1,6),rep(-1,5),rep(1,6),-1,rep(1,5),rep(-1,17),rep(1,6),rep(-1,6))
ervgew<-c(-48,rep(12,3),16,16,24,24,rep(-12,3),rep(-16,2),rep(4,6),-24,rep(6,3),rep(8,2),-24,rep(6,3),rep(8,2),rep(12,1),rep(-4,6),rep(-6,3),rep(-8,2),rep(2,6),rep(-6,3),rep(-8,2),rep(2,6),rep(-12,1),rep(3,3),rep(4,2),rep(-2,12),rep(-3,3),rep(-4,2),rep(1,6),rep(-1,6))


KMa<-rep(0,192)
KMa2<-rep(0,192)
KMa3<-rep(0,192)
KMa[c(2,5,6,7,8,9,10,11,12,15,18,21,24,27,28,29,30,31,32,33, 36,37,38,39,40,41,44,45,46,47,48,49,50,53,56,59,62,65,68,71,74,77,80,83,84,85,86,87,88,89,92,95,98,101,104,105,106,107,108,109,110,113,114,115,116,117,118,121,124,127,130,133,136,139,142,145,148,151,154,157,160,163,166,169,170,171,172,173,174,175,178,181,184,187,190)]<-ervgew
KMa2[c(3,5,6,7,8,9,10,11,13,16,19,22,25,27,28,29,30,31,32,34, 36,37,38,39,40,42,44,45,46,47,48,49,51,54,57,60,63,66,69,72,75,78,81,83,84,85,86,87,88,90,93,96,99,102,104,105,106,107,108,109,111,113,114,115,116,117,119,122,125,128,131,134,137,140,143,146,149,152,155,158,161,164,167,169,170,171,172,173,174,176,179,182,185,188,191)]<-ervgew
KMa3[c(4,5,6,7,8,9,10,11,14,17,20,23,26,27,28,29,30,31,32,35, 36,37,38,39,40,43,44,45,46,47,48,49,52,55,58,61,64,67,70,73,76,79,82,83,84,85,86,87,88,91,94,97,100,103,104,105,106,107,108,109,112,113,114,115,116,117,120,123,126,129,132,135,138,141,144,147,150,153,156,159,162,165,168,169,170,171,172,173,174,177,180,183,186,189,192)]<-ervgew

betfull<-RES$coefficients
betful2<-RES2$coefficients

Ergmat<-matrix(ncol=6,nrow=2,data=0)

Ergmat[1,1]<-KMa%*%betfull
Ergmat[2,1]<-linearHypothesis(RES,hypothesis.matrix=KMa)[2,6]
Ergmat[1,2]<-KMa2%*%betfull
Ergmat[2,2]<-linearHypothesis(RES,hypothesis.matrix=KMa2)[2,6]
Ergmat[1,3]<-KMa3%*%betfull
Ergmat[2,3]<-linearHypothesis(RES,hypothesis.matrix=KMa3)[2,6]
Ergmat[1,4]<-KMa%*%betful2
Ergmat[2,4]<-linearHypothesis(RES2,hypothesis.matrix=KMa)[2,6]
Ergmat[1,5]<-KMa2%*%betful2
Ergmat[2,5]<-linearHypothesis(RES2,hypothesis.matrix=KMa2)[2,6]
Ergmat[1,6]<-KMa3%*%betful2
Ergmat[2,6]<-linearHypothesis(RES2,hypothesis.matrix=KMa3)[2,6]





}else{
######## random functions


DATA2<-matrix(ncol=7,nrow=1,data=NA)
set.seed(1)
flag<-rep(0,8)
	for(rr in c(0.9,0.7,0.5)){
		RRR<-(-rr*10+11)/2
	#	exex<-RRR*4+6 #Optimal exex-Borders. Need to find an estimate based on current x-Values
		for(dve in 0:1){
			for(dvs in 0:1){
				
				colm<-dve+2*dvs+1

				for(seedvar in 1:10){
				dran<-sample(4:20,1)
				#finn<-sample(1:4,dran,replace=1)
				for(finnk in 1:4){
				finn<-rep(finnk,times=dran)
					for(tin in 1:dran){
					ting<-((tin-1)%%8)+1
					sigm<<-VLN[[finn[ting]]][[RRR]][ting,colm]
					}
				
				
				print(c(finn,colm,RRR))
					yvec<-omkv(finn,colm,0.5,1)
					v1<-c(yvec[1],yvec[5],"Quadallo",finnk,RRR,dve,dvs)
					v2<-c(yvec[c(2,6)],"Faust1",finnk,RRR,dve,dvs)
					v3<-c(yvec[c(3,7)],"Faust2",finnk,RRR,dve,dvs)
					v4<-c(yvec[c(4,8)],"Faust3",finnk,RRR,dve,dvs)
					DATA2<<-rbind(DATA2,v1)
					DATA2<<-rbind(DATA2,v2)
					DATA2<<-rbind(DATA2,v3)
					DATA2<<-rbind(DATA2,v4)
					}
				}
			}
		}
	}
dlen<-length(DATA2)/7
ANOMAT<-DATA2[2:dlen,]
yval1<-as.numeric(ANOMAT[,1])
yval2<-as.numeric(ANOMAT[,2])

if(norming==1){#leads to smaller R^2
	part=4
	inv=length(yval2)/part
	for(jump in 1:part){ 
		repl<-yval2[(1+inv*(jump-1)):(inv*jump)]
		repl<-repl/max(repl)
		yval2[(1+inv*(jump-1)):(inv*jump)]<-repl
	}
}
if(norming==2){#leads to extremely bad R^2
	inv=seedvar
	part=length(yval2)/inv
	for(jump in 1:part){ 
		repl<-yval2[(1+inv*(jump-1)):(inv*jump)]
		repl<-repl/max(repl)
		yval2[(1+inv*(jump-1)):(inv*jump)]<-repl
	}
}


yval1<-yval1/max(yval1)
yval2<-yval2/max(yval2)
algo<-ANOMAT[,3]
form<-as.numeric(ANOMAT[,4])
pert<-as.numeric(ANOMAT[,5])
elas<-as.numeric(ANOMAT[,6])
satu<-as.numeric(ANOMAT[,7])
ANODF1<-data.frame(optimality=yval1,Algorithm=relevel(factor(algo),ref="Quadallo"),Form=factor(form),Perturbations=factor(pert), Elasticities=factor(elas), Saturation=factor(satu))
ANODF2<-data.frame(Returns=yval2,Algorithm=relevel(factor(algo),ref="Quadallo"),Form=factor(form),Perturbations=factor(pert), Elasticities=factor(elas), Saturation=factor(satu))
RES<-lm(Returns~Algorithm*Form*Perturbations*Elasticities*Saturation, data=ANODF2)
RES2<-lm(optimality~Algorithm*Form*Perturbations*Elasticities*Saturation, data=ANODF1)
PURE1<-lm(Returns~Algorithm+Form+Perturbations+Elasticities+Saturation, data=ANODF2)
PURE2<-lm(optimality~Algorithm+Form+Perturbations+Elasticities+Saturation, data=ANODF1)
PAIR1<-lm(Returns~Algorithm+Form+Perturbations+Elasticities+Saturation+Algorithm:Form+Algorithm:Perturbations+Algorithm:Elasticities+Algorithm:Saturation+Form:Perturbations+Form:Elasticities+Form:Saturation+Perturbations:Elasticities+Perturbations:Saturation+Elasticities:Saturation, data=ANODF2)
PAIR2<-lm(optimality~Algorithm+Form+Perturbations+Elasticities+Saturation+Algorithm:Form+Algorithm:Perturbations+Algorithm:Elasticities+Algorithm:Saturation+Form:Perturbations+Form:Elasticities+Form:Saturation+Perturbations:Elasticities+Perturbations:Saturation+Elasticities:Saturation, data=ANODF1)
#PURE1<-lm(Returns~Algorithm+Form+Elasticities+Saturation, data=ANODF2)

#PAIRrc1<-lmres(residual_centering=1,Returns~Algorithm+Form+Perturbations+Elasticities+Saturation+Algorithm:Form+Algorithm:Perturbations+Algorithm:Elasticities+Algorithm:Saturation+Form:Perturbations+Form:Elasticities+Form:Saturation+Perturbations:Elasticities+Perturbations:Saturation+Elasticities:Saturation, data=ANODF2)
summary(PURE1)
RP1<-lm(Returns~Algorithm:Form+Algorithm:Perturbations+Algorithm:Elasticities+Algorithm:Saturation+Form:Perturbations+Form:Elasticities+Form:Saturation+Perturbations:Elasticities+Perturbations:Saturation+Elasticities:Saturation, data=ANODF2)
testreg<-lm(Returns~Algorithm:Form+Algorithm:Perturbations+Algorithm:Elasticities+Algorithm:Saturation+Form:Perturbations+Form:Elasticities+Form:Saturation+Perturbations:Elasticities+Perturbations:Saturation+Elasticities:Saturation -1, data=ANODF2)
PINT1<-lm(Returns~Algorithm+Form+Perturbations+Elasticities+Saturation+Algorithm:Form+Algorithm:Perturbations+Algorithm:Elasticities+Algorithm:Saturation, data=ANODF2)
PINT2<-lm(optimality~Algorithm+Form+Perturbations+Elasticities+Saturation+Algorithm:Form+Algorithm:Perturbations+Algorithm:Elasticities+Algorithm:Saturation, data=ANODF1)


ANODF2<-data.frame(Returns=yval2,Algorithm=relevel(factor(algo),ref="Faust1"),Form=factor(form),Perturbations=factor(pert), Elasticities=factor(elas), Saturation=factor(satu))
#####Contrasts/F-Tests
{
#ZM<-rep(0,32)
K1<-t(as.matrix(ZM))
K2<-t(as.matrix(ZM))
K3<-t(as.matrix(ZM))
KADB<-t(as.matrix(ZM))

##Unser Alg: 0+0+5+0+6+0+7+0+8+0+9+0+10+0+11
##Alg1:0+1+0+1+12+0+1+15+0+1+18+0+1+21+0+1+24+0+1+27+0+1+30
## Differenz: 5+6+7+8+9+10+11-12-15-18-21-24-27-30-8*1
K1[c(2,5,6,7,8,9,10,11,12,15,18,21,24,27,30)]<-c(-8,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1)
K2[c(3,5,6,7,8,9,10,11,13,16,19,22,25,28,31)]<-c(-8,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1)
K3[c(4,5,6,7,8,9,10,11,14,17,20,23,26,29,32)]<-c(-8,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1)
KADB[c(5,7,12,13,14,18,19,20)]<-c(1,-1,1,1,1,-1,-1,-1)
betpint1<-PINT1$coefficients
betpint2<-PINT1$coefficients
K1%*%betpint1
K2%*%betpint1
K3%*%betpint1
K1%*%betpint2
K2%*%betpint2
K3%*%betpint2

linearHypothesis(PINT1,hypothesis.matrix=K1)




##### Perttest, nur perturbations im Modell.
PERT1<-lm(Returns~Algorithm*Perturbations,data=ANODF2 )
PERT2<-lm(optimality~Algorithm*Perturbations,data=ANODF1 )
summary(PERT1)
K<-matrix(nrow=1,data=c(0,-4,0,0,0,0,2,0,0,2,0,0))
C<-matrix(nrow=1,data=c(0,-3,0,0,0,0,-1,0,0,-1,0,0))
t <- glht(PERT1, linfct = K)
summary(t)
betav<-PERT1$coefficients
REEEEE<-PERT1[[2]]
estvar<-sum(REEEEE^2)
DESMAT<-matrix(nrow=960,ncol=12,data=NA)
for(i in 1:960){
VV<-PERT1[[13]][i,]
DESMAT[i,1]<-1
DESMAT[i,2]<-if(VV[2]=="Faust1"){1}else{0}
DESMAT[i,3]<-if(VV[2]=="Faust2"){1}else{0}
DESMAT[i,4]<-if(VV[2]=="Faust3"){1}else{0}
DESMAT[i,5]<-if(VV[3]==2){1}else{0}
DESMAT[i,6]<-if(VV[3]==3){1}else{0}
DESMAT[i,7]<-DESMAT[i,2]*DESMAT[i,5]
DESMAT[i,8]<-DESMAT[i,3]*DESMAT[i,5]
DESMAT[i,9]<-DESMAT[i,4]*DESMAT[i,5]
DESMAT[i,10]<-DESMAT[i,2]*DESMAT[i,6]
DESMAT[i,11]<-DESMAT[i,3]*DESMAT[i,6]
DESMAT[i,12]<-DESMAT[i,4]*DESMAT[i,6]
}
FST<-t(C%*%betav)%*%solve(C%*%solve(t(DESMAT)%*%DESMAT)%*%t(C))%*%(C%*%betav)
pf(FST,1,960-11-1)
DATAAA<-linearHypothesis(PERT1,hypothesis.matrix=C)

betgigvec<-RES2$coefficients
KMGV<-c(rep(0,26),rep(1,6),rep(0,3),rep(1,5),rep(0,3),rep(1,149))
KMGVP<-rep(0,sum(KMGV));kt=1
for(i in 1:192){if(KMGV[i]==1){KMGVP[kt]<-i;kt=kt+1}}

KMG<-matrix(ncol=192,nrow=160,data=0)
for(i in 1:160){KMG[i,KMGVP[i]]<-1}
linearHypothesis(RES2,hypothesis.matrix=KMG)
##VGL Pair&ALGPAIR
KP<-matrix(nrow=17,ncol=49,data=0)
for(i in 1:17){KP[i,i+32]<-1}
linearHypothesis(PAIR2,hypothesis.matrix=KP)
###VGL Algo und Faust1 im RES Modell
#Diff:(falsch)
KMa<-rep(0,192)
KMa2<-rep(0,192)
KMa3<-rep(0,192)
KMa[c(2,5,6,7,8,9,10,11,12,15,18,21,24,27,28,29,30,31,32,33, 36,37,38,39,40,41,44,45,46,47,48,49,50,53,56,59,62,65,68,71,74,77,80,83,84,85,86,87,88,89,92,95,98,101,104,105,106,107,108,109,110,113,114,115,116,117,118,121,124,127,130,133,136,139,142,145,148,151,154,157,160,163,166,169,170,171,172,173,174,175,178,181,184,187,190)]<-ervgew
KMa2[c(3,5,6,7,8,9,10,11,13,16,19,22,25,27,28,29,30,31,32,34, 36,37,38,39,40,42,44,45,46,47,48,49,51,54,57,60,63,66,69,72,75,78,81,83,84,85,86,87,88,90,93,96,99,102,104,105,106,107,108,109,111,113,114,115,116,117,119,122,125,128,131,134,137,140,143,146,149,152,155,158,161,164,167,169,170,171,172,173,174,176,179,182,185,188,191)]<-ervgew
KMa3[c(4,5,6,7,8,9,10,11,14,17,20,23,26,27,28,29,30,31,32,35, 36,37,38,39,40,43,44,45,46,47,48,49,52,55,58,61,64,67,70,73,76,79,82,83,84,85,86,87,88,91,94,97,100,103,104,105,106,107,108,109,112,113,114,115,116,117,120,123,126,129,132,135,138,141,144,147,150,153,156,159,162,165,168,169,170,171,172,173,174,177,180,183,186,189,192)]<-ervgew

KMa%*%betfull
KMa2%*%betfull
KMa3%*%betfull
KMa%*%betful2
KMa2%*%betful2
KMa3%*%betful2


erv<-c(-1,rep(1,7),rep(-1,5),rep(1,6),-1,rep(1,5),-1,rep(1,6),rep(-1,11),rep(1,6),rep(-1,5),rep(1,6),-1,rep(1,5),rep(-1,17),rep(1,6),rep(-1,6))
ervgew<-c(-48,rep(12,3),16,16,24,24,rep(-12,3),rep(-16,2),rep(4,6),-24,rep(6,3),rep(8,2),-24,rep(6,3),rep(8,2),rep(12,1),rep(-4,6),rep(-6,3),rep(-8,2),rep(2,6),rep(-6,3),rep(-8,2),rep(2,6),rep(-12,1),rep(3,3),rep(4,2),rep(-2,12),rep(-3,3),rep(-4,2),rep(1,6),rep(-1,6))
}
#dlen<-length(DATA2)/6
#ANOMAT<-DATA2[2:dlen,]
#yval1<-as.numeric(ANOMAT[,1])
#yval2<-as.numeric(ANOMAT[,2])
#yval1<-yval1/max(yval1)
#yval2<-yval2/max(yval2)
#algo<-ANOMAT[,3]
#pert<-as.numeric(ANOMAT[,4])
#elas<-as.numeric(ANOMAT[,5])
#satu<-as.numeric(ANOMAT[,6])
#ANODF1<-data.frame(optimality=yval1,Algorithm=relevel(factor(algo),ref="Quadallo"),Perturbations=factor(pert), Elasticities=factor(elas), Saturation=factor(satu))
#ANODF2<-data.frame(Returns=yval2,Algorithm=relevel(factor(algo),ref="Quadallo"),Perturbations=factor(pert), Elasticities=factor(elas), Saturation=factor(satu))
#RES<-lm(Returns~Algorithm*Perturbations*Elasticities*Saturation, data=ANODF2)
#RES2<-lm(optimality~Algorithm*Perturbations*Elasticities*Saturation, data=ANODF1)
#PURE1<-lm(Returns~Algorithm+Perturbations+Elasticities+Saturation, data=ANODF2)
#PURE2<-lm(optimality~Algorithm+Perturbations+Elasticities+Saturation, data=ANODF1)
#PAIR1<-lm(Returns~Algorithm+Perturbations+Elasticities+Saturation+Algorithm:Perturbations+Algorithm:Elasticities+Algorithm:Saturation+Perturbations:Elasticities+Perturbations:Saturation+Elasticities:Saturation, data=ANODF2)
#PAIR2<-lm(optimality~Algorithm+Perturbations+Elasticities+Saturation+Algorithm:Perturbations+Algorithm:Elasticities+Algorithm:Saturation+Perturbations:Elasticities+Perturbations:Saturation+Elasticities:Saturation, data=ANODF1)
#PURE1<-lm(Returns~Algorithm+Form+Elasticities+Saturation, data=ANODF2)
#summary(PURE1)

}






#dv2<-c(abs(colm%%2-1),ceiling(colm/2)-1)




###### old index runners
#	for(i2 in form){
#		for(i3 in pert){
#			for(i4 in elas){
#				for(i5 in satu){
#					for(seedvar in 1:3){
#					sigm<<-
#					ci<-dv()
#					yvec<-omkv(fin,ci,0.5,1)
#					DATA<<-rbind(DATA,c(yvec[c(1,5)],"Quadallo",i2,i3,i4,i5))
#					DATA<<-rbind(DATA,c(yvec[c(2,6)],"Faust1",i2,i3,i4,i5))
#					DATA<<-rbind(DATA,c(yvec[c(3,7)],"Faust2",i2,i3,i4,i5))
#					
#					}
#				
#				}
#			
#			}
#		
#		
#		}
#	
#	
#	}
#
#}


 ### deleted things


#functest<-function(f){#function index:1,2,3,4
#if(min(deparse(f)==deparse(fmult))){return(1)}
#if(min(deparse(f)==deparse(fex))){return(2)}
#if(min(deparse(f)==deparse(fad))){return(3)}
#print("invalid function")
#return(NaN)
#}

#feval<- function(v) {
#	val<-0
#for( i in 1:d){
#	val<-val+fex(v[i],i)
#	}
#		return(val)
#}

#sigm<-VL[[1]][[1]][,1]
#omkv(1,1,0.5,1)

#coeff<-function(fin,detvec){# calculates all 18 coefficient sets, returns the necessary ones
#
#	#setting DB distribution
#	if(detvec[2]!=1){ust<-rep(2500000, times=d)}else{ust<-rep(c(1,1,4,4), times=ceiling(d/4))[1:d]*1000000}
#	#setting b/fc distribution
#	if(detvec[1]!=1){if(d==8){b<-c(0.26,0.27,0.28,0.29,0.31,0.32,0.33,0.34)}else{b<-seq(from=0.26, to=0.34, by=0.08/(d-1))}
#						}else{b<-seq(from=0.11, to=0.18,length=d);b[(floor(d/2)+1):d]<-seq(from=0.47, to=0.5, length=floor((d+1)/2))}
#	fc<<-b
#	b<<-b
#	fcS<<-fc+1
#	#setting Maximum Values 
#	if(detvec[4]!=1){if(d==8){Mex<<-c(6.1,6.2,6.3,6.4,6.6,6.7,6.8,6.9)}else{Mex<<-seq(from=6.1, to=6.9, by=0.8/(d-1))}
#			}else{Mex<<-rep(c(4.5,10),times=ceiling(d/2))[1:d]}
#	Mad<<-Mex
#	#calculating fmult coefficient a
#	aid<-rep(R/d, times=d)
#	a<<-ust/(aid^b)
#	#calculating fex coefficient h
#	h<<- log(1-ust/(Mex*1000000))/1000000
#	#calculating fad coefficient g
#	g<<-(Mad*1000000/ust -1 )*1000000^fc	
#	#calculating fadS coefficient gS
#	gS<<-(Mad*1000000/ust -1 )*1000000^fcS
#		#calculating starting solution
#	v<-quickopt(fin)
#	if(detvec[3]==0){x1<-rep(R/d, times=d)}else{if(detvec[3]==1){x1<-v[[2]]}else{x1<-v[[1]]}}
#		return(list(x1,v[[5]],v[[3]],v[[4]]))
#}

#VL1<-list(matrix(nrow=8,ncol=4,data=0),matrix(nrow=8,ncol=4,data=0),matrix(nrow=8,ncol=4,data=0))
#VL2<-list(matrix(nrow=8,ncol=4,data=0),matrix(nrow=8,ncol=4,data=0),matrix(nrow=8,ncol=4,data=0))
#VL3<-list(matrix(nrow=8,ncol=8,data=0),matrix(nrow=8,ncol=8,data=0),matrix(nrow=8,ncol=8,data=0))
#VL4<-list(matrix(nrow=8,ncol=8,data=0),matrix(nrow=8,ncol=8,data=0),matrix(nrow=8,ncol=8,data=0))
#VL<-list(VL1,VL2,VL3,VL4)

#varlist<-function(){ #produces variance list 
#for(fin in 1:4){ #1:4
#	for(rr in c(0.9,0.7,0.5)){ #c(0.9,0.7,0.5)
#			RRR<-(-rr*10+11)/2
#			if(fin<3){nc<-4}else{nc<-8}
#				for(colm in 1:nc){ #1:nc
#				if(fin==2&&colm>2){
#				omkv(fin,3*colm+6,rr,0)
#				VL[[fin]][[RRR]][,colm]<<-sigm
#				}else{
#			omkv(fin,3*colm,rr,0)
#			VL[[fin]][[RRR]][,colm]<<-sigm
#			}
#			dump("VL",file="Variance List.txt")
#			}
#		}
#	}	
#}
#varlist()


#reglist<-function(){ #Produces means of 200 omkv's with all Types
#for(fin in 1:4){
#	if(fin<3){nc<-12}else{nc<-24}
#	for(rr in c(0.5,0.7,0.9)){
#	RRR<-(-rr*10+11)/2
#		for(colm in 1:nc){
#		source("Variance List.txt")
#		if(fin!=2){
#		cind<-ceiling(colm/3)
#		v<-multreg(fin=fin,colm=colm,200,VL[[fin]][[RRR]][,cind])
#		}else{
#		if(colm<7){cind<-ceiling(colm/3);v<-multreg(fin=fin,colm=colm,200,VL[[fin]][[RRR]][,cind])
#		}else{cind<-ceiling(colm/3);v<-multreg(fin=fin,colm=colm+6,200,VL[[fin]][[RRR]][,cind])}
#		}
#		RL[[fin]][[RRR]][,colm]<<-v
#	}
#}
#}
#}
#reglist()

#R2list<-function(){#Produces a list of all R^2-Values
#source("Variance List.txt")
#for(fin in 1:4){
#	if(fin<3){nc<-12}else{nc<-24}
#	for(rr in c(0.9,0.7,0.5)){
#	RRR<-(-rr*10+11)/2 #RRR<-c(1,2,3)
#		for(cind in 1:(nc/3)){
#		colm<-cind*3
#		if(fin!=2){
#		v<-multR2(fin=fin,colm=colm,VL[[fin]][[RRR]][,cind])
#		}else{
#		if(colm<7){cind<-ceiling(colm/3);v<-multR2(fin=fin,colm=colm,VL[[fin]][[RRR]][,cind])
#		}else{cind<-ceiling(colm/3);v<-multR2(fin=fin,colm=colm+6,VL[[fin]][[RRR]][,cind])}
#		}
#		
#		R2L[[fin]][[RRR]][,cind]<<-v
#	}
#}
#}


#R2list()