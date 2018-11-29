d=1e1
M<-matrix(ncol=d,nrow=d,data=0)
count<-function(x,y) {M[x+1,y]+M[x+1,y+1]+M[x,y+1]+M[x-1,y+1]+M[x-1,y]+M[x-1,y-1]+M[x,y-1]+M[x+1,y-1]}

step=function(X){MEM<-X
for(xc in (2):(d-1)){
	for(yc in (2):(d-1)){
		numb<-count(xc,yc)
		if(numb==3){MEM[xc,yc]<-1}else{
			if(X[xc,yc]==1&numb!=2){MEM[xc,yc]<-0}}
	}
}

return(MEM)
}

 (M<-step(M))
#M<-matrix(ncol=d,nrow=d,data=sample(c(0,1),d^2,replace=1))

#M[4,4]=1;M[4,5]=1;M[4,6]=1 #Blinker
M[4,7]=1;M[4,8]=1;M[4,9]=1;M[3,7]=1;M[2,8]=1 #Gleiter