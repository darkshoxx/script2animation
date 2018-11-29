Dx=1:8
Dy=rep(0,8)
testy=c(4,6,8,2,7,1,3,5)
HitList=matrix(nrow=0,ncol=8)
colnames(HitList)=c("y1","y2","y3","y4","y5","y6","y7","y8")
addtolist=function(Dy){
	HitList<<-rbind(Dy,HitList)
}

showboard=function(Dy){
mat=diag(8)[Dy,]
colnames(mat)=1:8
rownames(mat)=c("A","B","C","D","E","F","G","H")
print(mat)
}

onboard=0

AddQueen=function(onboard){
flag=0
	if(onboard>8){print("too many Queens")}
	if(onboard<0){print("too few Queens")}
	if(onboard==8){
		addtolist(Dy)
		print(Dy)
		showboard(Dy)
		
		}else{
		for(i in 1:8){
		#print(Dy)
		if(tryY(onboard,i)){
		

			#onboard=onboard+1
			flag=Dy[onboard+1]
			Dy[onboard+1]<<-i
			#print(c(onboard,flag,Dy[onboard],i))
			AddQueen(onboard+1)
			#print(onboard)
			}
		}
	}
}


tryY=function(onboard,i){
Testresult=TRUE
if(onboard==0){return(Testresult)}
diag1=Dx+Dy
diag2=Dx-Dy
if(max(i==(Dy[1:onboard]))){Testresult=FALSE}
if(max((onboard+1+i)==(diag1[1:onboard]))){Testresult=FALSE}
if(max((onboard+1-i)==(diag2[1:onboard]))){Testresult=FALSE}
#print(c(onboard,Testresult))
return(Testresult)
}

#Dy=c(1,2,5,4,6,3,5,8)
#tryY(4,6)

AddQueen(0)
