Elx=function(t) {2.5*cos(t)}
Ely=function(t) {1.5+0.75*sin(t)}
ELC=function(t) {v=c(t,t,t+pi/2,t+pi/2,t+pi,t+pi,t+3*pi/2,t+3*pi/2)
xv<-Elx(v)
yv<-Ely(v)*c(1,-1,1,-1,1,-1,1,-1)
return(c(xv,yv))}

for(t in seq(0,pi*2,0.01)){
curve(1*x,-2.5,2.5,col="white")
#XV<-c(-2,-2,-1,-1,1,1,2,2)
#YV<-c(-2,2,-1,1,-1,1,-2,2)

V<-ELC(t)
XV<-V[1:8];YV<-V[9:16]
points(XV,YV)
segments(XV[1],YV[1],XV[2],YV[2])
segments(XV[3],YV[3],XV[4],YV[4])
segments(XV[5],YV[5],XV[6],YV[6])
segments(XV[7],YV[7],XV[8],YV[8])

segments(XV[1],YV[1],XV[7],YV[7])
segments(XV[3],YV[3],XV[5],YV[5])
segments(XV[4],YV[4],XV[6],YV[6])
segments(XV[2],YV[2],XV[8],YV[8])

segments(XV[1],YV[1],XV[3],YV[3])
segments(XV[2],YV[2],XV[4],YV[4])
segments(XV[5],YV[5],XV[7],YV[7])
segments(XV[6],YV[6],XV[8],YV[8])
}
