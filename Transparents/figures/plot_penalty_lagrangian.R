# plot lagrangian and external penalty function
f <- function(x){return((x-2)^2)}
g <- function(x){return(4-x)}
p <- function(x,r=1){
  res<-f(x)+r*(max(0,g(x)))^2
  return(res)
}
lag <- function(x,l=4){return(f(x)+l*g(x))}

xmin<--1
xmax<-7
x<-seq(from=xmin,to=xmax,length.out=1000)
fs<-f(x)
lz <- lag(x)
pz <- apply(X = as.matrix(x),FUN = p,MARGIN = 1)
pz10 <- apply(X = as.matrix(x),FUN = p,MARGIN = 1,r=10)
plot(x,fs,type="l",col="black",ylab = "")
lines(x = c(4,4),y = c(min(fs),max(fs)),col="black")
lines(x = x,y = lz,col="red")
lines(x = x,y = pz,col="green")
lines(x = x,y = pz10,col="darkgreen")
