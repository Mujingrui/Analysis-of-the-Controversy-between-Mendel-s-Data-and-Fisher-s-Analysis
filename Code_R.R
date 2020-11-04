W1<-read.csv("F:/Genetics Analysis/Week2.csv",header=TRUE,sep=",")
library(MASS)
library(stringr)
b <- as.data.frame(apply(str_split(W1$p0,"/",simplify = T),2,as.numeric))
W1$p0<-b[,1]/b[,2]
chi<-round((W1$n1-W1$n*W1$p0)/sqrt(W1$n*W1$p0*(1-W1$p0)),4)
p_value<-round((1-pchisq(chi^2,1)),4)
W1[,6]<-chi
W1[,7]<-p_value
names(W1) <- c("No.","n","n1","n-n1","p0","Chi","p_value")
library(stats)
y<-runif(100)
random_numbers<-rnorm(84,mean=0,sd=10^(-7))
x<-W1$p_value+random_numbers
# x2<-seq(0,1,by=0.01)
plot(ecdf(W1$p_value),  
     xlim = c(0,1),
     ylim = c(0,1),
     col = "blue",
     main = "Figure1",
     xlab = "x",
     ylab = "F(x)",
     pch=20,
     cex=0.5)
curve(x^2, 0, 1,add=TRUE, bty="l")
abline(coef = c(0, 1), col = "red")

library(dgof)
ks.test(x,'punif')

u2<-runif(100,min=0,max=1)
sample_f<-function(alpha,u){
  n<-length(u)
  X<-rep(0,n)
  for(i in 1:n){
    if(u[i]>=0 & u[i]<alpha){
      X[i]<-sqrt(u[i])
    }
    if(u[i]>alpha & u[i]<=1){
      X[i]<-(u[i]+alpha)/(1+alpha)
    }
  }
  return(X)
}
sample_f0<-sample_f(0.01,u2)
sample_f1<-sample_f(0.201,u2)
ks.test(x,sample_f0)
ks.test(x,sample_f1)

alpha_s<-seq(0,1,by=0.001)
samplef<-matrix(0,nrow=length(alpha_s),ncol=100)
for(i in 1:length(alpha_s)){
  samplef[i,]<-sample_f(alpha_s[i],u2)
}

ks_D<-rep(0,1001)
ks_p<-rep(0,1001)
# ks<-list(0,1001)
for(j in 1:1001){
  ks_p[j]<-ks.test(x,samplef[j,])$p.value
  ks_D[j]<-ks.test(x,samplef[j,])$statistic
  # ks_D[i]<-ks$p.value
}
data1<-data.frame(alpha_s,ks_p,ks_D)
p1<-ggplot(data1,aes(alpha_s,ks_p))+geom_point()+geom_smooth()
max(ks_p)
which.max(ks_p)

plot(ecdf(W1$p_value),  
     xlim = c(0,1),
     ylim = c(0,1),
     col = "blue",
     main = "Figure1",
     xlab = "x",
     ylab = "F(x)",
     pch=20,
     cex=0.5)

x11 <- seq(0, 1, 0.001)
f = function(x) {
  sapply(x, function(x) {
    y = numeric()
    if (x >= 0 & x <= alpha_s[99])
      y = x^2
    else if (x > alpha_s[99] & x <= 1)
      y = (1+alpha_s[99])*x-alpha_s[99]
    return(y)
  })
}
f(x11)
# y11 <- ifelse(x11 >= 0 & x11 <= alpha_s[99], x11^2, ifelse(x11 > alpha_s[99] & x11 <= 1, (1+alpha_s[99])*x11-alpha_s[99]))
plot(ecdf(W1$p_value),  
     xlim = c(0,1),
     ylim = c(0,1),
     col = "blue",
     main = "Figure2",
     xlab = "x",
     ylab = "F(x)",
     pch=20,
     cex=0.5)
curve(f(x),0,1,add=TRUE, bty="l")

