ca.pos.2000 <- read.csv("/home/mkim8/Data/RSS/BT-Cell/April-2013/Window/CA_Density/BT.2000.CA.positive", header=F, sep="\t")
ca.neg.2000 <- read.csv("/home/mkim8/Data/RSS/BT-Cell/April-2013/Window/CA_Density/BT.2000.CA.negative", header=F, sep="\t")
ca.pos.2000.mean <- colMeans(ca.pos.2000[,2:97])
ca.neg.2000.mean <- colMeans(ca.neg.2000[,2:97])
ca.pos.2000.plot <- cbind(seq(-4750,4750,100),ca.pos.2000.mean)
ca.neg.2000.plot <- cbind(seq(-4750,4750,100),ca.neg.2000.mean)

rss12.pos.2000 <- read.csv("/home/mkim8/Data/RSS/BT-Cell/April-2013/Window/RSS_Conc/BT.2000.12rss.positive",header=F, sep="\t")
rss12.neg.2000 <- read.csv("/home/mkim8/Data/RSS/BT-Cell/April-2013/Window/RSS_Conc/BT.2000.12rss.negative",header=F, sep="\t")
rss12.pos.2000.mean <- colMeans(rss12.pos.2000[,2:97])
rss12.neg.2000.mean <- colMeans(rss12.neg.2000[,2:97])

rss23.neg.2000 <- read.csv("/home/mkim8/Data/RSS/BT-Cell/April-2013/Window/RSS_Conc/BT.2000.23rss.negative",header=F, sep="\t")
rss23.pos.2000 <- read.csv("/home/mkim8/Data/RSS/BT-Cell/April-2013/Window/RSS_Conc/BT.2000.23rss.positive",header=F, sep="\t")
rss23.pos.2000.mean <- colMeans(rss23.pos.2000[,2:97])
rss23.neg.2000.mean <- colMeans(rss23.neg.2000[,2:97])


par(xpd = TRUE,mar=c(5,4,4,8))
plot(ca.pos.2000.plot,type="b",cex=0.5,col="blue",xlab="Distance from TSS", ylab="CA Density",ylim=c(0,0.1),bty="L") + points(ca.neg.2000.plot,col="red",cex=0.5) +axis(side=1, at=seq(-5000,5000,1000), las=1) 
text(2000,0.021,"12RSS")
text(2000,0.053,"23RSS")

par(new = TRUE)
plot(rss12.pos.2000.mean, type="l",axes=FALSE,xlab="",ylab="",ylim=c(0,0.012),col="blue",lwd=2) + lines(rss12.neg.2000.mean, col="red",lwd=2)
par(new = TRUE)
plot(rss23.pos.2000.mean, type="l",axes=FALSE,xlab="",ylab="",ylim=c(0,0.012),col="blue",lwd=2) + lines(rss23.neg.2000.mean, col="red",lwd=2)

z  <- seq(0,0.012,0.001)
axis(side=4, at = pretty(range(z)))
mtext("RSS_Concentration",side=4,col="black",line=3)
legend("topright",c("Rag1 Positive","Rag1 Negative","CA_Density","RSS_Conc"),inset=c(0.01,-0.1),lty=c(1,1,NA,1),lwd=c(2.5,2.5,NA,2.5),col=c("blue","red","black","black"),pch=c(NA,NA,1,NA))

