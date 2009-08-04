x<-seq(1,10)
xl<-seq(1,10) # seq(0,100,1)
e<-exp(1)

#postscript("den.ps", width = 5, height = 5, horizontal = FALSE, onefile = FALSE, paper = "special")

plot(xl, x^-x, type="p", xlab="x", ylab="y", xlim=range(1,10))
axis(1, at=seq(1,10,1))
axis(2, at=seq(0,1,0.1))
#lines(xl, e^-x, col=2, lty=3) 
points(xl, x^-(x+1), col=2, pch="x") 
#lines(xl, x^-e, col=3, lty=4) 
lines(xl, x^-3, col=3, lty=3) 
lines(xl, x^-2, col=4, lty=4) 
lines(xl, x^-1, col=5, lty=1) 
lines(xl, e^-(x-1), col=7, lty=2) 
#par(lty=1)
#lines(xl, x^-(0.5), col=7) 
#lines(xl, x^-(2.5), col=8) 
#lines(xl, x^-0, col=7) 
#plot(xl, x^-2, type="p", col=7) 

#legend("topright", c("x^-x", "e^-x", "x^-e", "x^-3", "x^-2", "x^-1", "x^-0"), col=c(1,2,3,4,5,6,7), lty=c(1,1,1,1,1,1,1))
legend("topright", c("x^-1", "x^-2", "x^-3", "x^-x", "x^-(x+1)", "10^-x"), col=c(5,4,3,1,2,7), lty=c(1,4,3,0,0,2), pch=c("","","","o","x",""))

#dev.off()

