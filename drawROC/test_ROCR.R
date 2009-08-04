library(ROCR)

v<-c(1,2,1,2,3,2,4,2,1,0)
l<-c(1,0,0,0,1,0,1,0,0,0)
pred<-prediction(v, l)
#perfROC<-performance(pred, "tpr", "fpr")
#perfPPV<-performance(pred, "ppv")
perfSens<-performance(pred, "sens")

#postscript("all_performance.eps", width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")

plot(perfROC, xlab="False Positive Rate / PPV / Sens", ylab="True Positive Rate / Score", type="p", lwd=2, col=2)
#plot(perfPPV, lwd=2, col=3) #, add=T)
#plot(perfSens, lwd=2, col=4, add=T)
savefont<-par(ps=10)
legend("bottomright", c("ROC", "PPV", "Sens"), lty=c(1,1,1), col=c(2,3,4)) 

#dev.off()

