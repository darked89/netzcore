library(ROCR)

v<-read.table("netscore2_predictions_n40_e0.001_R102_z_r_q.txt.sample25")
l<-read.table("netscore2_labels_n40_e0.001_R102_z_r_q.txt.sample25")
pred<-prediction(v, l)
perfNetZcore<-performance(pred, "tpr", "fpr")

v<-read.table("netscore2_predictions_n40_e0.001_R101_z_r.txt")
l<-read.table("netscore2_labels_n40_e0.001_R101_z_r.txt")
pred<-prediction(v, l)
perfNetZcore0<-performance(pred, "tpr", "fpr")

v<-read.table("netscore2_predictions_n2_e0.01_R102_z_a_q.txt")
l<-read.table("netscore2_labels_n2_e0.01_R102_z_a_q.txt")
pred<-prediction(v, l)
perfNetZcore1<-performance(pred, "tpr", "fpr")

v<-read.table("netscore2_predictions_n2_e0.01_R101_z.txt")
l<-read.table("netscore2_labels_n2_e0.01_R101_z.txt")
pred<-prediction(v, l)
perfNetZcore2<-performance(pred, "tpr", "fpr")

v<-read.table("functionalFlow_predictions_n2_e0.01.txt")
l<-read.table("functionalFlow_labels_n2_e0.01.txt")
pred<-prediction(v, l)
perfNetZcore3<-performance(pred, "tpr", "fpr")

v<-read.table("kMajority_predictions_n2_e0.01_e.txt")
l<-read.table("kMajority_labels_n2_e0.01_e.txt")
pred<-prediction(v, l)
perfNetZcore4<-performance(pred, "tpr", "fpr")

v<-read.table("netscore2_predictions_n40_e0.001_R1_N_r.txt")
l<-read.table("netscore2_labels_n40_e0.001_R1_N_r.txt")
pred<-prediction(v, l)
perfNetZcore5<-performance(pred, "tpr", "fpr")

#postscript("all_performance.eps", width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")

plot(perfNetZcore, xlab="False Positive Rate", ylab="True Positive Rate", plotCI.col=2, lwd=2,avg="vertical",spread.estimate="stddev", show.spread.at=seq(0,1,by=0.225), print.cutoffs.at=seq(-4,4, by=2), col=2)
plot(perfNetZcore0, plotCI.col=3, lwd=2,avg="vertical", spread.estimate="stddev", show.spread.at=seq(0,1,by=0.33), col=3, add=TRUE)
plot(perfNetZcore1, plotCI.col=4, lwd=2,avg="vertical", spread.estimate="stddev", show.spread.at=seq(0,1,by=0.15), col=4, add=TRUE)
plot(perfNetZcore2, plotCI.col=5, lwd=2,avg="vertical", spread.estimate="stddev", show.spread.at=seq(0,1,by=0.20), col=5, add=TRUE)
plot(perfNetZcore3, plotCI.col=6, lwd=2,avg="vertical", spread.estimate="stddev", show.spread.at=seq(0,1,by=0.25), col=6, add=TRUE)
plot(perfNetZcore4, plotCI.col=7, lwd=2,avg="vertical", spread.estimate="stddev", show.spread.at=seq(0,1,by=0.37), col=7, add=TRUE)
plot(perfNetZcore5, plotCI.col=8, lwd=2,avg="vertical", spread.estimate="stddev", show.spread.at=seq(0,1,by=0.43), col=8, add=TRUE)
savefont<-par(ps=10)
legend("bottomright", c("Norm1_R/d", "Norm1_TPw", "Norm1_R/d_Acc_NoRel (n=2)", "Norm1_TPw_NoRel (n=2)", "fFlow (n=2)", "kMajority (n=2)", "NoNorm"), lty=c(1,1,1,1,1,1,1), col=c(2,3,4,5,6,7,8)) 

#dev.off()

