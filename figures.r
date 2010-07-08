dirname<-"../data/summary/"

#postscript(paste(dirname, "compare.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special", title = "Comparison of scoring methods on aneurysm data")

postscript(paste(dirname, "2a.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
d<-read.table(paste(dirname, 'all_vs_all/auc_phenotypes.dat', sep=""))
boxplot(d, xlab="Prediction method", ylab="Average AUC over all diseases (%)", names=c("NetScore", "NetZcore", "NetShort", "Func.Flow", "T.Gene"), col=rainbow(5), notch=FALSE, varwidth=FALSE)
dev.off()

postscript(paste(dirname, "2b.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
d<-read.table(paste(dirname, 'all_vs_all/cov_phenotypes.dat', sep=""))
boxplot(d, xlab="Prediction method", ylab="Ratio of correctly predicted proteins among top 10% predictions (%)", names=c("NetScore", "NetZcore", "NetShort", "Func.Flow", "T.Gene"), col=rainbow(5), notch=FALSE, varwidth=FALSE)
dev.off()

postscript(paste(dirname, "3a.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
d<-read.table(paste(dirname, 'biana_no_tap-omim/auc_ppis.dat', sep=""))
e<-data.frame(d[,1],d[,2],d[,3])
barplot(as.matrix(t(e)), beside=TRUE, horiz=F, col = rainbow(3), ylab="AUC (%)", legend.text=c("NetScore", "NetZcore", "NetShort"))
text(seq(1.1, 0.1+length(d[,1])*4.1, 4.1), par("usr")[3] - 2, srt=45, adj=1, labels=names(as.data.frame(t(d))), xpd=T, cex=0.8)
dev.off()

postscript(paste(dirname, "3b.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
d<-read.table(paste(dirname, 'biana_no_tap-goh/auc_ppis.dat', sep=""))
e<-data.frame(d[,1],d[,2],d[,3])
barplot(as.matrix(t(e)), beside=TRUE, horiz=F, col = rainbow(3), ylab="AUC (%)", legend.text=c("NetScore", "NetZcore", "NetShort"))
text(seq(1.1, 0.1+length(d[,1])*4.1, 4.1), par("usr")[3] - 2, srt=45, adj=1, labels=names(as.data.frame(t(d))), xpd=T, cex=0.8)
dev.off()

postscript(paste(dirname, "3c.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
d<-read.table(paste(dirname, 'biana_no_tap-chen/auc_ppis.dat', sep=""))
e<-data.frame(d[,1],d[,2],d[,3])
barplot(as.matrix(t(e)), beside=TRUE, horiz=F, col = rainbow(3), ylab="AUC (%)", legend.text=c("NetScore", "NetZcore", "NetShort"))
text(seq(1.1, 0.1+length(d[,1])*4.1, 4.1), par("usr")[3] - 2, srt=45, adj=1, labels=names(as.data.frame(t(d))), xpd=T, cex=0.8)
dev.off()

postscript(paste(dirname, "3d.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
d<-read.table(paste(dirname, 'seed_biana_no_tap-omim/seeds.dat', sep=""))
e<-data.frame(d[,1],d[,2],d[,3])
barplot(as.matrix(t(e)), beside=TRUE, horiz=F, col = rainbow(3), ylab="value", legend.text=c("# of seeds", "# of neigh. seeds", "S. path bw/ seeds"))
text(seq(1.1, 0.1+length(d[,1])*4.1, 4.1), par("usr")[3] - 2, srt=45, adj=1, labels=names(as.data.frame(t(d))), xpd=T, cex=0.8)
dev.off()

postscript(paste(dirname, "4_.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
d<-read.table(paste(dirname, 'biana_no_tap-all_seed20below/auc_phenotypes.dat', sep=""))
e<-read.table(paste(dirname, 'biana_no_tap-all_seed20/auc_phenotypes.dat', sep=""))
f<-merge(d,e,all=T)
barplot(as.matrix(t(f)), beside=T, xlab="Diseases Groups (w.r.t. number of initially annotated proteins)", ylab="Average AUC over all diseases (%)", names.arg=c("Nseed < 20", "Nseed >= 20"), col=rainbow(5), legend.text=c("NetScore","NetZcore", "NetShort", "Func. Flow", "ToppGene"))
dev.off()

postscript(paste(dirname, "5_.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
d<-read.table(paste(dirname, 'all_vs_all/auc_phenotypes.dat', sep=""))
barplot(as.matrix(t(d)), beside=T, xlab="Prediction method", ylab="Average AUC over all diseases (%)", names.arg=c("Biana", "Bia. filtered", "Entrez", "Goh", "Bia. weighted"), col=rainbow(5), legend.text=c("NetScore","NetZcore", "NetShort", "Func. Flow", "ToppGene"))
dev.off()

