dir.name<-"../data/summary/"

#postscript(paste(dir.name, "compare.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special", title = "Comparison of scoring methods on aneurysm data")

postscript(paste(dir.name, "2a.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
#d<-read.table(paste(dir.name, 'all_vs_all/auc_phenotypes.dat', sep=""))
#boxplot(d, xlab="Prediction method", ylab="Average AUC over all diseases (%)", names=c("NetScore", "NetZcore", "NetShort", "Func.Flow", "T.Gene"), col=rainbow(5), notch=FALSE, varwidth=FALSE)
coords <- seq(0.4,4.8,by=1.1)
i = 1
cols <- rainbow(5) #2:6 #heat.colors(5) 
for(ppi in c("goh", "entrez", "biana", "biana_no_tap", "biana_no_tap_relevance")) {
    d<-read.table(paste(dir.name, ppi, '-all/auc_ppis.dat', sep=""))
    if(i == 1) {
	boxplot(d,col=cols[i],at=coords, boxwex=0.2, pars=list(xaxt="n"), xlab="Prediction method", ylab="Average AUC over all diseases (%)")
    } 
    else {
	boxplot(d,col=cols[i],at=coords, boxwex=0.2, pars=list(xaxt="n"), names=F, add=T)
    }
    if(i == 3) {
	axis(1,tick=F,labels=c("NetScore", "NetZcore", "NetShort", "Func.Flow", "T.Gene"),at=coords)
    }
    coords <- coords + 0.2
    i = i + 1
}
legend("bottomleft", c("GOH", "ENTREZ", "PPI", "bPPI", "bPPI weighted"), fill=rainbow(5))
dev.off()

postscript(paste(dir.name, "2b.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
#d<-read.table(paste(dir.name, 'all_vs_all/cov_phenotypes.dat', sep=""))
#boxplot(d, xlab="Prediction method", ylab="Ratio of correctly predicted proteins among top 10% predictions (%)", names=c("NetScore", "NetZcore", "NetShort", "Func.Flow", "T.Gene"), col=rainbow(5), notch=FALSE, varwidth=FALSE)
coords <- seq(0.4,4.8,by=1.1)
i = 1
cols <- rainbow(5)  
for(ppi in c("goh", "entrez", "biana", "biana_no_tap", "biana_no_tap_relevance")) {
    d<-read.table(paste(dir.name, ppi, '-all/cov_ppis.dat', sep=""))
    if(i == 1) {
	boxplot(d,col=cols[i],at=coords, boxwex=0.2, pars=list(xaxt="n"), xlab="Prediction method", ylab="Ratio of correctly predicted proteins among top 10% predictions (%)")
    } 
    else {
	boxplot(d,col=cols[i],at=coords, boxwex=0.2, pars=list(xaxt="n"), names=F, add=T)
    }
    if(i == 3) {
	axis(1,tick=F,labels=c("NetScore", "NetZcore", "NetShort", "Func.Flow", "T.Gene"),at=coords)
    }
    coords <- coords + 0.2
    i = i + 1
}
legend("bottomleft", c("GOH", "ENTREZ", "PPI", "bPPI", "bPPI weighted"), fill=rainbow(5))
dev.off()

postscript(paste(dir.name, "3a.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
e<-c()
f<-c()
for(p in seq(10,100,by=10)) {
    d <- read.table(paste(dir.name, 'biana_no_tap_permuted_p', p, '_10-omim/auc_ppis.dat', sep=""))
    m <- mean(d)
    n <- dim(d)[1]
    error <- qt(0.975, df=n-1) * sd(d) / sqrt(n)
    e <- cbind(e,m)
    f <- cbind(f,error)
}
m<-max(e)
b<-barplot(e, beside=TRUE, horiz=F, col = rainbow(3), width=0.1, space=c(0,0.4), names.arg=rep('',10), ylim=c(0,m+15), ylab="AUC (%)", xlab="Percentage of permuted interactions", legend.text=c("NetScore", "NetZcore", "NetShort"))
segments(b, e-f, b, e+f, col=2)
text(seq(0.2, 3.5, 0.34), par("usr")[3] - 2, srt=45, adj=1, labels=seq(10,100,by=10), xpd=T, cex=0.8)
dev.off()

postscript(paste(dir.name, "3b.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
e<-c()
for(p in seq(10,90,by=10)) {
    d<-read.table(paste(dir.name, 'biana_no_tap_pruned_p', p, '_10-omim/auc_ppis.dat', sep=""))
    e<-cbind(e,mean(d))
    barplot(e, beside=TRUE, horiz=F, col = rainbow(3), width=0.1, space=c(0,0.4), ylab="AUC (%)", xlab="Percentage of pruned interactions", legend.text=c("NetScore", "NetZcore", "NetShort"))
    text(seq(0.2, 3.2, 0.34), par("usr")[3] - 2, srt=45, adj=1, labels=seq(10,90,by=10), xpd=T, cex=0.8)
}
dev.off()

postscript(paste(dir.name, "5.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
e<-c()
for(p in seq(10,90,by=10)) {
    d<-read.table(paste(dir.name, 'biana_no_tap-omim_perturbed_', p, '_10/auc_phenotypes.dat', sep=""))
    e<-cbind(e,mean(d))
    barplot(e, beside=TRUE, horiz=F, col = rainbow(3), width=0.1, space=c(0,0.4), ylab="AUC (%)", xlab="Percentage of permuted seeds", legend.text=c("NetScore", "NetZcore", "NetShort"))
    text(seq(0.2, 3.2, 0.34), par("usr")[3] - 2, srt=45, adj=1, labels=seq(10,90,by=10), xpd=T, cex=0.8)
}
dev.off()

old # to stop running

# Old figures
postscript(paste(dir.name, "3a.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
d<-read.table(paste(dir.name, 'biana_no_tap-omim/auc_ppis.dat', sep=""))
e<-data.frame(d[,1],d[,2],d[,3])
barplot(as.matrix(t(e)), beside=TRUE, horiz=F, col = rainbow(3), ylab="AUC (%)", legend.text=c("NetScore", "NetZcore", "NetShort"))
text(seq(1.1, 0.1+length(d[,1])*4.1, 4.1), par("usr")[3] - 2, srt=45, adj=1, labels=names(as.data.frame(t(d))), xpd=T, cex=0.8)
dev.off()

postscript(paste(dir.name, "3b.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
d<-read.table(paste(dir.name, 'biana_no_tap-goh/auc_ppis.dat', sep=""))
e<-data.frame(d[,1],d[,2],d[,3])
barplot(as.matrix(t(e)), beside=TRUE, horiz=F, col = rainbow(3), ylab="AUC (%)", legend.text=c("NetScore", "NetZcore", "NetShort"))
text(seq(1.1, 0.1+length(d[,1])*4.1, 4.1), par("usr")[3] - 2, srt=45, adj=1, labels=names(as.data.frame(t(d))), xpd=T, cex=0.8)
dev.off()

postscript(paste(dir.name, "3c.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
d<-read.table(paste(dir.name, 'biana_no_tap-chen/auc_ppis.dat', sep=""))
e<-data.frame(d[,1],d[,2],d[,3])
barplot(as.matrix(t(e)), beside=TRUE, horiz=F, col = rainbow(3), ylab="AUC (%)", legend.text=c("NetScore", "NetZcore", "NetShort"))
text(seq(1.1, 0.1+length(d[,1])*4.1, 4.1), par("usr")[3] - 2, srt=45, adj=1, labels=names(as.data.frame(t(d))), xpd=T, cex=0.8)
dev.off()

postscript(paste(dir.name, "3d.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
d<-read.table(paste(dir.name, 'seed_biana_no_tap-omim/seeds.dat', sep=""))
e<-data.frame(d[,1],d[,2],d[,3])
barplot(as.matrix(t(e)), beside=TRUE, horiz=F, col = rainbow(3), ylab="value", legend.text=c("# of seeds", "# of neigh. seeds", "S. path bw/ seeds"))
text(seq(1.1, 0.1+length(d[,1])*4.1, 4.1), par("usr")[3] - 2, srt=45, adj=1, labels=names(as.data.frame(t(d))), xpd=T, cex=0.8)
dev.off()

postscript(paste(dir.name, "4_.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
d<-read.table(paste(dir.name, 'biana_no_tap-all_seed20below/auc_phenotypes.dat', sep=""))
e<-read.table(paste(dir.name, 'biana_no_tap-all_seed20/auc_phenotypes.dat', sep=""))
f<-merge(d,e,all=T)
barplot(as.matrix(t(f)), beside=T, xlab="Diseases Groups (w.r.t. number of initially annotated proteins)", ylab="Average AUC over all diseases (%)", names.arg=c("Nseed < 20", "Nseed >= 20"), col=rainbow(5), legend.text=c("NetScore","NetZcore", "NetShort", "Func. Flow", "ToppGene"))
dev.off()

postscript(paste(dir.name, "5_.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
d<-read.table(paste(dir.name, 'all_vs_all/auc_phenotypes.dat', sep=""))
barplot(as.matrix(t(d)), beside=T, xlab="Prediction method", ylab="Average AUC over all diseases (%)", names.arg=c("Biana", "Bia. filtered", "Entrez", "Goh", "Bia. weighted"), col=rainbow(5), legend.text=c("NetScore","NetZcore", "NetShort", "Func. Flow", "ToppGene"))
dev.off()

