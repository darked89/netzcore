dir.name<-"../data/summary/"

#postscript(paste(dir.name, "compare.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special", title = "Comparison of scoring methods on aneurysm data")

color5<-c(2,7,3,4,6) #rainbow(5)

postscript(paste(dir.name, "2a.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
d<-read.table(paste(dir.name, 'all_vs_all/auc_phenotypes.dat', sep=""))
# Average AUC over all diseases (%)
boxplot(d, xlab="Prediction method", ylab="AUC(%)", names=c("NetScore", "NetZcore", "NetShort", "Func.Flow", "T.Gene"), col=8, notch=FALSE, varwidth=FALSE, ylim=c(0,80))
dev.off()

postscript(paste(dir.name, "2b.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
d<-read.table(paste(dir.name, 'all_vs_all/cov_phenotypes.dat', sep=""))
#Ratio of correctly predicted proteins among top 10% predictions (%)
boxplot(d, xlab="Prediction method", ylab="PPV(%)", names=c("NetScore", "NetZcore", "NetShort", "Func.Flow", "T.Gene"), col=8, notch=FALSE, varwidth=FALSE, ylim=c(0,80))
dev.off()

#par(mar=c(5, 4, 4, 2) + 0.1) # to reset margins

postscript(paste(dir.name, "S2a.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
par(xpd=T, mar=par()$mar+c(0,0,2,0))
coords <- seq(0.4,4.8,by=1.1)
i = 1
cols <- color5 # rainbow(5) #2:6 #heat.colors(5)
for(ppi in c("goh", "entrez", "biana", "biana_no_tap", "biana_no_tap_relevance")) {
    d<-read.table(paste(dir.name, ppi, '-all/auc_ppis.dat', sep=""))
    if(i == 1) {
	boxplot(d,col=cols[i],at=coords, boxwex=0.2, pars=list(xaxt="n"), xlab="Prediction method", ylab="AUC(%)", ylim=c(0,100))
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
#legend("bottomright", c("GOH", "ENTREZ", "PPI", "bPPI", "bPPI weighted"), fill=color5, bty="n")
legend(0.5, 119, c("GOH", "ENTREZ", "PPI", "bPPI", "bPPI weighted"), fill=color5, bty="n", ncol=3) #horiz=T, x.intersp=0.5)
dev.off()

postscript(paste(dir.name, "S2b.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
par(xpd=T, mar=par()$mar+c(0,0,2,0))
coords <- seq(0.4,4.8,by=1.1)
i = 1
cols <- color5  
for(ppi in c("goh", "entrez", "biana", "biana_no_tap", "biana_no_tap_relevance")) {
    d<-read.table(paste(dir.name, ppi, '-all/cov_ppis.dat', sep=""))
    if(i == 1) {
	boxplot(d,col=cols[i],at=coords, boxwex=0.2, pars=list(xaxt="n"), xlab="Prediction method", ylab="PPV(%)", ylim=c(0,100))
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
legend(0.5, 119, c("GOH", "ENTREZ", "PPI", "bPPI", "bPPI weighted"), fill=color5, bty="n", ncol=3) 
dev.off()


postscript(paste(dir.name, "3a.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
e<-c()
f<-c()
percentages<-seq(0,70,by=10)
#methods<-c("NetScore", "NetZcore", "NetShort") 
methods<-c("NetScore", "NetZcore", "NetShort", "F.Flow", "T.Gene")
label<-"Percentage of permuted interactions"
for(p in percentages) {
    if(p == 0) {
	d <- read.table(paste(dir.name, 'biana_no_tap-omim/auc_phenotypes.dat', sep=""))
	#d <- d[,1:3]
	dir.name<-"../data/summary_gaudi/"
    }
    else {
	d <- read.table(paste(dir.name, 'biana_no_tap_permuted_p', p, '_10-omim/auc_phenotypes.dat', sep=""))
    }
    m <- mean(d)
    n <- dim(d)[1]
    error <- qt(0.975, df=n-1) * sd(d) / sqrt(n)
    e <- cbind(e,m)
    f <- cbind(f,error)
}

# ylim = c(0,max(e)+15)
b<-barplot(e, beside=TRUE, horiz=F, col = color5[1:length(methods)], width=0.1, space=c(0,0.4), names.arg=rep('', length(percentages)), ylim=c(0,100), ylab="AUC (%)", xlab=label) #, legend.text=methods)
segments(b, e-f, b, e+f, col=1, lty=1, lwd=1) #lty=3 dashed col=8 grey
points(b, e-f, col=1, lwd=2, pch='_', cex=0.9)
points(b, e+f, col=1, lwd=2, pch='_', cex=0.9)
#text(seq(0.2, 0.34*length(percentages)-0.1, 0.34), par("usr")[3] - 2, srt=45, adj=1, labels=percentages, xpd=T, cex=0.8)
text(seq(0.3, 0.55*length(percentages)-0.1, 0.55), par("usr")[3] - 2, srt=45, adj=1, labels=percentages, xpd=T, cex=0.8)
legend(0.5, 100, methods, fill=color5, bty="n", ncol=3)
dev.off()

dir.name<-"../data/summary/"
postscript(paste(dir.name, "3b.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
e<-c()
f<-c()
percentages<-seq(0,70,by=10)
methods<-c("NetScore", "NetZcore", "NetShort", "F.Flow", "T.Gene")
label<-"Percentage of pruned interactions"
for(p in percentages) {
    if(p == 0) {
	d <- read.table(paste(dir.name, 'biana_no_tap-omim/auc_phenotypes.dat', sep=""))
	dir.name<-"../data/summary_gaudi/"
    }
    else {
	d <- read.table(paste(dir.name, 'biana_no_tap_pruned_p', p, '_10-omim/auc_phenotypes.dat', sep=""))
    }
    m <- mean(d)
    n <- dim(d)[1]
    error <- qt(0.975, df=n-1) * sd(d) / sqrt(n)
    e <- cbind(e,m)
    f <- cbind(f,error)
}
b<-barplot(e, beside=TRUE, horiz=F, col = color5[1:length(methods)], width=0.1, space=c(0,0.4), names.arg=rep('', length(percentages)), ylim=c(0,100), ylab="AUC (%)", xlab=label)
segments(b, e-f, b, e+f, col=1, lty=1, lwd=1) 
points(b, e-f, col=1, pch='_', cex=0.9)
points(b, e+f, col=1, pch='_', cex=0.9)
text(seq(0.3, 0.55*length(percentages)-0.1, 0.55), par("usr")[3] - 2, srt=45, adj=1, labels=percentages, xpd=T, cex=0.8)
legend(0.5, 100, methods, fill=color5, bty="n", ncol=3)
dev.off()

dir.name<-"../data/summary/"
postscript(paste(dir.name, "3c.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
e<-c()
f<-c()
percentages<-seq(0,70,by=10)
methods<-c("NetScore", "NetZcore", "NetShort", "F.Flow", "T.Gene")
label<-"Percentage of removed non-seed interactions"
for(p in percentages) {
    if(p == 0) {
	d <- read.table(paste(dir.name, 'biana_no_tap-omim/auc_phenotypes.dat', sep=""))
	dir.name<-"../data/summary_gaudi/"
    }
    else {
	d <- read.table(paste(dir.name, 'biana_no_tap_pruned_non_seed_interactions_p', p, '_10-omim/auc_phenotypes.dat', sep=""))
    }
    m <- mean(d)
    n <- dim(d)[1]
    error <- qt(0.975, df=n-1) * sd(d) / sqrt(n)
    e <- cbind(e,m)
    f <- cbind(f,error)
}
b<-barplot(e, beside=TRUE, horiz=F, col = color5[1:length(methods)], width=0.1, space=c(0,0.4), names.arg=rep('', length(percentages)), ylim=c(0,100), ylab="AUC (%)", xlab=label)
segments(b, e-f, b, e+f, col=1, lty=1, lwd=1) 
points(b, e-f, col=1, pch='_', cex=0.9)
points(b, e+f, col=1, pch='_', cex=0.9)
text(seq(0.3, 0.55*length(percentages)-0.1, 0.55), par("usr")[3] - 2, srt=45, adj=1, labels=percentages, xpd=T, cex=0.8)
legend(0.5, 100, methods, fill=color5, bty="n", ncol=3)
dev.off()

dir.name<-"../data/summary/"
postscript(paste(dir.name, "5.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
e<-c()
f<-c()
percentages<-seq(0,70,by=10)
methods<-c("NetScore", "NetZcore", "NetShort", "F.Flow", "T.Gene")
label<-"Percentage of perturbed seeds"
for(p in percentages) {
    if(p == 0) {
	d <- read.table(paste(dir.name, 'biana_no_tap-omim_seedsabove19/auc_phenotypes.dat', sep=""))
	dir.name<-"../data/summary_gaudi/"
    }
    else {
	d<-read.table(paste(dir.name, 'biana_no_tap-omim_perturbed_p', p, '_10_seedsabove19/auc_perturbed.dat', sep=""))
    }
    m <- mean(d)
    n <- dim(d)[1]
    error <- qt(0.975, df=n-1) * sd(d) / sqrt(n)
    e <- cbind(e,m)
    f <- cbind(f,error)
}
b<-barplot(e, beside=TRUE, horiz=F, col = color5[1:length(methods)], width=0.1, space=c(0,0.4), names.arg=rep('', length(percentages)), ylim=c(0,100), ylab="AUC (%)", xlab=label)
segments(b, e-f, b, e+f, col=1, lty=1, lwd=1) 
points(b, e-f, col=1, pch='_', cex=0.9)
points(b, e+f, col=1, pch='_', cex=0.9)
text(seq(0.3, 0.55*length(percentages)-0.1, 0.55), par("usr")[3] - 2, srt=45, adj=1, labels=percentages, xpd=T, cex=0.8)
legend(0.5, 100, methods, fill=color5, bty="n", ncol=3)
dev.off()

dir.name<-"../data/summary/"
postscript(paste(dir.name, "5b.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
e<-c()
f<-c()
percentages<-seq(0,70,by=10)
methods<-c("NetScore", "NetZcore", "NetShort", "F.Flow", "T.Gene")
label<-"Percentage of perturbed seeds"
for(p in percentages) {
    if(p == 0) {
	d <- read.table(paste(dir.name, 'biana_no_tap-omim/auc_phenotypes.dat', sep=""))
	dir.name<-"../data/summary_gaudi/"
    }
    else {
	d<-read.table(paste(dir.name, 'biana_no_tap-omim_perturbed_p', p, '_10/auc_perturbed.dat', sep=""))
    }
    m <- mean(d)
    n <- dim(d)[1]
    error <- qt(0.975, df=n-1) * sd(d) / sqrt(n)
    e <- cbind(e,m)
    f <- cbind(f,error)
}
b<-barplot(e, beside=TRUE, horiz=F, col = color5[1:length(methods)], width=0.1, space=c(0,0.4), names.arg=rep('', length(percentages)), ylim=c(0,100), ylab="AUC (%)", xlab=label)
segments(b, e-f, b, e+f, col=1, lty=1, lwd=1) 
points(b, e-f, col=1, pch='_', cex=0.9)
points(b, e+f, col=1, pch='_', cex=0.9)
text(seq(0.3, 0.55*length(percentages)-0.1, 0.55), par("usr")[3] - 2, srt=45, adj=1, labels=percentages, xpd=T, cex=0.8)
legend(0.5, 100, methods, fill=color5, bty="n", ncol=3)
dev.off()

dir.name<-"../data/summary/"
postscript(paste(dir.name, "S4.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
d<-read.table(paste(dir.name, 'biana_no_tap-all_seed20below/auc_phenotypes.dat', sep=""))
e<-read.table(paste(dir.name, 'biana_no_tap-all_seed20/auc_phenotypes.dat', sep=""))
f<-merge(d,e,all=T)
#xlab="Diseases Groups (w.r.t. number of initially annotated proteins)", 
barplot(as.matrix(t(f)), beside=T, col=color5, ylim=c(0,100), ylab="AUC(%)", names.arg=c("Nseed < 20", "Nseed >= 20"))
legend(0.5, 100, methods, fill=color5, bty="n", horiz=T)
dev.off()



# Old figures
old_figures <- function() { 

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

postscript(paste(dir.name, "5_.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
d<-read.table(paste(dir.name, 'all_vs_all/auc_phenotypes.dat', sep=""))
barplot(as.matrix(t(d)), beside=T, xlab="Prediction method", ylab="Average AUC over all diseases (%)", names.arg=c("Biana", "Bia. filtered", "Entrez", "Goh", "Bia. weighted"), col=color5, legend.text=c("NetScore","NetZcore", "NetShort", "Func. Flow", "ToppGene"))
dev.off()

}

