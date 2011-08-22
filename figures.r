dir.name<-"../data/summary/"

color5<-c(2,7,3,4,6) #rainbow(5)
color2<-c(3,4) #rainbow(5)
methods<-c("NetScore", "NetZcore", "NetShort", "F.Flow", "T.Gene")

main <- function() {
    manuscript_figures()
    #navlakha_figures()
    #manuscript_tests()
}

### MANUSCRIPT FIGURES ###
manuscript_figures <- function() {

    # Average AUC (%) over all diseases and networks
    #postscript(paste(dir.name, "2a.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
    jpeg(paste(dir.name, "2a.jpeg", sep=""), width = 480, height = 480, quality = 100)
    #tiff(paste(dir.name, "2a.tiff", sep=""), width = 480, height = 480)
    d<-read.table(paste(dir.name, 'all_vs_all/auc_phenotypes.dat', sep=""))
    boxplot(d, xlab="Prediction method", ylab="AUC(%)", names=c("NetScore", "NetZcore", "NetShort", "Func.Flow", "T.Gene"), col=8, notch=FALSE, varwidth=FALSE, ylim=c(0,80))
    dev.off()

    # Ratio of correctly predicted proteins (%) among top 10% predictions over all diseases and networks
    #postscript(paste(dir.name, "2b.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
    jpeg(paste(dir.name, "2b.jpeg", sep=""), width = 480, height = 480, quality = 100)
    d<-read.table(paste(dir.name, 'all_vs_all/cov_phenotypes.dat', sep=""))
    boxplot(d, xlab="Prediction method", ylab="PPV(%)", names=c("NetScore", "NetZcore", "NetShort", "Func.Flow", "T.Gene"), col=8, notch=FALSE, varwidth=FALSE, ylim=c(0,80))
    dev.off()

    #par(mar=c(5, 4, 4, 2) + 0.1) # to reset margins

    # Average AUC (%) over all diseases for each network
    #postscript(paste(dir.name, "S2a.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
    jpeg(paste(dir.name, "S2a.jpeg", sep=""), width = 480, height = 480, quality = 100)
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
    legend(0.5, 119, c("GOH", "ENTREZ", "PPI", "bPPI", "bPPI weighted"), fill=color5, bty="n", ncol=3) #horiz=T, x.intersp=0.5)
    dev.off()

    # Ratio of correctly predicted proteins (%) among top 10% predictions over all diseases for each network
    #postscript(paste(dir.name, "S2b.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
    jpeg(paste(dir.name, "S2b.jpeg", sep=""), width = 480, height = 480, quality = 100)
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

    # Average AUC (%) at different levels of interaction permutation over OMIM disorders on bPPI network
    #postscript(paste(dir.name, "3a.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
    jpeg(paste(dir.name, "3a.jpeg", sep=""), width = 480, height = 480, quality = 100)
    e<-c()
    f<-c()
    percentages<-seq(0,70,by=10)
    label<-"Percentage of permuted interactions"
    for(p in percentages) {
	if(p == 0) {
	    d <- read.table(paste(dir.name, 'biana_no_tap-omim/auc_phenotypes.dat', sep=""))
	    #d <- d[,1:3]
	    dir.name<-"../data/summary_runs_on_random/"
	}
	else {
	    d <- read.table(paste(dir.name, 'biana_no_tap_permuted_p', p, '-omim/auc_phenotypes.dat', sep=""))
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
    segments(b-0.03, e-f, b+0.03, e-f, col=1, lty=1, lwd=0.5)
    segments(b-0.03, e+f, b+0.03, e+f, col=1, lty=1, lwd=0.5)
    #points(b, e-f, col=1, lwd=2, pch='_', cex=0.9)
    #points(b, e+f, col=1, lwd=2, pch='_', cex=0.9)
    text(seq(0.3, 0.55*length(percentages)-0.1, 0.55), par("usr")[3] - 2, srt=45, adj=1, labels=percentages, xpd=T, cex=0.8)
    legend(0.5, 100, methods, fill=color5, bty="n", ncol=3)
    dev.off()


    # Average AUC (%) at different levels of interaction permutation over OMIM disorders on GOH network
    dir.name<-"../data/summary/"
    jpeg(paste(dir.name, "S3a.jpeg", sep=""), width = 480, height = 480, quality = 100)
    e<-c()
    f<-c()
    percentages<-seq(0,70,by=10)
    label<-"Percentage of permuted interactions"
    for(p in percentages) {
	if(p == 0) {
	    d <- read.table(paste(dir.name, 'goh-omim/auc_phenotypes.dat', sep=""))
	    #d <- d[,1:3]
	    dir.name<-"../data/summary_runs_on_random/"
	}
	else {
	    d <- read.table(paste(dir.name, 'goh_permuted_p', p, '-omim/auc_phenotypes.dat', sep=""))
	}
	m <- mean(d)
	n <- dim(d)[1]
	error <- qt(0.975, df=n-1) * sd(d) / sqrt(n)
	e <- cbind(e,m)
	f <- cbind(f,error)
    }
    b<-barplot(e, beside=TRUE, horiz=F, col = color5[1:length(methods)], width=0.1, space=c(0,0.4), names.arg=rep('', length(percentages)), ylim=c(0,100), ylab="AUC (%)", xlab=label) 
    segments(b, e-f, b, e+f, col=1, lty=1, lwd=1) 
    segments(b-0.03, e-f, b+0.03, e-f, col=1, lty=1, lwd=0.5)
    segments(b-0.03, e+f, b+0.03, e+f, col=1, lty=1, lwd=0.5)
    text(seq(0.3, 0.55*length(percentages)-0.1, 0.55), par("usr")[3] - 2, srt=45, adj=1, labels=percentages, xpd=T, cex=0.8)
    legend(0.5, 100, methods, fill=color5, bty="n", ncol=3)
    dev.off()


    # Average AUC (%) at different levels of interaction pruning over OMIM disorders on bPPI network
    dir.name<-"../data/summary/"
    #postscript(paste(dir.name, "3b.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
    jpeg(paste(dir.name, "3b.jpeg", sep=""), width = 480, height = 480, quality = 100)
    e<-c()
    f<-c()
    percentages<-seq(0,70,by=10)
    label<-"Percentage of pruned interactions"
    for(p in percentages) {
	if(p == 0) {
	    d <- read.table(paste(dir.name, 'biana_no_tap-omim/auc_phenotypes.dat', sep=""))
	    dir.name<-"../data/summary_runs_on_random/"
	}
	else {
	    d <- read.table(paste(dir.name, 'biana_no_tap_pruned_p', p, '-omim/auc_phenotypes.dat', sep=""))
	}
	m <- mean(d)
	n <- dim(d)[1]
	error <- qt(0.975, df=n-1) * sd(d) / sqrt(n)
	e <- cbind(e,m)
	f <- cbind(f,error)
    }
    b<-barplot(e, beside=TRUE, horiz=F, col = color5[1:length(methods)], width=0.1, space=c(0,0.4), names.arg=rep('', length(percentages)), ylim=c(0,100), ylab="AUC (%)", xlab=label)
    segments(b, e-f, b, e+f, col=1, lty=1, lwd=1) 
    segments(b-0.03, e-f, b+0.03, e-f, col=1, lty=1, lwd=0.5)
    segments(b-0.03, e+f, b+0.03, e+f, col=1, lty=1, lwd=0.5)
    text(seq(0.3, 0.55*length(percentages)-0.1, 0.55), par("usr")[3] - 2, srt=45, adj=1, labels=percentages, xpd=T, cex=0.8)
    legend(0.5, 100, methods, fill=color5, bty="n", ncol=3)
    dev.off()

    # Average AUC (%) at different levels of interaction pruning over OMIM disorders on GOH network
    dir.name<-"../data/summary/"
    jpeg(paste(dir.name, "S3b.jpeg", sep=""), width = 480, height = 480, quality = 100)
    e<-c()
    f<-c()
    percentages<-seq(0,70,by=10)
    label<-"Percentage of pruned interactions"
    for(p in percentages) {
	if(p == 0) {
	    d <- read.table(paste(dir.name, 'goh-omim/auc_phenotypes.dat', sep=""))
	    dir.name<-"../data/summary_runs_on_random/"
	}
	else {
	    d <- read.table(paste(dir.name, 'goh_pruned_p', p, '-omim/auc_phenotypes.dat', sep=""))
	}
	m <- mean(d)
	n <- dim(d)[1]
	error <- qt(0.975, df=n-1) * sd(d) / sqrt(n)
	e <- cbind(e,m)
	f <- cbind(f,error)
    }
    b<-barplot(e, beside=TRUE, horiz=F, col = color5[1:length(methods)], width=0.1, space=c(0,0.4), names.arg=rep('', length(percentages)), ylim=c(0,100), ylab="AUC (%)", xlab=label)
    segments(b, e-f, b, e+f, col=1, lty=1, lwd=1) 
    segments(b-0.03, e-f, b+0.03, e-f, col=1, lty=1, lwd=0.5)
    segments(b-0.03, e+f, b+0.03, e+f, col=1, lty=1, lwd=0.5)
    text(seq(0.3, 0.55*length(percentages)-0.1, 0.55), par("usr")[3] - 2, srt=45, adj=1, labels=percentages, xpd=T, cex=0.8)
    legend(0.5, 100, methods, fill=color5, bty="n", ncol=3)
    dev.off()

    # Average AUC (%) at different levels of interaction pruning only between non-seeds over OMIM disorders on bPPI network
    #dir.name<-"../data/summary/"
    #postscript(paste(dir.name, "3c.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
    #e<-c()
    #f<-c()
    #percentages<-seq(0,70,by=10)
    #label<-"Percentage of removed non-seed interactions"
    #for(p in percentages) {
    #	if(p == 0) {
    #	    d <- read.table(paste(dir.name, 'biana_no_tap-omim/auc_phenotypes.dat', sep=""))
    #	    dir.name<-"../data/summary_runs_on_random/"
    #	}
    #	else {
    #	    d <- read.table(paste(dir.name, 'biana_no_tap_pruned_non_seed_interactions_p', p, '-omim/auc_phenotypes.dat', sep=""))
    #	}
    #	m <- mean(d)
    #	n <- dim(d)[1]
    #	error <- qt(0.975, df=n-1) * sd(d) / sqrt(n)
    #	e <- cbind(e,m)
    #	f <- cbind(f,error)
    #}
    #b<-barplot(e, beside=TRUE, horiz=F, col = color5[1:length(methods)], width=0.1, space=c(0,0.4), names.arg=rep('', length(percentages)), ylim=c(0,100), ylab="AUC (%)", xlab=label)
    #segments(b, e-f, b, e+f, col=1, lty=1, lwd=1) 
    #points(b, e-f, col=1, pch='_', cex=0.9)
    #points(b, e+f, col=1, pch='_', cex=0.9)
    #text(seq(0.3, 0.55*length(percentages)-0.1, 0.55), par("usr")[3] - 2, srt=45, adj=1, labels=percentages, xpd=T, cex=0.8)
    #legend(0.5, 100, methods, fill=color5, bty="n", ncol=3)
    #dev.off()

    # Average AUC (%) at different levels of seed permutation over OMIM disorders on bPPI network - with CI considering all randoms & phenotypes
    dir.name<-"../data/summary/"
    #postscript(paste(dir.name, "4c.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
    jpeg(paste(dir.name, "4b_ci_all.jpeg", sep=""), width = 480, height = 480, quality = 100)
    e<-c()
    f<-c()
    percentages<-seq(0,70,by=10)
    label<-"Percentage of perturbed seeds"
    for(p in percentages) {
	if(p == 0) {
	    #d <- read.table(paste(dir.name, 'biana_no_tap-omim_seedsabove19/auc_phenotypes.dat', sep="")) #!
	    d <- read.table(paste(dir.name, 'biana_no_tap-omim/auc_phenotypes.dat', sep="")) 
	    dir.name<-"../data/summary_runs_on_random/"
	}
	else {
	    d<-read.table(paste(dir.name, 'biana_no_tap-omim_perturbed_p', p, '/auc_ppis.dat', sep=""))
	}
	m <- mean(d)
	n <- dim(d)[1]
	error <- qt(0.975, df=n-1) * sd(d) / sqrt(n)
	e <- cbind(e,m)
	f <- cbind(f,error)
    }
    b<-barplot(e, beside=TRUE, horiz=F, col = color5[1:length(methods)], width=0.1, space=c(0,0.4), names.arg=rep('', length(percentages)), ylim=c(0,100), ylab="AUC (%)", xlab=label)
    segments(b, e-f, b, e+f, col=1, lty=1, lwd=1) 
    segments(b-0.03, e-f, b+0.03, e-f, col=1, lty=1, lwd=0.5)
    segments(b-0.03, e+f, b+0.03, e+f, col=1, lty=1, lwd=0.5)
    text(seq(0.3, 0.55*length(percentages)-0.1, 0.55), par("usr")[3] - 2, srt=45, adj=1, labels=percentages, xpd=T, cex=0.8)
    legend(0.5, 100, methods, fill=color5, bty="n", ncol=3)
    dev.off()

    # Average AUC (%) at different levels of seed permutation over OMIM disorders on bPPI network
    dir.name<-"../data/summary/"
    #postscript(paste(dir.name, "4b.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
    jpeg(paste(dir.name, "4b.jpeg", sep=""), width = 480, height = 480, quality = 100)
    e<-c()
    f<-c()
    percentages<-seq(0,70,by=10)
    label<-"Percentage of perturbed seeds"
    for(p in percentages) {
	if(p == 0) {
	    d <- read.table(paste(dir.name, 'biana_no_tap-omim/auc_phenotypes.dat', sep=""))
	    dir.name<-"../data/summary_runs_on_random/"
	}
	else {
	    d<-read.table(paste(dir.name, 'biana_no_tap-omim_perturbed_p', p, '/auc_perturbed.dat', sep=""))
	}
	m <- mean(d)
	n <- dim(d)[1]
	error <- qt(0.975, df=n-1) * sd(d) / sqrt(n)
	e <- cbind(e,m)
	f <- cbind(f,error)
    }
    b<-barplot(e, beside=TRUE, horiz=F, col = color5[1:length(methods)], width=0.1, space=c(0,0.4), names.arg=rep('', length(percentages)), ylim=c(0,100), ylab="AUC (%)", xlab=label)
    segments(b, e-f, b, e+f, col=1, lty=1, lwd=1) 
    segments(b-0.03, e-f, b+0.03, e-f, col=1, lty=1, lwd=0.5)
    segments(b-0.03, e+f, b+0.03, e+f, col=1, lty=1, lwd=0.5)
    text(seq(0.3, 0.55*length(percentages)-0.1, 0.55), par("usr")[3] - 2, srt=45, adj=1, labels=percentages, xpd=T, cex=0.8)
    legend(0.5, 100, methods, fill=color5, bty="n", ncol=3)
    dev.off()

    # Average AUC (%) at different levels of seed permutation over OMIM disorders on GOH network
    dir.name<-"../data/summary/"
    #postscript(paste(dir.name, "4b.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
    jpeg(paste(dir.name, "S4b.jpeg", sep=""), width = 480, height = 480, quality = 100)
    e<-c()
    f<-c()
    percentages<-seq(0,70,by=10)
    label<-"Percentage of perturbed seeds"
    for(p in percentages) {
	if(p == 0) {
	    d <- read.table(paste(dir.name, 'goh-omim/auc_phenotypes.dat', sep=""))
	    dir.name<-"../data/summary_runs_on_random/"
	}
	else {
	    d<-read.table(paste(dir.name, 'goh-omim_perturbed_p', p, '_10/auc_perturbed.dat', sep=""))
	}
	m <- mean(d)
	n <- dim(d)[1]
	error <- qt(0.975, df=n-1) * sd(d) / sqrt(n)
	e <- cbind(e,m)
	f <- cbind(f,error)
    }
    b<-barplot(e, beside=TRUE, horiz=F, col = color5[1:length(methods)], width=0.1, space=c(0,0.4), names.arg=rep('', length(percentages)), ylim=c(0,100), ylab="AUC (%)", xlab=label)
    segments(b, e-f, b, e+f, col=1, lty=1, lwd=1) 
    segments(b-0.03, e-f, b+0.03, e-f, col=1, lty=1, lwd=0.5)
    segments(b-0.03, e+f, b+0.03, e+f, col=1, lty=1, lwd=0.5)
    text(seq(0.3, 0.55*length(percentages)-0.1, 0.55), par("usr")[3] - 2, srt=45, adj=1, labels=percentages, xpd=T, cex=0.8)
    legend(0.5, 100, methods, fill=color5, bty="n", ncol=3)
    dev.off()

    # Average AUC (%) for two groups of OMIM disorders based on number of seeds on bPPI network
    dir.name<-"../data/summary/"
    #postscript(paste(dir.name, "4a.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
    jpeg(paste(dir.name, "4a.jpeg", sep=""), width = 480, height = 480, quality = 100)

    s<-read.table(paste(dir.name, 'biana_no_tap-all/seeds.dat', sep=""))
    d<-read.table(paste(dir.name, 'biana_no_tap-all/auc_ppis.dat', sep=""))

    x<-d[rownames(d) %in% rownames(s)[s$n_seed<=20],]
    y<-d[rownames(d) %in% rownames(s)[s$n_seed>20],]

    coords <- seq(0.6, 5.0, by=1.1)
    boxplot(x,col=3, at=coords, boxwex=0.4, pars=list(xaxt="n"), xlab="Prediction method", ylab="AUC(%)", ylim=c(0,100))
    axis(1,tick=F,labels=c("NetScore", "NetZcore", "NetShort", "Func.Flow", "T.Gene"),at=coords+0.2)
    coords <- coords + 0.4
    boxplot(y,col=4, at=coords, boxwex=0.4, pars=list(xaxt="n"), names=F, add=T)
    legend(1.5, 5, c("# of seeds < 20", "# of seeds >= 20"), fill=color2, bty="n", horiz=T)

    #d<-read.table(paste(dir.name, 'biana_no_tap-all_seed20below/auc_ppis.dat', sep=""))
    #e<-read.table(paste(dir.name, 'biana_no_tap-all_seed20/auc_ppis.dat', sep=""))
    #f<-cbind(mean(d),mean(e))
    #b<-barplot(f, beside=T, col=color5, ylim=c(0,100), ylab="AUC(%)", xlab="Number of seeds", names.arg=c("1-19", "20-292"))
    #n1 <- dim(d)[1]
    #n2 <- dim(e)[1]
    #e1 <- qt(0.975, df=n1-1) * sd(d) / sqrt(n1)
    #e2 <- qt(0.975, df=n2-1) * sd(e) / sqrt(n2)
    #e <- cbind(mean(d), mean(e))
    #f <- cbind(e1, e2)
    #segments(b, e-f, b, e+f, col=1, lty=1, lwd=1) 
    #points(b, e-f, col=1, pch='_', cex=0.9)
    #points(b, e+f, col=1, pch='_', cex=0.9)
    #legend(3, 100, methods, fill=color5, bty="n", ncol=3)
    dev.off()

    d<-read.table("module_summary_top5_mcl_go_n5union.dat", header=T)
    #d<-read.table("module_summary_top5_mcl_n5union.dat", header=T)
    methods<-c("no", "nn", "ff", "ns") 
    method.names<-c("PWAS-Neighborhood", "PWAS-Func.Flow", "PWAS-NetScore")
    #methods<-c("no", "nn", "nr", "ff", "nd", "nz", "ns") 
    #method.names<-c("GWAS", "PWAS-1", "PWAS-2")
    #method.names<-c("PWAS-Neighborhood", "PWAS-Func.Flow", "PWAS-NetScore")
    n<-length(methods)
    e<-c()
    f<-c()
    g<-c()
    for(i in 1:n) {
	print(methods[i])
	#f<-rbind(f, d[d$scoring==methods[i],]$n_seed_in_modules/d[d$scoring==methods[i],]$n_seed) 
	if(i != 1) {
	    f<-rbind(f, d[d$scoring==methods[i],]$ratio) 
	    e<-rbind(e, d[d$scoring==methods[i],]$n_module)
	    g<-rbind(g, d[d$scoring==methods[i],]$n_seed_in_modules/d[d$scoring==methods[i],]$n_all_in_modules)
	}
    }
    lambda<-function(x) { substring(x, 6) }
    postscript(paste(dir.name, "S5a.eps", sep=""), horizontal=F)
    barplot(e, beside=TRUE, horiz=F, col=(2:n)+1, ylab="Number of modules", xlab="Phenotypes", legend.text=c("PWAS-Neighborhood", "PWAS-ToppGene", "PWAS-Func.Flow", "PWAS-NetShort", "PWAS-NetZcore", "PWAS-NetScore")) #, width=0.1, space=c(0,1.3), ylim=c(0,30), axes=F)
    x1<-par("usr")[1]
    x2<-par("usr")[2]
    text(seq(x1+4.2, x2, by=(x2-x1-4)/(length(levels(d$phenotype)))), par("usr")[3]-0.02, srt=30, adj=1, labels=sapply(levels(d$phenotype), lambda), xpd=T, cex=0.7)
    dev.off()
    postscript(paste(dir.name, "S5b.eps", sep=""), horizontal=F) 
    #barplot(f, beside=TRUE, horiz=F, col=(1:n)+1, ylab="Ratio of covered seeds within modules (%)", xlab="Phenotypes", legend.text=c("GWAS", "PWAS-1", "PWAS-2")) 
    
    barplot(f, beside=TRUE, horiz=F, col=(2:n)+1, ylab="Average seed GO term enrichment of modules (%)", xlab="Phenotypes", legend.text=method.names)
    x1<-par("usr")[1]
    x2<-par("usr")[2]
    text(seq(x1+4.2, x2, by=(x2-x1-4)/(length(levels(d$phenotype)))), par("usr")[3]-0.02, srt=30, adj=1, labels=sapply(levels(d$phenotype), lambda), xpd=T, cex=0.7)
    dev.off()

    postscript(paste(dir.name, "S5c.eps", sep=""), horizontal=F)
    barplot(g, beside=TRUE, horiz=F, col=(2:n)+1, ylab="Coverage of seed go terms within modules (%)", xlab="Phenotypes", legend.text=method.names) 
    x1<-par("usr")[1]
    x2<-par("usr")[2]
    text(seq(x1+4.2, x2, by=(x2-x1-4)/(length(levels(d$phenotype)))), par("usr")[3], srt=30, adj=1, labels=sapply(levels(d$phenotype), lambda), xpd=T, cex=0.7)
    dev.off()

    postscript(paste(dir.name, "S5d.eps", sep=""), horizontal=F) #width = 480, height = 480, quality = 100)
    d<-read.table("module_summary_top5_mcl_n5union.dat", header=T)
    methods<-c("no", "nn", "nd") 
    n<-length(methods)
    x<-1:length(levels(d$phenotype))
    par(mar=c(5,4,4,5))
    e<-c()
    for(i in 1:n) {
	print(methods[i])
	if(i != 1) {
	    e<-rbind(e, d[d$scoring==methods[i],]$ratio)
	}
    }
    a<-barplot(e, beside=TRUE, horiz=F, col=(2:n)+1, ylim=c(0,1), legend.text=c("PWAS-1", "PWAS-2"), xlab="Phenotypes", ylab="Ratio of non-seed connections within modules (%)") #width=0.1, space=c(0,1.3), 
    #legend("topright", c("PWAS-1", "PWAS-2"), lty=rep(1, n), col=(2:n)+1)
    lambda<-function(x) { substring(x, 6) }
    text(1:length(levels(d$phenotype)), par("usr")[3]-0.02, srt=30, adj=1, labels=sapply(levels(d$phenotype), lambda), xpd=T, cex=0.7)
    for(i in 1:n) {
	print(methods[i])
	if(i != 1) {
	    lines(colMeans(a), d[d$scoring==methods[i],]$n_seed_in_modules/d[d$scoring==methods[i],]$n_all_in_modules, lty=4, col=i+1)
	}
    }
    #axis(4, 1:15, labels=1:15)
    #text(par("usr")[2]+1, 10, srt=90, adj=1, labels="Number of modules using MCL", xpd=T, cex=1)
    dev.off()

    postscript(paste(dir.name, "S5e.eps", sep=""), horizontal=F) #width = 480, height = 480, quality = 100)
    d<-read.table("module_summary_top5_mcl_n5union.dat", header=T)
    methods<-c("no", "nn", "nd") 
    n<-length(methods)
    x<-1:length(levels(d$phenotype))
    par(mar=c(5,4,4,5))
    e<-c()
    for(i in 1:n) {
	print(methods[i])
	if(i == 1) {
	    plot(x, d[d$scoring==methods[i],]$ratio, type="l", col=i+1, ylim=c(0,1), xlab="Phenotypes", ylab="Ratio of seeds within modules (%)", xaxt="n")
	    #plot(x, d[d$scoring==methods[i],]$n_seed_in_modules/d[d$scoring==methods[i],]$n_seed, type="l", col=i+1, ylim=c(0,1), xlab="Phenotypes", ylab="Ratio of seeds within modules (%)", xaxt="n")
	} 
	else {
	    lines(x, d[d$scoring==methods[i],]$ratio, lty=1, col=i+1)
	    #lines(x, d[d$scoring==methods[i],]$n_seed_in_modules/d[d$scoring==methods[i],]$n_seed, lty=1, col=i+1)
	    #lines(x, d[d$scoring==methods[i],]$n_seed_in_modules/d[d$scoring==methods[i],]$n_all_in_modules, lty=4, col=i+1)
	    lines(x, d[d$scoring==methods[i],]$n_seed/d[d$scoring==methods[i],]$n_all_in_modules, lty=4, col=i+1)
	}
	if(i != 1) {
	    e<-rbind(e, d[d$scoring==methods[i],]$n_module)
	}
    }
    legend("topright", c("GWAS", "PWAS-1", "PWAS-2"), lty=rep(1, n), col=(1:n)+1)
    lambda<-function(x) { substring(x, 6) }
    text(1:length(levels(d$phenotype)), par("usr")[3]-0.02, srt=30, adj=1, labels=sapply(levels(d$phenotype), lambda), xpd=T, cex=0.7)
    par(new=TRUE)
    e<-as.matrix(e)
    print(dim(e))
    barplot(e, beside=TRUE, horiz=F, col=(2:n)+1, width=0.1, space=c(0,1.3), ylim=c(0,30), axes=F)
    axis(4, 1:15, labels=1:15)
    text(par("usr")[2]+1, 10, srt=90, adj=1, labels="Number of modules using MCL", xpd=T, cex=1)
    dev.off()

}

### MANUSCRIPT TESTS ###
manuscript_tests <- function() {
    dirname<-"../data/summary/"


    # Significance of getting Alzheimer genes in Krauthammer data set
    #sum(dhyper(6:61,55,7224,61))

    # Correlation between seed connectivity measures and AUC
    d<-read.table(paste(dirname, "biana_no_tap-omim/auc_ppis.dat", sep=""))
    d<-d[order(rownames(d)),]
    e<-read.table(paste(dirname, "biana_no_tap-omim/seeds.dat", sep=""))
    e<-e[order(rownames(e)),]
    print(cor(d,e))

    networks<-c("biana-all", "biana_no_tap-all", "biana_no_tap_relevance-all", "goh-all", "entrez-all")
    for(i in 1:length(networks)) {
	print(networks[i])
	d<-read.table(paste(dirname, networks[i], "/auc_ppis.dat", sep=""))
	d<-d[order(rownames(d)),]
	e<-read.table(paste(dirname, networks[i], "/seeds.dat", sep=""))
	e<-e[order(rownames(e)),]
	print(cor(d,e))
    }

    #d<-read.table(paste(dirname, "biana-all/auc_ppis.dat", sep=""))
    #e<-read.table(paste(dirname, "biana_no_tap-all/auc_ppis.dat", sep=""))
    #print(wilcox.test(d[,4], e[,4]))
    #print(kruskal.test(list(d[,4],e[,4])))
    #print(kruskal.test(list(d[,1],e[,1])))

    networks<-c("biana-all", "biana_no_tap-all", "biana_no_tap_relevance-all", "goh-all", "entrez-all")
    methods<-c("ns", "nz", "nd", "ff", "nr")
    d<-read.table(paste(dirname, "all_vs_all/auc_phenotypes.dat", sep=""))
    #x<-0; y<-0; for(i in d) { x<-x+1; for(j in d) y<-y+1; e<-i-j; if(abs(e)>5) { print(c(x,y,e)) } }
    e<-abs(d[1,]-d[1,])>5
    print(e)



    y<-1
    for(k in 1:5) {
	print(c("---", methods[k], "---"))
	x<-data.frame(row.names=networks)
	x[1:5,]<-NA
	x[,1:5]<-NA
	for(i in 1:length(networks)) {
	    for(j in 1:length(networks)) {
		if(i > j) {
		    d<-read.table(paste(dirname, networks[i], "/auc_ppis.dat", sep=""))
		    e<-read.table(paste(dirname, networks[j], "/auc_ppis.dat", sep=""))
		    a<-wilcox.test(d[,k], e[,k])
		    x[i,j]<-a$p.value
		    if(a$p.value < 0.05) {
			#print(c(methods[k], networks[i], networks[j]))
			#print(c(networks[i], networks[j]))
			#print(a$p.value)
			#x[y]<-a$p.value
			y<-y+1
			#print(kruskal.test(list(d[,k],e[,k])))
		    }
		}
	    }
	}
	names(x)<-networks
	print(x)
    }
    #print(x)
    #print(min(x))

    # Significance of AUC changes among networks
    networks1<-c("goh-all", "entrez-all")
    networks2<-c("entrez-all", "biana-all", "biana_no_tap-all")
    for(k in 1:5) {
	print(c("---", methods[k], "---"))
	x<-data.frame(row.names=networks1)
	for(i in 1:length(networks1)) {
	    for(j in 1:length(networks2)) {
		d<-read.table(paste(dirname, networks1[i], "/auc_ppis.dat", sep=""))
		e<-read.table(paste(dirname, networks2[j], "/auc_ppis.dat", sep=""))
		d<-d[order(row.names(d)),]
		e<-e[order(row.names(e)),] 
		a<-wilcox.test(d[,k], e[,k], alternative="less", paired=T)
		x[i,j]<-a$p.value
	    }
	}
	names(x)<-networks2
	print(x)
    }

    # Significance of connectivity change between PPI - bPPI
    d<-read.table(paste(dirname, "biana-all", "/seeds.dat", sep=""))
    e<-read.table(paste(dirname, "biana_no_tap-all", "/seeds.dat", sep=""))
    wilcox.test(e[,2], d[,2], alternative="less", paired=T)
    wilcox.test(e[,3], d[,3], alternative="greater", paired=T)

    # Significance of AUC change between PPI - bPPI
    networks1<-c("biana-all")
    networks2<-c("biana_no_tap-all")
    for(k in 1:5) {
	print(c("---", methods[k], "---"))
	x<-data.frame(row.names=networks1)
	for(i in 1:length(networks1)) {
	    for(j in 1:length(networks2)) {
		d<-read.table(paste(dirname, networks1[i], "/auc_ppis.dat", sep=""))
		e<-read.table(paste(dirname, networks2[j], "/auc_ppis.dat", sep=""))
		d<-d[order(row.names(d)),]
		e<-e[order(row.names(e)),] 
		a<-wilcox.test(d[,k], e[,k], alternative="less", paired=T)
		x[i,j]<-a$p.value
	    }
	}
	names(x)<-networks2
	print(x)
    }


    d<-read.table(paste(dirname, "biana-all", "/auc_ppis.dat", sep=""))
    e<-read.table(paste(dirname, "biana_no_tap-all", "/auc_ppis.dat", sep=""))
    wilcox.test(e[,2], d[,2], alternative="less", paired=T)
    wilcox.test(e[,3], d[,3], alternative="greater", paired=T)


    x<-c()
    y<-1
    for(k in 1:5) {
	print(c("---", methods[k], "---"))
	d<-read.table(paste(dirname, "biana_no_tap-all_seed20/auc_ppis.dat", sep=""))
	e<-read.table(paste(dirname, "biana_no_tap-all_seed20below/auc_ppis.dat", sep=""))
	a<-wilcox.test(d[,k], e[,k])
	#if(a$p.value < 0.05) {
	    print(a$p.value)
	    x[y]<-a$p.value
	    y<-y+1
	    #print(kruskal.test(list(d[,k],e[,k])))
    }
    print(x)
    print(min(x))

    # Significance of AUC change between bPPI n_seed <= 20 and > 20
    s<-read.table(paste(dirname, 'biana_no_tap-all/seeds.dat', sep=""))
    d<-read.table(paste(dirname, 'biana_no_tap-all/auc_ppis.dat', sep=""))

    x<-d[rownames(d) %in% rownames(s)[s$n_seed<=20],]
    y<-d[rownames(d) %in% rownames(s)[s$n_seed>20],]

    wilcox.test(x$ns, x$nd, alternative="less", paired=T)
    wilcox.test(x$nd, y$nd, alternative="greater", paired=F)

}

### NAVLAKHA ###
navlakha_figures <- function() {
    dir.name<-"../data/compare/navlakha/"
    postscript(paste(dir.name, "results.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
    d <- read.table(paste(dir.name, "results.dat", sep=""))
    par(mar=c(6,6,0.2,0))
    methods <- c("ns", "nz", "nd", "ff", "nr")
    colors <- c(2,6,4,5,3)
    chars <- 21:25
    for(i in 1:length(methods)) {
	e<-subset(d, substr(rownames(d),1,2)==methods[i])
	if(i == 1) {
	    plot(e$sens, e$ppv, col=colors[i], pch=chars[i], cex=3, xlim=c(0,1), ylim=c(0,1))
	} else {
	    points(e$sens, e$ppv, col=colors[i], pch=chars[i], cex=3)
	}   
	#lines(e$sens, e$ppv, col=colors[i])
	text(e$sens[length(e$ppv)], e$ppv[length(e$ppv)]-0.05, col=colors[i], toupper(methods[i]))
    }
    dev.off()
}

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

    postscript(paste(dir.name, "4.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
    d<-read.table(paste(dir.name, 'biana_no_tap-all_seed20below/auc_phenotypes.dat', sep=""))
    e<-read.table(paste(dir.name, 'biana_no_tap-all_seed20/auc_phenotypes.dat', sep=""))
    f<-merge(d,e,all=T)
    #xlab="Diseases Groups (w.r.t. number of initially annotated proteins)", 
    barplot(as.matrix(t(f)), beside=T, col=color5, ylim=c(0,100), ylab="AUC(%)", xlab="Number of seeds", names.arg=c("1-19", "20-292"))
    legend(3, 100, methods, fill=color5, bty="n", ncol=3)
    dev.off()

    postscript(paste(dir.name, "S5_a.eps", sep=""), horizontal=F) #width = 480, height = 480, quality = 100)
    d<-read.table("module_summary_top1_mcl_n3.dat", header=T)
    #a<-read.table("module_summary_top1_mcl.dat", header=T)
    #d<-rbind(d,a[a$scoring=="nr",])
    methods<-c("no", "nn", "n3") 
    n<-length(methods)
    x<-1:length(levels(d$phenotype))
    par(mar=c(5,4,4,5))
    e<-c()
    for(i in 1:n) {
	print(methods[i])
	if(i == 1) {
	    plot(x, d[d$scoring==methods[i],]$n_seed_in_modules/d[d$scoring==methods[i],]$n_seed, type="l", col=i+1, ylim=c(0,1), xlab="Phenotypes", ylab="Ratio of covered seeds within modules", xaxt="n")
	} 
	else {
	    lines(x, d[d$scoring==methods[i],]$n_seed_in_modules/d[d$scoring==methods[i],]$n_seed, lty=1, col=i+1)
	}
	if(i != 1) {
	    e<-rbind(e, d[d$scoring==methods[i],]$n_module)
	}
    }
    legend("topright", c("GWAS", "PWAS-1", "PWAS-2"), lty=rep(1, n), col=(1:n)+1)
    lambda<-function(x) { substring(x, 6) }
    text(1:length(levels(d$phenotype)), par("usr")[3]-0.02, srt=30, adj=1, labels=sapply(levels(d$phenotype), lambda), xpd=T, cex=0.7)
    par(new=TRUE)
    e<-as.matrix(e)
    print(dim(e))
    barplot(e, beside=TRUE, horiz=F, col=(2:n)+1, width=0.1, space=c(0,1.3), ylim=c(0,30), axes=F)
    axis(4, 1:15, labels=1:15)
    text(par("usr")[2]+1, 10, srt=90, adj=1, labels="Number of modules using MCL", xpd=T, cex=1)
    dev.off()

    postscript(paste(dir.name, "S5_b.eps", sep=""), horizontal=F) #width = 480, height = 480, quality = 100)
    d<-read.table("module_summary_top1_mcl_n5.dat", header=T)
    methods<-c("no", "nn", "n5") 
    n<-length(methods)
    x<-1:length(levels(d$phenotype))
    par(mar=c(5,4,4,5))
    e<-c()
    for(i in 1:n) {
	print(methods[i])
	if(i == 1) {
	    plot(x, d[d$scoring==methods[i],]$n_seed_in_modules/d[d$scoring==methods[i],]$n_all_in_modules, type="n", col=i+1, ylim=c(0,1), xlab="Phenotypes", ylab="Ratio of seeds within modules", xaxt="n")
	} 
	else {
	    lines(x, d[d$scoring==methods[i],]$n_seed_in_modules/d[d$scoring==methods[i],]$n_all_in_modules, lty=1, col=i+1)
	}
	if(i != 1) {
	    e<-rbind(e, d[d$scoring==methods[i],]$n_module)
	}
    }
    legend("topright", c("GWAS", "PWAS-1", "PWAS-2"), lty=rep(1, n), col=(1:n)+1)
    lambda<-function(x) { substring(x, 6) }
    text(1:length(levels(d$phenotype)), par("usr")[3]-0.02, srt=30, adj=1, labels=sapply(levels(d$phenotype), lambda), xpd=T, cex=0.7)
    par(new=TRUE)
    e<-as.matrix(e)
    print(dim(e))
    barplot(e, beside=TRUE, horiz=F, col=(2:n)+1, width=0.1, space=c(0,1.3), ylim=c(0,30), axes=F)
    axis(4, 1:15, labels=1:15)
    text(par("usr")[2]+1, 10, srt=90, adj=1, labels="Number of modules using MCL", xpd=T, cex=1)
    dev.off()

    postscript(paste(dir.name, "S5_1.eps", sep=""), horizontal=F)#width = 480, height = 480, quality = 100)
    d<-read.table("module_summary_top1_connected.dat", header=T)
    methods<-c("no", "nn", "nr") 
    #methods<-c("no", "nn", "nr", "ff", "nz", "nd", "ns") #levels(d$scoring)
    n<-length(methods)
    x<-1:length(levels(d$phenotype))
    par(mar=c(5,4,4,5))
    e<-c()
    for(i in 1:n) {
	print(methods[i])
	if(i == 1) {
	    plot(x, d[d$scoring==methods[i],]$n_seed, type="l", col=1, ylim=c(0,100), xlab="Phenotypes", ylab="Number of seeds", xaxt="n")
	} 
	e<-rbind(e, d[d$scoring==methods[i],]$n_seed_in_modules)
    }
    legend("topright", c("# of seeds", "GWAS coverage", "PWAS-1 coverage", "PWAS-2 coverage"), lty=rep(1, n+1), col=(0:n)+1)
    lambda<-function(x) { substring(x, 6) }
    text(1:length(levels(d$phenotype)), par("usr")[3]-0.02, srt=30, adj=1, labels=sapply(levels(d$phenotype), lambda), xpd=T, cex=0.7)
    par(new=TRUE)
    e<-as.matrix(e)
    print(dim(e))
    barplot(e, beside=TRUE, horiz=F, col=(1:n)+1, width=0.1, space=c(0,1.3), ylim=c(0,100), axes=F)
    axis(4, 1:10, labels=1:10)
    text(par("usr")[2]+1, 70, srt=90, adj=1, labels="Number of covered seeds within modules using CC", xpd=T, cex=1)
    dev.off()

    postscript(paste(dir.name, "S5_2.eps", sep=""), horizontal=F) #width = 480, height = 480, quality = 100)
    d<-read.table("module_summary_top1_mcl.dat", header=T)
    methods<-c("no", "nn", "nr") 
    #methods<-c("no", "nn", "nr", "ff", "nz", "nd", "ns") #levels(d$scoring)
    n<-length(methods)
    x<-1:length(levels(d$phenotype))
    par(mar=c(5,4,4,5))
    e<-c()
    for(i in 1:n) {
	print(methods[i])
	if(i == 1) {
	    plot(x, d[d$scoring==methods[i],]$n_seed, type="l", col=1, ylim=c(0,100), xlab="Phenotypes", ylab="Number of seeds", xaxt="n")
	} 
	e<-rbind(e, d[d$scoring==methods[i],]$n_seed_in_modules)
    }
    #legend("topright", c("# of seeds", "Connected seeds", "Neighboring seeds", "High scoring seeds"), lty=rep(1, n+1), col=(0:n)+1)
    legend("topright", c("# of seeds", "GWAS coverage", "PWAS-1 coverage", "PWAS-2 coverage"), lty=rep(1, n+1), col=(0:n)+1)
    lambda<-function(x) { substring(x, 6) }
    text(1:length(levels(d$phenotype)), par("usr")[3]-0.02, srt=30, adj=1, labels=sapply(levels(d$phenotype), lambda), xpd=T, cex=0.7)
    par(new=TRUE)
    e<-as.matrix(e)
    print(dim(e))
    barplot(e, beside=TRUE, horiz=F, col=(1:n)+1, width=0.1, space=c(0,1.3), ylim=c(0,100), axes=F)
    axis(4, 1:100, labels=1:100)
    text(par("usr")[2]+1, 70, srt=90, adj=1, labels="Number of covered seeds within modules using MCL", xpd=T, cex=1)
    dev.off()
 
    postscript(paste(dir.name, "S5_3.eps", sep=""), horizontal=F)#width = 480, height = 480, quality = 100)
    d<-read.table("module_summary_top1_connected.dat", header=T)
    methods<-c("no", "nn", "nr") 
    #methods<-c("no", "nn", "nr", "ff", "nz", "nd", "ns") #levels(d$scoring)
    n<-length(methods)
    x<-1:length(levels(d$phenotype))
    par(mar=c(5,4,4,5))
    e<-c()
    for(i in 1:n) {
	print(methods[i])
	if(i == 1) {
	    plot(x, d[d$scoring==methods[i],]$n_seed_in_modules/d[d$scoring==methods[i],]$n_seed, type="l", col=i+1, ylim=c(0,1), xlab="Phenotypes", ylab="Ratio of covered seeds within modules", xaxt="n")
	} 
	else {
	    lines(x, d[d$scoring==methods[i],]$n_seed_in_modules/d[d$scoring==methods[i],]$n_seed, lty=1, col=i+1)
	}
	e<-rbind(e, d[d$scoring==methods[i],]$n_module)
    }
    #legend("topright", methods, lty=rep(1, n), col=(1:n)+1)
    legend("topright", c("GWAS", "PWAS-1", "PWAS-2"), lty=rep(1, n), col=(1:n)+1)
    lambda<-function(x) { substring(x, 6) }
    text(1:length(levels(d$phenotype)), par("usr")[3]-0.02, srt=30, adj=1, labels=sapply(levels(d$phenotype), lambda), xpd=T, cex=0.7)
    par(new=TRUE)
    e<-as.matrix(e)
    print(dim(e))
    barplot(e, beside=TRUE, horiz=F, col=(1:n)+1, width=0.1, space=c(0,1.3), ylim=c(0,10), axes=F)
    axis(4, 1:10, labels=1:10)
    text(par("usr")[2]+1, 7, srt=90, adj=1, labels="Number of modules using CC", xpd=T, cex=1)
    dev.off()

    postscript(paste(dir.name, "S5_4.eps", sep=""), horizontal=F) #width = 480, height = 480, quality = 100)
    d<-read.table("module_summary_top1_mcl.dat", header=T)
    methods<-c("no", "nn", "nr") 
    #methods<-c("no", "nn", "nr", "ff", "nz", "nd", "ns") #levels(d$scoring)
    n<-length(methods)
    x<-1:length(levels(d$phenotype))
    par(mar=c(5,4,4,5))
    e<-c()
    for(i in 1:n) {
	print(methods[i])
	if(i == 1) {
	    plot(x, d[d$scoring==methods[i],]$n_seed_in_modules/d[d$scoring==methods[i],]$n_seed, type="l", col=i+1, ylim=c(0,1), xlab="Phenotypes", ylab="Ratio of covered seeds within modules", xaxt="n")
	} 
	else {
	    lines(x, d[d$scoring==methods[i],]$n_seed_in_modules/d[d$scoring==methods[i],]$n_seed, lty=1, col=i+1)
	}
	e<-rbind(e, d[d$scoring==methods[i],]$n_module)
    }
    #legend("topright", methods, lty=rep(1, n), col=(1:n)+1)
    legend("topright", c("GWAS", "PWAS-1", "PWAS-2"), lty=rep(1, n), col=(1:n)+1)
    lambda<-function(x) { substring(x, 6) }
    text(1:length(levels(d$phenotype)), par("usr")[3]-0.02, srt=30, adj=1, labels=sapply(levels(d$phenotype), lambda), xpd=T, cex=0.7)
    par(new=TRUE)
    e<-as.matrix(e)
    print(dim(e))
    barplot(e, beside=TRUE, horiz=F, col=(1:n)+1, width=0.1, space=c(0,1.3), ylim=c(0,15), axes=F)
    axis(4, 1:15, labels=1:15)
    text(par("usr")[2]+1, 10, srt=90, adj=1, labels="Number of modules using MCL", xpd=T, cex=1)
    dev.off()

    postscript(paste(dir.name, "S5_5.eps", sep=""), horizontal=F)#width = 480, height = 480, quality = 100)
    d<-read.table("module_summary_top1_connected.dat", header=T)
    methods<-c("no", "nn", "nr") 
    #methods<-c("no", "nn", "nr", "ff", "nz", "nd", "ns") #levels(d$scoring)
    n<-length(methods)
    x<-1:length(levels(d$phenotype))
    par(mar=c(5,4,4,5))
    e<-c()
    for(i in 1:n) {
	print(methods[i])
	if(i == 1) {
	    plot(x, d[d$scoring==methods[i],]$n_seed_in_modules/d[d$scoring==methods[i],]$n_all_in_modules, type="l", col=i+1, ylim=c(0,1), xlab="Phenotypes", ylab="Ratio of seeds within modules", xaxt="n")
	} 
	else {
	    lines(x, d[d$scoring==methods[i],]$n_seed_in_modules/d[d$scoring==methods[i],]$n_all_in_modules, lty=1, col=i+1)
	}
	e<-rbind(e, d[d$scoring==methods[i],]$n_module)
    }
    #legend("topright", methods, lty=rep(1, n), col=(1:n)+1)
    legend("topright", c("GWAS", "PWAS-1", "PWAS-2"), lty=rep(1, n), col=(1:n)+1)
    lambda<-function(x) { substring(x, 6) }
    text(1:length(levels(d$phenotype)), par("usr")[3]-0.02, srt=30, adj=1, labels=sapply(levels(d$phenotype), lambda), xpd=T, cex=0.7)
    par(new=TRUE)
    e<-as.matrix(e)
    print(dim(e))
    barplot(e, beside=TRUE, horiz=F, col=(1:n)+1, width=0.1, space=c(0,1.3), ylim=c(0,10), axes=F)
    axis(4, 1:10, labels=1:10)
    text(par("usr")[2]+1, 7, srt=90, adj=1, labels="Number of modules using CC", xpd=T, cex=1)
    dev.off()

    postscript(paste(dir.name, "S5_6.eps", sep=""), horizontal=F) #width = 480, height = 480, quality = 100)
    d<-read.table("module_summary_top1_mcl.dat", header=T)
    methods<-c("no", "nn", "nr") 
    #methods<-c("no", "nn", "nr", "ff", "nz", "nd", "ns") #levels(d$scoring)
    n<-length(methods)
    x<-1:length(levels(d$phenotype))
    par(mar=c(5,4,4,5))
    e<-c()
    for(i in 1:n) {
	print(methods[i])
	if(i == 1) {
	    plot(x, d[d$scoring==methods[i],]$n_seed_in_modules/d[d$scoring==methods[i],]$n_all_in_modules, type="l", col=i+1, ylim=c(0,1), xlab="Phenotypes", ylab="Ratio of seeds within modules", xaxt="n")
	} 
	else {
	    lines(x, d[d$scoring==methods[i],]$n_seed_in_modules/d[d$scoring==methods[i],]$n_all_in_modules, lty=1, col=i+1)
	}
	e<-rbind(e, d[d$scoring==methods[i],]$n_module)
    }
    #legend("topright", methods, lty=rep(1, n), col=(1:n)+1)
    legend("topright", c("GWAS", "PWAS-1", "PWAS-2"), lty=rep(1, n), col=(1:n)+1)
    lambda<-function(x) { substring(x, 6) }
    text(1:length(levels(d$phenotype)), par("usr")[3]-0.02, srt=30, adj=1, labels=sapply(levels(d$phenotype), lambda), xpd=T, cex=0.7)
    par(new=TRUE)
    e<-as.matrix(e)
    print(dim(e))
    barplot(e, beside=TRUE, horiz=F, col=(1:n)+1, width=0.1, space=c(0,1.3), ylim=c(0,15), axes=F)
    axis(4, 1:15, labels=1:15)
    text(par("usr")[2]+1, 10, srt=90, adj=1, labels="Number of modules using MCL", xpd=T, cex=1)
    dev.off()
}

main()
