dir.name<-"../data/summary/"
cols<-c(2,7,3,4,6,5,8) #rainbow(7) #2:8 #heat.colors(7)
color5<-cols[1:5]
color2<-c(3,4)
scoring.methods<-c("N.Score", "N.Zcore", "N.Short", "N.Combo", "F.Flow", "P.Rank", "R.Walk", "N.Prop")
scoring.methods.full<-c("NetScore", "NetZcore", "NetShort", "NetCombo", "Functional\nFlow", "Page\nRank", "Random\nWalk", "Network\n Propagation")
#scoring.methods<-c("N.Score", "N.Zcore", "N.Short", "F.Flow", "P.Rank")
scoring.method.ids<-c("ns", "nz", "nd", "nc3", "ff", "nr", "rw", "np")

main <- function() {
    #manuscript()
    #manuscript2()

    disease_category_figures()

    #case_study_figures()
    #navlakha_figures()
    #aneurysm_figures()
    #yeast_figures()
}

manuscript<-function() {
    manuscript_figures()
    manuscript_tests()
}

manuscript2<-function() {
    manuscript2_figures()
    #manuscript2_tests()
}


###### MANUSCRIPT2 FIGURES ######
manuscript2_figures<-function() {
    neighborhood_figures()
    robustness_figures()
    disease_category_figures()
    module_figures()
    omim_similarity_figures()
}

lambda<-function(x) { x<-substring(x, 6); words<-unlist(strsplit(x, "_")); x<-paste(words, collapse=" "); return(x) }


###### NEIGHBORHOOD FIGURES ######
neighborhood_figures <- function() {

    # Number of genes in neighborhood of alzheimer over randomly permuted bPPI network 
    cairo_ps(paste(dir.name, "Figure P1a.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    par(family = "Arial") 
    par(mar=c(5, 4, 4, 5) + 0.1)
    d<-read.table(paste(dir.name, '../compare/biana_no_tap-omim_alzheimer-nn/permuted/results.dat', sep=""))
    a<-barplot(rbind(d$picked, d$picked_good),beside=T,xlab="Percentage of permuted interactions (%)",ylab="Number of genes", legend.text=c("All genes in n.hood", "AD genes in n.hood"),names.arg=seq(0,80,by=10),ylim=c(0,300))
    par(new=T)
    plot(colMeans(a),100*d$picked_good/d$picked,col=2,xaxt="n",yaxt="n",xlab="",ylab="",type='l',bty="n",ylim=c(0,10))
    axis(4, xpd=T, col=2, col.axis=2)
    mtext("Ratio (%)", side=4, line=3, col=2)
    dev.off()

    # Number of genes in neighborhood of alzheimer over randomly pruned bPPI network 
    cairo_ps(paste(dir.name, "Figure P1b.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    par(family = "Arial") 
    par(mar=c(5, 4, 4, 5) + 0.1)
    d<-read.table(paste(dir.name, '../compare/biana_no_tap-omim_alzheimer-nn/pruned/results.dat', sep=""))
    a<-barplot(rbind(d$picked, d$picked_good),beside=T,xlab="Percentage of pruned interactions (%)",ylab="Number of genes", legend.text=c("All genes in n.hood", "AD genes in n.hood"),names.arg=seq(0,80,by=10),ylim=c(0,300))
    par(new=T)
    plot(colMeans(a),100*d$picked_good/d$picked,col=2,xaxt="n",yaxt="n",xlab="",ylab="",type='l',bty="n",ylim=c(0,10))
    axis(4, xpd=T, col=2, col.axis=2)
    mtext("Ratio (%)", side=4, line=3, col=2)
    dev.off()
}




###### ROBUSTNESS FIGURES ######
robustness_figures <- function() {

    dir.name<-"../data/summary_draft_before_revision/"
    percentages<-seq(0,80,by=10) # c(0, 10, 30, 50, 70)
    scoring.methods<-c("N.Score", "N.Zcore", "N.Short", "F.Flow", "P.Rank")

    # Average AUC (%) at different levels of interaction permutation over OMIM disorders on bPPI network
    cairo_ps(paste(dir.name, "Figure P2a.eps", sep=""), width = 6, height = 6, onefile = TRUE) # ps: horizontal = FALSE, onefile = FALSE, paper = "special")
    par(family = "Arial") 
    e<-c()
    f<-c()
    label<-"Percentage of permuted interactions"
    for(p in percentages) {
	if(p == 0) {
	    d <- read.table(paste(dir.name, 'biana_no_tap-omim/auc_phenotypes.dat', sep=""))
	    dir.name<-"../data/summary_runs_on_random/"
	}
	else {
	    d <- read.table(paste(dir.name, 'biana_no_tap_permuted_p', p, '-omim/auc_phenotypes.dat', sep=""))
	    if(!all(dim(d) == c(100, 5))) {
		stop("Object dimension is different from expected!")
	    }
	}
	m <- mean(d)
	n <- dim(d)[1]
	error <- qt(0.975, df=n-1) * sd(d) / sqrt(n)
	e <- cbind(e,m)
	f <- cbind(f,error)
    }

    b<-barplot(e, beside=TRUE, horiz=F, col = color5[1:length(scoring.methods)], width=0.1, space=c(0,0.4), names.arg=rep('', length(percentages)), ylim=c(0,100), ylab="AUC (%)", xlab=label) 
    segments(b, e-f, b, e+f, col=1, lty=1, lwd=1) 
    segments(b-0.03, e-f, b+0.03, e-f, col=1, lty=1, lwd=0.5)
    segments(b-0.03, e+f, b+0.03, e+f, col=1, lty=1, lwd=0.5)
    text(seq(0.3, 0.55*length(percentages)-0.1, 0.55), par("usr")[3] - 2, srt=45, adj=1, labels=percentages, xpd=T, cex=0.8)
    legend(0.5, 100, scoring.methods, fill=color5, bty="n", ncol=3)
    dev.off()

    # Average AUC (%) at different levels of interaction permutation over OMIM disorders on GOH network
    dir.name<-"../data/summary_draft_before_revision/"
    cairo_ps(paste(dir.name, "Figure PS1a.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    par(family = "Arial") 
    e<-c()
    f<-c()
    label<-"Percentage of permuted interactions"
    for(p in percentages) {
	if(p == 0) {
	    d <- read.table(paste(dir.name, 'goh-omim/auc_phenotypes.dat', sep=""))
	    #d <- d[,1:3]
	    dir.name<-"../data/summary_runs_on_random/"
	}
	else {
	    d <- read.table(paste(dir.name, 'goh_permuted_p', p, '-omim/auc_phenotypes.dat', sep=""))
	    if(!all(dim(d) == c(100, 5))) {
		stop("Object dimension is different from expected!")
	    }
	}
	m <- mean(d)
	n <- dim(d)[1]
	error <- qt(0.975, df=n-1) * sd(d) / sqrt(n)
	e <- cbind(e,m)
	f <- cbind(f,error)
    }
    b<-barplot(e, beside=TRUE, horiz=F, col = color5[1:length(scoring.methods)], width=0.1, space=c(0,0.4), names.arg=rep('', length(percentages)), ylim=c(0,100), ylab="AUC (%)", xlab=label) 
    segments(b, e-f, b, e+f, col=1, lty=1, lwd=1) 
    segments(b-0.03, e-f, b+0.03, e-f, col=1, lty=1, lwd=0.5)
    segments(b-0.03, e+f, b+0.03, e+f, col=1, lty=1, lwd=0.5)
    text(seq(0.3, 0.55*length(percentages)-0.1, 0.55), par("usr")[3] - 2, srt=45, adj=1, labels=percentages, xpd=T, cex=0.8)
    legend(0.5, 100, scoring.methods, fill=color5, bty="n", ncol=3)
    dev.off()

    # Average AUC (%) at different levels of interaction pruning over OMIM disorders on bPPI network
    dir.name<-"../data/summary_draft_before_revision/"
    cairo_ps(paste(dir.name, "Figure P2b.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    par(family = "Arial") 
    e<-c()
    f<-c()
    label<-"Percentage of pruned interactions"
    for(p in percentages) {
	if(p == 0) {
	    d <- read.table(paste(dir.name, 'biana_no_tap-omim/auc_phenotypes.dat', sep=""))
	    dir.name<-"../data/summary_runs_on_random/"
	}
	else {
	    d <- read.table(paste(dir.name, 'biana_no_tap_pruned_p', p, '-omim/auc_phenotypes.dat', sep=""))
	    # Below for Average AUC (%) at different levels of interaction pruning only between non-seeds over OMIM disorders on bPPI network
	    #d <- read.table(paste(dir.name, 'biana_no_tap_pruned_non_seed_interactions_p', p, '-omim/auc_phenotypes.dat', sep=""))
	    if(!all(dim(d) == c(100, 5))) {
		stop("Object dimension is different from expected!")
	    }
	}
	m <- mean(d)
	n <- dim(d)[1]
	error <- qt(0.975, df=n-1) * sd(d) / sqrt(n)
	e <- cbind(e,m)
	f <- cbind(f,error)
    }
    b<-barplot(e, beside=TRUE, horiz=F, col = color5[1:length(scoring.methods)], width=0.1, space=c(0,0.4), names.arg=rep('', length(percentages)), ylim=c(0,100), ylab="AUC (%)", xlab=label)
    segments(b, e-f, b, e+f, col=1, lty=1, lwd=1) 
    segments(b-0.03, e-f, b+0.03, e-f, col=1, lty=1, lwd=0.5)
    segments(b-0.03, e+f, b+0.03, e+f, col=1, lty=1, lwd=0.5)
    text(seq(0.3, 0.55*length(percentages)-0.1, 0.55), par("usr")[3] - 2, srt=45, adj=1, labels=percentages, xpd=T, cex=0.8)
    legend(0.5, 100, scoring.methods, fill=color5, bty="n", ncol=3)
    dev.off()

    # Average AUC (%) at different levels of interaction pruning over OMIM disorders on GOH network
    dir.name<-"../data/summary_draft_before_revision/"
    cairo_ps(paste(dir.name, "Figure PS1b.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    par(family = "Arial") 
    e<-c()
    f<-c()
    label<-"Percentage of pruned interactions"
    for(p in percentages) {
	if(p == 0) {
	    d <- read.table(paste(dir.name, 'goh-omim/auc_phenotypes.dat', sep=""))
	    dir.name<-"../data/summary_runs_on_random/"
	}
	else {
	    d <- read.table(paste(dir.name, 'goh_pruned_p', p, '-omim/auc_phenotypes.dat', sep=""))
	    if(!all(dim(d) == c(100, 5))) {
		stop("Object dimension is different from expected!")
	    }
	}
	m <- mean(d)
	n <- dim(d)[1]
	error <- qt(0.975, df=n-1) * sd(d) / sqrt(n)
	e <- cbind(e,m)
	f <- cbind(f,error)
    }
    b<-barplot(e, beside=TRUE, horiz=F, col = color5[1:length(scoring.methods)], width=0.1, space=c(0,0.4), names.arg=rep('', length(percentages)), ylim=c(0,100), ylab="AUC (%)", xlab=label)
    segments(b, e-f, b, e+f, col=1, lty=1, lwd=1) 
    segments(b-0.03, e-f, b+0.03, e-f, col=1, lty=1, lwd=0.5)
    segments(b-0.03, e+f, b+0.03, e+f, col=1, lty=1, lwd=0.5)
    text(seq(0.3, 0.55*length(percentages)-0.1, 0.55), par("usr")[3] - 2, srt=45, adj=1, labels=percentages, xpd=T, cex=0.8)
    legend(0.5, 100, scoring.methods, fill=color5, bty="n", ncol=3)
    dev.off()

    # Average AUC (%) for two groups of OMIM disorders based on number of seeds on bPPI network 
    dir.name<-"../data/summary_draft_before_revision/"
    cairo_ps(paste(dir.name, "Figure P3a.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    par(family = "Arial") 

    s<-read.table(paste(dir.name, 'biana_no_tap-omim/seeds.dat', sep=""))
    d<-read.table(paste(dir.name, 'biana_no_tap-omim/auc_ppis.dat', sep=""))
    cutoff<-median(s$n_seed)
    #cutoff<-median(s$n_path)
    x<-d[rownames(s)[s$n_seed<=cutoff],]
    y<-d[rownames(s)[s$n_seed>cutoff],]

    coords <- seq(0.6, 5.0, by=1.1)
    boxplot(x, col=color2[1], at=coords, boxwex=0.4, pars=list(xaxt="n"), xlab="Prediction method", ylab="AUC(%)", ylim=c(0,100))
    axis(1,tick=F, labels=scoring.methods[1:5], at=coords+0.2)
    coords <- coords + 0.4
    boxplot(y, col=color2[2], at=coords, boxwex=0.4, pars=list(xaxt="n"), names=F, add=T)
    legend(1.5, 5, c(paste("# of seeds < ", cutoff, sep=""), paste("# of seeds >= ", cutoff, sep="")), fill=color2, bty="n", horiz=T)
    #legend(1.0, 15, c("Group with shorter seed-connecting paths", "Group with longer seed-connecting paths"), fill=color2, bty="n") #1.5, 5, horiz=T, cex=0.9)
    dev.off()

    # Average AUC (%) for two groups of OMIM disorders based on number of seeds on GOH network
    dir.name<-"../data/summary_draft_before_revision/"
    cairo_ps(paste(dir.name, "Figure PS2a.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    par(family = "Arial") 

    s<-read.table(paste(dir.name, 'goh-omim/seeds.dat', sep=""))
    d<-read.table(paste(dir.name, 'goh-omim/auc_ppis.dat', sep=""))

    cutoff<-median(s$n_seed)
    x<-d[rownames(s)[s$n_seed<=cutoff],]
    y<-d[rownames(s)[s$n_seed>cutoff],]

    coords <- seq(0.6, 5.0, by=1.1)
    boxplot(x,col=3, at=coords, boxwex=0.4, pars=list(xaxt="n"), xlab="Prediction method", ylab="AUC(%)", ylim=c(0,100))
    axis(1,tick=F,labels=scoring.methods,at=coords+0.2)
    coords <- coords + 0.4
    boxplot(y,col=4, at=coords, boxwex=0.4, pars=list(xaxt="n"), names=F, add=T)
    legend(1.5, 5, c(paste("# of seeds < ", cutoff, sep=""), paste("# of seeds >= ", cutoff, sep="")), fill=color2, bty="n", horiz=T)
    dev.off()

    # Average AUC (%) at different levels of seed permutation over OMIM disorders on bPPI network
    dir.name<-"../data/summary_draft_before_revision/"
    cairo_ps(paste(dir.name, "Figure P3b.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    par(family = "Arial") 
    e<-c()
    f<-c()
    label<-"Percentage of perturbed seeds"
    for(p in percentages) {
	if(p == 0) {
	    d <- read.table(paste(dir.name, 'biana_no_tap-omim/auc_phenotypes.dat', sep=""))
	    dir.name<-"../data/summary_runs_on_random/"
	}
	else {
	    d<-read.table(paste(dir.name, 'biana_no_tap-omim_perturbed_p', p, '/auc_perturbed.dat', sep=""))
	    if(!all(dim(d) == c(100, 5))) {
		stop("Object dimension is different from expected!")
	    }
	}
	m <- mean(d)
	n <- dim(d)[1]
	error <- qt(0.975, df=n-1) * sd(d) / sqrt(n)
	e <- cbind(e,m)
	f <- cbind(f,error)
    }
    b<-barplot(e, beside=TRUE, horiz=F, col = color5[1:length(scoring.methods)], width=0.1, space=c(0,0.4), names.arg=rep('', length(percentages)), ylim=c(0,100), ylab="AUC (%)", xlab=label)
    segments(b, e-f, b, e+f, col=1, lty=1, lwd=1) 
    segments(b-0.03, e-f, b+0.03, e-f, col=1, lty=1, lwd=0.5)
    segments(b-0.03, e+f, b+0.03, e+f, col=1, lty=1, lwd=0.5)
    text(seq(0.3, 0.55*length(percentages)-0.1, 0.55), par("usr")[3] - 2, srt=45, adj=1, labels=percentages, xpd=T, cex=0.8)
    legend(0.5, 100, scoring.methods, fill=color5, bty="n", ncol=3)
    dev.off()

    # Average AUC (%) at different levels of seed permutation over OMIM disorders on GOH network
    dir.name<-"../data/summary_draft_before_revision/"
    cairo_ps(paste(dir.name, "Figure PS2b.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    par(family = "Arial") 
    e<-c()
    f<-c()
    label<-"Percentage of perturbed seeds"
    for(p in percentages) {
	if(p == 0) {
	    d <- read.table(paste(dir.name, 'goh-omim/auc_phenotypes.dat', sep=""))
	    dir.name<-"../data/summary_runs_on_random/"
	}
	else {
	    d<-read.table(paste(dir.name, 'goh-omim_perturbed_p', p, '/auc_perturbed.dat', sep=""))
	    if(!all(dim(d) == c(100, 5))) {
		stop("Object dimension is different from expected!")
	    }
	}
	m <- mean(d)
	n <- dim(d)[1]
	error <- qt(0.975, df=n-1) * sd(d) / sqrt(n)
	e <- cbind(e,m)
	f <- cbind(f,error)
    }
    b<-barplot(e, beside=TRUE, horiz=F, col = color5[1:length(scoring.methods)], width=0.1, space=c(0,0.4), names.arg=rep('', length(percentages)), ylim=c(0,100), ylab="AUC (%)", xlab=label)
    segments(b, e-f, b, e+f, col=1, lty=1, lwd=1) 
    segments(b-0.03, e-f, b+0.03, e-f, col=1, lty=1, lwd=0.5)
    segments(b-0.03, e+f, b+0.03, e+f, col=1, lty=1, lwd=0.5)
    text(seq(0.3, 0.55*length(percentages)-0.1, 0.55), par("usr")[3] - 2, srt=45, adj=1, labels=percentages, xpd=T, cex=0.8)
    legend(0.5, 100, scoring.methods, fill=color5, bty="n", ncol=3)
    dev.off()
}




###### DISEASE CATEGORY FIGURES ######
disease_category_figures<-function() {

    #common.up<-c('omim_anemia', 'omim_breast_cancer', 'omim_leukemia', 'omim_lymphoma', 'omim_systemic_lupus_erythematosus')
    #common.down<-c('omim_asthma', 'omim_ataxia', 'omim_cataract', 'omim_neuropathy', 'omim_schizophrenia', 'omim_spastic_paraplegia')
    #cols<-2:23
    method<-"ns"
    percentages<-seq(0,80,by=10)
    dir.summary<-"../data/summary_draft_before_revision/"

    # Change on AUC of diseases over permuted and pruned interactions
    dir.name<-dir.summary
    label<-"Percentage of perturbated interactions"
    container.permuted<-data.frame()
    container.pruned<-data.frame()
    for(p in percentages) {
	if(p == 0) {
	    d<-read.table(paste(dir.name, 'biana_no_tap-omim/auc_ppis.dat', sep=""))
	    dir.name<-"../data/summary_runs_on_random/"
	    container.permuted<-d[order(rownames(d)), method]
	    container.pruned<-d[order(rownames(d)), method]
	    phenotypes<-rownames(d)[order(rownames(d))]
	}
	else {
	    d<-read.table(paste(dir.name, 'biana_no_tap_permuted_p', p, '-omim/auc_ppis.dat', sep=""))
	    container.permuted<-cbind(container.permuted, d[order(rownames(d)), method])
	    d<-read.table(paste(dir.name, 'biana_no_tap_pruned_p', p, '-omim/auc_ppis.dat', sep=""))
	    container.pruned<-cbind(container.pruned, d[order(rownames(d)), method])
	}
    }	
    colnames(container.permuted)<-percentages
    colnames(container.pruned)<-percentages
    rownames(container.permuted)<-phenotypes
    rownames(container.pruned)<-phenotypes

    # Take common to both permuted and pruned

    # Robustness cutoff: auc-(auc-50)/2
    common.up<-c()
    common.down<-c()
    for(pheno in phenotypes) {
	d<-container.permuted[pheno,] 
	if(d["0"] <= 50) { next }
	cutoff<-d["0"]-(d["0"]-50)/2
	idx<-match(F, d>cutoff)
	if(is.na(idx)) { idx<-length(d) }
	idx.permuted<-idx
	d<-container.pruned[pheno,] 
	idx<-match(F, d>cutoff)
	if(is.na(idx)) { idx<-length(d) }
	idx.pruned<-idx
	if(idx.permuted < 6 & idx.pruned < 6) {
	    common.down<-c(common.down, pheno)
	} else if(idx.permuted > 6 & idx.pruned > 6) {
	    common.up<-c(common.up, pheno)
	}
    }
    container<-container.permuted
    container.permuted.up<-container[common.up,]
    container.permuted.down<-container[common.down,]
    container<-container.pruned
    container.pruned.up<-container[common.up,]
    container.pruned.down<-container[common.down,]

    write.table(cbind(container.permuted[,1], container.permuted[,1]-(container.permuted[,1]-50)/2, container.permuted[,6], container.pruned[,6]), paste(dir.summary, 'test.dat', sep=""), sep="\t", col.names=c("AUC", "critical AUC", "50% swap AUC", "50% deletion AUC"))

    # Before robustness cutoff: median(aucs)
    container<-container.permuted
    a<-container[, "0"]-container[, "50"]
    cutoff<-median(a) 
    x<-a[a<cutoff] # up
    #container.permuted.up<-container[names(x),] # for median cutoff
    y<-a[a>=cutoff] # down
    #container.permuted.down<-container[names(y),]

    container<-container.pruned
    a<-container[, "0"]-container[, "50"]
    cutoff<-median(a)
    x<-a[a<cutoff] # up
    #container.pruned.up<-container[names(x),]
    y<-a[a>=cutoff] # down
    #container.pruned.down<-container[names(y),]

    common.up<-intersect(rownames(container.permuted.up), rownames(container.pruned.up))
    #common.up<-rownames(container.pruned.up)
    #common.up<-rownames(container.permuted.up)
    common.up<-common.up[order(common.up)]
    common.down<-intersect(rownames(container.permuted.down), rownames(container.pruned.down))
    #common.down<-rownames(container.pruned.down)
    #common.down<-rownames(container.permuted.down)
    common.down<-common.down[order(common.down)]

    common<-union(common.up, common.down)
    non.common<-phenotypes[!(phenotypes %in% common)]
    # Now two categories: robust vs non-robust
    common.down<-union(common.down, non.common)
    container<-container.permuted
    container.permuted.down<-container[common.down,]
    container<-container.pruned
    container.pruned.down<-container[common.down,]

    dir.name<-dir.summary
    cairo_ps(paste(dir.name, "Figure P4a.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    par(family = "Arial") 

    for(i in 1:length(common.up)) {
	if(i==1) {
	    container<-container.permuted.up[common.up,]
	    plot(percentages, container[i,], type='l', col=cols[i], xaxt="n", xlab=label, ylab="AUC(%)", ylim=c(0,100))
	    container<-container.pruned.up[common.up,]
	    #plot(percentages, container[i,], type='l', col=cols[i], xaxt="n", xlab=label, ylab="AUC(%)", ylim=c(0,100))
	    lines(percentages, container[i,], col=cols[i], lty=2)
	} else {
	    container<-container.permuted.up[common.up,]
	    lines(percentages, container[i,], col=cols[i])
	    container<-container.pruned.up[common.up,]
	    lines(percentages, container[i,], col=cols[i], lty=2)
	}
    }
    axis(1,tick=F, labels=percentages, at=percentages)
    legend(0.1, 40, sapply(rownames(container), lambda), lty=rep(1, dim(container)[1]), col=cols, bty="n", ncol=1)
    legend(40, 20, c("Interaction permutation", "Interaction pruning"), lty=c(1,2), col=c(1,1), bty="n")
    dev.off()

    dir.name<-dir.summary
    cairo_ps(paste(dir.name, "Figure P4b.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    par(family = "Arial") 

    for(i in 1:length(common.down)) {
	if(i==1) {
	    container<-container.permuted.down[common.down,]
	    plot(percentages, container[i,], type='l', col=cols[i], xaxt="n", xlab=label, ylab="AUC(%)", ylim=c(0,100))
	    container<-container.pruned.down[common.down,]
	    #plot(percentages, container[i,], type='l', col=cols[i], xaxt="n", xlab=label, ylab="AUC(%)", ylim=c(0,100))
	    lines(percentages, container[i,], col=cols[i], lty=2)
	} else {
	    container<-container.permuted.down[common.down,]
	    lines(percentages, container[i,], col=cols[i])
	    container<-container.pruned.down[common.down,]
	    lines(percentages, container[i,], col=cols[i], lty=2)
	}
    }
    axis(1,tick=F, labels=percentages, at=percentages)
    legend(0.1, 40, sapply(rownames(container), lambda), lty=rep(1, dim(container)[1]), col=cols, bty="n", ncol=1)
    legend(40, 20, c("Interaction permutation", "Interaction pruning"), lty=c(1,2), col=c(1,1), bty="n")
    dev.off()

    # Functional enrichment of modules in robust vs non-robust diseases
    #method<-"nn"
    dir.name<-"../data/module/"
    d<-read.table(paste(dir.name, "biana_no_tap-omim/module_summary.dat", sep=""), header=T)
    e<-d[d$scoring==method & d$phenotype %in% common.up,]
    f<-d[d$scoring==method & d$phenotype %in% common.down,]
    #g<-d[d$scoring==method & d$phenotype %in% non.common,]

    dir.name<-dir.summary
    #labels<-c("Robust", "Uncharacterized", "Non-robust")
    labels<-c("Robust", "Non-robust")
    cairo_ps(paste(dir.name, "Figure P5a.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    par(family = "Arial") 
    #boxplot(e$n_module, g$n_module, f$n_module, col=8, names=labels, xlab="Disease category", ylab="Number of modules", ylim=c(0,20))
    boxplot(e$n_module, f$n_module, col=8, names=labels, xlab="Disease category", ylab="Number of modules", ylim=c(0,15))
    #boxplot(e$n_module, f$n_module, col=8, boxwex=0.4, at=c(0.75, 1.25), xlim=c(0.5,1.5), names=labels, xlab="Disease category", ylab="Number of modules", ylim=c(0,15))
    dev.off()
    a<-wilcox.test(e$n_module, f$n_module)
    print(c("module:", a$p.value))

    cairo_ps(paste(dir.name, "Figure P5b.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    par(family = "Arial") 
    #boxplot(100*e$n_seed_go_in_modules/e$n_seed_go, 100*g$n_seed_go_in_modules/g$n_seed_go, 100*f$n_seed_go_in_modules/f$n_seed_go, col=8, names=labels, xlab="Disease category", ylab="Average seed GO term enrichment among the modules (%)", ylim=c(0,100))
    boxplot(100*e$n_seed_go_in_modules/e$n_seed_go, 100*f$n_seed_go_in_modules/f$n_seed_go, col=8, names=labels, xlab="Disease category", ylab="Average seed GO term enrichment among the modules (%)", ylim=c(0,100))
    dev.off()
    a<-wilcox.test(e$n_seed_go_in_modules/e$n_seed_go, f$n_seed_go_in_modules/f$n_seed_go)
    print(c("enrichment: ", a$p.value)) #, mean(e$n_seed_go_in_modules/e$n_seed_go, na.rm=T), mean(f$n_seed_go_in_modules/f$n_seed_go, na.rm=T)))

    cairo_ps(paste(dir.name, "Figure P5c.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    par(family = "Arial") 
    #boxplot(100*e$n_seed_go_in_modules/e$n_go_in_modules, 100*g$n_seed_go_in_modules/g$n_go_in_modules, 100*f$n_seed_go_in_modules/f$n_seed_go, col=8, names=labels, xlab="Disease category", ylab="Average coverage of seed GO terms among the modules (%)", ylim=c(0,100))
    boxplot(100*e$n_seed_go_in_modules/e$n_go_in_modules, 100*f$n_seed_go_in_modules/f$n_seed_go, col=8, names=labels, xlab="Disease category", ylab="Average coverage of seed GO terms among the modules (%)", ylim=c(0,100))
    dev.off()
    a<-wilcox.test(e$n_seed_go_in_modules/e$n_go_in_modules, f$n_seed_go_in_modules/f$n_go_in_modules)
    print(c("ratio: ", a$p.value))

    cairo_ps(paste(dir.name, "Figure P5d.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    par(family = "Arial") 
    #boxplot(e$n_seed_go, g$n_seed_go, f$n_seed_go, col=8, names=labels, xlab="Disease category", ylab="Number of seed GO terms", ylim=c(0,150))
    boxplot(e$n_seed_go, f$n_seed_go, col=8, names=labels, xlab="Disease category", ylab="Number of seed GO terms", ylim=c(0,150))
    dev.off()
    a<-wilcox.test(e$n_seed_go, f$n_seed_go)
    print(c("seed go: ", a$p.value))

    cairo_ps(paste(dir.name, "Figure P5e.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    par(family = "Arial") 
    #boxplot(e$n_go_in_modules, g$n_go_in_modules, f$n_go_in_modules, col=8, names=labels, xlab="Disease category", ylab="Number of GO terms in the modules", ylim=c(0,250))
    boxplot(e$n_go_in_modules, f$n_go_in_modules, col=8, names=labels, xlab="Disease category", ylab="Number of GO terms in the modules", ylim=c(0,350))
    dev.off()
    a<-wilcox.test(e$n_go_in_modules, f$n_go_in_modules)
    print(c("go: ", a$p.value))

    cairo_ps(paste(dir.name, "Figure P5f.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    par(family = "Arial") 
    #boxplot(e$n_go_in_modules/e$n_module, g$n_go_in_modules/g$n_module, f$n_go_in_modules/f$n_module, col=8, names=labels, xlab="Disease category", ylab="Average number of GO terms per module", ylim=c(0,100))
    boxplot(e$n_go_in_modules/e$n_module, f$n_go_in_modules/f$n_module, col=8, names=labels, xlab="Disease category", ylab="Average number of GO terms per module", ylim=c(0,100))
    dev.off()
    a<-wilcox.test(e$n_go_in_modules/e$n_module, f$n_go_in_modules/f$n_module)
    print(c("go per module: ", a$p.value))

    cairo_ps(paste(dir.name, "Figure P5g.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    par(family = "Arial") 
    #boxplot(e$n_seed_go_in_modules, g$n_seed_go_in_modules, f$n_seed_go_in_modules, col=8, names=labels, xlab="Disease category", ylab="Number of seed GO terms in the modules", ylim=c(0,100))
    boxplot(e$n_seed_go_in_modules, f$n_seed_go_in_modules, col=8, names=labels, xlab="Disease category", ylab="Number of seed GO terms in the modules", ylim=c(0,150))
    dev.off()
    a<-wilcox.test(e$n_seed_go_in_modules, f$n_seed_go_in_modules)
    print(c("seed go in modules: ", a$p.value)) 


    # n_seed and n_path in robust vs non-robust diseases
    dir.name<-dir.summary
    s<-read.table(paste(dir.name, 'biana_no_tap-omim/seeds.dat', sep=""))
    cairo_ps(paste(dir.name, "Figure PS3a.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    par(family = "Arial") 
    #boxplot(s[common.up,"n_seed"], s[non.common,"n_seed"], s[common.down,"n_seed"], col=8, names=labels, xlab="Disease category", ylab="Number of seeds", ylim=c(0,120))
    boxplot(s[common.up,"n_seed"], s[common.down,"n_seed"], col=8, names=labels, xlab="Disease category", ylab="Number of seeds", ylim=c(0,120))
    dev.off()
    a<-wilcox.test(s[common.up,"n_seed"], s[common.down,"n_seed"])
    print(c("seed:", a$p.value))
    cairo_ps(paste(dir.name, "Figure PS3b.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    par(family = "Arial") 
    #boxplot(s[common.up,"n_path"], s[non.common,"n_path"], s[common.down,"n_path"], col=8, names=labels, xlab="Disease category", ylab="Average length of seed connecting paths", ylim=c(0,5))
    boxplot(s[common.up,"n_path"], s[common.down,"n_path"], col=8, names=labels, xlab="Disease category", ylab="Average length of seed connecting paths", ylim=c(0,5))
    dev.off()
    a<-wilcox.test(s[common.up,"n_path"], s[common.down,"n_path"])
    print(c("path:", a$p.value))


    # comparison of the age of the genes
    cairo_ps(paste(dir.name, "Figure PS6.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    d<-read.table("../data/omim/2009_Aug_27/extended/age_category.dat")
    d<-d*100
    coords <- seq(0.7, 3.7, by=1)
    boxplot(d[common.up,], col="green", at=coords, boxwex=0.3, pars=list(xaxt="n"), xlab="Phylogenetic category", ylab="Ratio of genes (%)", ylim=c(0,100))
    coords <- coords + 0.32
    axis(1,tick=F,labels=c("Eukarya", "Metazoans", "Vertebrates", "Mammals"), at=coords)
    #boxplot(d[non.common,], col="blue", at=coords, boxwex=0.3, pars=list(xaxt="n"),add=T)
    #coords <- coords + 0.32
    boxplot(d[common.down,], col="red", at=coords, boxwex=0.3, pars=list(xaxt="n"),add=T)
    legend("topright", c("robust", "uncharacterized", "non-robust"), fill=c("green", "blue", "red"), bty="n", horiz=F)
    dev.off()
    a<-wilcox.test(d[common.up,1], d[common.down,1])
    print(c("eukarya:", a$p.value))
    a<-wilcox.test(d[common.up,2], d[common.down,2])
    print(c("metzoans:", a$p.value))
    a<-wilcox.test(d[common.up,3], d[common.down,3])
    print(c("vertebrates:", a$p.value))
    a<-wilcox.test(d[common.up,4], d[common.down,4])
    print(c("mammals:", a$p.value))

    return() 
    
    # Below gives error in batch execution due to margins of the heatmaps but images can be generated in interactive mode

    library(gplots)
    library(RColorBrewer)
    #cols<-c(rep("red", length(common.up)), rep("grey", length(non.common)), rep("green", length(common.down)))
    cols<-c(rep("green", length(common.up)), rep("lightgrey", length(common.down)))
    val.cols <- brewer.pal(9,"Blues") #greenred
    # Common functions in pairwise  enrichment
    dir.name<-dir.summary
    cairo_ps(paste(dir.name, "Figure PS5a.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    dir.name<-"../data/module/"
    d<-read.table(paste(dir.name, "biana_no_tap-omim/functional_similarity_matrix.dat", sep=""), header=T)
    #e<-d[c(common.up, non.common, common.down),c(common.up, non.common, common.down)]
    e<-d[c(common.up, common.down),c(common.up, common.down)]
    e[upper.tri(e)]<-NA
    heatmap.2(as.matrix(e), revC=F, trace="none", dendrogram="none", col=val.cols, density.info="none", margins=c(12,12), RowSideColors=cols, Rowv=NA, Colv=NA, labRow=sapply(rownames(e), lambda), labCol=sapply(rownames(e), lambda), keysize=0.1)
    par(fig=c(0.9, 1.0, 0.0, 0.1), new = T)
    par(xaxt="n", yaxt="n")
    image(as.matrix(1:5),col=greenred(5), mar=c(2,2))
    text(c(0, 0.5, 1), rep(0, 3), c(0, 0.5, 1), col="white")
    #heatmap.2(as.matrix(e), revC=F, trace="none", dendrogram="none", col=greenred, density.info="none", margins=c(10,10), keysize=0.8, RowSideColors=cols, Rowv=NA, Colv=NA, labRow=NA, labCol=NA) #labCol=sapply(rownames(e), lambda))
    x1<-0.11 #par("usr")[1]
    x2<-0.8 #par("usr")[2]
    y1<-0.1
    y2<-0.94
    n<-dim(e)[1]#+8
    interval<-(x2-x1)/(n-1)
    x<-seq(x1, x2, by=interval)
    interval<-(y2-y1)/(n-1)
    y<-rev(seq(y1, y2, by=interval))
    text(x,y,sapply(rownames(e), lambda), cex=0.8, srt=30, pos=4)
    dev.off()

    dir.name<-dir.summary
    cairo_ps(paste(dir.name, "Figure PS5b.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    dir.name<-"../data/module/"
    d<-read.table(paste(dir.name, "biana_no_tap-omim/gene_similarity_matrix.dat", sep=""), header=T)
    #e<-d[c(common.up, non.common, common.down),c(common.up, non.common, common.down)]
    e<-d[c(common.up, common.down),c(common.up, common.down)]
    e[upper.tri(e)]<-NA
    heatmap.2(as.matrix(e), revC=F, trace="none", dendrogram="none", col=val.cols, density.info="none", margins=c(12,12), RowSideColors=cols, Rowv=NA, Colv=NA, labRow=sapply(rownames(e), lambda), labCol=sapply(rownames(e), lambda), keysize=0.1)
    dev.off()

    dir.name<-"../data/module/"
    d<-read.table(paste(dir.name, "biana_no_tap-omim/phenotype_vs_functions.dat", sep=""), header=T)
    dir.name<-dir.summary
    cairo_ps(paste(dir.name, "Figure PS5c.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    #selected<-c("omim_breast_cancer", "omim_lung_cancer", "omim_prostate_cancer", "omim_leukemia", "omim_diabetes", "omim_obesity", "omim_insulin")
    selected<-common.up
    e<-d[,selected]
    e<-e[rowSums(e)>2,]
    e<-e[,colSums(e)>0]
    lambda2<-function(x) { words<-unlist(strsplit(x, "_")); x<-paste(words, collapse=" "); return(x) }
    heatmap.2(as.matrix(e), revC=F, trace="none", dendrogram="none", col=val.cols, density.info="none", margins=c(7,30), keysize=0.1, labCol=sapply(colnames(e), lambda), labRow=sapply(rownames(e), lambda2), cexCol=0.9)
    dev.off()

    cairo_ps(paste(dir.name, "Figure PS5d.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    selected<-common.down
    e<-d[,selected]
    e<-e[rowSums(e)>2,]
    e<-e[,colSums(e)>0]
    lambda2<-function(x) { words<-unlist(strsplit(x, "_")); x<-paste(words, collapse=" "); return(x) }
    heatmap.2(as.matrix(e), revC=F, trace="none", dendrogram="none", col=val.cols, density.info="none", margins=c(7,30), keysize=0.1, labCol=sapply(colnames(e), lambda), labRow=sapply(rownames(e), lambda2), cexCol=0.7, cexRow=0.7)
    dev.off()

    dir.name<-"../data/module/"
    d<-read.table(paste(dir.name, "biana_no_tap-omim/phenotype_vs_genes.dat", sep=""), header=T)
    dir.name<-dir.summary
    cairo_ps(paste(dir.name, "Figure PS5e.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    selected<-common.up
    e<-d[,selected]
    e<-e[rowSums(e)>1,]
    e<-e[,colSums(e)>0]
    print(c("common.up", rownames(e)))
    heatmap.2(as.matrix(e), revC=F, trace="none", dendrogram="none", col=val.cols, density.info="none", margins=c(7,30), keysize=0.1, labCol=sapply(colnames(e), lambda), labRow=rownames(e), cexCol=0.9)
    dev.off()

    cairo_ps(paste(dir.name, "Figure PS5f.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    selected<-common.down
    e<-d[,selected]
    e<-e[rowSums(e)>1,]
    e<-e[,colSums(e)>0]
    print(c("common.down", rownames(e)))
    heatmap.2(as.matrix(e), revC=F, trace="none", dendrogram="none", col=val.cols, density.info="none", margins=c(7,30), keysize=0.1, labCol=sapply(colnames(e), lambda), labRow=rownames(e), cexCol=0.9)
    dev.off()

    return()

    # Was trying an alternative definition
    # Robustness cutoff: median of sum of auc_i-(auc_0/2) for each disease where i = 0 - 80 (%)
    n<-9

    d<-container.permuted[,]-container.permuted[,1]/2
    d<-rowSums(d[,1:n])/n
    cutoff.d<-median(d)

    e<-container.pruned[,]-container.pruned[,1]/2
    e<-rowSums(e[,1:n])/n
    cutoff.e<-median(e)

    common.up<-names(which(d>cutoff.d & e>cutoff.e))
    common.down<-names(which(d<cutoff.d & e<cutoff.e))

    common<-union(common.up, common.down)
    non.common<-phenotypes[!(phenotypes %in% common)]

    container<-container.permuted
    container.permuted.up<-container[common.up,]
    container.permuted.down<-container[common.down,]
    container<-container.pruned
    container.pruned.up<-container[common.up,]
    container.pruned.down<-container[common.down,]

}



###### MODULE FIGURES ######
module_figures<-function() {
    dir.name<-"../data/module/"
    d<-read.table(paste(dir.name, "biana_no_tap-omim/module_summary.dat", sep=""), header=T)
    #d2<-read.table(paste(dir.name, "module_summary_top5_mcl_n5union.dat", sep=""), header=T)
    method.ids<-c("nn", "ns", "nz", "nd", "ff", "nr") #, "nc") 
    method.names<-c("N.hood", scoring.methods) #"NetScore", "NetZcore", "NetShort", "FunctionalFlow", "ToppGene") #scoring.methods) #, "NetComb.")
    #method.names<-c("Neighborhood", "Func.Flow", "ToppGene", "NetScore", "NetZcore", "NetShort")
    n<-length(method.ids)
    e.ratio<-c()
    e.go<-c()
    #e.genes<-c()
    e.modules<-c()
    #e.seeds<-c()
    for(i in 1:n) {
	print(method.ids[i])
	e.modules<-rbind(e.modules, d[d$scoring==method.ids[i],]$n_module)
	#e.ratio<-rbind(e.ratio, 100*d[d$scoring==method.ids[i],]$ratio) 
	e.ratio<-rbind(e.ratio, 100*d[d$scoring==method.ids[i],]$n_seed_go_in_modules/d[d$scoring==method.ids[i],]$n_go_in_modules)
	e.go<-rbind(e.go, 100*d[d$scoring==method.ids[i],]$n_seed_go_in_modules/d[d$scoring==method.ids[i],]$n_seed_go)
	#e.genes<-rbind(e.genes, d2[d2$scoring==method.ids[i],]$n_all_in_modules)
	#e.seeds<-rbind(e.seeds, d2[d2$scoring==method.ids[i],]$n_seed_in_modules)
    }
    dir.name<-"../data/summary/"
    cairo_ps(paste(dir.name, "Figure PS4a.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    par(family = "Arial") 
    #col = c(8, cols), names=method.names
    boxplot(t(e.modules), col=8, names=NA, xlab="Prediction method", ylab="Number of modules", ylim=c(0,20), pars=list(xaxt="n"))
    range<-par("usr")[2] - par("usr")[1]
    text(par("usr")[1]-0.23+(1:6)*(range/6.5), par("usr")[3]-0.5, labels=method.names, cex=0.9, xpd=NA)
    #axis(1, par("usr")[1]-0.23+(1:6)*(range/6.5), labels=method.names, cex=0.9)
    dev.off()
    cairo_ps(paste(dir.name, "Figure PS4b.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    par(family = "Arial") 
    boxplot(t(e.ratio),col=8, names=NA, xlab="Prediction method", ylab="Average seed GO term enrichment of the modules (%)", ylim=c(0,70), pars=list(xaxt="n")) 
    range<-par("usr")[2] - par("usr")[1]
    text(par("usr")[1]-0.23+(1:6)*(range/6.5), par("usr")[3]-2, labels=method.names, cex=0.9, xpd=NA)
    dev.off()
    cairo_ps(paste(dir.name, "Figure PS4c.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    par(family = "Arial") 
    boxplot(t(e.go),col=8, names=NA, xlab="Prediction method", ylab="Coverage of seed go terms within the modules (%)", ylim=c(0,100), pars=list(xaxt="n")) 
    range<-par("usr")[2] - par("usr")[1]
    text(par("usr")[1]-0.23+(1:6)*(range/6.5), par("usr")[3]-3, labels=method.names, cex=0.9, xpd=NA)
    dev.off()
}

###### OMIM DISEASE SIMILARITY FIGURES #####
omim_similarity_figures<-function() {
    library(RColorBrewer)
    val.cols <- brewer.pal(9,"Blues") 
    tiff("omim.tif", width=2000, height=2000, res=300)
    d<-read.table("../data/omim/2009_Aug_27/similarity.dat")
    heatmap(as.matrix(d), revC=T, col=val.cols, margins=c(9,9), Rowv=NA, Colv=NA, scale="none") # labRow=scoring.methods.full, labCol=scoring.methods.full) 
    dev.off()
    tiff("omim_in_ppi.tif", width=2000, height=2000, res=300)
    d<-read.table("../data/omim/2009_Aug_27/similarity_in_ppi.dat")
    heatmap(as.matrix(d), revC=T, col=val.cols, margins=c(9,9), Rowv=NA, Colv=NA, scale="none") 
    dev.off()
    tiff("omim_extended.tif", width=2000, height=2000, res=300)
    d<-read.table("../data/omim/2009_Aug_27/extended/similarity.dat")
    heatmap(as.matrix(d), revC=T, col=val.cols, margins=c(9,9), Rowv=NA, Colv=NA, scale="none") 
    dev.off()
}



###### MANUSCRIPT2 TESTS ######
manuscript2_tests<-function() {

    dir.summary<-"../data/summary_draft_before_revision/"

    # Significance of modules detected by methods vs nn
    d<-read.table(paste(dir.name, "../module/", "biana_no_tap-omim/", "module_summary.dat", sep=""), header=T)
    for(k in 1:length(scoring.methods)) {
	print(c("k", scoring.methods[k], "---"))
	x<-d[d$scoring=="nn",]$n_module
	y<-d[d$scoring==scoring.methods[k],]$n_module
	a<-wilcox.test(x, y, alternative="greater", paired=T)
	print(a$p.value)
	x<-d[d$scoring=="nn",]$ratio
	y<-d[d$scoring==scoring.methods[k],]$ratio
	a<-wilcox.test(x, y, alternative="less", paired=T)
	print(a$p.value)
	#for(i in 1:length(scoring.methods)) {
	#    print(c("i", scoring.methods[i], "---"))
	#    x<-d[d$scoring==scoring.methods[i],]$n_module
	#    y<-d[d$scoring==scoring.methods[k],]$n_module
	#    a<-wilcox.test(x, y, paired=T)
	#    print(a$p.value)
	#}
    }

    # Significance of AUC difference of methods between two groups based on number of seeds
    dir.name<-dir.summary
    s<-read.table(paste(dir.name, 'biana_no_tap-omim/seeds.dat', sep=""))
    d<-read.table(paste(dir.name, 'biana_no_tap-omim/auc_ppis.dat', sep=""))
    cutoff<-median(s$n_seed)
    x<-d[rownames(s)[s$n_seed<cutoff],]
    y<-d[rownames(s)[s$n_seed>=cutoff],]
    for(k in 1:length(scoring.methods)) {
	print(c("k", scoring.methods[k], "---"))
	a<-wilcox.test(x[,k], y[,k])
	print(a$p.value)
    }
    cor(d[order(rownames(d)),], s[order(rownames(s)),])

    # Linear regression using n_path & ns 
    y<-d[rownames(s),"ns"]
    x<-s[rownames(s),"n_path"]
    cor(x, y)
    fit<-lm(y~x)
    summary(fit)
    cor(fit$fitted.values, y)^2 # fit$coefficients[[2]]*x+fit$coefficients[[1]]
    1 - sum(fit$residuals^2) / sum((y-mean(y))^2)
    plot(x, y, ylim=c(0,100), xlim=c(0,100))
    abline(fit)

    #y<-d[rownames(s),]
    #x<-s[rownames(s),]
    #fit <- lm(y$ns ~ n_seed + n_linker + n_path, data=x)
    #fit.other <- lm(y$ns ~ n_path, data=x)
    #anova(fit, fit.other)

}




###### MANUSCRIPT FIGURES ######
manuscript_figures <- function() {

    #par(mar=c(5, 4, 4, 2) + 0.1) # to reset margins

    # Average AUC (%) over all diseases for each network
    #cairo_ps(paste(dir.name, "Figure 1.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    #par(family = "Arial") 
    tiff(paste(dir.name, "Figure 1.tif", sep=""), width=2000, height=2000, res=300)
    #par(xpd=T, mar=par()$mar+c(0,0,0,0))
    #coords <- seq(0.4,4.8,by=1.1)
    #coords <- seq(0,7.8,by=1.3)
    coords <- seq(0.05,7.75,by=1.1)
    i = 1
    for(ppi in c("goh", "entrez", "biana", "biana_no_tap", "biana_no_tap_relevance")) {
	d<-read.table(paste(dir.name, ppi, '_vs_all-wo_LI/auc_ppis.dat', sep=""))
	d<-d[,scoring.method.ids]
	if(i == 1) {
	    boxplot(d,col=cols[i],at=coords, boxwex=0.2, pars=list(xaxt="n"), xlab="Prediction method", ylab="AUC (%)", ylim=c(0,100), xlim=c(0.18,8.42))
	} 
	else {
	    boxplot(d,col=cols[i],at=coords, boxwex=0.2, pars=list(xaxt="n"), names=F, add=T)
	}
	if(i == 4) {
	    text.coords<-coords
	    #axis(1,tick=F,labels=scoring.methods,at=coords)
	}
	coords <- coords + 0.2
	i = i + 1
    }
    legend(0.5, 119, c("Goh", "Entrez", "PPI", "bPPI", "weighted bPPI"), fill=color5, bty="n", ncol=3, cex=0.9) #horiz=T, x.intersp=0.5)
    text(text.coords+0.15, par("usr")[3]-3, srt=40, adj=1, labels=scoring.methods.full, xpd=T, cex=0.9)
    dev.off()

    # Ratio of correctly predicted proteins (%) among top 10% predictions over all diseases for each network
    #cairo_ps(paste(dir.name, "Supplementary Figure 2.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    #par(family = "Arial") 
    tiff(paste(dir.name, "Supplementary Figure 2.tif", sep=""), width=2000, height=2000, res=300)
    par(xpd=T, mar=par()$mar+c(0,0,0,0))
    #coords <- seq(0.4,4.8,by=1.1)
    #coords <- seq(0,7.8,by=1.3)
    coords <- seq(0.05,7.75,by=1.1)
    i = 1
    for(ppi in c("goh", "entrez", "biana", "biana_no_tap", "biana_no_tap_relevance")) {
	d<-read.table(paste(dir.name, ppi, '_vs_all-wo_LI/cov_ppis.dat', sep=""))
	d<-d[,scoring.method.ids]
	if(i == 1) {
	    boxplot(d,col=cols[i],at=coords, boxwex=0.2, pars=list(xaxt="n"), xlab="Prediction method", ylab="Sensitivity at top 1% (%)", ylim=c(0,100), xlim=c(0.18,8.42))
	} 
	else {
	    boxplot(d,col=cols[i],at=coords, boxwex=0.2, pars=list(xaxt="n"), names=F, add=T)
	}
	if(i == 4) {
	    text.coords<-coords
	    #axis(1,tick=F,labels=scoring.methods,at=coords)
	}
	coords <- coords + 0.2
	i = i + 1
    }
    legend(0.5, 119, c("Goh", "Entrez", "PPI", "bPPI", "weighted bPPI"), fill=color5, bty="n", ncol=3, cex=0.9) 
    text(text.coords+0.15, par("usr")[3]-4, srt=40, adj=1, labels=scoring.methods.full, xpd=T, cex=0.9)
    dev.off()

    # Average AUC (%) for two groups of OMIM disorders based on number of seeds on bPPI network
    dir.name<-"../data/summary/"
    #postscript(paste(dir.name, "Figure 4a.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
    #cairo_ps(paste(dir.name, "Figure 2.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    #par(family = "Arial") 
    tiff(paste(dir.name, "Figure 2.tif", sep=""), width=2000, height=2000, res=300)
    #par(xpd=T, mar=par()$mar+c(0,0,0,0))

    #prefix<-"biana_no_tap_vs_omim-wo_LI"
    prefix<-"biana_no_tap_vs_all-wo_LI"
    s<-read.table(paste(dir.name, prefix, '/seeds.dat', sep=""))
    d<-read.table(paste(dir.name, prefix, '/auc_ppis.dat', sep=""))

    d<-d[,scoring.method.ids]

    cutoff<-median(s$n_seed)
    x<-d[rownames(s[s$n_seed<cutoff,]),]
    y<-d[rownames(s[s$n_seed>=cutoff,]),]

    #coords <- seq(0.6, 5.0, by=1.1)
    #coords <- seq(0.6, 7.2, by=1.1)
    coords <- seq(0.6, 8.3, by=1.1)
    boxplot(x,col="lightgray", at=coords, boxwex=0.4, pars=list(xaxt="n"), xlab="Prediction method", ylab="AUC(%)", ylim=c(0,100), xlim=c(0.6,8.7))
    #axis(1,tick=F,labels=scoring.methods,at=coords+0.2)
    text(coords+0.55, par("usr")[3]-3, srt=40, adj=1, labels=scoring.methods.full, xpd=T, cex=0.9)
    coords <- coords + 0.4
    boxplot(y,col="darkgray", at=coords, boxwex=0.4, pars=list(xaxt="n"), names=F, add=T)
    legend(1.5, 5, c(paste("# of seeds < ", cutoff, sep=""), paste("# of seeds >= ", cutoff, sep="")), fill=c("lightgray","darkgray"), bty="n", horiz=T, cex=0.9)
    dev.off()

    #cairo_ps(paste(dir.name, "Supplementary Figure 3.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    tiff(paste(dir.name, "Supplementary Figure 3.tif", sep=""), width=2000, height=2000, res=300)
    prefix<-"all_new_vs_all_new"
    # Significance of AUC & COV difference between methods on bPPI
    n<-length(scoring.method.ids)
    pvals<-matrix(1, n, n)
    print("AUC & COV differences between methods on bPPI")
    d<-read.table(paste(dir.name, prefix, "/auc_ppis.dat", sep=""))
    e<-read.table(paste(dir.name, prefix, "/cov_ppis.dat", sep=""))
    d<-d[order(row.names(d)),scoring.method.ids]
    e<-e[order(row.names(e)),scoring.method.ids]
    for(i in 1:length(scoring.method.ids)) {
	for(j in 1:length(scoring.method.ids)) {
		if(i<j) {
		    a<-wilcox.test(d[,i], d[,j], alternative="greater", paired=T)
		    b<-wilcox.test(e[,i], e[,j], alternative="greater", paired=T)
		    if(a$p.value <= 0.05 || b$p.value <= 0.05) {
			print(c(scoring.method.ids[i], scoring.method.ids[j], a$p.value, b$p.value))
		    }
		}
		if(i==j) {
		    val<-1
		} else {
		    val<-wilcox.test(d[,i], d[,j], alternative="greater", paired=T)$p.value
		}
		pvals[i,j]<-val 
	}
    }
    library(gplots)
    library(RColorBrewer)
    pvals<-ifelse((pvals<=0.05) == TRUE, 1, 0)
    val.cols <- brewer.pal(9,"Blues") 
    #heatmap.2(as.matrix(pvals), revC=F, trace="none", dendrogram="none", col=val.cols, density.info="none", margins=c(9,9), Rowv=NA, Colv=NA, labRow=scoring.methods.full, labCol=scoring.methods.full, key=F) #keysize=0.1) #, key=F)
    heatmap(as.matrix(pvals), revC=T, col=val.cols, margins=c(9,9), Rowv=NA, Colv=NA, labRow=scoring.methods.full, labCol=scoring.methods.full, scale="none") 
    dev.off()

}




###### MANUSCRIPT TESTS ######
manuscript_tests <- function() {
    dir.name<-"../data/summary/"

    #prefix<-"all_new_vs_omim-w_LI"
    #prefix<-"all_new_vs_omim-wo_LI"
    #prefix<-"all_new_vs_goh-wo_LI"
    #prefix<-"all_new_vs_chen-wo_LI"
    prefix<-"biana_no_tap_vs_all-wo_LI"

    # Significance of AUC & COV difference between methods on bPPI
    n<-length(scoring.method.ids)
    print("AUC & COV differences between methods on bPPI")
    d<-read.table(paste(dir.name, prefix, "/auc_ppis.dat", sep=""))
    e<-read.table(paste(dir.name, prefix, "/cov_ppis.dat", sep=""))
    d<-d[order(row.names(d)),scoring.method.ids]
    e<-e[order(row.names(e)),scoring.method.ids]
    for(i in 1:length(scoring.method.ids)) {
	for(j in 1:length(scoring.method.ids)) {
		if(i<j) {
		    a<-wilcox.test(d[,i], d[,j], alternative="greater", paired=T)
		    b<-wilcox.test(e[,i], e[,j], alternative="greater", paired=T)
		    if(a$p.value <= 0.05 || b$p.value <= 0.05) {
			print(c(scoring.method.ids[i], scoring.method.ids[j], a$p.value, b$p.value))
		    }
		}
	}
    }

    prefices<-c("biana_vs_all-wo_LI", "biana_no_tap_vs_all-wo_LI", "biana_no_tap_relevance_vs_all-wo_LI", "goh_vs_all-wo_LI", "entrez_vs_all-wo_LI")
    #prefices<-c("biana_vs_omim-wo_LI", "biana_no_tap_vs_omim-wo_LI", "biana_no_tap_relevance_vs_omim-wo_LI", "goh_vs_omim-wo_LI", "entrez_vs_omim-wo_LI")
    #prefices<-c("all_new_vs_omim-w_LI", "all_new_vs_omim-wo_LI", "all_new_vs_goh-wo_LI", "all_new_vs_chen-wo_LI")
    for(prefix in prefices) {
	# Correlation between seed connectivity measures and AUC on bPPI
	print(c("Correlations on ", prefix))
	d<-read.table(paste(dir.name, prefix, "/auc_ppis.dat", sep=""))
	d<-d[order(rownames(d)),]
	s<-read.table(paste(dir.name, prefix, "/seeds.dat", sep=""))
	s<-s[order(rownames(s)),]
	print(cor(d,s)) 
    }

    prefix<-"biana_no_tap_vs_all-wo_LI"
    # Significance of AUC change between bPPI n_seed < cutoff and >= cutoff
    d<-read.table(paste(dir.name, prefix, "/auc_ppis.dat", sep=""))
    s<-read.table(paste(dir.name, prefix, "/seeds.dat", sep=""))

    cutoff<-median(s$n_seed)
    x<-d[rownames(d) %in% rownames(s)[s$n_seed<cutoff],scoring.method.ids]
    y<-d[rownames(d) %in% rownames(s)[s$n_seed>=cutoff],scoring.method.ids]
    print(c(dim(x), dim(y)))

    print("AUC difference between two groups based on n_seed using bPPI")
    #wilcox.test(x$ns, x$nd, alternative="less", paired=T)
    #wilcox.test(x$nd, y$nd, alternative="greater", paired=F)
    #wilcox.test(x$np, y$np, alternative="greater", paired=F)
    for(i in 1:length(scoring.method.ids)) {
	a<-wilcox.test(x[,i], y[,i], alternative="greater", paired=F)
	print(c(scoring.method.ids[i], a$p.value))
    }

    prefix<-"biana_no_tap_vs_all-wo_LI"
    prefix2<-"biana_no_tap_relevance_vs_all-wo_LI"
    # Significance between bPPI and weighted bPPI
    print("AUC difference between bPPI and weighted bPPI")
    d<-read.table(paste(dir.name, prefix, "/auc_ppis.dat", sep=""))
    e<-read.table(paste(dir.name, prefix2, "/auc_ppis.dat", sep=""))
    d<-d[order(row.names(d)),scoring.method.ids]
    e<-e[order(row.names(e)),scoring.method.ids]
    for(i in 1:length(scoring.method.ids)) {
	a<-wilcox.test(e[,i], d[,i], alternative="greater", paired=T)
	print(c(scoring.method.ids[i], a$p.value))
    }

    prefix<-"all_new_vs_omim-w_LI"
    prefix2<-"all_new_vs_omim-wo_LI"
    # Significance between bPPI and weighted bPPI
    print("AUC difference between OMIM-LI and OMIM")
    d<-read.table(paste(dir.name, prefix, "/auc_ppis.dat", sep=""))
    e<-read.table(paste(dir.name, prefix2, "/auc_ppis.dat", sep=""))
    d<-d[order(row.names(d)),scoring.method.ids]
    e<-e[order(row.names(e)),scoring.method.ids]
    for(i in 1:length(scoring.method.ids)) {
	a<-wilcox.test(e[,i], d[,i], alternative="greater", paired=T)
	print(c(scoring.method.ids[i], a$p.value))
    }


    prefix<-"goh_vs_all-wo_LI"
    prefix2<-"biana_no_tap_vs_all-wo_LI"
    # Significance between Goh and bPPI 
    print("AUC difference between Goh and bPPI")
    d<-read.table(paste(dir.name, prefix, "/auc_ppis.dat", sep=""))
    e<-read.table(paste(dir.name, prefix2, "/auc_ppis.dat", sep=""))
    d<-d[order(row.names(d)),scoring.method.ids]
    e<-e[order(row.names(e)),scoring.method.ids]
    for(i in 1:length(scoring.method.ids)) {
	a<-wilcox.test(e[,i], d[,i], alternative="greater", paired=T)
	print(c(scoring.method.ids[i], a$p.value))
    }


}





###### YEAST FIGURES ######
yeast_figures<-function() {


    # AUC over networks
    dir.name<-"../data/summary/"
    cairo_ps(paste(dir.name, "Figure P6a.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    par(family = "Arial")
    par(xpd=T, mar=par()$mar+c(0,0,2,0))

    dir.name<-"/sbi/users/emre/data/netzcore/summary_rob/"
    coords <- seq(0.5,4.9,by=1.1)
    i = 1
    #for(ppi in c("rob_biogrid_no_tap", "rob_biogrid", "rob_biogrid_only_genetic", "rob_biogrid_with_genetic", "rob_yeastnet2")) {
    for(ppi in c("rob_biogrid", "rob_biogrid_only_genetic", "rob_biogrid_with_genetic", "rob_yeastnet2")) {  
	d<-read.table(paste(dir.name, ppi, '/auc_ppis.dat', sep=""))
        if(i == 1) {
            boxplot(d,col=cols[i],at=coords, boxwex=0.2, pars=list(xaxt="n"), xlab="Prediction method", ylab="AUC(%)", ylim=c(0,100)) #, title="5-fold cross-validation AUCs over 17 gene sets for the prediction methods")
        }
        else {
            boxplot(d,col=cols[i],at=coords, boxwex=0.2, pars=list(xaxt="n"), names=F, add=T)
        }
        if(i == 3) {
            axis(1, tick=F, labels=scoring.methods, at=coords)
        }
        coords <- coords + 0.2
        i = i + 1
    }
    #legend(0.5, 119, c("BG-no_tap", "BG", "BG_genetic", "BG_w/genetic", "YeastNet2"), fill=color5, bty="n", ncol=3) #horiz=T, x.intersp=0.5)
    legend(0.5, 119, c("PPI", "GI", "PPI+GI", "FI"), fill=color5, bty="n", horiz=T) #, x.intersp=0.5) #ncol=3) 
    dev.off()

    # Top 5% seed coverage over networks
    dir.name<-"../data/summary/"
    cairo_ps(paste(dir.name, "Figure P6b.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    par(family = "Arial")
    par(xpd=T, mar=par()$mar+c(0,0,2,0))

    dir.name<-"/sbi/users/emre/data/netzcore/summary_rob/"
    coords <- seq(0.5,4.9,by=1.1)
    i = 1
    #for(ppi in c("rob_biogrid_no_tap", "rob_biogrid", "rob_biogrid_only_genetic", "rob_biogrid_with_genetic", "rob_yeastnet2")) {
    for(ppi in c("rob_biogrid", "rob_biogrid_only_genetic", "rob_biogrid_with_genetic", "rob_yeastnet2")) {  
	d<-read.table(paste(dir.name, ppi, '/cov_ppis.dat', sep=""))
        if(i == 1) {
            boxplot(d,col=cols[i],at=coords, boxwex=0.2, pars=list(xaxt="n"), xlab="Prediction method", ylab="Seeds covered at top 5% (%)", ylim=c(0,100)) #, title="5-fold cross-validation AUCs over 17 gene sets for the prediction methods")
        }
        else {
            boxplot(d,col=cols[i],at=coords, boxwex=0.2, pars=list(xaxt="n"), names=F, add=T)
        }
        if(i == 3) {
            axis(1, tick=F, labels=scoring.methods, at=coords)
        }
        coords <- coords + 0.2
        i = i + 1
    }
    #legend(0.5, 119, c("BG-no_tap", "BG", "BG_genetic", "BG_w/genetic", "YeastNet2"), fill=color5, bty="n", ncol=3) #horiz=T, x.intersp=0.5)
    legend(0.5, 119, c("PPI", "GI", "PPI+GI", "FI"), fill=color5, bty="n", horiz=T) #, x.intersp=0.5) #ncol=3) 
    dev.off()
}




###### AD CASE STUDY FIGURES ######
case_study_figures <- function() {
    library(VennDiagram)
    A=251
    A.B=14
    A.B.C=1
    A.C=11
    B.C=3
    B=32
    C=80
    x1 = list(
		'Predicted' = c(1:A, (A+1):(A+A.B), (A+A.B+1):(A+A.B+A.B.C), (A+A.B+A.B.C+1):(A+A.B+A.B.C+A.C)),
		'AD related' = c((A+A.B+A.B.C+A.C+1):(A+A.B+A.B.C+A.C+B), (A+1):(A+A.B), (A+A.B+1):(A+A.B+A.B.C), (A+A.B+A.B.C+A.C+B+1):(A+A.B+A.B.C+A.C+B+B.C)),
		'Aging related' = c((A+A.B+A.B.C+A.C+B+B.C+1):(A+A.B+A.B.C+A.C+B+B.C+C), (A+A.B+1):(A+A.B+A.B.C), (A+A.B+A.B.C+1):(A+A.B+A.B.C+A.C), (A+A.B+A.B.C+A.C+B+1):(A+A.B+A.B.C+A.C+B+B.C))
		    )
    A=458
    A.B=23
    A.B.C=1
    A.C=28
    B.C=2
    B=24
    C=64
    x2 = list(
		'Predicted' = c(1:A, (A+1):(A+A.B), (A+A.B+1):(A+A.B+A.B.C), (A+A.B+A.B.C+1):(A+A.B+A.B.C+A.C)),
		'AD related' = c((A+A.B+A.B.C+A.C+1):(A+A.B+A.B.C+A.C+B), (A+1):(A+A.B), (A+A.B+1):(A+A.B+A.B.C), (A+A.B+A.B.C+A.C+B+1):(A+A.B+A.B.C+A.C+B+B.C)),
		'Aging related' = c((A+A.B+A.B.C+A.C+B+B.C+1):(A+A.B+A.B.C+A.C+B+B.C+C), (A+A.B+1):(A+A.B+A.B.C), (A+A.B+A.B.C+1):(A+A.B+A.B.C+A.C), (A+A.B+A.B.C+A.C+B+1):(A+A.B+A.B.C+A.C+B+B.C))
		    )
    venn.diagram( x = x1, filename = paste(dir.name, "Figure 6c.tif", sep=""), col = "transparent", fill = c("red", "blue", "green"), alpha = 0.5, label.col = c("darkred", "white", "darkblue", "white", "white", "white", "darkgreen"), cex = 2.5, fontfamily = "arial", fontface = "bold", cat.default.pos = "text", cat.col = c("darkred", "darkblue", "darkgreen"), cat.cex = 2.5, cat.fontfamily = "arial", cat.dist = c(0.08, 0.08, -0.08), cat.pos = 0, width=6, height=6, units="in", resolution=300) # helvetica
    venn.diagram( x = x2, filename = paste(dir.name, "Figure 6d.tif", sep=""), col = "transparent", fill = c("red", "blue", "green"), alpha = 0.5, label.col = c("darkred", "white", "darkblue", "white", "white", "white", "darkgreen"), cex = 2.5, fontfamily = "arial", fontface = "bold", cat.default.pos = "text", cat.col = c("darkred", "darkblue", "darkgreen"), cat.cex = 2.5, cat.fontfamily = "arial", cat.dist = c(0.08, 0.08, -0.08), cat.pos = 0, width=6, height=6, units="in", resolution=300)
}





###### NAVLAKHA ######
navlakha_figures <- function() {
    library(RColorBrewer)
    dir.name<-"../data/compare/navlakha/"
    postscript(paste(dir.name, "results.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
    d <- read.table(paste(dir.name, "results.dat", sep=""))
    par(mar=c(6,6,0.2,0))
    #scoring.methods <- c("ns", "nz", "nd", "ff", "nr", "mcl")
    #scoring.methods <- c("mcl")
    #scoring.methods <- c("ns", 'np')
    scoring.methods <- c("ns", "nz", "nd", "ff", "nr", "rw", "np", "nc", "nc3")
    #colors <- c(2,6,4,5,3,8,1,7,9)
    colors<-c(brewer.pal(9, "Set1"), brewer.pal(8, "Accent"))
    chars <- c(0:2,5:6,9:14)
    for(i in 1:length(scoring.methods)) {
	if(scoring.methods[i] %in% c("mcl", "nc3")) {
	    e<-subset(d, substr(rownames(d),1,3)==scoring.methods[i])
	}
	else {
	    e<-subset(d, substr(rownames(d),1,2)==scoring.methods[i])
	}
	e<-e[order(e$sens),]
	n<-dim(e)[1]
	x<-c(e$sens[1])
	y<-c(e$ppv[1])
	for(j in seq(i+1, n, by=5)) {
	    x<-c(x, e$sens[j])
	    y<-c(y, e$ppv[j])
	}
	x<-c(x, e$sens[n])
	y<-c(y, e$ppv[n])
	if(i == 1) {
	    plot(e$sens, e$ppv, col=colors[i], type="l", xlim=c(0,1), ylim=c(0,1))
	    points(x, y, col=colors[i], pch=chars[i], cex=2)
	} else {
	    lines(e$sens, e$ppv, col=colors[i], lwd=1)
	    points(x, y, col=colors[i], pch=chars[i], cex=2)
	}   
	#lines(e$sens, e$ppv, col=colors[i])
	text(e$sens[1], e$ppv[1]+0.05, col=colors[i], toupper(scoring.methods[i]))
    }
    legend("topright", scoring.methods, col=colors, pch=chars)
    dev.off()
}




###### ANEURYSM FIGURES ######
aneurysm_figures <- function() {

    # Average AUC (%) over all diseases for each network
    #postscript(paste(dir.name, "Figure S1a.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
    cairo_ps(paste(dir.name, "Figure.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    par(family = "Arial") 
    par(xpd=T, mar=par()$mar+c(-1,0,1,0))
    d<-read.table(paste(dir.name, 'bppi_coexp/auc_phenotypes.dat', sep=""))
    e<-as.matrix(d[c(3,1,6,4,5,2),])
    barplot(e,col=cols[1:6],beside=T,names.arg=scoring.methods[1:5],ylab="AUC (%)",xlab="Prediction method")
    legend(0.5, 90, c("bPPI", "bPPI_coexp", "bPPI_Dcoexp", "bPPI_loc", "bPPI_coexp_loc", "bPPI_Dcoexp_loc"), fill=cols[1:6], bty="n", ncol=3, cex=0.9) 
    dev.off()
}





###### OBSOLETE FIGURES ######
old_figures <- function() { 

    # Average AUC (%) over all diseases and networks
    #postscript(paste(dir.name, "Figure 2a.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special", family="arial", bg="white")
    #svg(paste(dir.name, "Figure 2a.svg", sep=""), width = 6, height = 6, onefile = TRUE)
    cairo_ps(paste(dir.name, "Figure S1a.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    par(family = "Arial") 
    d<-read.table(paste(dir.name, 'all_new_vs_all_new-wo_LI/auc_phenotypes.dat', sep=""))
    d<-d[,scoring.method.ids]
    boxplot(d, xlab="Prediction method", ylab="AUC(%)", names=scoring.methods, col=8, notch=FALSE, varwidth=FALSE, ylim=c(0,80))
    dev.off()

    # Ratio of correctly predicted proteins (%) among top 10% predictions over all diseases and networks
    #postscript(paste(dir.name, "Figure 2b.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
    cairo_ps(paste(dir.name, "Figure S1b.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    par(family = "Arial") 
    d<-read.table(paste(dir.name, 'all_new_vs_all_new-wo_LI/cov_phenotypes.dat', sep=""))
    d<-d[,scoring.method.ids]
    boxplot(d, xlab="Prediction method", ylab="Seeds covered at top 1% (%)", names=scoring.methods, col=8, notch=FALSE, varwidth=FALSE, ylim=c(0,80))
    dev.off()


    postscript(paste(dir.name, "3a.eps", sep=""), width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
    d<-read.table(paste(dir.name, 'biana_no_tap-omim/auc_ppis.dat', sep=""))
    e<-data.frame(d[,1],d[,2],d[,3])
    barplot(as.matrix(t(e)), beside=TRUE, horiz=F, col = rainbow(3), ylab="AUC (%)", legend.text=c("NetScore", "NetZcore", "NetShort"))
    text(seq(1.1, 0.1+length(d[,1])*4.1, 4.1), par("usr")[3] - 2, srt=45, adj=1, labels=names(as.data.frame(t(d))), xpd=T, cex=0.8)
    dev.off()

    postscript(paste(dir.name, "S5_a.eps", sep=""), horizontal=F) #width = 480, height = 480, quality = 100)
    d<-read.table("module_summary_top1_mcl_n3.dat", header=T)
    scoring.methods<-c("no", "nn", "n3") 
    n<-length(scoring.methods)
    x<-1:length(levels(d$phenotype))
    par(mar=c(5,4,4,5))
    e<-c()
    for(i in 1:n) {
	print(scoring.methods[i])
	if(i == 1) {
	    plot(x, d[d$scoring==scoring.methods[i],]$n_seed_in_modules/d[d$scoring==scoring.methods[i],]$n_seed, type="l", col=i+1, ylim=c(0,1), xlab="Phenotypes", ylab="Ratio of covered seeds within modules", xaxt="n")
	} 
	else {
	    lines(x, d[d$scoring==scoring.methods[i],]$n_seed_in_modules/d[d$scoring==scoring.methods[i],]$n_seed, lty=1, col=i+1)
	}
	if(i != 1) {
	    e<-rbind(e, d[d$scoring==scoring.methods[i],]$n_module)
	}
    }
    legend("topright", c("GWAS", "PWAS-1", "PWAS-2"), lty=rep(1, n), col=(1:n)+1)
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
    scoring.methods<-c("no", "nn", "nr") 
    #scoring.methods<-c("no", "nn", "nr", "ff", "nz", "nd", "ns") #levels(d$scoring)
    n<-length(scoring.methods)
    x<-1:length(levels(d$phenotype))
    par(mar=c(5,4,4,5))
    e<-c()
    for(i in 1:n) {
	print(scoring.methods[i])
	if(i == 1) {
	    plot(x, d[d$scoring==scoring.methods[i],]$n_seed, type="l", col=1, ylim=c(0,100), xlab="Phenotypes", ylab="Number of seeds", xaxt="n")
	} 
	e<-rbind(e, d[d$scoring==scoring.methods[i],]$n_seed_in_modules)
    }
    legend("topright", c("# of seeds", "GWAS coverage", "PWAS-1 coverage", "PWAS-2 coverage"), lty=rep(1, n+1), col=(0:n)+1)
    text(1:length(levels(d$phenotype)), par("usr")[3]-0.02, srt=30, adj=1, labels=sapply(levels(d$phenotype), lambda), xpd=T, cex=0.7)
    par(new=TRUE)
    e<-as.matrix(e)
    print(dim(e))
    barplot(e, beside=TRUE, horiz=F, col=(1:n)+1, width=0.1, space=c(0,1.3), ylim=c(0,100), axes=F)
    axis(4, 1:10, labels=1:10)
    text(par("usr")[2]+1, 70, srt=90, adj=1, labels="Number of covered seeds within modules using CC", xpd=T, cex=1)
    dev.off()

    dir.name<-"../data/summary/"
    d<-read.table(paste(dir.name, "module/biana_no_tap_relevance-omim/module_summary.dat", sep=""), header=T)
    #d2<-read.table(paste(dir.name, "module/module_summary_top5_mcl_n5union.dat", sep=""), header=T)
    scoring.methods<-c("no", "nn", "nr", "ns") 
    method.names<-c("Neighborhood", "Func.Flow", "NetScore")
    #scoring.methods<-c("no", "nn", "nr", "ff", "nd", "nz", "ns") 
    #method.names<-c("GWAS", "PWAS-1", "PWAS-2")
    #method.names<-c("PWAS-Neighborhood", "PWAS-ToppGene", "PWAS-Func.Flow", "PWAS-NetShort", "PWAS-NetZcore", "PWAS-NetScore")
    n<-length(scoring.methods)
    e<-c()
    f<-c()
    g<-c()
    #e.genes<-c()
    e.modules<-c()
    #e.seeds<-c()
    for(i in 1:n) {
	print(scoring.methods[i])
	#f<-rbind(f, d[d$scoring==scoring.methods[i],]$n_seed_in_modules/d[d$scoring==scoring.methods[i],]$n_seed) 
	if(scoring.methods[i] != "no") {
	    e<-rbind(e, d[d$scoring==scoring.methods[i],]$n_module)
	    f<-rbind(f, 100*d[d$scoring==scoring.methods[i],]$ratio) 
	    #g<-rbind(g, 100*d[d$scoring==scoring.methods[i],]$n_seed_in_modules/d[d$scoring==scoring.methods[i],]$n_all_in_modules)
	    #g<-rbind(g, 100*d[d$scoring==scoring.methods[i],]$n_seed_go_in_modules/d[d$scoring==scoring.methods[i],]$n_go_in_modules)
	    g<-rbind(g, 100*d[d$scoring==scoring.methods[i],]$n_seed_go_in_modules/d[d$scoring==scoring.methods[i],]$n_seed_go)
	    #e.genes<-rbind(e.genes, d2[d2$scoring==scoring.methods[i],]$n_all_in_modules)
	}
	if(scoring.methods[i] == "no" | scoring.methods[i] == "ns") {
	    e.modules<-rbind(e.modules, d[d$scoring==scoring.methods[i],]$n_module)
	    #e.seeds<-rbind(e.seeds, d2[d2$scoring==scoring.methods[i],]$n_seed_in_modules)
	}
    }
    postscript(paste(dir.name, "S5a.eps", sep=""), horizontal=F)
    barplot(e, beside=TRUE, horiz=F, col=(2:n), ylab="Number of modules", ylim=c(0,16), xlab="Phenotypes", legend.text=method.names) #, width=0.1, space=c(0,1.3), ylim=c(0,30), axes=F)
    x1<-par("usr")[1]
    x2<-par("usr")[2]
    text(seq(x1+4.2, x2, by=(x2-x1-4)/(length(levels(d$phenotype)))), par("usr")[3]-0.1, srt=30, adj=1, labels=sapply(levels(d$phenotype), lambda), xpd=T, cex=0.7)
    dev.off()

   postscript(paste(dir.name, "S5f.eps", sep=""), horizontal=F) #width = 480, height = 480, quality = 100)
    d<-read.table(paste(dir.name, "module_summary_top5_mcl_n5union.dat", sep=""), header=T)
    scoring.methods<-c("no", "nn", "nd") 
    n<-length(scoring.methods)
    x<-1:length(levels(d$phenotype))
    par(mar=c(5,4,4,5))
    e<-c()
    for(i in 1:n) {
	print(scoring.methods[i])
	if(i != 1) {
	    e<-rbind(e, d[d$scoring==scoring.methods[i],]$ratio)
	}
    }
    a<-barplot(e, beside=TRUE, horiz=F, col=(2:n)+1, ylim=c(0,1), legend.text=c("PWAS-1", "PWAS-2"), xlab="Phenotypes", ylab="Ratio of non-seed connections within modules (%)") #width=0.1, space=c(0,1.3), 
    #legend("topright", c("PWAS-1", "PWAS-2"), lty=rep(1, n), col=(2:n)+1)
    text(1:length(levels(d$phenotype)), par("usr")[3]-0.02, srt=30, adj=1, labels=sapply(levels(d$phenotype), lambda), xpd=T, cex=0.7)
    for(i in 1:n) {
	print(scoring.methods[i])
	if(i != 1) {
	    lines(colMeans(a), d[d$scoring==scoring.methods[i],]$n_seed_in_modules/d[d$scoring==scoring.methods[i],]$n_all_in_modules, lty=4, col=i+1)
	}
    }
    #axis(4, 1:15, labels=1:15)
    #text(par("usr")[2]+1, 10, srt=90, adj=1, labels="Number of modules using MCL", xpd=T, cex=1)
    dev.off()



    # Average AUC (%) on different perturbations for two groups of disorders based on number of seeds on bPPI network
    dir.name<-"../data/summary/"
    s<-read.table(paste(dir.name, 'biana_no_tap-omim/seeds.dat', sep=""))
    cutoff<-median(s$n_seed)
    selected.names<-rownames(s)[s$n_seed<=cutoff]
    not.selected.names<-rownames(s)[s$n_seed>cutoff]

    percentages<-c(0,10,30,50,70) 
    method<-"ns"

    # Pruned interactions
    dir.name<-"../data/summary/"
    cairo_ps(paste(dir.name, "Figure PS3b.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    par(family = "Arial") 

    label<-"Percentage of pruned interactions"
    container.x<-data.frame()
    container.y<-data.frame()
    for(p in percentages) {
	if(p == 0) {
	    d<-read.table(paste(dir.name, 'biana_no_tap-omim/auc_ppis.dat', sep=""))
	    dir.name<-"../data/summary_runs_on_random/"
	    x<-d[selected.names, method]
	    y<-d[not.selected.names, method]
	    container.x<-x
	    container.y<-y
	}
	else {
	    d<-read.table(paste(dir.name, 'biana_no_tap_pruned_p', p, '-omim/auc_ppis.dat', sep=""))
	    if(!all(dim(d) == c(23, 5))) {
		stop("Object dimension is different from expected!")
	    }
	    x<-d[selected.names, method]
	    y<-d[not.selected.names, method]
	    container.x<-cbind(container.x,x)
	    container.y<-cbind(container.y,y)
	}
    }
    colnames(container.x)<-percentages
    colnames(container.y)<-percentages

    coords <- seq(0.6, 5.0, by=1.1)
    boxplot(container.x, col=3, at=coords, boxwex=0.4, pars=list(xaxt="n"), xlab=label, ylab="AUC(%)", ylim=c(0,100))
    axis(1,tick=F, labels=percentages, at=coords+0.2)
    coords <- coords + 0.4
    boxplot(container.y, col=4, at=coords, boxwex=0.4, pars=list(xaxt="n"), names=F, add=T)
    legend(1.0, 15, c("Group with shorter seed-connecting paths", "Group with longer seed-connecting paths"), fill=color2, bty="n")
    dev.off()

    # Perturbed seeds
    dir.name<-"../data/summary/"
    cairo_ps(paste(dir.name, "Figure P3a.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    par(family = "Arial") 

    #lambda<-function(x) { paste("perturbed_", x, sep="") }
    #selected.names<-apply(as.matrix(selected.names), 1, lambda)

    label<-"Percentage of perturbed seeds"
    container.x<-data.frame()
    container.y<-data.frame()
    for(p in percentages) {
	if(p == 0) {
	    d<-read.table(paste(dir.name, 'biana_no_tap-omim/auc_ppis.dat', sep=""))
	    dir.name<-"../data/summary_runs_on_random/"
	    x<-d[selected.names, method]
	    y<-d[not.selected.names, method]
	    container.x<-x
	    container.y<-y
	}
	else {
	    d<-read.table(paste(dir.name, 'biana_no_tap-omim_perturbed_p', p, '/auc_ppis.dat', sep=""))
	    if(!all(dim(d) == c(2300, 5))) {
		stop("Object dimension is different from expected!")
	    }
	    f<-data.frame()
	    for(current in rownames(s)) {
		e<-data.frame()
		for(i in 1:100) {
		    e<-rbind(e,d[paste("perturbed_",current,"_p",p,"_",i,sep=""),])
		} 
		f<-rbind(f,mean(e))
	    }
	    rownames(f)<-rownames(s)
	    colnames(f)<-colnames(d)
	    #selected.names.inner<-c()
	    #for(row in rownames(d)) {
	    #	words<-unlist(strsplit(row, "_"))
	    #	i<-grep(paste("p", p, sep=""), words)
	    #	current<-paste(words[2:(i-1)], collapse="_")
	    #	if(current %in% selected.names) {
	    #	    selected.names.inner<-c(selected.names.inner, row)
	    #	}
	    #}
	    x<-f[selected.names, method]
	    y<-f[not.selected.names, method]
	    container.x<-cbind(container.x,x)
	    container.y<-cbind(container.y,y)
	}
    }
    colnames(container.x)<-percentages
    colnames(container.y)<-percentages

    coords <- seq(0.6, 5.0, by=1.1)
    boxplot(container.x, col=3, at=coords, boxwex=0.4, pars=list(xaxt="n"), xlab=label, ylab="AUC(%)", ylim=c(0,100))
    axis(1,tick=F, labels=percentages, at=coords+0.2)
    coords <- coords + 0.4
    boxplot(container.y, col=4, at=coords, boxwex=0.4, pars=list(xaxt="n"), names=F, add=T)
    legend(1.0, 15, c("Group with shorter seed-connecting paths", "Group with longer seed-connecting paths"), fill=color2, bty="n")
    dev.off()

    # Permuted interactions
    method<-"ns"
    dir.name<-"../data/summary/"
    percentages<-seq(0,80,by=10)
    cols<-1:23

    cairo_ps(paste(dir.name, "Figure P4a_permuted.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    #cairo_ps(paste(dir.name, "Figure P4a_pruned.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    par(family = "Arial") 

    label<-"Percentage of permuted interactions"
    #label<-"Percentage of pruned interactions"
    container<-data.frame()
    phenotypes<-c()
    for(p in percentages) {
	if(p == 0) {
	    d<-read.table(paste(dir.name, 'biana_no_tap-omim/auc_ppis.dat', sep=""))
	    dir.name<-"../data/summary_runs_on_random/"
	    container<-d[order(rownames(d)),method]
	    phenotypes<-rownames(d)[order(rownames(d))]
	}
	else {
	    d<-read.table(paste(dir.name, 'biana_no_tap_permuted_p', p, '-omim/auc_ppis.dat', sep=""))
	    if(!all(dim(d) == c(23, 5))) {
		stop("Object dimension is different from expected!")
	    }
	    container<-cbind(container, d[order(rownames(d)),method])
	}
    }
    colnames(container)<-percentages
    rownames(container)<-phenotypes

    #a<-container[,1]-container[,2]
    #cutoff<-median(abs(a))
    #a<-a[abs(a)>cutoff]
    #a<-a[a<=(-cutoff)] # up
    #a<-a[a>=(cutoff)] # down
    #container<-container[names(a),]

    for(i in 1:dim(container)[1]) {
	if(i==1) {
	    plot(percentages, container[i,], type='l', col=cols[i], xaxt="n", xlab=label, ylab="AUC(%)", ylim=c(0,100))
	} else {
	    lines(percentages, container[i,], col=cols[i])
	}
    }
    axis(1,tick=F, labels=percentages, at=percentages)
    #legend("bottomleft", sapply(rownames(container), lambda), lty=rep(1, dim(container)[1]), col=cols, bty="n", ncol=2)
    dev.off()

    dir.name<-"../data/summary/"
    percentages<-seq(0,80,by=10) # c(0, 10, 30, 50, 70)

    # Average AUC (%) at different levels of interaction permutation over OMIM disorders on bPPI network
    cairo_ps(paste(dir.name, "Figure P2a.eps", sep=""), width = 6, height = 6, onefile = TRUE) # ps: horizontal = FALSE, onefile = FALSE, paper = "special")
    par(family = "Arial") 
    e<-c()
    f<-c()
    label<-"Percentage of permuted interactions"
    for(p in percentages) {
	if(p == 0) {
	    d <- read.table(paste(dir.name, 'biana_no_tap-omim/auc_phenotypes.dat', sep=""))
	    dir.name<-"../data/summary_runs_on_random/"
	}
	else {
	    d <- read.table(paste(dir.name, 'biana_no_tap_permuted_p', p, '-omim/auc_phenotypes.dat', sep=""))
	    if(!all(dim(d) == c(100, 5))) {
		stop("Object dimension is different from expected!")
	    }
	}
	m <- mean(d)
	n <- dim(d)[1]
	error <- qt(0.975, df=n-1) * sd(d) / sqrt(n)
	e <- cbind(e,m)
	f <- cbind(f,error)
    }

    for(i in 1:length(scoring.methods)) {
	if(i==1) {
	    plot(percentages, e[i,], col=cols[i], type="l", ylim=c(0,100), ylab="AUC (%)", xlab=label)
	}
	else {
	    lines(percentages, e[i,], col=cols[i])
	}
    }
    segments(percentages, e-f, percentages, e+f, col=cols[1:5], lty=1, lwd=1) 
    segments(percentages-2, e-f, percentages+2, e-f, col=cols[1:5], lty=1, lwd=0.5)
    segments(percentages-2, e+f, percentages+2, e+f, col=cols[1:5], lty=1, lwd=0.5)

    legend(0.5, 100, scoring.methods, fill=color5, bty="n", ncol=3)
    dev.off()

}


main()
