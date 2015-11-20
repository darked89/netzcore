cols<-c(2,7,3,4,6,5,8) #rainbow(7) #2:8 #heat.colors(7)
color5<-cols[1:5]
color2<-c(3,4)
#scoring.methods<-c("N.Score", "N.Zcore", "N.Short", "N.Combo", "F.Flow", "P.Rank", "R.Walk", "N.Prop")
#scoring.methods.full<-c("NetScore", "NetZcore", "NetShort", "NetCombo", "Functional\nFlow", "Page\nRank", "Random\nWalk", "Network\n Propagation")
scoring.methods<-c("N.Score", "N.Zcore", "N.Short", "F.Flow", "P.Rank")
#scoring.method.ids<-c("ns", "nz", "nd", "nc3", "ff", "nr", "rw", "np")
scoring.method.ids<-c("ns", "nz", "nd", "ff", "nr")

dir.name<-"../data/summary/"
out.dir<-"../doc/draft/plos_one_plasticity/img/"

main <- function() {
    #manuscript()

    disease_category_figures() 
    #! manuscript2()

    #omim_similarity_figures() 
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
    module_figures() 
    disease_category_figures() 
    bc_case_study_figures()
}

lambda<-function(x) { x<-substring(x, 6); words<-unlist(strsplit(x, "_")); x<-paste(words, collapse=" "); return(x) }


###### NEIGHBORHOOD FIGURES ######
neighborhood_figures <- function() {

    dir.name<-"../data/summary_draft_before_revision/"
    # Number of genes in neighborhood of alzheimer over randomly permuted bPPI network 
    #cairo_ps(paste(dir.name, "Figure P1a.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    svg(paste(out.dir, "AD_a.svg", sep=""))
    par(family = "Arial") 
    par(mar=c(5, 4, 4, 5) + 0.1)
    d<-read.table(paste(dir.name, '../compare/biana_no_tap-omim_alzheimer-nn/permuted/results.dat', sep=""))
    a<-barplot(rbind(d$picked, d$picked_good),beside=T,xlab="Percentage of permuted interactions (%)",ylab="Number of genes", names.arg=seq(0,80,by=10),ylim=c(0,300)) # legend.text=c("All genes in n.hood", "AD genes in n.hood"),
    par(new=T)
    plot(colMeans(a),100*d$picked_good/d$picked,col=2,xaxt="n",yaxt="n",xlab="",ylab="",type='l',bty="n",ylim=c(0,10))
    axis(4, xpd=T, col=2, col.axis=2)
    mtext("Ratio (%)", side=4, line=3, col=2)
    d<-read.table(paste(dir.name, '../compare/biana_no_tap-omim_alzheimer-nn/permuted/results_random.dat', sep=""))
    ycoords<-c()
    for(i in seq(0, 80, by=10)) {
	e<-d[d$percentage==i,]
	ycoords<-c(ycoords, mean(100*e$picked_good/e$picked))
    }
    lines(colMeans(a), ycoords, lty=2, col=2)
    legend("topright", c("All genes in n.hood", "AD genes in n.hood", "Ratio (observed)", "Ratio (random)"), pch=c(15,15, NA, NA), pt.cex=2, col=c("grey30", "grey70",2,2), lty=c(0,0,1,2), bty="n")
    dev.off()

    # Number of genes in neighborhood of alzheimer over randomly pruned bPPI network 
    #cairo_ps(paste(dir.name, "Figure P1b.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    svg(paste(out.dir, "AD_b.svg", sep=""))
    par(family = "Arial") 
    par(mar=c(5, 4, 4, 5) + 0.1)
    d<-read.table(paste(dir.name, '../compare/biana_no_tap-omim_alzheimer-nn/pruned/results.dat', sep=""))
    a<-barplot(rbind(d$picked, d$picked_good),beside=T,xlab="Percentage of pruned interactions (%)",ylab="Number of genes", names.arg=seq(0,80,by=10),ylim=c(0,300)) #legend.text=c("All genes in n.hood", "AD genes in n.hood"),
    par(new=T)
    plot(colMeans(a),100*d$picked_good/d$picked,col=2,xaxt="n",yaxt="n",xlab="",ylab="",type='l',bty="n",ylim=c(0,10))
    axis(4, xpd=T, col=2, col.axis=2)
    mtext("Ratio (%)", side=4, line=3, col=2)
    d<-read.table(paste(dir.name, '../compare/biana_no_tap-omim_alzheimer-nn/pruned/results_random.dat', sep=""))
    ycoords<-c()
    for(i in seq(0, 80, by=10)) {
	e<-d[d$percentage==i,]
	ycoords<-c(ycoords, mean(100*e$picked_good/e$picked))
    }
    lines(colMeans(a), ycoords, lty=2, col=2)
    legend("topright", c("All genes in n.hood", "AD genes in n.hood", "Ratio (observed)", "Ratio (random)"), pch=c(15,15, NA, NA), pt.cex=2, col=c("grey30", "grey70",2,2), lty=c(0,0,1,2), bty="n")
    dev.off()
}




###### ROBUSTNESS FIGURES ######
robustness_figures <- function() {

    #dir.name<-"../data/summary_draft_before_revision/"
    dir.name<-"../data/summary_runs_on_random/"

    percentages<-seq(0,80,by=10) # c(0, 10, 30, 50, 70)
    scoring.methods<-c("N.Score", "N.Zcore", "N.Short", "F.Flow", "P.Rank")

    svg(paste(out.dir, "permuted2.svg", sep=""))
    coords <- seq(0.05,7.75,by=0.95)
    i = 1
    for(scoring.method.id in scoring.method.ids) {
	e<-matrix(nrow=100, ncol=length(percentages))
	j<-1
	for(p in percentages) {
	    if(p == 0) {
		d <- read.table(paste(dir.name, 'biana_no_tap-omim/auc_phenotypes.dat', sep=""))
	    }
	    else {
		d <- read.table(paste(dir.name, 'biana_no_tap_permuted_p', p, '-omim/auc_phenotypes.dat', sep=""))
		if(!all(dim(d) == c(100, 5))) {
		    stop("Object dimension is different from expected!")
		}
	    }
	    d<-d[,scoring.method.id]
	    e[,j]<-as.vector(d)
	    j<-j + 1
	}
	if(i == 1) {
	    boxplot(e,col=cols[i],at=coords, boxwex=0.15, pars=list(xaxt="n"), xlab="Perturbation level", ylab="AUC (%)", ylim=c(40,80), xlim=c(0.18,8.42))
	} 
	else {
	    boxplot(e,col=cols[i],at=coords, boxwex=0.15, pars=list(xaxt="n"), names=F, add=T)
	}
	if(i == 4) {
	    text.coords<-coords
	    #axis(1,tick=F,labels=scoring.methods,at=coords)
	}
	coords <- coords + 0.2
	i = i + 1
    }
    #legend(0.5, 119, c("Goh", "Entrez", "PPI", "bPPI", "weighted bPPI"), fill=color5, bty="n", ncol=3, cex=0.9) #horiz=T, x.intersp=0.5)
    text(text.coords+0.15, par("usr")[3]-3, srt=40, adj=1, labels=percentages, xpd=T, cex=0.9)
    dev.off()

    # Average AUC (%) at different levels of interaction permutation over OMIM disorders on bPPI network
    #cairo_ps(paste(dir.name, "Figure P2a.eps", sep=""), width = 6, height = 6, onefile = TRUE) # ps: horizontal = FALSE, onefile = FALSE, paper = "special")
    svg(paste(out.dir, "permuted.svg", sep=""))
    par(family = "Arial") 
    e<-c()
    f<-c()
    label<-"Percentage of permuted interactions"
    for(p in percentages) {
	if(p == 0) {
	    d <- read.table(paste(dir.name, 'biana_no_tap-omim/auc_phenotypes.dat', sep=""))
	}
	else {
	    d <- read.table(paste(dir.name, 'biana_no_tap_permuted_p', p, '-omim/auc_phenotypes.dat', sep=""))
	    if(!all(dim(d) == c(100, 5))) {
		stop("Object dimension is different from expected!")
	    }
	}
	d<-d[,scoring.method.ids]
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
    svg(paste(out.dir, "permuted_goh.svg", sep=""))
    par(family = "Arial") 
    e<-c()
    f<-c()
    label<-"Percentage of permuted interactions"
    for(p in percentages) {
	if(p == 0) {
	    d <- read.table(paste(dir.name, 'goh-omim/auc_phenotypes.dat', sep=""))
	}
	else {
	    d <- read.table(paste(dir.name, 'goh_permuted_p', p, '-omim/auc_phenotypes.dat', sep=""))
	    if(!all(dim(d) == c(100, 5))) {
		stop("Object dimension is different from expected!")
	    }
	}
	d<-d[,scoring.method.ids]
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
    svg(paste(out.dir, "pruned.svg", sep=""))
    par(family = "Arial") 
    e<-c()
    f<-c()
    label<-"Percentage of pruned interactions"
    for(p in percentages) {
	if(p == 0) {
	    d <- read.table(paste(dir.name, 'biana_no_tap-omim/auc_phenotypes.dat', sep=""))
	}
	else {
	    d <- read.table(paste(dir.name, 'biana_no_tap_pruned_p', p, '-omim/auc_phenotypes.dat', sep=""))
	    # Below for Average AUC (%) at different levels of interaction pruning only between non-seeds over OMIM disorders on bPPI network
	    #d <- read.table(paste(dir.name, 'biana_no_tap_pruned_non_seed_interactions_p', p, '-omim/auc_phenotypes.dat', sep=""))
	    if(!all(dim(d) == c(100, 5))) {
		stop("Object dimension is different from expected!")
	    }
	}
	d<-d[,scoring.method.ids]
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
    svg(paste(out.dir, "pruned_goh.svg", sep=""))
    par(family = "Arial") 
    e<-c()
    f<-c()
    label<-"Percentage of pruned interactions"
    for(p in percentages) {
	if(p == 0) {
	    d <- read.table(paste(dir.name, 'goh-omim/auc_phenotypes.dat', sep=""))
	}
	else {
	    d <- read.table(paste(dir.name, 'goh_pruned_p', p, '-omim/auc_phenotypes.dat', sep=""))
	    if(!all(dim(d) == c(100, 5))) {
		stop("Object dimension is different from expected!")
	    }
	}
	d<-d[,scoring.method.ids]
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

    # Average AUC (%) at different levels of seed permutation over OMIM disorders on bPPI network
    svg(paste(out.dir, "perturbed.svg", sep=""))
    par(family = "Arial") 
    e<-c()
    f<-c()
    label<-"Percentage of perturbed seeds"
    for(p in percentages) {
	if(p == 0) {
	    d <- read.table(paste(dir.name, 'biana_no_tap-omim/auc_phenotypes.dat', sep=""))
	}
	else {
	    d<-read.table(paste(dir.name, 'biana_no_tap-omim_perturbed_p', p, '/auc_perturbed.dat', sep=""))
	    if(!all(dim(d) == c(100, 5))) {
		stop("Object dimension is different from expected!")
	    }
	}
	d<-d[,scoring.method.ids]
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
    svg(paste(out.dir, "perturbed_goh.svg", sep=""))
    par(family = "Arial") 
    e<-c()
    f<-c()
    label<-"Percentage of perturbed seeds"
    for(p in percentages) {
	if(p == 0) {
	    d <- read.table(paste(dir.name, 'goh-omim/auc_phenotypes.dat', sep=""))
	}
	else {
	    d<-read.table(paste(dir.name, 'goh-omim_perturbed_p', p, '/auc_perturbed.dat', sep=""))
	    if(!all(dim(d) == c(100, 5))) {
		stop("Object dimension is different from expected!")
	    }
	}
	d<-d[,scoring.method.ids]
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

remove_go_term_redundancy<-function(input.file, out.file) {
    # Remove redundancy among seed go terms
    library(GOSemSim)
    d<-read.table(input.file, header=T)
    d2<-d[0,1:2]
    k<-1
    phenotypes<-as.vector(levels(d$phenotype))
    for(pheno in phenotypes) {
    	s <- d[d$phenotype == pheno, ]
	to.remove<-c()
	for(i in 1:nrow(s)) {
	    for(j in 1:nrow(s)) {
		if(i<j) {
		    go1<-as.character(s[i,"go"])
		    go2<-as.character(s[j,"go"])
		    a<-goSim(go1, go2, ont="BP", measure="Wang")
		    if(!is.na(a) & a>=0.9) {
			#print(c(go1, go2, a))
			if(s[i,"level"] < s[j,"level"]) {
			    to.remove<-c(to.remove,go2)
			} else {
			    to.remove<-c(to.remove,go1)
			}
		    }
		}
	    }
	}
	print(c(pheno, union(to.remove, to.remove))) # self union to convert it to set
	for(go in setdiff(d[d$phenotype == pheno,"go"],to.remove)) {
	    #d2<-rbind(d2, c(pheno, go))
	    d2[k,1]<-pheno
	    d2[k,2]<-go
	    k<-k+1
	}
    }
    names(d2)<-names(d)[1:2]
    write.table(d2, out.file, row.names=F)
    return()
}

remove_redundant_go_terms<-function(go.term.levels, similarity) {
    library(GOSemSim)
    d<-go.term.levels
    to.remove<-c()
    for(i in 1:nrow(d)) {
	for(j in 1:nrow(d)) {
	    if(i<j) {
		go1<-as.character(d[i,"go"])
		go2<-as.character(d[j,"go"])
		a<-goSim(go1, go2, ont="BP", measure="Wang")
		if(!is.na(a) & a>=similarity) {
		    #print(c(go1, go2, a))
		    if(d[i,"level"] < d[j,"level"]) {
			to.remove<-c(to.remove,go2)
		    } else {
			to.remove<-c(to.remove,go1)
		    }
		}
	    }
	}
    }
    print(union(to.remove, to.remove)) 
    go.terms.nonredundant<-setdiff(d[,"go"],to.remove)
    return(go.terms.nonredundant)
}


disease_category_figures<-function() {
    #common.up<-c("omim_breast_cancer", "omim_cardiomyopathy", "omim_diabetes", "omim_leukemia", "omim_obesity", "omim_parkinson_disease")
    #common.down<-c("omim_anemia", "omim_ataxia", "omim_cataract", "omim_epilepsy", "omim_hypertension", "omim_mental_retardation", "omim_schizophrenia", 
    #		    "omim_alzheimer", "omim_lung_cancer", "omim_lymphoma", "omim_myopathy", "omim_prostate_cancer", "omim_systemic_lupus_erythematosus")
    #print(common.up);
    #print(common.down);
    method<-"ns" #! 
    categories<-get_disease_categories(method)
    common.up<-categories$up
    common.down<-categories$down
    print(c("up", common.up))
    print(c("down", common.down))

    #return(); 
    # below uses "ns" specific file for go, 
    disease_category_comparison(method, common.up, common.down)
    # Below gives error in batch execution due to margins of the heatmaps but images can be generated in interactive mode
    disease_category_functional_comparison(common.up, common.down)
}


get_disease_categories<-function(method) {
    #common.up.old<-c('omim_anemia', 'omim_breast_cancer', 'omim_leukemia', 'omim_lymphoma', 'omim_systemic_lupus_erythematosus')
    #common.down.old<-c('omim_asthma', 'omim_ataxia', 'omim_cataract', 'omim_neuropathy', 'omim_schizophrenia', 'omim_spastic_paraplegia')
    #cols<-2:23
    #to.remove<-c("omim_insulin", "omim_neuropathy", "omim_asthma", "omim_spastic_paraplegia") # they are already out in the new analysis files
    #method<-"ns"
    percentages<-seq(0,80,by=10)
    dir.name<-"../data/summary_runs_on_random/"

    # Change on AUC of diseases over permuted and pruned interactions
    label<-"Percentage of perturbated interactions"
    container.permuted<-data.frame()
    container.pruned<-data.frame()
    for(p in percentages) {
	if(p == 0) {
	    d<-read.table(paste(dir.name, 'biana_no_tap-omim/auc_ppis.dat', sep=""))
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

    if(method == "ns") {
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
	    } else if(idx.permuted >= 6 & idx.pruned >= 6) { #! if >= several diseases are added to up
		common.up<-c(common.up, pheno)
	    }
	}
	write.table(cbind(container.permuted[,1], container.permuted[,1]-(container.permuted[,1]-50)/2, container.permuted[,6], container.pruned[,6]), paste(dir.name, 'critical_auc.dat', sep=""), sep="\t", col.names=c("AUC", "critical AUC", "50% swap AUC", "50% deletion AUC"))
    } else if(method == "nd") {
	# Robustness cutoff: 100-(100-auc)/2
	# new way of categorization - based on only pruning 
	common.up<-c()
	common.down<-c()
	for(pheno in phenotypes) {
	    #d<-container.permuted[pheno,] 
	    #print(c(pheno, d))
	    #if(d["0"] <= 50) { next }
	    #cutoff<-d["0"]/2+25 # same as above when auc/2 + 50 is used none of them is tolerant
	    #idx<-match(F, d>cutoff)
	    #if(is.na(idx)) { idx<-length(d) }
	    #idx.permuted<-idx
	    d<-container.pruned[pheno,] 
	    #print(c(pheno, d))
	    if(d["0"] <= 50) { next }
	    cutoff<-d["0"]/2+25 # same as above when auc/2 + 50 is used none of them is tolerant in case of pruned/permuted
	    #print(c("cutoff", cutoff)) 
	    idx<-match(F, d>cutoff)
	    if(is.na(idx)) { idx<-length(d) }
	    idx.pruned<-idx
	    #if(idx.permuted < 6 & idx.pruned < 6) {
	    #	common.down<-c(common.down, pheno)
	    #} else if(idx.permuted >= 6 & idx.pruned >= 6) {
	    #	common.up<-c(common.up, pheno)
	    #}
	    #print(c("idx", idx.pruned))
	    if(idx.pruned < 6) {
	    	common.down<-c(common.down, pheno)
	    } else if(idx.pruned > 6) { # to make it compatible with above
	    	common.up<-c(common.up, pheno)
	    }
	}
	write.table(cbind(container.permuted[,1], container.permuted[,1]/2+25, container.permuted[,6], container.pruned[,6]), paste(dir.name, 'critical_auc_nd.dat', sep=""), sep="\t", col.names=c("AUC", "critical AUC", "50% swap AUC", "50% deletion AUC"))
    }

    container<-container.permuted
    container.permuted.up<-container[common.up,]
    container.permuted.down<-container[common.down,]
    container<-container.pruned
    container.pruned.up<-container[common.up,]
    container.pruned.down<-container[common.down,]

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

    svg(paste(out.dir, "tolerant.svg", sep=""))
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

    svg(paste(out.dir, "intolerant.svg", sep=""))
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

    s<-read.table(paste(dir.name, 'biana_no_tap-omim/seeds.dat', sep=""))
    e<-intersect(rownames(s[s$n_seed<50,]),common.up)
    f<-intersect(rownames(s[s$n_seed<50,]),common.down)
    a<-wilcox.test(container.permuted[e,"50"], container.permuted[f,"50"])
    print(c("auc-permuted-seed50:", a$p.value))
    a<-wilcox.test(container.pruned[e,"50"], container.pruned[f,"50"])
    print(c("auc-pruned-seed50:", a$p.value))

    return(list(up=common.up, down=common.down));

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

disease_category_comparison<-function(method, common.up, common.down) {
    # Functional enrichment of modules in robust vs non-robust diseases
    #method<-"ns" # "nn"
    module.dir<-"../data/module/"
    d<-read.table(paste(module.dir, "biana_no_tap-omim/module_summary.dat", sep=""), header=T) 
    e2<-d[d$scoring==method & d$phenotype %in% common.up, "n_module"]
    f2<-d[d$scoring==method & d$phenotype %in% common.down, "n_module"]
    d<-read.table(paste(module.dir, "biana_no_tap-omim/module_summary_ns-only_bp.dat", sep=""), header=T) # before it was no_parent
    e<-d[d$scoring==method & d$phenotype %in% common.up,]
    f<-d[d$scoring==method & d$phenotype %in% common.down,]
    #g<-d[d$scoring==method & d$phenotype %in% non.common,]
    e<-cbind(e, n_module=e2)
    f<-cbind(f, n_module=f2)
    #i<-which(e$n_module == 0)
    #e.m<-e[-i,]
    #j<-which(f$n_module == 0)
    #f.m<-f[-j,]
    print(e)
    print("---")
    print(f)

    #labels<-c("Robust", "Uncharacterized", "Non-robust")
    labels<-c("Tolerant", "Non-tolerant")
    #svg(paste(out.dir, "module.svg", sep=""))
    #par(family = "Arial") 
    ##boxplot(e$n_module, g$n_module, f$n_module, col=8, names=labels, xlab="Disease category", ylab="Number of modules", ylim=c(0,20))
    #boxplot(e$n_module, f$n_module, col=8, names=labels, xlab="Disease category", ylab="Number of modules", ylim=c(0,15))
    #dev.off()
    a<-wilcox.test(e$n_module, f$n_module)
    print(c("module:", a$p.value))

    #svg(paste(out.dir, "seed_go_per_module.svg", sep=""))
    #par(family = "Arial") 
    #boxplot(100*e$n_seed_go_in_modules/e$n_seed_go, 100*f$n_seed_go_in_modules/f$n_seed_go, col=8, names=labels, xlab="Disease category", ylab="Average seed GO term enrichment among the modules (%)", ylim=c(0,100))
    #dev.off()
    a<-wilcox.test(e$n_seed_go_in_modules/e$n_seed_go, f$n_seed_go_in_modules/f$n_seed_go)
    print(c("seed go enrichment in modules: ", a$p.value)) #, mean(e$n_seed_go_in_modules/e$n_seed_go, na.rm=T), mean(f$n_seed_go_in_modules/f$n_seed_go, na.rm=T)))

    a<-wilcox.test(e$n_seed_go_in_modules/e$n_module, f$n_seed_go_in_modules/f$n_module)
    print(c("seed go in modules per module: ", a$p.value)) 

    #svg(paste(out.dir, "seed_go_to_all_go_in_module.svg", sep=""))
    #par(family = "Arial") 
    #boxplot(100*e$n_seed_go_in_modules/e$n_go_in_modules, 100*f$n_seed_go_in_modules/f$n_seed_go, col=8, names=labels, xlab="Disease category", ylab="Average coverage of seed GO terms among the modules (%)", ylim=c(0,100))
    #dev.off()
    a<-wilcox.test(e$n_seed_go_in_modules/e$n_go_in_modules, f$n_seed_go_in_modules/f$n_go_in_modules)
    print(c("ratio: ", a$p.value))

    #svg(paste(out.dir, "seed_go.svg", sep=""))
    #par(family = "Arial") 
    #boxplot(e$n_seed_go, f$n_seed_go, col=8, names=labels, xlab="Disease category", ylab="Number of seed GO terms", ylim=c(0,150))
    #dev.off()
    a<-wilcox.test(e$n_seed_go, f$n_seed_go)
    print(c("seed go: ", a$p.value))

    svg(paste(out.dir, "go_at_top.svg", sep=""))
    par(family = "Arial") 
    boxplot(e$n_go, f$n_go, col=8, names=labels, xlab="Disease category", ylab="Number of GO terms among top ranking genes") #, ylim=c(0,3))
    dev.off()
    a<-wilcox.test(e$n_go, f$n_go)
    print(c("go: ", a$p.value))

    #svg(paste(out.dir, "go.svg", sep=""))
    #par(family = "Arial") 
    #boxplot(e$n_go_in_modules, f$n_go_in_modules, col=8, names=labels, xlab="Disease category", ylab="Number of GO terms in the modules", ylim=c(0,350))
    #dev.off()
    a<-wilcox.test(e$n_go_in_modules, f$n_go_in_modules)
    print(c("go in modules: ", a$p.value))

    #svg(paste(out.dir, "go_per_module.svg", sep=""))
    #par(family = "Arial") 
    #boxplot(e$n_go_in_modules/e$n_module, f$n_go_in_modules/f$n_module, col=8, names=labels, xlab="Disease category", ylab="Average number of GO terms per module", ylim=c(0,100))
    #dev.off()
    a<-wilcox.test(e$n_go_in_modules/e$n_module, f$n_go_in_modules/f$n_module)
    print(c("go per module: ", a$p.value))

    #svg(paste(out.dir, "seed_go_in_module.svg", sep=""))
    #par(family = "Arial") 
    #boxplot(e$n_seed_go_in_modules, f$n_seed_go_in_modules, col=8, names=labels, xlab="Disease category", ylab="Number of seed GO terms in the modules", ylim=c(0,100))
    #dev.off()
    a<-wilcox.test(e$n_seed_go_in_modules, f$n_seed_go_in_modules)
    print(c("seed go in modules: ", a$p.value)) 

    s<-read.table(paste(dir.name, 'biana_no_tap-omim/seeds.dat', sep=""))
    e<-merge(e, s, by.x="phenotype", by.y="row.names")
    f<-merge(f, s, by.x="phenotype", by.y="row.names")
    svg(paste(out.dir, "seed_go_per_seed.svg", sep=""))
    par(family = "Arial") 
    boxplot(e$n_seed_go/e$n_seed, f$n_seed_go/f$n_seed, col=8, names=labels, xlab="Disease category", ylab="Number of seed GO terms per seed", ylim=c(0,3))
    dev.off()
    a<-wilcox.test(e$n_seed_go/e$n_seed, f$n_seed_go/f$n_seed)
    print(c("seed go per seed: ", a$p.value)) 

    a<-wilcox.test(e$n_go_in_modules/e$n_seed, f$n_go_in_modules/f$n_seed)
    print(c("go per seed in modules: ", a$p.value))

    # Comparison of n_seed in robust vs non-robust diseases
    svg(paste(out.dir, "seed.svg", sep=""))
    par(family = "Arial") 
    boxplot(s[common.up,"n_seed"], s[common.down,"n_seed"], col=8, names=labels, xlab="Disease category", ylab="Number of seeds", ylim=c(0,120))
    dev.off()
    a<-wilcox.test(s[common.up,"n_seed"], s[common.down,"n_seed"])
    print(c("seed:", a$p.value))
    # Comparison of n_path (path length)
    svg(paste(out.dir, "path.svg", sep=""))
    par(family = "Arial") 
    boxplot(s[common.up,"n_path"], s[common.down,"n_path"], col=8, names=labels, xlab="Disease category", ylab="Average length of seed connecting paths", ylim=c(0,5))
    dev.off()
    a<-wilcox.test(s[common.up,"n_path"], s[common.down,"n_path"])
    print(c("path:", a$p.value))

    a<-wilcox.test(s[common.up,"n_linker"], s[common.down,"n_linker"])
    print(c("n_linker:", a$p.value))

    a<-wilcox.test(s[common.up,"n_degree"], s[common.down,"n_degree"])
    print(c("degree:", a$p.value))

    # Comparison of number of alternative paths
    s<-read.table(paste(dir.name, 'biana_no_tap-omim/path_counts.dat', sep=""))

    svg(paste(out.dir, "alternative_path.svg", sep=""))
    par(family = "Arial") 
    boxplot(s[common.up,"n_path"], s[common.down,"n_path"], col=8, names=labels, xlab="Disease category", ylab="Average length of seed connecting paths") #, ylim=c(0,5))
    dev.off()
    a<-wilcox.test(s[common.up,"n_path"], s[common.down,"n_path"])
    print(c("alt. path:", a$p.value))
    a<-wilcox.test(s[common.up,"n_path_per_seed"], s[common.down,"n_path_per_seed"])
    print(c("alt. path per seed:", a$p.value))
    a<-wilcox.test(s[common.up,"path_length"], s[common.down,"path_length"])
    print(c("path length:", a$p.value))

    # Comparison with non-redundant number of go terms
    s<-read.table(paste(dir.name, 'biana_no_tap-omim/seeds.dat', sep=""))
    out.file<-paste(module.dir, "biana_no_tap-omim/module_summary_ns-seed_go_non_redundant.dat", sep="")
    #remove_go_term_redundancy(paste(module.dir, "biana_no_tap-omim/module_summary_ns-seed_go.dat", sep=""), out.file)
    d2<-read.table(out.file, header=T)
    #d2.up<-as.vector(table(as.vector(d2[d2$phenotype %in% common.up,]$phenotype)))
    #d2.down<-as.vector(table(as.vector(d2[d2$phenotype %in% common.down,]$phenotype)))
    d2.up<-c()
    for(pheno in common.up) {
	d2.up<-c(d2.up, nrow(d2[d2$phenotype==pheno,]))
    }
    d2.down<-c()
    for(pheno in common.down) {
	d2.down<-c(d2.down, nrow(d2[d2$phenotype==pheno,]))
    }

    a<-wilcox.test(d2.up, d2.down)
    print(c("seed go non-redundant:", a$p.value))
    a<-wilcox.test(d2.up/s[common.up,]$n_seed, d2.down/s[common.down,]$n_seed)
    print(c("seed go non-redundant per seed:", a$p.value))
    print(c("seed go non-redundant per seed - ratio:", a$p.value))

    out.file<-paste(module.dir, "biana_no_tap-omim/module_summary_ns-all_go_non_redundant.dat", sep="")
    #remove_go_term_redundancy(paste(module.dir, "biana_no_tap-omim/module_summary_ns-all_go.dat", sep=""), out.file)
    d2<-read.table(out.file, header=T)
    d2.up<-c()
    for(pheno in common.up) {
	d2.up<-c(d2.up, nrow(d2[d2$phenotype==pheno,]))
    }
    d2.down<-c()
    for(pheno in common.down) {
	d2.down<-c(d2.down, nrow(d2[d2$phenotype==pheno,]))
    }    
    a<-wilcox.test(d2.up, d2.down)
    print(c("all go non-redundant:", a$p.value))

    # Comparison of the age of the genes
    svg(paste(out.dir, "age.svg", sep=""))
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
}
    
disease_category_functional_comparison<-function(common.up, common.down) {
    library(gplots)
    library(RColorBrewer)
    #cols<-c(rep("red", length(common.up)), rep("grey", length(non.common)), rep("green", length(common.down)))
    cols<-c(rep("green", length(common.up)), rep("orange", length(common.down)))
    #val.cols <- brewer.pal(9,"Blues") #greenred
    # Common functions in pairwise  enrichment
    svg(paste(out.dir, "functional_similarity.svg", sep=""))
    dir.name<-"../data/module/"
    d<-read.table(paste(dir.name, "biana_no_tap-omim/functional_similarity_matrix.dat", sep=""), header=T)
    #e<-d[c(common.up, non.common, common.down),c(common.up, non.common, common.down)]
    e<-d[c(common.up, common.down),c(common.up, common.down)]
    e[e==0]<-NA
    e[upper.tri(e)]<-0
    palette.breaks <- seq(0, 1, by=0.0001)
    val.cols  <- colorRampPalette(c("#FFFFFF", "#6D6DFF", "#0000FF"))(length(palette.breaks) - 1)
    heatmap.2(as.matrix(e), revC=F, trace="none", dendrogram="none", col=val.cols, density.info="none", margins=c(12,12), RowSideColors=cols, Rowv=NA, Colv=NA, labRow=sapply(rownames(e), lambda), labCol=sapply(rownames(e), lambda), keysize=0.1,  breaks = palette.breaks, na.color="grey")
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

    #svg(paste(out.dir, "gene_similarity.svg", sep=""))
    #d<-read.table(paste(dir.name, "biana_no_tap-omim/gene_similarity_matrix.dat", sep=""), header=T)
    ##e<-d[c(common.up, non.common, common.down),c(common.up, non.common, common.down)]
    #e<-d[c(common.up, common.down),c(common.up, common.down)]
    #e[upper.tri(e)]<-NA
    #heatmap.2(as.matrix(e), revC=F, trace="none", dendrogram="none", col=val.cols, density.info="none", margins=c(12,12), RowSideColors=cols, Rowv=NA, Colv=NA, labRow=sapply(rownames(e), lambda), labCol=sapply(rownames(e), lambda), keysize=0.1)
    #dev.off()

    d<-read.table(paste(dir.name, "biana_no_tap-omim/phenotype_vs_functions.dat", sep=""), header=T)
    svg(paste(out.dir, "functions.svg", sep=""))
    selected<-c(common.up, common.down)
    cols<-c(rep("green", length(common.up)), rep("orange", length(common.down)))
    e<-d[,selected]
    e<-e[rowSums(e)>3,]
    e<-e[order(rownames(e)),]
    #e<-e[,colSums(e)>0]
    #val.cols<-4
    #Rowv <- rowMeans(e, na.rm = F)
    #a <- as.dendrogram(hclust(as.dist(e)))
    #a <- reorder(a, Rowv)
    a <- F
    #e[e==0]<-NA
    lambda2<-function(x) { words<-unlist(strsplit(x, "_")); x<-paste(words, collapse=" "); return(x) }
    labels.col<-as.vector(sapply(colnames(e), lambda))
    labels.col<-c(labels.col[-19], "systemic lup. ery.")
    heatmap.2(as.matrix(e), revC=F, trace="none", dendrogram="none", col=val.cols, margins=c(7,34), ColSideColors=cols, Colv=NA, Rowv=a, keysize=0.1, labCol=labels.col, labRow=sapply(rownames(e), lambda2), cexCol=0.9, cexRow=0.9, na.color="grey")
    dev.off()

    return()

    d<-read.table(paste(dir.name, "biana_no_tap-omim/phenotype_vs_functions.dat", sep=""), header=T)
    svg(paste(out.dir, "functions_tolerant.svg", sep=""))
    #selected<-c("omim_breast_cancer", "omim_lung_cancer", "omim_prostate_cancer", "omim_leukemia", "omim_diabetes", "omim_obesity", "omim_insulin")
    selected<-common.up
    e<-d[,selected]
    e<-e[rowSums(e)>1,]
    e<-e[,colSums(e)>0]
    lambda2<-function(x) { words<-unlist(strsplit(x, "_")); x<-paste(words, collapse=" "); return(x) }
    heatmap.2(as.matrix(e), revC=F, trace="none", dendrogram="none", col=val.cols, density.info="none", margins=c(7,30), keysize=0.1, labCol=sapply(colnames(e), lambda), labRow=sapply(rownames(e), lambda2), cexCol=0.9)
    dev.off()

    svg(paste(out.dir, "functions_intolerant.svg", sep=""))
    selected<-common.down
    e<-d[,selected]
    e<-e[rowSums(e)>1,]
    e<-e[,colSums(e)>0]
    lambda2<-function(x) { words<-unlist(strsplit(x, "_")); x<-paste(words, collapse=" "); return(x) }
    heatmap.2(as.matrix(e), revC=F, trace="none", dendrogram="none", col=val.cols, density.info="none", margins=c(7,30), keysize=0.1, labCol=sapply(colnames(e), lambda), labRow=sapply(rownames(e), lambda2), cexCol=0.7, cexRow=0.7)
    dev.off()

    d<-read.table(paste(dir.name, "biana_no_tap-omim/phenotype_vs_genes.dat", sep=""), header=T)
    svg(paste(out.dir, "genes_tolerant.svg", sep=""))
    selected<-common.up
    e<-d[,selected]
    e<-e[rowSums(e)>1,]
    e<-e[,colSums(e)>0]
    print(c("common.up", rownames(e)))
    heatmap.2(as.matrix(e), revC=F, trace="none", dendrogram="none", col=val.cols, density.info="none", margins=c(7,30), keysize=0.1, labCol=sapply(colnames(e), lambda), labRow=rownames(e), cexCol=0.9)
    dev.off()

    svg(paste(out.dir, "genes_intolerant.svg", sep=""))
    selected<-common.down
    e<-d[,selected]
    e<-e[rowSums(e)>1,]
    e<-e[,colSums(e)>0]
    print(c("common.down", rownames(e)))
    heatmap.2(as.matrix(e), revC=F, trace="none", dendrogram="none", col=val.cols, density.info="none", margins=c(7,30), keysize=0.1, labCol=sapply(colnames(e), lambda), labRow=rownames(e), cexCol=0.9)
    dev.off()
}


###### BC CASE STUDY FIGURES ######
bc_case_study_figures<-function() {
    data.dir <- "../data/summary_runs_on_random/breast_cancer_pruned/"
    # Remove redundancy among GO terms enriched in the component + seed GO terms
    #d<-read.table(paste(data.dir, "functional_comparison_goids.dat", sep=""), header=T)
    #print(length(d))
    #out.file<-paste(data.dir, "functional_comparison_goids-non_redundant.dat", sep="")
    #d2<-remove_redundant_go_terms(d, 0.8)
    #print(length(d2))
    #write(d2, out.file)
    #d2<-as.vector(read.table(out.file,header=F)[,1])
    
    # Draw heatmap of the functions among GO terms enriched in the component + seed GO terms
    library(gplots)
    d<-read.table(paste(data.dir, "functional_comparison.dat", sep=""), header=T, row.names=1, sep="\t", check.names=F)
    #d<-d[d2,]
    e<-d[order(d[,1]),2:3]
    e[e[,1]==1,1]<-2
    #e[e==0]<-NA
    rownames(e)<-d[order(d[,1]),1]
    svg(paste(out.dir, "functional_comparison.svg", sep=""))
    #lambda2<-function(x) { words<-unlist(strsplit(x, ".")); x<-paste(words, collapse=" "); return(x) }
    heatmap.2(as.matrix(e), revC=F, trace="none", dendrogram="none", col=c("lightgrey","blue","red"), margins=c(7,40), Colv=NA, Rowv=F, keysize=0.1, cexCol=0.9, cexRow=0.9, na.color="grey")
    dev.off()
}


###### MODULE FIGURES ######
module_figures<-function() {
    to.remove<-c("omim_insulin", "omim_neuropathy", "omim_asthma", "omim_spastic_paraplegia") 
    dir.name<-"../data/module/"
    d<-read.table(paste(dir.name, "biana_no_tap-omim/module_summary.dat", sep=""), header=T)
    i<-which(d$phenotype %in% to.remove)
    d<-d[-i,]
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
    svg(paste(out.dir, "modules.svg", sep=""))
    par(family = "Arial") 
    #col = c(8, cols), names=method.names
    boxplot(t(e.modules), col=8, names=NA, xlab="Prediction method", ylab="Number of modules", ylim=c(0,20), pars=list(xaxt="n"))
    range<-par("usr")[2] - par("usr")[1]
    text(par("usr")[1]-0.23+(1:6)*(range/6.5), par("usr")[3]-0.5, labels=method.names, cex=0.9, xpd=NA)
    #axis(1, par("usr")[1]-0.23+(1:6)*(range/6.5), labels=method.names, cex=0.9)
    dev.off()

    svg(paste(out.dir, "seed_go_ratio_in_modules.svg", sep=""))
    par(family = "Arial") 
    boxplot(t(e.ratio),col=8, names=NA, xlab="Prediction method", ylab="Average seed GO term enrichment of the modules (%)", ylim=c(0,70), pars=list(xaxt="n")) 
    range<-par("usr")[2] - par("usr")[1]
    text(par("usr")[1]-0.23+(1:6)*(range/6.5), par("usr")[3]-2, labels=method.names, cex=0.9, xpd=NA)
    dev.off()

    svg(paste(out.dir, "coverage.svg", sep=""), width = 6, height = 6, onefile = TRUE)
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
    tiff("omim.tif", width=2000, height=2000, res=300, compression="lzw")
    d<-read.table("../data/omim/2009_Aug_27/similarity.dat")
    heatmap(as.matrix(d), revC=T, col=val.cols, margins=c(9,9), Rowv=NA, Colv=NA, scale="none") # labRow=scoring.methods.full, labCol=scoring.methods.full) 
    dev.off()
    tiff("omim_in_ppi.tif", width=2000, height=2000, res=300, compression="lzw")
    d<-read.table("../data/omim/2009_Aug_27/similarity_in_ppi.dat")
    heatmap(as.matrix(d), revC=T, col=val.cols, margins=c(9,9), Rowv=NA, Colv=NA, scale="none") 
    dev.off()
    tiff("omim_extended.tif", width=2000, height=2000, res=300, compression="lzw")
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
    tiff(paste(dir.name, "Figure 1.tif", sep=""), width=2000, height=2000, res=300, compression="lzw")
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
    tiff(paste(dir.name, "Figure S2.tif", sep=""), width=2000, height=2000, res=300, compression="lzw")
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
    tiff(paste(dir.name, "Figure 2.tif", sep=""), width=2000, height=2000, res=300, compression="lzw")
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
    boxplot(x,col="lightgrey", at=coords, boxwex=0.4, pars=list(xaxt="n"), xlab="Prediction method", ylab="AUC(%)", ylim=c(0,100), xlim=c(0.6,8.7))
    #axis(1,tick=F,labels=scoring.methods,at=coords+0.2)
    text(coords+0.55, par("usr")[3]-3, srt=40, adj=1, labels=scoring.methods.full, xpd=T, cex=0.9)
    coords <- coords + 0.4
    boxplot(y,col="darkgrey", at=coords, boxwex=0.4, pars=list(xaxt="n"), names=F, add=T)
    legend(1.5, 5, c(paste("# of seeds < ", cutoff, sep=""), paste("# of seeds >= ", cutoff, sep="")), fill=c("lightgrey","darkgrey"), bty="n", horiz=T, cex=0.9)
    dev.off()

    #cairo_ps(paste(dir.name, "Supplementary Figure 3.eps", sep=""), width = 6, height = 6, onefile = TRUE)
    tiff(paste(dir.name, "Figure S1.tif", sep=""), width=2000, height=2000, res=300, compression="lzw")
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

    case_study_score_distributions()
}

case_study_score_distributions <- function() {

    ctd.dir.name<-"../data/ctd/"
    compare.dir.name<-"../data/compare/"

    #labels<-c("Alzheimer's Disease", "Diabetes", "AIDS")
    i = 1

    is.wholenumber <-function(x, tol = .Machine$double.eps^0.5)  { abs(x - round(x)) < tol }

    tiff(paste(dir.name, "Figure 3.tif", sep=""), width=2000, height=2000, res=300, compression="lzw")
    op<-par()
    par(mfrow=c(3,1))
    for(p in c("alzheimer", "diabetes", "aids")) {
	print(p)
	e<-read.table(paste(ctd.dir.name, p, "_nc.dat", sep=""), row.names=1, header=T)
	seeds<-read.table(paste(compare.dir.name, "biana_no_tap_relevance-new_omim_", p, "-gad/seeds.txt", sep=""))

	# direct vs rest
	print("#direct vs rest")
	d<-e[!rownames(e) %in% seeds$V1,]
	x<-d[which(d$is_direct==1), "nc_score"]
	y<-d[which(d$is_direct==0), "nc_score"]
	a<-t.test(x, y, alternative="greater")
	#print(c(a$p.value, format(mean(x), 2), format(mean(y), 2)))
	print(a$p.value)
	print(c(format(mean(x), 2), format(sd(x), 2)))
	print(c(format(mean(y), 2), format(sd(y), 2)))

	direct<-x
	rest<-y

	# direct vs non-ctd
	print("# direct vs non-ctd")
	d<-e[!rownames(e) %in% seeds$V1,]
	x<-d[which(d$is_direct==1), "nc_score"]
	y<-d[which(d$in_ctd==0), "nc_score"]
	a<-t.test(x, y, alternative="greater")
	print(a$p.value)
	print(c(format(mean(x), 2), format(sd(x), 2)))
	print(c(format(mean(y), 2), format(sd(y), 2)))
	
	nonctd<-y

	# ctd vs rest
	print("# ctd vs rest")
	d<-e[!rownames(e) %in% seeds$V1,]
	x<-d[which(d$in_ctd==1), "nc_score"]
	y<-d[which(d$in_ctd==0), "nc_score"]
	a<-t.test(x, y, alternative="greater")
	print(a$p.value)
	print(c(format(mean(x), 2), format(sd(x), 2)))
	print(c(format(mean(y), 2), format(sd(y), 2)))

	# direct vs indirect
	print("# direct vs indirect")
	d<-e[!rownames(e) %in% seeds$V1,]
	x<-d[which(d$is_direct==1), "nc_score"]
	y<-d[which(d$is_direct==0 & d$in_ctd==1), "nc_score"]
	a<-t.test(x, y, alternative="greater")
	print(a$p.value)
	print(c(format(mean(x), 2), format(sd(x), 2)))
	print(c(format(mean(y), 2), format(sd(y), 2)))

	indirect<-y

	print("")

	#f<-density(direct)
	#f<-cumsum(f$y*c(0,diff(f$x)))
	#f<-cumsum(f$y*c(0,diff(f$x)))

	breaks = seq(0, 0.4, by=0.1) 
	k = length(breaks) - 1

	d.cut = cut(direct, breaks, right=FALSE) 
	d.freq = table(d.cut) 
	d.cumfreq = cumsum(rev(d.freq))
	f = 100*d.cumfreq/cumsum(d.freq)[k]
	#f = spline(f)
	print(cbind(d.freq))
	print(cbind(f))
	
	plot(c(1,k), c(0, 100), type='n', xaxt="n", xlim=c(0.99,k+0.001), axes=F, xlab=paste("", "NetCombo score", sep=""), ylab="Cumulative % of genes") # , xlog=T, alpha = 0.5
	axis(1, 1:k, labels=c(">0.3", ">0.2", ">0.1", ">0.0"), pos=0, tcl=-0.3)
	#mtext(c(">0.3", ">0.2", ">0.1", ">0.0"), 1, at=1:k)
	axis(2, seq(0,100,by=20), col = "grey30", tcl = -0.3)
	axis(4, seq(0,100,by=20), col = "grey30", tcl = -0.3, pos=k)
	polygon(c(0, 1:k, k+0.001), c(0, f, 0), col="grey30", border=1)
	lines(f, col="black") # spline(f)
	#polygon(c(f$x,max(f$x)+0.001), c(f$y, 0), col="grey30", border=1)
	a<-f[3]
	points(3, a, pch=21, bg="white")
	rect(3+0.02, a+1, 3+0.13, a+11, col="white", border="grey")
	b<-round(a, digits=1)
	if(is.wholenumber(b)) # 50%
	    b<-paste(b, ".0", sep="")
	text(3+0.07, a+6, b, col="black")

	#d.cut = cut(rest, breaks, right=FALSE) 
	d.cut = cut(nonctd, breaks, right=FALSE) 
	d.freq = table(d.cut) 
	d.cumfreq = cumsum(rev(d.freq))
	f = 100*d.cumfreq/cumsum(d.freq)[k]
	#f = spline(f)
	print(cbind(d.freq))
	print(cbind(f))

	polygon(c(0, 1:k, k+0.001), c(0, f, 0), col="grey70", border=1)
	lines(f, col="black") 
	#polygon(c(f$x,max(f$x)+0.001), c(f$y, 0), col="grey70", border=1)
	a<-f[3]
	points(3, a, pch=21, bg="white")
	rect(3+0.02, a-2, 3+0.13, a+8, col="white", border="grey")
	text(3+0.07, a+3, round(a, digits=1), col="black")

	#rect(3+0.1, a-7, 3+0.3, a-1, col="white")
	#text(3+0.2, a-0.2, round(a, digits=2))

	#legend("topleft", c("direct CTD disease-gene associations", "rest (indirect or no association)"), fill=c("grey30", "grey70"), bty='n')
	legend("topleft", c("direct CTD disease-gene associations", "no association"), fill=c("grey30", "grey70"), bty='n')

	i = i+1
    }
    par(op)
    #box()
    dev.off()

    #require(splines)
    #yy <-predict(interpSpline(x, y))
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

