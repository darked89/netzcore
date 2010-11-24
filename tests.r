
dirname<-"../data/summary/"

# Significance of getting Alzheimer genes in Krauthammer data set
sum(dhyper(6:61,55,7224,61))


d<-read.table(paste(dirname, "biana_no_tap-omim/auc_ppis.dat", sep=""))
e<-read.table(paste(dirname, "biana_no_tap-omim/seeds.dat", sep=""))
print(cor(d,e))

#d<-read.table(paste(dirname, "biana-all/auc_ppis.dat", sep=""))
#e<-read.table(paste(dirname, "biana_no_tap-all/auc_ppis.dat", sep=""))
#print(wilcox.test(d[,4], e[,4]))
#print(kruskal.test(list(d[,4],e[,4])))
#print(kruskal.test(list(d[,1],e[,1])))

d<-read.table(paste(dirname, "all_vs_all/auc_phenotypes.dat", sep=""))
#x<-0; y<-0; for(i in d) { x<-x+1; for(j in d) y<-y+1; e<-i-j; if(abs(e)>5) { print(c(x,y,e)) } }
e<-abs(d[1,]-d[1,])>5
print(e)


networks<-c("biana-all", "biana_no_tap-all", "biana_no_tap_relevance-all", "goh-all", "entrez-all")
methods<-c("ns", "nz", "nd", "ff", "nr")

x<-c()
y<-1
for(k in 1:5) {
    print(c("---", methods[k], "---"))
    for(i in 1:length(networks)) {
	for(j in 1:length(networks)) {
	    if(i > j) {
		d<-read.table(paste(dirname, networks[i], "/auc_ppis.dat", sep=""))
		e<-read.table(paste(dirname, networks[j], "/auc_ppis.dat", sep=""))
		a<-wilcox.test(d[,k], e[,k])
		if(a$p.value < 0.05) {
		    #print(c(methods[k], networks[i], networks[j]))
		    print(c(networks[i], networks[j]))
		    print(a$p.value)
		    x[y]<-a$p.value
		    y<-y+1
		    #print(kruskal.test(list(d[,k],e[,k])))
		}
	    }
	}
    }
}

print(x)
print(min(x))

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

