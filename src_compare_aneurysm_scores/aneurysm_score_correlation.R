
df<-read.table("aneurysm_seed_scores.txt")
cor(df[,1], 1/df[,2])
#plot(df[,1], 1/df[,2], xlab="tm score", ylab="1 / p value")

