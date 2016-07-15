setwd("C:/Users/Debha/Documents/UCSD/Gage Lab/test/OntologyBasedClustering")

pdf("data_fig.pdf")

sum_gene_muts <- read.csv("sum_gene_muts.csv", header=FALSE)
summary(sum_gene_muts$V2)
hist(sum_gene_muts$V2,main= "Mutated Subjects per Gene")
boxplot(sum_gene_muts$V2, horizontal = TRUE, main = "Mutated Subjects per Gene")
sum_gene_muts <- sum_gene_muts[with(sum_gene_muts, order(-V2)),]
plot(sum_gene_muts$V3, sum_gene_muts$V2, main="Mutations in Gene over Gene Length", xlab="Gene Length", ylab="Mutations")
top10 = round(.1*dim(sum_gene_muts)[1])
write.csv(sum_gene_muts[1:top10,], file = "top_genes.csv",row.names=FALSE, na="")

sum_sub_muts <- read.csv("sum_sub_muts.csv", header=FALSE)
summary(sum_sub_muts$V2)
hist(sum_sub_muts$V2,main= "Mutated Genes per Subject")
boxplot(sum_sub_muts$V2, horizontal = TRUE, main = "Mutated Genes per Subject")

dev.off()

pairwise <- read.csv("pairwise_distance.txt", header=FALSE)
d <- density(pairwise$V1)
plot(d)
summary(pairwise)
hist(pairwise$V1)

sd(pairwise)

jacD <- read.csv("jac_g11.txt", header=FALSE)
summary(jacD)
d <- density(jacD$V1)
plot(d)