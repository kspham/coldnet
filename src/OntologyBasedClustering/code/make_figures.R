setwd("~/coldnet/src/OntologyBasedClustering/outputs")

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

summary(pairwise)
hist(pairwise$V1)

sd(pairwise)

jacD <- read.csv("jac_g11.txt", header=FALSE)
summary(jacD)
d <- density(jacD$V1)
plot(d)

library(gplots)
pdf("dendrogram.pdf")
sim_mat <- as.matrix(read.table("similarity_matrix.csv", header = FALSE, sep=","))
sim_mat[sim_mat == 0] <- NA
row.names(sim_mat) <- 1:345
col.names(sim_mat) <- 1:345
heatmap.2(sim_mat, dendrogram=c("none"), na.rm = TRUE, key = TRUE, keysize = 1.2, cexRow = 0.1, cexCol = 0.1, labCol = rownames(sim_mat), trace = "none")
title(main = "Similarity Matrix")
dev.off()

pdf("test.pdf")
test <- sim_mat[1:50,1:50]
heatmap.2(test, dendrogram=c("none"), na.rm = TRUE, key = TRUE, trace = 'none')
dev.off()