library(stats)
library(ade4)
library(ape)
library(adegenet)
library(phangorn)



dna <- fasta2DNAbin(file="../data/clustal/COI-align.fasta")
D <- dist.dna(dna, model="TN93")
temp <- as.data.frame(as.matrix(D))
table.paint(temp, cleg=0, clabel.row=.5, clabel.col=.5)

temp <- t(as.matrix(D))
temp <- temp[,ncol(temp):1]
par(mar=c(1,5,5,1))
image(x=1:6, y=1:6, temp, col=rev(heat.colors(100)), xaxt="n", yaxt="n",
      xlab="",ylab="")
axis(side=2, at=1:6, lab=rev(rownames(dna)), las=2, cex.axis=.5)
axis(side=3, at=1:6, lab=rownames(dna), las=3, cex.axis=.5)

tre <- nj(D)
pdf('../data/clustal/COI-NJ.pdf')
plot(tre, cex=.6)
title("A simple NJ tree - COI")
dev.off()

dna2 <- as.phyDat(dna)
#tre.ini <- nj(dist.dna(dna,model="raw"))
#tre.pars <- optim.parsimony(tre.ini, dna2)
trees <- bab(dna2)
pdf('../data/clustal/COI-MP.pdf')
plot(trees[[1]], cex=.6)
title("MP tree - COI")
dev.off()


dna <- fasta2DNAbin(file="../data/clustal/cytb-align.fasta")
D <- dist.dna(dna, model="TN93")
tre <- nj(D)
pdf('../data/clustal/cytb-NJ.pdf')
plot(tre, cex=.6)
title("A simple NJ tree - cyt b")
dev.off()

dna2 <- as.phyDat(dna)
trees <- bab(dna2)
pdf('../data/clustal/cytb-MP-1.pdf')
plot(trees[[1]], cex=.6)
title("MP tree - cyt b - 1")
dev.off()

pdf('../data/clustal/cytb-MP-2.pdf')
plot(trees[[2]], cex=.6)
title("MP tree - cyt b - 2")
dev.off()

dna <- fasta2DNAbin(file="../data/clustal/mito-align.fasta")
D <- dist.dna(dna, model="TN93")
tre <- nj(D)
pdf('../data/clustal/mito-NJ.pdf')
plot(tre, cex=.6)
title("A simple NJ tree - mitochondria")
dev.off()

dna2 <- as.phyDat(dna)
trees <- bab(dna2)
pdf('../data/clustal/mito-MP.pdf')
plot(trees[[1]], cex=.6)
title("MP tree - mitochondria")
dev.off()