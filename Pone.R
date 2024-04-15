setwd("/Users/khristinafediv/Documents/bioinformatics_course/pone")
fname = "pone.fas"
InstallPackages = TRUE
if (InstallPackages) {
  if (!requireNamespace("BiocManager", quietly=TRUE)) 
    install.packages("BiocManager")
    install.packages("Biostrings")
    install.packages("BiocGenerics")
  BiocManager::install("msa")
  BiocManager::install("ggtree")
  
  install.packages("ade4")
  install.packages("adegenet")
  install.packages("ape")
  install.packages("ggtree")
  install.packages("ggplot2")
  install.packages("ips")
  install.packages("bios2mds")
  install.packages("haplotypes")
  install.packages("pegas")
  install.packages("phytools")
  install.packages("stats")
  install.packages("treeio")
  install.packages("amap")
  install.packages("e1071")
  install.packages("scales")
  install.packages("cluster")
  install.packages("rgl")
  install.packages("pegas")
  install.packages("phytools")
}

library(BiocManager)
library(Biostrings)
library(BiocGenerics)
library(ade4)
library(adegenet)
library(ape)
library(ggtree)
library(ggplot2)
library(stats)
library(ips)
library(msa)

library(amap)
library(e1071)
library(scales)
library(cluster)
library(rgl)

AlignNeeded = TRUE
if (AlignNeeded){
  file <- readDNAStringSet(fname)
  cb<- msa(file)
  cv<-msaConvert(cb, type=c("bios2mds::align"))
  library(bios2mds)
  export.fasta(cv, outfile = "outfile.fas", ncol(cb), open = "w")
}

nbin<-as.DNAbin(cb)

TRIM = TRUE
if (TRIM) {
  nbin<-trimEnds(nbin)
}

an<-as.alignment(nbin)
nm<-as.matrix(an)
nbinmat<-as.matrix(labels(nbin))
class(nbin)
dnbin<-dist.dna(nbin, model = "K80")
tree<-nj(dnbin)
ggt<-ggtree(tree, cex = 0.8, aes(color=branch.length))+scale_color_continuous(high='lightskyblue1', low='coral4')+geom_tiplab(align=TRUE, size=2)+geom_treescale(y = 5, color = "coral4", fontsize = 4)
njmsaplot<-msaplot(ggt, nbin, offset = 0.009, width=1, height = 0.5, color = c(rep("rosybrown",1), rep("sienna1", 1), rep("lightgoldenrod1",1), rep("lightskyblue1", 1)))
njmsaplot

dev.new()
njmsaplot

pdf("njmsaplot.pdf", width = 11, height = 9)
njmsaplot
dev.off()

nrow(nm)
ncol(nm)

sat2 <- NULL
for (i in 1:nrow(nm)) {
  sat2[i] <- paste(nm[i, ], collapse="")
}

sat2 <- toupper(sat2)
sat3 <- unique(sat2)
sat3
hfreq <- NULL
for (i in 1:length(sat3)) {
  hcount = 0 
  s3 <-sat3[i]
  for (j in 1:length(sat2)) {
    s2 <- sat2[j]
    if (s3==s2) {
      hcount <- (hcount+1)
    }
  }
  hname<-(paste("H",i, sep=""))
  hfreq[i] <- hcount
}

len <- nchar(sat3[1])
cnt <- 1
sat4 = list()
for (j in 1:len) {
  same <- TRUE
  first <- substr(sat3[1], j, j)
  for (i in 2:length(sat3)) {
    ch1 <- substr (sat3[i], j, j)
    if (first !=ch1) {
      str <- paste(j, first, ch1)
      print(str)
      same <- FALSE
      break
    }
  }
  if (!same) {
    ss <- NULL
    for (i in 1:length(sat3)) {
      ss <- paste(ss, substr(sat3[i], j, j), sep="")
    }
    sat4[cnt] <- ss
    cnt <- cnt +1
  }
}

len <- nchar(sat3[1])
cnt <- 1
sat5 = list()
for (j in 1:len) {
  same <- TRUE
  first <- substr(sat3[1], j, j)
  scol <- first
  for (i in 2:length(sat3)) {
    ch1 <- substr(sat3[i], j, j)
    scol <- paste(scol, ch1, sep="")
    if (first != ch1) {
      str <-paste(j, first, ch1)
      same <- FALSE
    
    }
  }
  if (!same) {
    scol <- paste("V_", cnt, " ", scol, sep="")
    ss <- NULL
    for (i in 1:length(sat3)) {
      ss <- paste(ss, substr(sat3[i], j, j), sep="")
    }
    sat5[cnt] <- ss
    cnt <- cnt + 1
  }
}

sat6 <- as.matrix(sat5)
mat6 = matrix(nrow=nrow(sat6), ncol=nchar(sat6[1]))
for (i in 1:nrow(mat6)) {
  s <- as.vector(strsplit(as.character(sat5[i]), ""))
  for (j in 1:ncol(mat6)) {
  mat6[i, j] <- as.character(s[[1]] [j])
  }
}

mat7 <- t(mat6)
write.table(mat7, file="mat7.txt", quote=FALSE, sep="\t")
hname<-paste("H", 1:nrow(mat7), sep = "")
rownames(mat7)=hname
write.table(mat7,file="mat7.txt", quote=FALSE, sep="\t") 

str4 <- NULL
str4[1] <- paste(mat7[1, ], collapse="")
for (i in 2:nrow(mat7)) {
  tmp <- NULL
  for (j in 1:ncol(mat7)) {
    chr = "."
    if(mat7[i, j] != mat7[1, j]) chr = mat7[i, j]
    tmp <- paste(tmp, chr, sep="")
  }
  str4[i] <- paste(tmp, collapse="")
}
nchar(str4[1]) 
mstr4 <- as.matrix(str4)
rownames(mstr4) <- hname
colnames(mstr4) <- paste("sequences length","(", ncol(mat7), "base pairs", ")")
pct <- round((as.matrix(hfreq)*100/colSums(as.matrix(hfreq))), 2)
colnames(pct) <- c("pct")
cmstr4 <- as.data.frame(cbind(mstr4, hfreq, pct))
cmstr4
write.table(cmstr4, file="cmstr4.txt", quote=FALSE, sep="\t") 

library(haplotypes)

kn<-as.dna(nbin)
kh<-haplotypes::haplotype(kn)

ncb <- as.matrix(labels(nbin))
n2 <- NULL
for (i in 1:nrow(ncb)) {
  n2[i] <- strsplit(ncb[i], '_')[[1]][1]
}
n2

for (i in 1:nrow(ncb)) {
  n2[i] <- strsplit(n2[i], ' ')[[1]][1] 
}
n2

hf<-grouping(kh, factors=n2)
hf[["hapvec"]] <- NULL
dhf<-as.data.frame(hf$hapmat) 
rownames(dhf)<-paste("H", 1:nrow(mat7), sep = "")
dhf
write.table(dhf,file="dhf.txt", quote=FALSE, sep="\t")

library(pegas)
dh<-dist.hamming(mat7) 
dh
dhm<-as.matrix(dh)
write.table(dhm,file="dhm.txt", quote=FALSE, sep="\t")

sat2 <- NULL
for (i in 1:nrow(nm)) {
  sat2[i] <- paste(nm[i, ], collapse="")
}
sat2
sat2 <- toupper(sat2)
sat3 <- unique(sat2)
comat = matrix(nrow=length(sat3), ncol=length(sat3))
for (i in 1:length(sat3)) { 
  si <- sat3[i]
  for (j in 1:length(sat3)) { 
    sj <- sat3[j]
    difcnt = 0
    s1 = as.vector(strsplit(as.character(si), ""))
    s2 = as.vector(strsplit(as.character(sj), ""))
    for (k in 1:length(s1[[1]])) {
      if (s1[[1]][k] != s2[[1]][k]) {
        difcnt = difcnt + 1
      }
      comat[i, j] = difcnt
    }
  }
}
comat	
colnames(comat)<-paste("H", 1:nrow(comat), sep = "")
rownames(comat)<-paste("H", 1:nrow(comat), sep = "")
heatmap(comat,scale="none",col=heat.colors(100),keep.dendro=TRUE, symm=TRUE) #stats package

dev.new()
heatmap(comat,scale="none",col=heat.colors(100),keep.dendro=TRUE, symm=TRUE)


pdf("heatmap.pdf", width = 7, height = 7)
heatmap(comat,scale="none",col=heat.colors(100),keep.dendro=TRUE, symm=TRUE)
dev.off()


h<-pegas::haplotype(nbin, strict = FALSE, trailingGapsAsN = TRUE)
hname<-paste("H", 1:nrow(h), sep = "")
rownames(h)= paste(hname)
net<-haploNet(h, d = NULL, getProb = TRUE)
net
ind.hap<-with(
  utils::stack(setNames(attr(h, "index"), rownames(h))),
  table(hap=ind, individuals=rownames(nbin))
)

par(mar=c(0.01,0.01,0.01,15))
plot(net, size=attr(net, "freq"), scale.ratio = 2, cex = 0.6, labels=TRUE, pie = ind.hap, show.mutation=1, font=2, fast=TRUE)
legend(x= 57,y=15, colnames(ind.hap), fill=rainbow(ncol(ind.hap)), cex=0.52, ncol=6, x.intersp=0.2, text.width=11)

dev.new()
plot(net, size=attr(net, "freq"), scale.ratio = 2, cex = 0.6, labels=TRUE, pie = ind.hap, show.mutation=1, font=2, fast=TRUE)
legend(x= 57,y=30, colnames(ind.hap), fill=rainbow(ncol(ind.hap)), cex=0.5, ncol=6, x.intersp=0,003, text.width=9)

pdf("hapind.pdf", width = 11, height = 5)
par(mar=c(0.01,0.01,0.01,15))
plot(net, size=attr(net, "freq"), scale.ratio = 2, cex = 0.6, labels=TRUE, pie = ind.hap, show.mutation=1, font=2, fast=TRUE)
legend(x= 57,y=15, colnames(ind.hap), fill=rainbow(ncol(ind.hap)), cex=0.52, ncol=6, x.intersp=0.2, text.width=11)
dev.off()

h<-pegas::haplotype(nbin, strict = FALSE, trailingGapsAsN = TRUE)#extracts haplotypes from DNAbin object
hname<-paste("H", 1:nrow(h), sep = "")
rownames(h)= paste(hname)
net<-haploNet(h, d = NULL, getProb = TRUE) #constructs the haplotype network
net
ind.hap<-with(
  utils::stack(setNames(attr(h, "index"), rownames(h))),
  table(hap=ind, individuals=rownames(nbin)[values])
)
bg<-c(rep("dodgerblue4", 15), rep("olivedrab4",15), rep("royalblue2", 15), rep("red",15), rep("olivedrab3",15), 
      rep("skyblue1", 15), rep("olivedrab1", 15),  rep("darkseagreen1", 15))
hapcol<-c("Aksu", "Demre", "Kumluca", "Firm", "Bayatbadem", "Geyikbayir", "Phaselis", "Termessos")
ubg<-c(rep("dodgerblue4",1), rep("royalblue2",1), rep("skyblue1",1), rep("red",1), rep("olivedrab4",1), 
       rep("olivedrab3",1), rep("olivedrab1",1), rep("darkseagreen1",1))

par(mar=c(0.001,0.001,0.001,0.001))
plot(net, size=attr(net, "freq"), bg = bg, scale.ratio = 2, cex = 0.7, labels=TRUE, pie = ind.hap, show.mutation=1, font=2, fast=TRUE)
legend(x=-36,y=53, hapcol, fill=ubg, cex=0.8, ncol=1, bty="n", x.intersp = 0.2)

dev.new()
par(mar=c(0.001,0.001,0.001,0.001))
plot(net, size=attr(net, "freq"), bg = bg, scale.ratio = 2, cex = 0.7, labels=TRUE, pie = ind.hap, show.mutation=1, font=2, fast=TRUE)
legend(x=-36,y=53, hapcol, fill=ubg, cex=0.8, ncol=1, bty="n", x.intersp = 0.2)

pdf("happop.pdf", width = 7, height = 4)
par(mar=c(0.001,0.001,0.001,0.001))
plot(net, size=attr(net, "freq"), bg = bg, scale.ratio = 2, cex = 0.7, labels=TRUE, pie = ind.hap, show.mutation=1, font=2, fast=TRUE)
legend(x=-36,y=53, hapcol, fill=ubg, cex=0.8, ncol=1, bty="n", x.intersp = 0.2)
dev.off()


nmname<-c(rep("Nature", 3), rep("Nature", 6), rep("Nature", 6), rep("Nature", 4), rep("Firm", 4), 
          rep("Nature", 11), rep("Greenhouse", 3), rep("Nature", 3), rep("Greenhouse", 2), rep("Greenhouse", 12), 
          rep("Nature", 6), rep("Firm", 4),  rep("Nature", 7), rep("Greenhouse", 13), rep("Nature", 6), rep("Greenhouse", 15), 
          rep("Firm", 3), rep("Nature", 6), rep("Nature", 2), rep("Firm", 4))

ng<-nbin
rownames(ng)<-nmname
hg<-pegas::haplotype(ng, strict = FALSE, labels = hname, trailingGapsAsN = TRUE) #extracts haplotypes from DNAbin object
hg
netg<-haploNet(hg, d = NULL, getProb = TRUE) #constructs the haplotype network
netg
ind.hapg<-with(
  utils::stack(setNames(attr(hg, "index"), rownames(hg))),
  table(hap=ind, individuals=rownames(ng)[values])
)
gbg<-c(rep("red"), rep("blue"), rep("green"))

par(mar=c(0.001,0.001,0.001, 0.001))
plot(netg, size=attr(netg, "freq"), bg = gbg, scale.ratio = 2, cex = 0.7, labels=TRUE, pie = ind.hapg, show.mutation=1, font=2, fast=TRUE)
legend(x=-35,y=45, colnames(ind.hapg), fill=c("red","blue", "green"), cex=0.8, ncol=1, bty = "n", x.intersp = 0.2)

dev.new()
par(mar=c(0.001,0.001,0.001, 0.001))
plot(netg, size=attr(netg, "freq"), bg = gbg, scale.ratio = 2, cex = 0.7, labels=TRUE, pie = ind.hapg, show.mutation=1, font=2, fast=TRUE)
legend(x=-35,y=45, colnames(ind.hapg), fill=c("red","blue", "green"), cex=0.8, ncol=1, bty = "n", x.intersp = 0.2)


pdf("hapgrp.pdf", width = 8, height = 5) #save as pdf file
par(mar=c(0.001,0.001,0.001, 0.001))
plot(netg, size=attr(netg, "freq"), bg = gbg, scale.ratio = 2, cex = 0.7, labels=TRUE, pie = ind.hapg, show.mutation=1, font=2, fast=TRUE)
legend(x=-35,y=45, colnames(ind.hapg), fill=c("red","blue", "green"), cex=0.8, ncol=1, bty = "n", x.intersp = 0.2)
dev.off()

library(pegas)
library(ggtree)
library(ggplot2)
library(phytools)
library(stats)

nname<-as.matrix(sort(labels(nbin)))
krp<-list(Aksu = nname[1:15,],
          Bayatbadem = nname[16:30,],
          Demre = nname[31:45,],
          Firm = nname[46:60,],
          Geyikbayir = nname[61:75,],
          Phaselis = nname[76:90,],
          Kumluca = nname[91:105,],
          Termessos = nname[106:120,])

nbin
class(nbin)
dnbin<-dist.dna(nbin, model = "K80") 
tree<-nj(dnbin) 

emos<-ggtree(tree, layout = 'circular', branch.length='branch.length', lwd = 0.5)+xlim(-0.1, NA)
groupOTU(emos, krp, 'Populations') + aes(color=Populations)+theme(legend.position="right")+geom_tiplab(names(nbin), cex = 1.7, offset=0.002)+guides(color = guide_legend(override.aes = list(size = 2.5)))+geom_treescale(x=-0.1, color = "coral4", fontsize = 3, offset = 9)

dev.new()
emos<-ggtree(tree, layout = 'circular', branch.length='branch.length', lwd = 0.5)+xlim(-0.1, NA)
groupOTU(emos, krp, 'Populations') + aes(color=Populations)+theme(legend.position="right")+geom_tiplab(names(nbin), cex = 1.7, offset=0.002)+guides(color = guide_legend(override.aes = list(size = 2.5)))+geom_treescale(x=-0.1, color = "coral4", fontsize = 3, offset = 9)


pdf("njtreepop.pdf", width = 6.5, height = 5)
emos<-ggtree(tree, layout = 'circular', branch.length='branch.length', lwd = 0.5)+xlim(-0.1, NA)
groupOTU(emos, krp, 'Populations') + aes(color=Populations)+theme(legend.position="right")+geom_tiplab(names(nbin), cex = 1.7, offset=0.002)+guides(color = guide_legend(override.aes = list(size = 2.5)))+geom_treescale(x=-0.1, color = "coral4", fontsize = 3, offset = 9)
dev.off()


nbin
class(nbin)
dnbin<-dist.dna(nbin, model = "K80")
tree<-nj(dnbin)

njdistree<-ggtree(tree,layout = 'circular', branch.length='branch.length', aes(color=branch.length), lwd = 0.5)+xlim(-0.1, NA)+
  geom_tiplab(names(nbin), size = 1.7, offset=0.002)+scale_color_continuous(high='lightskyblue1',low='coral4')+
  geom_treescale(x=-0.1, color = "coral4", fontsize = 3, offset = 9) 

njdistree

dev.new()
njdistree


pdf("njdistree.pdf", width = 6.5, height = 5)
njdistree
dev.off()


library(treeio)
D <- dist.hamming(mat7)
class(D)
htre<-nj(D)#This function performs the neighbor-joining tree estimation of Saitou and Nei (1987).
bp <- boot.phylo(htre, mat7, B=100, function(x) nj(dist.hamming(x)))
bp2 <- data.frame(node=1:Nnode(htre) + Ntip(htre), bootstrap = bp)
htree <- full_join(htre, bp2, by="node")
boothap<-ggtree(htree, size=1, branch.length='branch.length')+geom_tiplab(size=4)+
  geom_nodepoint(aes(fill=cut(bootstrap, c(0,50,70,85,100))), shape=21, size=4)+
  theme_tree(legend.position=c(0.85, 0.2))+ 
  scale_fill_manual(values=c("black", "red", "pink1", "white"), guide='legend', name='Bootstrap Percentage (BP)',breaks=c('(85,100]', '(70,85]', '(50,70]', '(0,50]'), labels=expression(BP>=85, 70<=BP*"<85",50<=BP*"<70", BP<50))

boothap

dev.new()
boothap


pdf("boothap.pdf", width = 11, height = 5)
boothap
dev.off()

