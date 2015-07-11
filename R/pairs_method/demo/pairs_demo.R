require(ggplot2)
require(grid)
load("../core/pairs_functions.RData")





#The package comes already with a list of marker pairs selected from our training dataset. 
#So, if you want, you can use it straight away with your data!
#
#If a new list of markers needs to be found (e.g., with a different training dataset), please use the commands below



#FIND THE MARKER PAIRS #########


#training.data = matrix with training data (genes are in rows and cells in columns)
#id.G1, id.S, id.G2M = names of G1, S and G2M samples in the training data
#genes.cc = genes that must be used for training. By default, these are cell-cycle annotated genes.


#Filter out genes that are not in the training dataset or are not expressed in any cell
genes.training<-intersect(row.names(training.data), genes.cc)
genes.training<-genes.training[rowSums(training.data[genes.training,])!=0]


#Find marker pairs, i.e., pairs of genes that change their relative ranking in phase 1 compared to phase 2 and 3.
#id1 = names of the samples in phase 1
#id2, id3 = names of the samples in phase 2 and 3
#genes.list = list of genes to consider
#thr.frac = fraction of mismatches allowed (default: 0.5)
G1.marker.pairs<-find.markers(id1=id.G1, 
                              id2=id.S, id3=id.G2M, 
                              genes.list=genes.training, thr.frac=0.5)

S.marker.pairs<-find.markers(id1=id.S, 
                             id2=id.G1, id3=id.G2M, 
                             genes.list=genes.training,thr.frac=0.5)

G2M.marker.pairs<-find.markers(id1=id.G2M, 
                               id2=id.G1, id3=id.S, 
                               genes.list=genes.training,thr.frac=0.5)

#The output is a data frame that looks like this:
#Gene.1             Gene.2 Sign
#ENSMUSG00000000001 ENSMUSG00000001785    1
#ENSMUSG00000000001 ENSMUSG00000005470    1
#ENSMUSG00000000001 ENSMUSG00000012443    1
#...
#Gene.1 expression level is greater than Gene.2 in at least a fraction "thr.frac" of cells in the phase "id1".
#The opposite is true (i.e., Gene.1 < Gene.2) in at least a fraction "thr.frac" of cells in the phases "id2" and "id3"



#ALLOCATE CELLS TO CELL-CYCLE STAGE #########

#load data: a matrix with data must be loaded. Genes are in rows, cells in columns.
require(rhdf5)
qseq.data<-h5read("../../data/normCounts_mESCquartz.h5f","LogNcountsQuartz")
gene.names<-h5read("../../data/normCounts_mESCquartz.h5f","gene_names_all")
cell.names<-h5read("../../data/normCounts_mESCquartz.h5f","cell_names")
row.names(qseq.data)<-gene.names
colnames(qseq.data)<-cell.names
#qseq.data =This is the dataset published in Sasagawa et al, Genome Biology, 14:R31, 2013


#Label samples according to their phase
labels.qseq<-factor(c(rep("G1", length(grep(colnames(qseq.data), pattern="ES_G1"))),
         rep("S", length(grep(colnames(qseq.data), pattern="ES_S"))),
         rep("G2M", length(grep(colnames(qseq.data), pattern="ES_G2M")))),
       levels=c("G1","S", "G2M"))






#Allocate cells to cell-cycle phase
#
#data = matrix with data
#N = Number of times random pairs are sampled
#Nmin.couples = In each random sample, pairs of genes with equal expression levels (ties) are excluded. Nmin.couples is the minimum number of pairs (exc. ties) that must be used for allocation. 
#Nmin = A random sample fails when it generates less than Nmin.couples pairs with no ties. Nmin is the minimum number of random samples that must not fail in order to carry out the allocation.
#G1.marker.pairs, S.marker.pairs, G2M.marker.pairs = Marker pairs for G1, S and G2M phase
#genes.training = list of genes used for training
results.qseq<-predictor(data = qseq.data, 
                        N = 1e3,  Nmin.couples = 50, Nmin = 100,
                        G1.marker.pairs = G1.marker.pairs, 
                        S.marker.pairs = S.marker.pairs, 
                        G2M.marker.pairs = G2M.marker.pairs, 
                        genes.training = genes.training)
#The output "results.qseq" is a list with the following elements:
#G1.marker.pairs, S.marker.pairs, G2M.marker.pairs = Marker pairs used for allocation
#Scores = Unnormalised G1, S and G2M scores. 
#Normalised.scores = Normalised G1,S and G2M scores (i.e., score of a phase divided by the sum of all scores). 



#Plot the results with the unnormalised scores (we found this to be more reliable than the allocation based on the normalised scores).
#Cells with a G1 or a G2M score greater than 0.5 are allocated to the phase with the highest score between G1 and G2M. 
#If both the G1 and the G2M scores are less than 0.5, the cell is allocated to the S phase.
df.scores<-data.frame(x=results.qseq$Scores[,"score.G1"], 
                      y=results.qseq$Scores[,"score.G2M"], 
                      z=labels.qseq )

plot.scores(df.scores,legend.title="True cell-cycle phase",plot.title="Cell cycle allocation")
#df.scores = data frame with the G1 and G2M unnormalised scores (x and y) and data labels (z)






#Plot the results with the normalised scores. 
#Each cell is allocated to the phase with the highest normalised score
df.norm.scores<-data.frame(x=results.qseq$Normalised.scores[,"score.G1"], 
                           y=results.qseq$Normalised.scores[,"score.G2M"], 
                           z=labels.qseq )

plot.normalised.scores(df=df.norm.scores,legend.title="True cell-cycle phase",plot.title="Cell cycle allocation")
#df.norm.scores = data frame with the G1 and G2M normalised scores (x and y) and data labels (z)




