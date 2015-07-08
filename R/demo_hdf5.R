#generate hdf5 files for python
source('~/Code/ccClassify/R/write_h5.R')
#read in text file with normalised reads
data = read.table('../data/liver.reads.Deng.txt')
#data = read.table('../data/GSE46980_ID_CombinedMoleculeCounts_clean.tab.txt')
#look at data
head(data)

#write out hdf5 file for ccClassify

write_h5(data, filename = './data_hipsci.h5',len_str=4)

#genes_mmus=inpIDMapper(rownames(data),'HOMSA','MUSMU',srcIDType='ENSEMBL',destIDType='ENSEMBL')
#genes_mmus_vec = unlist(lapply(genes_mmus,function(x)x[1]))
#genes_ss_vec = names(lapply(genes_mmus,function(x)names(x)))

#data = data[rownames(data) %in% genes_ss_vec,]
#rownames(data) = genes_mmus_vec


# #load and write new gene lists
# dataCB=read.table(file='~/Downloads/genes_cycle_base2.txt', header=T)
# hu2musAll=inpIDMapper(dataCB[1:50,3],'HOMSA','MUSMU',srcIDType='ENSEMBL',destIDType='ENSEMBL')
# Hists1 <- read.table('~/Downloads/histones_1.txt', sep = '\t')
# Hists2 <- read.table('~/Downloads/histones_2.txt', sep = '\t')
# 
# 
# gene_names = c(as.character(Hists1[2:30,2]), as.character(Hists2[2:88,2]))
# 
# 
# xxsymeg <- as.list(org.Hs.egSYMBOL2EG)
# gene_names_eg<-unlist(xxsymeg[gene_names])
# 
# x <- org.Hs.egENSEMBL
# # Get the entrez gene IDs that are mapped to an Ensembl ID
# mapped_genes <- mappedkeys(x)
# # Convert to a list
# xx <- as.list(x[mapped_genes])
# ens_hist = unlist(xx[as.character(gene_names_eg)])
# ensMmus=(inpIDMapper(as.character(ens_hist),'HOMSA','MUSMU',srcIDType='ENSEMBL',destIDType='ENSEMBL'))
# 
# gene_set = c(as.character(unlist(hu2musAll)), as.character(ensMmus))
# h5save(gene_set, file = './small_gene_set.h5')

