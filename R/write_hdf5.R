write_h5 <- function(dataFrame, filename = 'data.hf5', ens = TRUE, write_labels = TRUE, len_str=3){
  require(org.Mm.eg.db)
  require(rhdf5)
  if(ens==T){
    ensIds = rownames(dataFrame)
    counts = as.matrix(dataFrame)
    
    gene_names = ensIds
    x <- org.Mm.egSYMBOL
    # Get the gene symbol that are mapped to an entrez gene identifiers
    mapped_genes <- mappedkeys(x)
    # Convert to a list
    xx <- as.list(x[mapped_genes])
    xxenseg <- as.list(org.Mm.egENSEMBL2EG)
    gene_syms=unlist(xx[unlist(xxenseg[gene_names])])
    gene_names_list<-(lapply(xxenseg[gene_names],function(x){if(is.null(x)){x=NA}else{x=x[1]}}))
    sym_names=unlist(lapply(xx[unlist(gene_names_list)],function(x){if(is.null(x)){x=NA}else{x=x[1]}}))
    sym_names[is.na(sym_names)]=gene_names[is.na(sym_names)]
    if(write_labels==T)
      {labels=as.numeric(as.factor(sapply(colnames(dataFrame),function(x)substr(x,1,len_str))))
       class_labels = sort(unique(sapply(colnames(data),function(x)substr(x,1,len_str))))
       h5save(counts,ensIds, labels,sym_names,class_labels,file=filename)
       #h5save(counts,ensIds, labels,sym_names,class_labels,file=filename)
    }else{
    h5save(counts,ensIds,sym_names,file=filename)
    #h5save(counts,ensIds,sym_names,file=filename)
    }
  }else{
    ensIds = rownames(dataFrame)
    counts = as.matrix(dataFrame)
    if(write_labels==T)
      {labels=as.numeric(as.factor(sapply(colnames(dataFrame),function(x)substr(x,1,len_str))))
       class_labels = sort(unique(sapply(colnames(data),function(x)substr(x,1,len_str))))
       h5save(counts,ensIds, labels,class_labels,file=filename)
       }else{
       #h5save(counts,ensIds, labels,class_labels,file=filename)
         h5save(counts,ensIds,file=filename)
         #h5save(counts,ensIds,file=filename)
                        }
  }
}