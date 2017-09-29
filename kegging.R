#!/usr/bin/env Rscript

.libPaths(c("/data/sw1/Dropbox/themetagenomics/packrat/lib/x86_64-pc-linux-gnu/3.3.3",
            "/data/sw1/Dropbox/themetagenomics/packrat/lib-ext",
            "/data/sw1/Dropbox/themetagenomics/packrat/lib-R"))

library(KEGGREST)
library(stringr)
library(Biostrings)

query <- keggGet('K01977')
genes <- query[[1]]$GENES

queries <- unlist(sapply(genes, function(gene){
  
  matches <- str_match_all(gene,'^([A-Z]+)\\: (.*)$')[[1]]
  org <- str_to_lower(matches[2])
  ids <- gsub('\\(.*\\)','',str_split(matches[3],' ')[[1]])
  
  paste0(org,':',ids)
  
},simplify=TRUE,USE.NAMES=FALSE))


cats <- c('query','entry','name','organism','position')
mat <- matrix('',length(queries),length(cats),dimnames=list(NULL,cats))
seqs <- DNAStringSet()
for (i in seq_along(queries)){
  
  query_info <- keggGet(queries[i])[[1]]
  mat[i,'query'] <- queries[i]
  
  query_update <- names(query_info)[names(query_info) %in% c('ENTRY','NAME','ORGANISM','POSITION')]
  for (q in query_update) mat[i,str_to_lower(q)] <- query_info[[q]]
  
  links <- gsub('^(.*)\\: (.*)$','\\2',query_info$DBLINKS)
  names(links) <- str_to_lower(gsub('^(.*)\\: (.*)$','\\1',query_info$DBLINKS))
  for (j in seq_along(links)){
    if (names(links)[j] %in% colnames(mat)){
      mat[i,names(links)[j]] <- links[j]
    }else{
      mat <- cbind(mat,'')
      colnames(mat)[ncol(mat)] <- names(links)[j]
      mat[i,names(links)[j]] <- links[j]
    }
    
  }
  
  seqs <- append(seqs,query_info$NTSEQ)
  names(seqs)[i] <- queries[i]

  if ((i-1) %% 50 == 0) cat(i,' ')
  if ((i-1) %% 500 == 0){
    out <- list(metadata=mat,sequences=seqs)
    attr(out,'iteration') <- i
    attr(out,'query') <- queries[i]
    saveRDS(out,file='~/kegg_out_tmp.rds')
  }
  
}

cnn <- table(gsub('^(.*)\\:.*$','\\1',mat[,'query']))
org <- names(cnn)

ko <- keggList('ko')
names(ko) <- gsub('^ko\\:(K[0-9]+)$','\\1',names(ko))

metagenome <- matrix(0,length(org),length(ko),dimnames=list(org,names(ko)))

for (o in org){
  cat('.')
  kegg <- keggLink('ko',o)
  kegg <- table(gsub('^ko\\:(K[0-9]+)$','\\1',kegg))
  
  overlap <- intersect(names(kegg),colnames(metagenome))
  metagenome[o,overlap] <- metagenome[o,overlap] + kegg[overlap]
}

out <- list(metadata=mat,sequences=seqs,cnn=cnn,ko=ko,metagenome=metagenome)
saveRDS(out,file='~/kegg_out.rds')
