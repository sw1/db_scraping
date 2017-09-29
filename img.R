#!/usr/bin/env Rscript

library(XML)
library(RCurl)
library(stringr)
library(tidyverse)
library(Biostrings)

add_item <- function(x,item,key,n,mode,keys){
  if (is.null(x[[item]][key])){
    update <- vector(mode=mode,length=n)
    names(update) <- keys
    x[[item]] <- update
    attr(x,'update') <- TRUE
  }else{
    attr(x,'update') <- FALSE
  }
  return(x)
}

gg_to_img <- read_delim('~/rcrust/gg_13_5_img.txt.gz','\t') %>%
  dplyr::rename(gg_id=`#gg_id`) %>%
  mutate_all(as.character)

img_ids <- unique(gg_to_img$img_genome_id)
n_img <- length(img_ids)

IMG_META <- matrix(0.0,n_img,2,dimnames=list(img_ids,c('cn','gc')))
KEGG_META <- vector(mode='list')
COG_META <- vector(mode='list')
METACYC_META <- vector(mode='list')
KEGG <- vector(mode='list')
COG <- vector(mode='list')
METACYC <- vector(mode='list')
SEQS <- vector(mode='list',length=n_img)
names(SEQS) <- img_ids

url <- 'https://img.jgi.doe.gov/cgi-bin/m/'

for (idx in seq_along(img_ids)){

  cat('.')

  img <- img_ids[idx]

  genome_url <- sprintf('%smain.cgi?section=TaxonDetail&page=taxonDetail&taxon_oid=%s',url,img)
  genome_html <- htmlParse(readLines(genome_url))
  genome_q <- xpathSApply(genome_html,"//body//tr[./td[@class='img']]//td//a",saveXML)
  cn <- as.integer(gsub('^.*onclick=\\"\\">([0-9+])<\\/a>$','\\1',genome_q[grepl('16S',genome_q)]))
  genome_q <- xpathSApply(genome_html,"//body//tr[./td[@class='img']]",saveXML)
  gc <- as.numeric(gsub('^.*align=\\"right\\">(.*)%.*$','\\1',genome_q[grepl('G\\+C',genome_q)]))/100
  IMG_META[img,] <- c(cn,gc)

  agent <- 'Mozilla/5.0'
  curl <- getCurlHandle()
  curlSetOpt(cookiejar='',useragent=agent,followlocation=TRUE,curl=curl)
  s16_url <- getURL(sprintf('%smain.cgi?section=TaxonDetail&page=rnas&taxon_oid=%s&locus_type=rRNA&gene_symbol=16S',
                            url,img),curl=curl)
  s16_url <- gsub('^.*var myDataSource = new YAHOO.util.DataSource\\("(json.*&)"\\);.*$','\\1',s16_url)
  s16_html <- htmlParse(getURL(sprintf('%s%s',url,s16_url),curl=curl),asText=TRUE)
  s16_q <- xpathSApply(s16_html,"//input",saveXML)
  s16 <- unique(gsub('^.*value=\\"([0-9]+).*$','\\1',s16_q))

  seq <- vector(mode='character',length=length(s16))
  names(seq) <- s16
  for (s in seq_along(s16)){
    seq_url <- sprintf('%smain.cgi?exportGenes=1&exportType=nucleic&gene_oid=%s&up_stream=0&down_stream=0',url,s16[s])
    seq_html <- htmlParse(readLines((seq_url)))
    seq_q <- xpathSApply(seq_html,"//pre",saveXML)
    seq[s] <- gsub('^.*\\/>([A-Z]+).*$','\\1',gsub('\n','',seq_q))
  }
  SEQS[[img]] <- DNAStringSet(seq)


  agent <- 'Mozilla/5.0'
  curl <- getCurlHandle()
  curlSetOpt(cookiejar='',useragent=agent,followlocation=TRUE,curl=curl)
  ko_url <- getURL(sprintf('%smain.cgi?section=TaxonDetail&page=ko&taxon_oid=%s',
                           url,img),curl=curl)
  ko_url <- gsub('^.*var myDataSource = new YAHOO.util.DataSource\\("(json.*&)"\\);.*$','\\1',ko_url)
  ko_html <- htmlParse(getURL(sprintf('%s%s',url,ko_url),curl=curl),asText=TRUE)
  ko_q <- strsplit(as(ko_html,'character'),"<a href=\"main.cgi?section=TaxonDetail&amp;page=koGenes&amp;",fixed=TRUE)[[1]]
  ko_q <- ko_q[-length(ko_q)]


  agent <- 'Mozilla/5.0'
  curl <- getCurlHandle()
  curlSetOpt(cookiejar='',useragent=agent,followlocation=TRUE,curl=curl)
  cog_url <- getURL(sprintf('%smain.cgi?section=TaxonDetail&page=cogs&taxon_oid=%s',
                           url,img),curl=curl)
  cog_url <- gsub('^.*var myDataSource = new YAHOO.util.DataSource\\("(json.*&)"\\);.*$','\\1',cog_url)
  cog_html <- htmlParse(getURL(sprintf('%s%s',url,cog_url),curl=curl),asText=TRUE)
  cog_q <- strsplit(as(cog_html,'character'),"<a href=\"main.cgi?section=TaxonDetail&amp;page=cogGeneList&amp;",fixed=TRUE)[[1]]
  cog_q <- cog_q[-1]


  agent <- 'Mozilla/5.0'
  curl <- getCurlHandle()
  curlSetOpt(cookiejar='',useragent=agent,followlocation=TRUE,curl=curl)
  metacyc_url <- getURL(sprintf('%smain.cgi?section=TaxonDetail&page=metacyc&taxon_oid=%s',
                            url,img),curl=curl)
  metacyc_url <- gsub('^.*var myDataSource = new YAHOO.util.DataSource\\("(json.*&)"\\);.*$','\\1',metacyc_url)
  metacyc_html <- htmlParse(getURL(sprintf('%s%s',url,metacyc_url),curl=curl),asText=TRUE)
  metacyc_q <- strsplit(as(metacyc_html,'character'),'<a href=\"main.cgi?section=TaxonDetail&amp;page=metaCycGenes&amp;',fixed=TRUE)[[1]]
  metacyc_q <- metacyc_q[-1]


  for (i in seq_along(ko_q)){

    item <- gsub('^.*KO\\:(K[0-9]+).*$','\\1',ko_q[[i]])
    KEGG <- add_item(KEGG,item,img,n_img,'integer',img_ids)
    KEGG[[item]][img] <- as.integer(gsub('^.*GeneCount\\":\\"([0-9]+)\\".*$','\\1',ko_q[[i]]))

    if (attr(KEGG,'update')){
      KEGG_META <- add_item(KEGG_META,item,NULL,2,'character',c('def','pw'))
      KEGG_META[[item]][1] <- gsub('^.*Definition\\":\\"(.*)\\",\\"KOID.*$','\\1',ko_q[[i]])
      pw <- gsub('^.*(\\[EC.*\\]).*$','\\1',ko_q[[i]])
      KEGG_META[[item]][2] <- ifelse(nchar(pw) == nchar(ko_q[[i]]),'',pw)
    }

  }

  for (i in seq_along(cog_q)){

    item <- gsub('^.*(COG[0-9]+).*$','\\1',cog_q[[i]])
    COG <- add_item(COG,item,img,n_img,'integer',img_ids)
    if(is.na(as.integer(gsub('^.*GeneCount\\":\\"([0-9]+)\\".*$','\\1',cog_q[[i]])))) print(i)
    COG[[item]][img] <- as.integer(gsub('^.*GeneCount\\":\\"([0-9]+)\\".*$','\\1',cog_q[[i]]))

    if (attr(COG,'update')){
      COG_META <- add_item(COG_META,item,NULL,1,'character',c('name'))
      cog_name <- trimws(gsub('^.*\\COGName\\":\\"(.*)\\",\\"COGNameDisp.*$','\\1',cog_q[[i]]))
      COG_META[[item]][1] <- ifelse(nchar(cog_name) == nchar(cog_q[[i]]),'',cog_name)
    }

  }

  for (i in seq_along(metacyc_q)){

    item <- gsub('^.*unique_id=(.*?)&amp;.*$','\\1',metacyc_q[[i]])
    METACYC <- add_item(METACYC,item,img,n_img,'integer',img_ids)
    METACYC[[item]][img] <- as.integer(gsub('^.*GeneCount\\":\\"([0-9]+)\\".*$','\\1',metacyc_q[[i]]))

    if (attr(COG,'update')){
      METACYC_META <- add_item(COG_META,item,NULL,1,'character',c('pw'))
      metacyc_name <- trimws(gsub('^.*MetaCycPathway\\":\\"(.*)\\",\\"GeneCount\\".*$','\\1',metacyc_q[[i]]))
      metacyc_name <- gsub('<i>|<\\/i>','',metacyc_name)
      METACYC_META[[item]][1] <- ifelse(nchar(metacyc_name) == nchar(metacyc_q[[i]]),'',metacyc_name)
    }

  }

  if ((idx-1) %% 100 == 0){
    cat(sprintf('\n\nidx:%s\n',idx))
    OUT <- list(gg_to_img=gg_to_img,
                img_meta=IMG_META,
                kegg_meta=KEGG_META,
                cog_meta=COG_META,
                metacyc_meta=METACYC_META,
                kegg=KEGG,
                cog=COG,
                seqs=SEQS,
                metacyc=METACYC)
    attr(OUT,'tmp_idx') <- idx
    saveRDS(OUT,file='~/kegging/img_out_tmp.rds')
  }

}

OUT <- list(gg_to_img=gg_to_img,
            img_meta=IMG_META,
            kegg_meta=KEGG_META,
            cog_meta=COG_META,
            metacyc_meta=METACYC_META,
            kegg=KEGG,
            cog=COG,
            seqs=SEQS,
            metacyc=METACYC)
saveRDS(OUT,file='~/kegging/img_out.rds')
