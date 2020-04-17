#!/usr/bin/env Rscript

library("dplyr");
library("data.table");


args = commandArgs(trailingOnly=TRUE);
if(length(args) < 4) {
  stop("requires chr, strand, sample-dir, and annotation-dir");
}
chr <- args[1];
strand <- args[2];
sample.dir <- args[3];
annotation.dir <- args[4];


out.dir <- paste( sample.dir,  "sj_exon_map", sep="/");
if(!dir.exists( out.dir ) ) {
  dir.create( out.dir, recursive = T);
}

sj.dir <- paste( sample.dir,  "sj", sep="/");

measure.closest.delta.with.id <- function( a.pos, a.id, b.pos, b.id) {
  delta.list <- list( a.id=c(), b.id=c(), delta=c() );
  k <- 1;
  range.l <- 0;
  range.h <- 0;
  len.a <- length(a.pos);
  len.b <- length(b.pos);
  
  j <- 1;
  for( i in 1:len.a ) {
    query <- a.pos[[i]];
    min.delta <- Inf;
    range.h <- len.b;
    while( j <= len.b ) {
      target <- b.pos[[j]];
      delta <- target - query;
      if( abs(delta) < abs(min.delta)) {
        min.delta <- delta;
        range.l <- j;
        j <- j + 1;
      } else if( abs(delta) == abs(min.delta) ) {
        j <- j + 1;
      } else {
        range.h <- j - 1;
        for(m in range.l:range.h) {
          delta.list$a.id[[k]] <- a.id[[i]];
          delta.list$b.id[[k]] <- b.id[[m]];
          delta.list$delta[[k]] <- min.delta;
          k <- k + 1;
        }
        j <- range.l;
        break;
      }
    }
  }
  as.data.frame(delta.list);
}  
  

id <- paste0( chr, strand );
exon.list.file <- paste0( annotation.dir, "/exon_list_", id, ".tsv" );
sj.file <- paste0( sj.dir, "/",  id);

exon.df <- fread(file=exon.list.file, data.table = F, col.names=c("left","right","gene", "exon"));
sj.df <- fread(file=sj.file, data.table = F, 
               col.names=c("chr","first_intron","last_intron","strand","intron_motif","annotation","unique_map","multi_map","max_overhang")); 

exon.df <- exon.df[order(exon.df$right),];
sj.df <- sj.df[order(sj.df$first_intron),];
left.sj.side <- measure.closest.delta.with.id( sj.df$first_intron - 1, sj.df$first_intron, exon.df$right, exon.df$exon  );

exon.df <- exon.df[order(exon.df$left),];
sj.df <- sj.df[order(sj.df$last_intron),];
right.sj.side <- measure.closest.delta.with.id( sj.df$last_intron + 1, sj.df$last_intron, exon.df$left, exon.df$exon  );

sj.input <- data.frame( first_intron=sj.df$first_intron, last_intron=sj.df$last_intron, count=sj.df$unique_map);

temp.df <- full_join( x=sj.input, y=left.sj.side, by=c("first_intron" = "a.id"));
setnames( temp.df, c("b.id", "delta"), c("left_exon", "left_delta") );

sj.exon.df <- full_join( x=temp.df, y=right.sj.side, by=c("last_intron" = "a.id"));
setnames( sj.exon.df, c("b.id", "delta"), c("right_exon", "right_delta") );

sj.exon.df <- sj.exon.df[order(sj.exon.df$first_intron, sj.exon.df$last_intron),];

out.file <- paste0(out.dir, "/", id,  ".tsv");
write.table(format(distinct(sj.exon.df), scientific = F, trim = T), out.file, col.names = T, row.names = F, sep="\t", quote=F);

