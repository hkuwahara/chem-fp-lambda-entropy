#!/usr/bin/env Rscript


toplevel.dir <- ".";
drugbank.dir <- "../drug-bank-data-analysis";

get.tanimoto.sim <- function( i, a, b ) {
  x <- a[i,];
  y <- b[i,];
  intsxn <- sum( x == 1 & y == 1 );
  intsxn / (1e-10 + sum( x == 1) + sum( y == 1) - intsxn);
}
fp.types <- c("maccs", "pubchem");
for( fp.type in fp.types ) {
  used.fp.file <- paste( drugbank.dir, fp.type, paste0("drug-reduced-0.3-", fp.type, "-fp.list"), sep="/");
  used.fp.list <- read.delim( used.fp.file, header = F );
  selected.fp.list <- as.character(unlist(used.fp.list));
  
  list.a <-  paste(toplevel.dir, paste0( "A-", fp.type, ".tsv"), sep="/");
  list.b <-  paste(toplevel.dir, paste0( "B-", fp.type, ".tsv"), sep="/");
  a.df <- read.delim( list.a, header = F, row.names = 1);
  b.df <- read.delim( list.b, header = F, row.names = 1);
  colnames(a.df) <- paste0( "F", 1:dim(a.df)[2]);
  colnames(b.df) <- paste0( "F", 1:dim(b.df)[2]);
  
  orginal.s <- sapply(1:dim(a.df)[1], get.tanimoto.sim, a.df, b.df );

  selected.a.df <- a.df[, selected.fp.list];  
  selected.b.df <- b.df[, selected.fp.list];  
  selected.s <- sapply(1:dim(selected.a.df)[1], get.tanimoto.sim, selected.a.df, selected.b.df );
  rec <- data.frame( id=rownames(a.df), v_0=orginal.s, v_0.3= selected.s);
  
  colnames(rec) <- c("id", paste0(fp.type, "_", c("0", "0.3")));
  out.file <- paste( toplevel.dir, paste0( "sim-", fp.type,  ".tsv"), sep="/");
  write.table(rec, out.file, col.names = T, row.names = F, sep="\t", quote=F);
}
