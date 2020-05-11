#!/usr/bin/env Rscript


toplevel.dir <- ".";
get.tanimoto.sim <- function( i, a, b ) {
  x <- a[i,];
  y <- b[i,];
  intsxn <- sum( x == 1 & y == 1 );
  intsxn / (1e-10 + sum( x == 1) + sum( y == 1) - intsxn);
}

fp.types <- c("maccs", "pubchem");
for( fp.type in fp.types ) {
  list.a <-  paste(toplevel.dir, paste0( "A", "-", fp.type, ".tsv"), sep="/");
  list.b <-  paste(toplevel.dir, paste0( "B", "-", fp.type, ".tsv"), sep="/");
  
  a.df <- read.delim( list.a, header = F, row.names = 1);
  b.df <- read.delim( list.b, header = F, row.names = 1);
  colnames(a.df) <- paste0( "F", 1:dim(a.df)[2]);
  colnames(b.df) <- paste0( "F", 1:dim(b.df)[2]);
  
  
  sim.val <- sapply(1:dim(a.df)[1], get.tanimoto.sim, a.df, b.df );
  rec <- data.frame( id=rownames(a.df), level0=sim.val );
  
  reduced.levels <- c("0.1", "0.2", "0.3");
  reduced.level <- "0.1";
  for( reduced.level in reduced.levels ) {
    fp.list.file <- paste(toplevel.dir, paste0( "blood", "-", fp.type, "-", reduced.level, "-fp-list.tsv"), sep="/");
    fp.list.df <-  read.delim( fp.list.file, header = F);
    fp.list <- as.character(unlist(fp.list.df[1,]));
    reduced.a.df <- a.df[, fp.list];
    reduced.b.df <- b.df[, fp.list];
    s <- sapply(1:dim(a.df)[1], get.tanimoto.sim, reduced.a.df, reduced.b.df );
    rec <- cbind( rec, data.frame(x=s) );
  }
  
  colnames(rec) <- c("id", paste0(fp.type, "-", c("0", reduced.levels)));
  out.file <- paste( toplevel.dir, paste0( "sim-", fp.type,  ".tsv"), sep="/");
  write.table(rec, out.file, col.names = T, row.names = F, sep="\t", quote=F);
}
