#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE);
if(length(args) < 2) {
  stop("requires analysis.dir, id");
}
analysis.dir <- args[1];
id <- args[2];


input.file <- paste( analysis.dir, paste0(id, ".tsv"), sep="/");
# dat <- read.delim( input.file, header = F, col.names = c("hmdb_id", paste0("F", 1:166)), row.names = 1);
input.dat <- read.delim( input.file, header = T, row.names = 1);

compounds.num <- length(input.dat[,1]);
compound.ids <- rownames(input.dat);
sim.matrix <- matrix( NA, nrow=compounds.num, ncol=compounds.num, dimnames = list(row.names(input.dat), row.names(input.dat)) );


get.tanimoto.sim <- function( x, y ) {
  intsxn <- sum( x == 1 & y == 1 );
  intsxn / ( sum( x == 1) + sum( y == 1) - intsxn);
}

apply.tanimoto.sim.on.dat2 <- function( ref.compound, target.compounds ) {
  apply( target.compounds, 1, get.tanimoto.sim, ref.compound  );
}


r <- input.dat[1,];
v <- apply.tanimoto.sim.on.dat2( r, input.dat );
a <- rep( rownames(r), compounds.num );
b <- compound.ids;
rec <- data.frame( a=a, b=b, sim=v );

for( i in 1:(compounds.num-1) ) {
  r <- input.dat[i+1,];
  v <- apply.tanimoto.sim.on.dat2( r, input.dat[-(1:i),] );
  a <- rep( rownames(r), compounds.num - i );
  b <- compound.ids[-(1:i)];  
  rec <- rbind( rec, data.frame( a=a, b=b, sim=v ) );
}

out.file <- paste( analysis.dir, paste0("sim-", id, "-melted.tsv"), sep="/");
write.table(rec, out.file, col.names = T, row.names = F, sep="\t", quote=F);

