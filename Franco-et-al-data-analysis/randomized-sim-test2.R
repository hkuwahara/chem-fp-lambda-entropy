#!/usr/bin/env Rscript

toplevel.dir <- ".";

args = commandArgs(trailingOnly=TRUE);
# fp.types <- c("maccs", "pubchem");
fp.type <- args[1];
num.keys <- as.integer(args[2]);
num.samples <- as.integer(args[3]);

# toplevel.dir <- "/home/hiro/research/kaust/chem-fp/revision/chem-fp-lambda-entropy/Franco-et-al-data-analysis";
# fp.type <- "maccs";
# reduced.level <- "0.3";
# tmp.df <- read.delim(paste(toplevel.dir, fp.type, paste0("A-reduced-",reduced.level, "-", fp.type, ".tsv"), sep="/"), header = T, row.names = 1);
# used.fp <- colnames(tmp.df);
# num.keys <- dim(tmp.df)[2];
# num.samples <- 10;

# pairs.df <- read.delim( paste(toplevel.dir, "sim_testset.tsv", sep="/"), header = T);
# target.pairs <- pairs.df[,1];
# target.pairs <- target.pairs[pairs.df[,4] >= 80];

list.a <- paste(toplevel.dir, fp.type, paste0("A-reduced-0-", fp.type, ".tsv"), sep="/");
a.df <- read.delim( list.a, header = T, row.names = 1);
list.b <- paste(toplevel.dir, fp.type, paste0("B-reduced-0-", fp.type, ".tsv"), sep="/");
b.df <- read.delim( list.b, header = T, row.names = 1);

fp.size <- dim(a.df)[2];

results.mat <- matrix( ncol=dim(a.df)[1], nrow=num.samples );

get.tanimoto.sim <- function( i, a, b ) {
  x <- a[i,];
  y <- b[i,];
  intsxn <- sum( x == 1 & y == 1 );
  intsxn / (1e-10 + sum( x == 1) + sum( y == 1) - intsxn);
}

# selected.indices <- sample( fp.size, num.keys );
# selected.a <- a.df[, selected.indices];
# selected.b <- b.df[, selected.indices];
# pair.size <- dim(a.df)[1];
# sim <- sapply( 1:pair.size, get.tanimoto.sim, selected.a, selected.b );

apply.tanimoto.sim <- function( index ) {
  selected.indices <- sample( fp.size, num.keys );
  selected.a <- a.df[, selected.indices];
  selected.b <- b.df[, selected.indices];
  pair.size <- dim(a.df)[1];
  sim <- sapply( 1:pair.size, get.tanimoto.sim, selected.a, selected.b );
  # results.mat[index, ] <- sim;
}

results.mat <- sapply( 1:num.samples, apply.tanimoto.sim);

rec <- as.data.frame( results.mat );
rec <- cbind( rownames(a.df), rec );
colnames(rec) <- c( "id", paste0("run-", seq( 1, num.samples)) );

out.file <- paste( toplevel.dir, paste0( "random-results-", fp.type,  "-keys-", num.keys, ".tsv"), sep="/");
write.table(rec, out.file, col.names = T, row.names = F, sep="\t", quote=F);

