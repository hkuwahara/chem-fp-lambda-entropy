#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE);
if(length(args) < 2) {
  stop("requires fp.scheme, query-id");
}
analysis.dir <- ".";
fp.scheme <- args[1];
query.drug.id <- args[2];

# analysis.dir <- "/home/hiro/research/kaust/chem-fp/revision/chem-fp-lambda-entropy/drug-bank-data-analysis";
# fp.scheme <- "maccs";
# query.drug.id <- "DB00760";

used.fp.file <- paste( analysis.dir, fp.scheme, paste0("drug-reduced-0.3-", fp.scheme, "-fp.list"), sep="/");
used.fp.list <- read.delim( used.fp.file, header = F );
selected.fp.list <- as.character(unlist(used.fp.list));

candidate.file <- paste( analysis.dir, paste0("approved-drug-", fp.scheme, "-candidate.tsv"), sep="/");
# dat <- read.delim( input.file, header = F, col.names = c("hmdb_id", paste0("F", 1:166)), row.names = 1);
candidate.dat <- read.delim( candidate.file, header = F, row.names = 1);

fp.size <- dim(candidate.dat)[2];
colnames(candidate.dat) <- paste0( "F", 1:fp.size );

query.file <- paste( analysis.dir, paste0("approved-drug-", fp.scheme, "-query.tsv"), sep="/");
query.dat <- read.delim( query.file, header = F, row.names = 1);
query.fp <- query.dat[query.drug.id,];
colnames(query.fp) <- paste0( "F", 1:fp.size );


gen.get.tanimoto.sim <- function( x ) {
  function(y) {
    intsxn <- sum( x == 1 & y == 1 );
    intsxn / ( sum( x == 1) + sum( y == 1) - intsxn + 1e-10);
  }
}

gen.get.contrib.numerator <- function(x, selected.fp.list) {
  selected.x <- x[selected.fp.list];
  function(y) {
    selected.y <- y[selected.fp.list];
    intsxn <- sum( x == 1 & y == 1 ) + 1e-10;
    selected.intsxn <- sum( selected.x == 1 & selected.y == 1 );
    selected.intsxn / intsxn;
  }
}

gen.get.contrib.denominator <- function(x, selected.fp.list) {
  selected.x <- x[selected.fp.list];
  function(y) {
    selected.y <- y[selected.fp.list];
    union.sum <- sum( x == 1 | y == 1 ) + 1e-10;
    selected.union.sum <- sum( selected.x == 1 | selected.y == 1 );
    selected.union.sum / union.sum;
  }
}



get.tanimoto.sim <- gen.get.tanimoto.sim( query.fp );
get.contrib.numerator <- gen.get.contrib.numerator( query.fp, selected.fp.list );
get.contrib.denominator <- gen.get.contrib.denominator( query.fp, selected.fp.list );


sim.results <- apply( candidate.dat, 1, get.tanimoto.sim);
numerator.results <- apply( candidate.dat, 1, get.contrib.numerator);
denominator.results <- apply( candidate.dat, 1, get.contrib.denominator);

rec <- data.frame( query=rep(query.drug.id, dim(candidate.dat)[1] ), 
                   candidate=rownames(candidate.dat), 
                   tanimoto.sim=sim.results,
                   related.contrib.numerator=ifelse(numerator.results < 1e-10,  0, (1.0 - numerator.results)),
                   selected.contrib.numerator=numerator.results,
                   related.contrib.denominator=ifelse(denominator.results < 1e-10,  0, (1.0 - denominator.results)),
                   selected.contrib.denominator=denominator.results
                   );

out.file <- paste( analysis.dir, fp.scheme, paste0(fp.scheme, "-", query.drug.id, "-results.tsv"), sep="/");
write.table(rec, out.file, col.names = T, row.names = F, sep="\t", quote=F);

