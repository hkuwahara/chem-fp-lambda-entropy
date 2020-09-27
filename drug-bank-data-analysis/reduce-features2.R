#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE);
if(length(args) < 1) {
  stop("requires fp.scheme");
}
analysis.dir <- ".";
fp.scheme <- args[1];

results.dir <- paste( analysis.dir, fp.scheme, sep="/");
if( !dir.exists(results.dir )) {
  dir.create( results.dir, recursive = T);
}

input.file <- paste( analysis.dir, paste0("approved-drug-", fp.scheme, ".tsv"), sep="/");
original.dat <- read.delim( input.file, header = F, row.names = 1);

entry.size <- dim(original.dat)[1];

original.num.features <- length(original.dat[1,]);
colnames(original.dat) <-  paste0("F", 1:original.num.features);

num.ones <- sapply( 1:original.num.features, function(i, dat){sum(dat[,i])}, original.dat );
valid.features <- num.ones > 0 & num.ones < length(original.dat[,1]);
dat <- original.dat[, valid.features];
num.ones <- num.ones[valid.features];
original.rank <- qr(dat)$rank;

get.svd.entropy <- function(M) {
  ei.val <- svd(M, 0, 0)$d^2;
  normed.ei.val <- ei.val/sum(ei.val);
  -sum(normed.ei.val * log10(normed.ei.val + 1e-10))/log10(length(M[1,]));
}

ref.value <- get.svd.entropy(original.dat);
feature.values <- sapply( 1:length(original.dat[1,]), function(i, dat){ t <- dat; t[,i] <- 0; get.svd.entropy(t) }, original.dat );
score <- ref.value - feature.values;

out.file <- paste( results.dir, paste0("svd-entropy-", fp.scheme, ".tsv"), sep="/");
write.table( data.frame( id=names(original.dat), score=score, entropy=feature.values ), out.file, col.names = T, row.names = F, sep="\t", quote=F);

disp <- sqrt(mean(score * score));
mean.score <- mean(score);
sd.score <- sd(score);
rec <- data.frame(scale=c(0.0), rank=c(original.rank), feature.count=c(length(original.dat[1,])));
# scale.vals <- c(0.25, 0.5, 0.75, 1.0, 1.5, 2.0);
# scale.vals <- seq(0.1, 0.3, 0.1)
scale.vals <- c(0.3);
for( scale.val in scale.vals ) {
  selected.indices <- (feature.values >= (ref.value + scale.val * disp) | feature.values <= (ref.value - scale.val * disp)); 
  dat <- original.dat[, selected.indices];
  dat.rank <- sum( svd(dat)$d >= 1e-5);
  rec <- rbind(rec, data.frame(scale=c(scale.val), rank=c(dat.rank), feature.count=c(length(dat[1,]))));
  
  out.file <- paste( results.dir, paste0("drug-reduced-", scale.val, "-", fp.scheme, ".tsv"), sep="/");
  write.table(dat, out.file, col.names = T, row.names = rownames(original.dat), sep="\t", quote=F);
}

out.file <- paste( results.dir, paste0("ce-stat-", fp.scheme, ".tsv"), sep="/");
write.table(rec, out.file, col.names = T, row.names = F, sep="\t", quote=F);


# rank.results <- sapply( 1:length(dat[1,]), function(i, dat){qr(dat[,-i])$rank;}, dat );


