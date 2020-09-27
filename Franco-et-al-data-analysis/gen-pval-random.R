#!/usr/bin/env Rscript

library("dplyr");
library("data.table");
library("reshape2");
library("ggplot2");
library("scales");
library("ggsci");
library("gridExtra");
library("ggpubr");

fp.type <- "maccs";
num.keys <- 132;

toplevel.dir <- ".";

testset.file <- paste(toplevel.dir, "sim_testset.tsv", sep="/");
testset.df <- read.delim( testset.file, header = T);
selected.pairs <- testset.df[testset.df[,4] >= 80,1];

selected.pairs <- as.character(selected.pairs);

gen.results.df <- function( selected.pairs, fp.type, num.keys ) {
  sim.results.file <- paste0(toplevel.dir, "/",  "sim-", fp.type, ".tsv");
  sim.results.df <- read.delim( sim.results.file, header = T);
  colnames(sim.results.df) <- c( "id", paste0("level_", c("0", "0.3") ) );
  
  selected.sim.results.df <- subset(sim.results.df, subset = id %in% selected.pairs );
  ref.results <- selected.sim.results.df$level_0;
  
  rel.change.sim <- (selected.sim.results.df$level_0.3 - ref.results)/ ref.results;
  
  random.results.file <- paste0(toplevel.dir, "/",  "random-results-", fp.type, "-keys-", num.keys, "-t.tsv");
  
  random.results.df <- read.delim( random.results.file, header = T);   
  colnames(random.results.df) <- c(id, paste0(1:100, "a"));
  selected.random.results.df <- subset(random.results.df, select = selected.pairs);
  rel.change.random.sim <- apply( selected.random.results.df, 1, function(v){ (v-ref.results)/ref.results;} );
  
  valid.pairs <- sapply( 1:length(rel.change.sim), function(i){
      s <- rel.change.sim[i]; 
      r <- quantile( rel.change.random.sim[i,], c(0.25, 0.75) );
      (max(abs(rel.change.random.sim[i,])) > 1e-6) && (s < r[1] || s > r[2]);
    } 
    );
  
  
  # valid.pairs <- apply( rel.change.random.sim, 1, function(v){ max(abs(v)) > 1e-6; } );
  # valid.pairs <- rel.change.sim > 1e-6;
  
  dist.df <- cbind(rel.change.sim[valid.pairs], rel.change.random.sim[valid.pairs,]);
  p.vals <- apply( dist.df, 1, function(v){ ifelse( v[1] > 0, sum(v[-1] >= v[1])/length(v[-1]),  sum(v[-1] <= v[1])/length(v[-1])); } );
  adj.p <- p.adjust( p.vals, "fdr" );

  rec.df <- data.frame( id=selected.sim.results.df$id[valid.pairs],
                        rel.change.sim=rel.change.sim[valid.pairs],
                        p.val=p.vals,
                        adj.p=adj.p);
  rec.df;  
}

maccs.rec.df <- gen.results.df(selected.pairs, "maccs", 132); 
pubchem.rec.df <- gen.results.df(selected.pairs, "pubchem", 411); 

#
# in MACCS
# id rel.change.sim  p.val  adj.p
# 33a     0.07210031 0.0331 0.1988182
# 60a     0.07236269 0.0269 0.1988182

# in pubchem
# id rel.change.sim  p.val  adj.p
# 44a    0.051707505 0.0006 0.0048
# 60a    0.061965080 0.0012 0.0064
# 66a    0.008676236 0.0193 0.0772
# 77a    0.085128205 0.0001 0.0016

