#!/usr/bin/env Rscript

results.file <- "./sim-results.tsv";
df <- read.delim( results.file, header = T);

ref.cutoff.val <- 0.8;
maccs.cutoff.val <- 0.8;
maccs.cutoff.low <- 0.5;
pubchem.cutoff.val <- 0.8;
pubchem.cutoff.low <- 0.5;

subset.df <- df[df$ref >= ref.cutoff.val, ];
dim(subset.df);
sum(subset.df$maccs_0 >= maccs.cutoff.val );
sum(subset.df$maccs_0.3 >= maccs.cutoff.val );
sum(subset.df$maccs_0 < maccs.cutoff.low );
sum(subset.df$maccs_0.3 < maccs.cutoff.low );
mean(subset.df$maccs_0);
mean(subset.df$maccs_0.3);

sum(subset.df$pubchem_0 >= pubchem.cutoff.val );
sum(subset.df$pubchem_0.3 >= pubchem.cutoff.val );
sum(subset.df$pubchem_0 < pubchem.cutoff.low );
sum(subset.df$pubchem_0.3 < pubchem.cutoff.low );
mean(subset.df$pubchem_0);
mean(subset.df$pubchem_0.3);

# total.mean <- sapply( 3:10, function(i){mean(df[,i]);} );

summary.df <- data.frame(
  fp=rep( c("MACCS", "Pubchem" ), each=2 ),
  level=rep( c("0.0", "0.3"), 2 ),
  mean=sapply( 3:6, function(i){mean(subset.df[,i]);} ),
  sd=apply( subset.df[,c(-1,-2),], 2, sd ),
  ge.80=sapply( 3:6, function(i){sum(subset.df[,i] >= 0.8);}),
  ge.75=sapply( 3:6, function(i){sum(subset.df[,i] >= 0.75);}),
  ge.70=sapply( 3:6, function(i){sum(subset.df[,i] >= 0.70);}),
  # lt.50=sapply( 3:6, function(i){sum(subset.df[,i] < 0.5);}),
  min.sim=sapply( 3:6, function(i){min(subset.df[,i]);})
);
summary.df;
