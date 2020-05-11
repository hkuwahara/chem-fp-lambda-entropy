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
sum(subset.df$maccs.0 >= maccs.cutoff.val );
sum(subset.df$maccs.0.1 >= maccs.cutoff.val );
sum(subset.df$maccs.0.2 >= maccs.cutoff.val );
sum(subset.df$maccs.0.3 >= maccs.cutoff.val );
sum(subset.df$maccs.0 < maccs.cutoff.low );
sum(subset.df$maccs.0.1 < maccs.cutoff.low );
sum(subset.df$maccs.0.2 < maccs.cutoff.low );
sum(subset.df$maccs.0.3 < maccs.cutoff.low );
mean(subset.df$maccs.0);
mean(subset.df$maccs.0.1);
mean(subset.df$maccs.0.2);
mean(subset.df$maccs.0.3);

sum(subset.df$pubchem.0 >= pubchem.cutoff.val );
sum(subset.df$pubchem.0.1 >= pubchem.cutoff.val );
sum(subset.df$pubchem.0.2 >= pubchem.cutoff.val );
sum(subset.df$pubchem.0.3 >= pubchem.cutoff.val );
sum(subset.df$pubchem.0 < pubchem.cutoff.low );
sum(subset.df$pubchem.0.1 < pubchem.cutoff.low );
sum(subset.df$pubchem.0.2 < pubchem.cutoff.low );
sum(subset.df$pubchem.0.3 < pubchem.cutoff.low );
mean(subset.df$pubchem.0);
mean(subset.df$pubchem.0.1);
mean(subset.df$pubchem.0.2);
mean(subset.df$pubchem.0.3);

total.mean <- sapply( 3:10, function(i){mean(df[,i]);} );

summary.df <- data.frame(
  fp=rep( c("MACCS", "Pubchem" ), each=4 ),
  level=rep( c("0.0", "0.1", "0.2", "0.3"), 2 ),
  mean=sapply( 3:10, function(i){mean(subset.df[,i]);} ),
  ge.80=sapply( 3:10, function(i){sum(subset.df[,i] >= 0.8);}),
  lt.50=sapply( 3:10, function(i){sum(subset.df[,i] < 0.5);})
);
summary.df;
