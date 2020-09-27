#!/usr/bin/env Rscript

library("dplyr");
library("data.table");
library("reshape2");
library("ggplot2");
library("scales");
library("ggsci");
library("gridExtra");
library("ggpubr");

library("RColorBrewer");
n <- 10;
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',];
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)));
                    

args = commandArgs(trailingOnly=TRUE);
# if(length(args) < 3) {
#   stop("requires analysis.dir, fp.scheme, query-id");
# }
# analysis.dir <- args[1];
# fp.scheme <- args[2];
# query.drug.id <- args[3];

# plot.colors <- pal_npg()(10);
plot.colors <- col_vector[1:10];
show_col(plot.colors);

analysis.dir <- ".";
fp.scheme <- "maccs";

query.drug.file <- paste( analysis.dir, "query_id.txt", sep="/" );
query.ids <- as.character(read.delim( query.drug.file, header = F )$V1);

num.entries <- 50;

gen.rec <- function(fp.scheme) {
  query.id <- query.ids[1];
  input.file <- paste( analysis.dir, fp.scheme, paste0(fp.scheme, "-", query.id, "-results.tsv" ), sep="/" );
  dat <- read.delim( input.file, header = T );
  dat <- dat[order(dat$tanimoto.sim, decreasing = T), ];
  top.dat <- dat[1:num.entries, ];
  rec <- data.frame( query=top.dat$query, 
                     contrib.num=top.dat$related.contrib.numerator, 
                     contrib.den=top.dat$related.contrib.denominator, 
                     rank=1:num.entries, 
                     col=plot.colors[1] );
  
  for( i in 2:length(query.ids) ) {
    query.id <- query.ids[i];
    input.file <- paste( analysis.dir, fp.scheme, paste0(fp.scheme, "-", query.id, "-results.tsv" ), sep="/" );
    dat <- read.delim( input.file, header = T );
    dat <- dat[order(dat$tanimoto.sim, decreasing = T), ];
    top.dat <- dat[1:num.entries, ];
    rec <- rbind( rec, data.frame( query=top.dat$query, 
                          contrib.num=top.dat$related.contrib.numerator, 
                          contrib.den=top.dat$related.contrib.denominator, 
                          rank=1:num.entries, 
                          col=plot.colors[i] )
                  );
    
  }
  rec;
}

maccs.rec <- gen.rec("maccs");
pubchem.rec <- gen.rec("pubchem");

rec <- maccs.rec;
gen.plot <- function( rec, lower.lim=0.0, upper.lim=0.26 ) {
  plt <- ggplot(rec, aes(x = contrib.den, y = contrib.num));
  plt <- plt + geom_abline(intercept = 0, slope = 1, linetype = "dashed");
  plt <- plt + geom_point(aes(color=query), size=2.2, alpha=0.7, shape=21, stroke=1.4);
  plt <- plt + scale_color_manual(breaks=unique(rec$query), values=plot.colors, name="query drug");
  plt <- plt + labs(x="contribution to the union", y="contribution to the intersection");
  plt <- plt + xlim(lower.lim,upper.lim);
  plt <- plt + ylim(lower.lim,upper.lim);
  plt <- plt + coord_fixed();
  plt <- plt + theme_bw();
  plt <- plt + theme(
    axis.text.x = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  );
}

maccs.plt <- gen.plot( maccs.rec, 0.0, 0.27 );
pubchem.plt <- gen.plot( pubchem.rec, 0.0, 0.27 );

print(pubchem.plt);

plt <- ggarrange(maccs.plt, 
                 pubchem.plt, 
                 ncol = 2, nrow = 1, labels="AUTO");
print(plt);
ggsave(filename=paste(analysis.dir, "related-fp-contrib-scatterplot.svg", sep="/"), device="svg", plot=plt, width=28, height= 14, units="cm");


