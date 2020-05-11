#!/usr/bin/env Rscript

library("dplyr");
library("data.table");
analysis.home <- ".";
outdir <- analysis.home;
reduce.levels <- c("0", "0.1", "0.2", "0.3");
compounds.classes <- c("drug", "microbial", "plant", "endogenous");
sample.type <- "blood";

rec.df <- data.frame( scheme=c(), level=c(), drug=c(), microbial=c(), plant=c(), endogenous=c(), all=c() );
fp.schemes <- c("maccs", "pubchem");
all.compound.list.file <- paste( analysis.home,  paste0( "filtered-", sample.type, "-compounds.list"), sep="/"); 
all.compound.df <- fread(file=all.compound.list.file, header=F, data.table=F, col.names=c("id"));

for( fp.scheme in fp.schemes ) {
  sim.dir <- paste0(analysis.home, "/", sample.type, "-", fp.scheme);
  for( reduce.level in reduce.levels ) {
    t <- c();
    t["scheme"] <- fp.scheme;
    t["level"] <- reduce.level;
    sim.file <- paste(sim.dir, paste0("sim-ce-reduced-", reduce.level, "-", sample.type, "-", fp.scheme, "-melted.tsv"), sep="/" );
    sim.df <- subset( fread(file=sim.file, header=T, data.table=F), (a != b) );
    sim.df <- subset( sim.df, (a %in% all.compound.df$id) & (b %in% all.compound.df$id) );
    sim.df$sim[is.na(sim.df$sim)] <- 0;
    t["all"] <- mean(sim.df$sim);
    for( comp.class in compounds.classes ) {
      class.compound.list.file <- paste( analysis.home,  paste0( sample.type, "-", comp.class, "-compounds.list"), sep="/"); 
      class.compound.df <- fread(file=class.compound.list.file, header=F, data.table=F, col.names=c("id"));
      class.sim.df <- subset( sim.df, (a %in% class.compound.df$id) & (b %in% class.compound.df$id) );
      t[comp.class] <- mean(class.sim.df$sim);
    }
    rec.df <- rbind( rec.df, data.frame( scheme=t["scheme"], level=t["level"], drug=t["drug"], 
                                         microbial=t["microbial"], plant=t["plant"], toxin=t["toxin"], 
                                         cosmetic=t["cosmetic"], endogenous=t["endogenous"], all=t["all"] ), make.row.names = F);
  }
}
out.file <- paste( outdir, "subclass-sim-results.tsv", sep="/");
write.table(rec.df, out.file, col.names = T, row.names = F, sep="\t", quote=F);
