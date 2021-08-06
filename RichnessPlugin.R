library(microbiome)
library(ggplot2)
#library(phyloseq)
library(ape)
library(psadd)

dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")


input <- function(inputfile) {
  pfix = prefix()
  parameters <<- read.table(inputfile, as.is=T);
  rownames(parameters) <<- parameters[,1]; 
   # Need to get the three files
   otu.path <<- paste(pfix, parameters["otufile", 2], sep="/")
   tree.path <<- paste(pfix, parameters["tree", 2], sep="/")
   map.path <<- paste(pfix, parameters["mapping", 2], sep="/")
   column <<- parameters["column", 2]
   #HMP <<- import_qiime(otu.path, map.path, tree.path, parseFunction = parse_taxonomy_qiime)
}
run <- function() {
   #samples.to.keep <<- sample_sums(HMP) >= 1000
   #HMP <<- prune_samples(samples.to.keep, HMP)
   #HMP <<- filter_taxa(HMP, function(x) sum(x >3) > (0.01*length(x)), TRUE)
   physeq <<- read_csv2phyloseq(otu.file=otu.path, taxonomy.file=tree.path, metadata.file=map.path)
   mytree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
   physeq <<- merge_phyloseq(physeq, mytree)
}
output <- function(outputfile) {
  pdf(paste(outputfile,"pdf",sep="."))#,  width = 10*300,        # 5 x 300 pixels
  #height = 10*300); #,)
  print("Generating plot...")
  #result <<- PCoA(physeq)
  y <- plot_richness(physeq, color=column)
  #y <- plot_sparsity(p0)
  #print(str(y))
  print("Generating CSV...")
  rich <- richness(physeq)
  #print(str(y$data))
  write.csv(rich, paste(outputfile,"csv",sep="."))
  print(y)#plot_bar(HMP, x="Description", fill=diffcol))
  dev.off()
}
#input("plugins/Bar/example/parameters.txt")
#run()
#output("plugins/Bar/example/yes.pdf")

