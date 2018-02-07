#!/usr/bin/env Rscript

# dependencies
errors <- 0
if( !(library(phangorn, verbose = F, warn.conflicts	= F, quietly = T, logical.return = T)) ){
  sprintf(" - ERROR: phangorn not installed")
  errors <- errors + 1 
}
if( !(library(dplyr, verbose = F, warn.conflicts	= F, quietly = T, logical.return = T)) ){
  sprintf(" - ERROR: dplyr not installed")
  errors <- errors + 1 
}
if( !(library(treeWAS, verbose = F, warn.conflicts	= F, quietly = T, logical.return = T)) ){
  sprintf(" - ERROR: treeWAS not installed")
  errors <- errors + 1 
}
if( !(library(data.table, verbose = F, warn.conflicts	= F, quietly = T, logical.return = T)) ){
  sprintf(" - ERROR: data.table not installed")
  errors <- errors + 1 
}

# check dependencies are installed
if( errors > 0 ){
  stop("dependencies missing -  see above")
}

# data files
args = commandArgs(trailingOnly=TRUE)

# check number of command line arguements
no_args <- length(args)
if ( no_args<4  ){
  stre <- sprintf( "not enough commands supplied script.")
  str2 <- sprintf( " run_treeWAS.R [snps/variants] [metadata] [tree] [output directory] [optional args]" )
  str3 <- sprintf( "    --r    rank phenotype data [default:off]")
  
  text <- writeLines( c("", str2, "", " Options:", str3))
  stop( stre )
}

# check for optional arguements
ranked <- 0
if(no_args>4){
  for(i in 5:no_args){
    if(args[i] == "--r"){
      ranked <- 1
    }
  }
}

# set input output files
snps_file <- args[1]
meta_file <- args[2]
tree_file <- args[3]
output_dir <-args[4]

# check output directory exists - make if possible
dir.create(output_dir, showWarnings = F )
if (!dir.exists(output_dir)){
  stop("output directory does not exist or could not be created")
}

# load data
raw_metadata <- read.delim (meta_file, sep="\t")
tree_raw <- midpoint(read.tree(tree_file))

# fread 
snps_org <- fread(snps_file, sep="\t", header=T)
no_sites <- length(snps_org[1,])

# convert to matrix
snps_raw <- as.matrix(snps_org[,2:no_sites], header=TRUE)
rownames(snps_raw) <- snps_org$Samples

# alt - conversion
#snps_raw <- data.matrix(snps_org[,2:no_sites])
#rownames(snps_raw) <- snps_org$Samples

# get phenotypes
no_cols <- length(colnames(raw_metadata))
myCols <- colnames(raw_metadata)[2:no_cols]

for (field in myCols){
  
  writeLines(c(sprintf(" - processing: %s",field)))

  # subset metadata for treewas
  metadata <- raw_metadata %>%
        select( id, match(field,names(.) ) )
  colnames(metadata) <- c("id", "pheno")
  metadata <- metadata[complete.cases(metadata$pheno),]  # only complete cases
  phen <- as.factor(metadata$pheno)
  names(phen) <-  metadata$id
  
  # remove samples not in metadata, tips and tree
  tips <- tree_raw$tip.label 
  meta_ids <- metadata$id
  snp_ids <- rownames(snps_raw)
  
  # intersection between all lists
  include <- intersect(intersect(meta_ids,snp_ids), tips)
  
  # drop any tips not in metadata
  drop.tips <- tree_raw$tip.label[!(tree_raw$tip.label %in% include)]
  tree <- midpoint(drop.tip(tree_raw, drop.tips))
  
  # drop metadata not in tips 
  metadata <- metadata[metadata$id %in% include,] 
  
  # remove any snps that are not in phenotype 
  snps <- snps_raw[ rownames(snps_raw) %in% include, ] 
  rownames(snps) <- rownames(snps_raw)[rownames(snps_raw) %in% include]
  
  # make phenotype factor
  phen <- as.vector(metadata$pheno)
  names(phen) <-  metadata$id
  
  # make snps binary and filter multiple reciprocal freq snps 
  suffixes <- keepLastN(colnames(snps), n = 2)
  suffixes <- unique(suffixes)
  
  # if input is snp data then remove extraneous data
  if(all(suffixes %in% c(".a", ".c", ".g", ".t"))){
    snps <- get.binary.snps(snps)
  }
  
  # Check that samples remain to test
  no_samples <- length(metadata[,1])
  if ( no_samples == 0 ){
    (sprintf (" - WARNING: 0 samples intersect between variants/metadata/tree"))
  }else{
    
    # feedback
    writeLines(c(sprintf (" - %s samples intersect between variants/metadata/tree", no_samples)))
  
    # Examine distribution of ranks:
    phen.rank <- rank(phen, ties.method = "average")
    
    # Check distribution of phenotypes as histograms
    out_hist_name = sprintf("%s/%s.histograms.pdf", output_dir, field)
    pdf(out_hist_name)
    hist(phen.rank)
    hist(phen)
    dev.off()
    
    # [optional] replace phenotype with ranks if selected.
    phen.in <- phen
    if( ranked  == 1 ){
      phen.in <- phen.rank
    }
    
    # ensure the tree does not have negative branch lengths and is binary
    #tree$edge.length[tree$edge.length < 0] <- abs(min(tree$edge.length))
    
    # check if tree is binary 
    if ( !is.binary.tree(tree) ){
      (sprintf( " - WARNING: the tree supplied was not bifurcating. It has been converted but the input tree should be checked." ))
      tree <- multi2di(tree)
    }
    
    # output filename
    out_plot_name <- sprintf("%s/%s.treeWAS_plots.pdf", output_dir, field)
    
    # run treewas
    out <- treeWAS(snps = snps,
                   phen = phen.in,
                   tree = tree,
                   n.subs = NULL,
                   n.snps.sim = ncol(snps)*10,
                   chunk.size = ncol(snps),
                   test = c("terminal", "simultaneous", "subsequent"),
                   snps.reconstruction = "parsimony",
                   snps.sim.reconstruction = "parsimony",
                   phen.reconstruction = "parsimony",
                   phen.type = NULL,
                   na.rm = TRUE,
                   p.value = 0.01,
                   p.value.correct = "bonf",
                   p.value.by = "count",
                   dist.dna.model = "JC69", ### ignored as tree supplied
                   plot.tree = TRUE,
                   plot.manhattan = TRUE,
                   plot.null.dist = TRUE,
                   plot.dist = FALSE,
                   snps.assoc = NULL, 
                   filename.plot = out_plot_name,
                   seed = 1)
    
    # Print treewas to files as tables 
    min_p <- out$terminal$min.p.value
    
    if ( !is.character(out$terminal$sig.snps) ){
      
      # replace 0 with min_p
      out$terminal$sig.snps$p.value[out$terminal$sig.snps$p.value == 0] <- sprintf("<%s",min_p)
      
      out_term_name = sprintf("%s/%s.terminal.significant.tab", output_dir, field)
      out_term <- out$terminal$sig.snps[order(out$terminal$sig.snps$p.value), ]
      write.table(out_term, file = out_term_name, append = F, quote = F, sep = "\t",
                  eol = "\n", na = "NA", dec = ".", row.names = TRUE,
                  col.names = TRUE)
    }
    
    min_p <- out$simultaneous$min.p.value
    
    if ( !is.character(out$simultaneous$sig.snps) ){
      
      # replace 0 with min_p
      out$simultaneous$sig.snps$p.value[out$simultaneous$sig.snps$p.value == 0] <- sprintf("<%s",min_p)
      
      out_simul_name = sprintf("%s/%s.simultaneous.significant.tab", output_dir, field)
      out_simul <- out$simultaneous$sig.snps[order(out$simultaneous$sig.snps$p.value), ]
      write.table(out_simul, file = out_simul_name, append = F, quote = F, sep = "\t",
                  eol = "\n", na = "NA", dec = ".", row.names = TRUE,
                  col.names = TRUE)
    }
    
    min_p <- out$subsequent$min.p.value
    
    if ( !is.character(out$subsequent$sig.snps) ){
      
      # replace 0 with min_p
      out$subsequent$sig.snps$p.value[out$subsequent$sig.snps$p.value == 0] <- sprintf("<%s",min_p)
      
      out_subs_name = sprintf("%s/%s.subsequent.significant.tab", output_dir, field)
      out_subs <- out$subsequent$sig.snps[order(out$subsequent$sig.snps$p.value), ]
      write.table(out_subs, file = out_subs_name, append = F, quote = F, sep = "\t",
                eol = "\n", na = "NA", dec = ".", row.names = TRUE,
                col.names = TRUE)
    }
    
    # make data.frame of all correlations
    
    # terminal
    out_all <- sprintf("%s/%s.terminal.all_loci.tab", output_dir, field)
    header <- sprintf("# significance threshold: %f, minimum p-value: %f",  out$terminal$sig.thresh,  out$terminal$min.p.value)
    writeLines(header, out_all)
    
    df_out <- data.frame( loci = as.character(names(out$terminal$corr.dat)),
            correlation = out$terminal$corr.dat,
            pvalue = out$terminal$p.vals)
    
    # replace zero by min p-value and sort
    df_out$pvalue[df_out$pvalue == 0] <- out$terminal$min.p.value
    df_out <- df_out[order(df_out$pvalue),] 
    
    write.table(df_out, file = out_all, append = T, quote = F, sep = "\t",
                eol = "\n", na = "NA", dec = ".", row.names = F,
                col.names = TRUE)
    
    # simultaneous
    out_all <- sprintf("%s/%s.simultaneous.all_loci.tab", output_dir, field)
    header <- sprintf("# significance threshold: %f, minimum p-value: %f",  out$simultaneous$sig.thresh,  out$simultaneous$min.p.value)
    writeLines(header, out_all)
    
    df_out <- data.frame( loci = as.character(names(out$simultaneous$corr.dat)),
                          correlation = out$simultaneous$corr.dat,
                          pvalue = out$simultaneous$p.vals)
    
    # replace zero by min p-value and sort
    df_out$pvalue[df_out$pvalue == 0] <- out$simultaneous$min.p.value
    df_out <- df_out[order(df_out$pvalue),] 
    
    write.table(df_out, file = out_all, append = T, quote = F, sep = "\t",
                eol = "\n", na = "NA", dec = ".", row.names = F,
                col.names = TRUE)
    
    # subsequent
    out_all <- sprintf("%s/%s.subsequent.all_loci.tab", output_dir, field)
    header <- sprintf("# significance threshold: %f, minimum p-value: %f",  out$subsequent$sig.thresh,  out$subsequent$min.p.value)
    writeLines(header, out_all)
    
    df_out <- data.frame( loci = as.character(names(out$subsequent$corr.dat)),
                          correlation = out$subsequent$corr.dat,
                          pvalue = out$subsequent$p.vals)
    
    # replace zero and sort
    df_out$pvalue[df_out$pvalue == 0] <- out$subsequent$min.p.value
    df_out <- df_out[order(df_out$pvalue),] 
    
    write.table(df_out, file = out_all, append = T, quote = F, sep = "\t",
                eol = "\n", na = "NA", dec = ".", row.names = F,
                col.names = TRUE)
    
    # store treewas variable as R object
    out_rds <- sprintf("%s/%s.treewas.rds", output_dir, field)
    saveRDS(out, out_rds)
    
  }

}
