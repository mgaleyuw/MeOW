plotRibbonGenes <- function(cpg.df, testsample, exclusions="", minmax=TRUE,nrow=NA,cpgs="NA",genes="NA") {
  if(is.null(match.call()$cpgs) && is.null(match.call()$genes)){
    print("Supply either a list of genes or a list of CPG island names")
    return()
  } else if(is.null(match.call()$cpgs)){
    plot.df <- cpg.df %>% filter(gene %in% genes) %>% dplyr::select(-any_of(exclusions)) %>% reshape2::melt(id.vars=c("pos", "chr", "cpg", "gene"))
  } else {
    plot.df <- cpg.df %>% filter(cpg %in% cpgs) %>% dplyr::select(-any_of(exclusions)) %>% reshape2::melt(id.vars=c("pos", "chr", "cpg", "gene"))
  }
  colnames(plot.df) <- c("Position", "Chromosome", "CpgIsland", "Gene", "Sample", "Methylation")
  plot.df$PlotGene <- paste(plot.df$Chromosome, plot.df$Gene, sep=":")
  ctrl.avgs <- plot.df %>% filter(Sample!=testsample) %>% group_by(Chromosome, PlotGene, CpgIsland, Position) %>% dplyr::summarize(mnMethylation=mean(Methylation), maxMethylation=max(Methylation), minMethylation=min(Methylation), sdMethylation=sd(Methylation), n=n()) %>% mutate(CI=1.95*sdMethylation/(sqrt(n)))
  
  if(minmax){
    g <- ggplot() + geom_point(data=plot.df %>% filter(Sample==testsample), aes(x=Position, y=Methylation), color='red', shape=4, alpha=0.75) + geom_ribbon(data=ctrl.avgs, aes(x=Position, ymin=minMethylation, ymax=maxMethylation), fill="gray", alpha=0.5) + ylab("Mean Methylation Frequency") + scale_x_continuous(label=comma) + theme(axis.text.x = element_text(angle=45, hjust=1)) + xlab("Position within CpG Island")
  } else {
    g <- ggplot() + geom_point(data=plot.df %>% filter(Sample==testsample), aes(x=Position, y=Methylation), color='red', shape=4, alpha=0.75) + geom_ribbon(data=ctrl.avgs, aes(x=Position, ymin=mnMethylation - CI, ymax=mnMethylation + CI), fill="gray", alpha=0.5) + ylab("Mean Methylation Frequency") + scale_x_continuous(label=comma) + theme(axis.text.x = element_text(angle=45, hjust=1)) + xlab("Position within CpG Island")
  }
  if(is.na(nrow)){
    return(g+facet_wrap(~ PlotGene + CpgIsland, scales="free_x"))
  } else {
    return(g+facet_wrap(~ PlotGene + CpgIsland, scales="free_x", nrow=nrow))
  }
  return(g)
}

plotRibbonRegion <- function(cpg.df, region, testsample, exclusions="", minmax=TRUE,nrow=NA){
  chr_in=strsplit(region,":")[[1]][1]
  coords=strsplit(strsplit(region,":")[[1]][2],"-")[[1]]
  start=as.numeric(coords[1])
  stop=as.numeric(coords[2])
  print(c(chr_in, start, stop))
  plot.df <- cpg.df %>% filter(chr==chr_in) %>% filter(between(pos,as.numeric(start),as.numeric(stop))) %>% dplyr::select(-any_of(exclusions)) %>% reshape2::melt(id.vars=c("pos", "chr", "cpg", "gene"))
 
  colnames(plot.df) <- c("Position", "Chromosome", "CpgIsland", "Gene", "Sample", "Methylation")
  plot.df$PlotGene <- paste(plot.df$Chromosome, plot.df$Gene, sep=":")
  
  ctrl.avgs <- plot.df %>% filter(Sample!=testsample) %>% group_by(Chromosome, PlotGene, CpgIsland, Position) %>% dplyr::summarize(mnMethylation=mean(Methylation), maxMethylation=max(Methylation), minMethylation=min(Methylation), sdMethylation=sd(Methylation), n=n()) %>% mutate(CI=1.95*sdMethylation/(sqrt(n)))
  
  if(minmax){
    g <- ggplot() + geom_point(data=plot.df %>% filter(Sample==testsample), aes(x=Position, y=Methylation), color='red', shape=4, alpha=0.75) + geom_ribbon(data=ctrl.avgs, aes(x=Position, ymin=minMethylation, ymax=maxMethylation), fill="gray", alpha=0.5) + ylab("Mean Methylation Frequency") + scale_x_continuous(label=comma) + theme(axis.text.x = element_text(angle=45, hjust=1)) + xlab("Position within CpG Island")
  } else {
    g <- ggplot() + geom_point(data=plot.df %>% filter(Sample==testsample), aes(x=Position, y=Methylation), color='red', shape=4, alpha=0.75) + geom_ribbon(data=ctrl.avgs, aes(x=Position, ymin=mnMethylation - CI, ymax=mnMethylation + CI), fill="gray", alpha=0.5) + ylab("Mean Methylation Frequency") + scale_x_continuous(label=comma) + theme(axis.text.x = element_text(angle=45, hjust=1)) + xlab("Position within CpG Island")
  }
  if(is.na(nrow)){
    return(g+facet_wrap(~ PlotGene + CpgIsland, scales="free_x"))
  } else {
    return(g+facet_wrap(~ PlotGene + CpgIsland, scales="free_x", nrow=nrow))
  }
  return(g)
}

compareRibbons <- function(cpg.df, genes, group1, group2, exclusions="", minmax=TRUE, pcloud2=FALSE){
  plot.df <- cpg.df %>% filter(gene %in% genes) %>% dplyr::select(-any_of(exclusions)) %>% reshape2::melt(id.vars=c("pos", "chr", "cpg", "gene"))
  colnames(plot.df) <- c("Position", "Chromosome", "CpgIsland", "Gene", "Sample", "Methylation")
  
  grp1.avgs <- plot.df %>% filter(Sample %in% group1) %>% group_by(Chromosome, Gene, CpgIsland, Position) %>% dplyr::summarize(mnMethylation=mean(Methylation), maxMethylation=max(Methylation), minMethylation=min(Methylation), sdMethylation=sd(Methylation), n=n()) %>% mutate(CI=1.95*sdMethylation/(sqrt(n)))
  
  grp2.avgs <- plot.df %>% filter(Sample %in% group2) %>% group_by(Chromosome, Gene, CpgIsland, Position) %>% dplyr::summarize(mnMethylation=mean(Methylation), maxMethylation=max(Methylation), minMethylation=min(Methylation), sdMethylation=sd(Methylation), n=n()) %>% mutate(CI=1.95*sdMethylation/(sqrt(n)))
  
  if(minmax){
    if(pcloud2){
    g <- ggplot() + geom_point(data=plot.df %>% filter(Sample %in% group2), aes(x=Position, y=Methylation), color='red', shape=4, alpha=0.75) + geom_ribbon(data=grp1.avgs, aes(x=Position, ymin=minMethylation, ymax=maxMethylation), fill="gray", alpha=0.5) + facet_wrap(~ Gene+CpgIsland, scales="free_x") + ylab("Mean Methylation Frequency") + scale_x_continuous(label=comma) + theme(axis.text.x = element_text(angle=45, hjust=1)) + xlab("Position within CpG Island")
    } else {
      g <- ggplot() + geom_ribbon(data=grp1.avgs, aes(x=Position, ymin=minMethylation, ymax=maxMethylation), fill="blue", alpha=0.5) +  geom_ribbon(data=grp2.avgs, aes(x=Position, ymin=minMethylation, ymax=maxMethylation), fill="red", alpha=0.5) + facet_wrap(~ Gene+CpgIsland, scales="free_x") + ylab("Mean Methylation Frequency") + scale_x_continuous(label=comma) + theme(axis.text.x = element_text(angle=45, hjust=1)) + xlab("Position within CpG Island")
    }
  } else {
    if(pcloud2){
     g <- ggplot() + geom_point(data=plot.df %>% filter(Sample %in% group2), aes(x=Position, y=Methylation), color='red', shape=4, alpha=0.75) + geom_ribbon(data=grp1.avgs, aes(x=Position, ymin=mnMethylation - CI, ymax=mnMethylation + CI), fill="gray", alpha=0.5) + facet_wrap(~ Gene+CpgIsland, scales="free_x") + ylab("Mean Methylation Frequency") + scale_x_continuous(label=comma) + theme(axis.text.x = element_text(angle=45, hjust=1)) + xlab("Position within CpG Island")
    } else {
      g <- ggplot() + geom_ribbon(data=grp1.avgs, aes(x=Position, ymin=mnMethylation - CI, ymax=mnMethylation + CI), fill="blue", alpha=0.5) +  geom_ribbon(data=grp2.avgs, aes(x=Position, ymin=mnMethylation - CI, ymax=mnMethylation +CI), fill="red", alpha=0.5) + facet_wrap(~ Gene+CpgIsland, scales="free_x") + ylab("Mean Methylation Frequency") + scale_x_continuous(label=comma) + theme(axis.text.x = element_text(angle=45, hjust=1)) + xlab("Position within CpG Island")
    }
  }
  return(g)
}

cpgBoxPlots <- function(cpg.df, testsample, exclusions="", minmax=TRUE,nrow=NA, cpgs="NA", genes="NA"){
  if(is.null(match.call()$cpgs) && is.null(match.call()$genes)){
    print("Supply either a list of genes or a list of CPG island names")
    return()
  } else if(is.null(match.call()$cpgs)){
    print("trying genes")
    plot.df <- cpg.df %>% filter(gene %in% genes) %>% dplyr::select(-any_of(exclusions)) %>% reshape2::melt(id.vars=c("pos", "chr", "cpg", "gene"))
  } else {
    print("trying cpgs")
    plot.df <- cpg.df %>% filter(cpg %in% cpgs) %>% dplyr::select(-any_of(exclusions)) %>% reshape2::melt(id.vars=c("pos", "chr", "cpg", "gene"))
  }
  colnames(plot.df) <- c("Position", "Chromosome", "CpgIsland", "Gene", "Sample", "Methylation")
  plot.df$SampleType <- factor(plot.df$Sample==testsample)
  levels(plot.df$SampleType) <- c("Reference", testsample)
  plot.df$PlotGene <- paste(plot.df$Chromosome, plot.df$Gene, sep=":")
  
  g <- ggplot(plot.df, aes(x=SampleType,y=Methylation,fill=SampleType)) + geom_boxplot(name="Sample") + geom_hline(yintercept=80,linetype='dashed',color='gray', alpha=70) + geom_hline(yintercept=40,linetype='dotted',color='gray', alpha=70) + scale_fill_manual(values=c("gray", "red")) + theme(axis.text.x = element_text(angle=45,hjust=1)) + xlab("") + facet_wrap(~ PlotGene+CpgIsland, scales="free_x") + ylab("Mean Methylation Frequency")

  return(g)  
}