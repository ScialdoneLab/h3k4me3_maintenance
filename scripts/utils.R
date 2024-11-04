packages <- c("here", "ggplot2", "ChIPpeakAnno", "reshape2", "AnnotationDbi", "eulerr", "PerformanceAnalytics", "factoextra", "chromstaR", "BSgenome",
"clusterProfiler", "org.Xl.eg.db", "simplifyEnrichment", "gridExtra", "grid", "plyr", "ggalluvial", "tidyr", "scales", "corrplot", "RColorBrewer", "ComplexUpset",
"dplyr", "ggpubr", "Biostrings", "EnrichedHeatmap", "rtracklayer", "SummarizedExperiment", "Seurat", "circlize", "ComplexHeatmap", "reshape2")
lapply(packages, require, character.only = TRUE)

stages <- list("Spermatid", "Sperm", "Pre-ZGA", "ZGA")

go_enrichment <- function(data) {
  Xl_ENTREZ <- read.delim("data/raw/genome/GenePageLaevisEntrezGeneUnigeneMapping.txt", header = FALSE)
  entrez <- Xl_ENTREZ[Xl_ENTREZ$V2 %in% data, ] 
  colnames(entrez) <- c("Xenbase ID", "Gene symbol", "Entrez ID", "Xl ID unknown")
  GO <- enrichGO(entrez$`Entrez ID`, OrgDb = org.Xl.eg.db, ont = "BP", readable = TRUE)
  return(GO)
}
plot_go <- function(enrichGO, name){
  dp <- dotplot(enrichGO)
  cn <- cnetplot(enrichGO)
  grid.arrange(dp, cn, top=name, nrow=2)
}

readGranges <- function(filename) {
  granges <- readRDS(filename)
  suppressWarnings({
    granges <- resize(granges, 200)
  })
  granges <- trim(granges)
  granges_chr <- granges[seqnames(granges) %in% chr]
}

enrichment <- function(coverage, title, mode="ranges", value_name="H3K4me3 enrichment", enriched_score=FALSE){
  if(mode=="ranges"){
    coverage <- normalizeToMatrix(coverage, promoters_chr, extend = 1000, mean_mode = "coverage")
  } else if(mode=="coverage"){
    coverage <- normalizeToMatrix(coverage, promoters_chr, extend = 1000, mean_mode = "coverage", value_column = "score")
  }
  coverage <- log2(coverage+1)
  hm <- EnrichedHeatmap(coverage, 
                        name = value_name,
                        column_title = title,
                        col = c("blue", "red"),
                        axis_name = c("-1kb", "TSS", "+1kb"))
  draw(hm)
  if(enriched_score){
    return(enriched_score(coverage))
  }else{

  }
  return(coverage)
}

label_dyn <- function(gene, dyns){for(i in 1:length(dyns)){if(any(gene %in% dyns[[i]])){return(names(dyns)[i])}}; return(NA)}
plot_peak <- function(stage, dyns, cols, onlypeaks=FALSE, remove_NA=TRUE){
  stage_raw <- data.frame(chip_norm_mats[[paste(stage)]])
  if(onlypeaks){stage_raw <- stage_raw[chip_dfs[[paste(stage)]]$peak, ]}
  stage_raw$dyn <- sapply(rownames(stage_raw), label_dyn, dyns = dyns)
  stage_raw$gene <- rownames(stage_raw)
  if(remove_NA){stage_raw <- stage_raw[!is.na(stage_raw$dyn),]}
  melted <- melt(stage_raw)
  colnames(melted) <- c("dynamics", "Gene", "Genomic location", "Coverage")
  plot <- ggline(melted, x = "Genomic location", y = "Coverage",
   add = "mean_se", conf.int = TRUE, color="dynamics") + scale_color_manual(values=cols) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + labs(y=expression(atop(Mean~H3K4me3~ChIP-Seq~intensity,log[2](Coverage + 1))), x = "+/-1kb of TSS") + ggtitle(stage) + guides(color = "none")
}

plot_stage_dynamics <- function(data, dyns, cols, timepoint, type, onlypeaks=FALSE, legend=FALSE, remove_NA=TRUE, grouped=FALSE, levels=NULL, only_violin=FALSE){
  plot_df <- data[[timepoint]]
  if(onlypeaks){plot_df <- plot_df[plot_df$peak, ]}
  plot_df$genes <- rownames(plot_df)
  
  plot_df$dyn <- factor(sapply(plot_df$genes, label_dyn, dyns = dyns))
  if(remove_NA){plot_df <- plot_df[!is.na(plot_df$dyn),]}
  #plot_df <- na.omit(plot_df)
  
  if(!is.null(levels)){
    plot_df <- plot_df[as.vector(plot_df$dyn) %in% levels,]
    plot_df$dyn <- factor(plot_df$dyn)
  }
  if(grouped){
    plot_df$dyn_group <- factor(sapply(plot_df$genes, label_dyn, dyns = chip_dyns))
    cluster_order <- levels(reorder(plot_df$dyn_group, plot_df[,paste(type)], median, decreasing=TRUE))
    pos <- seq(1, length(cluster_order)*2, 2)
    get_sub_order <- function(cluster){
      sub_order <-  levels(droplevels(reorder(plot_df[plot_df$dyn_group==cluster, "dyn"], 
                                   plot_df[plot_df$dyn_group == cluster,paste(type)], 
                                   median, decreasing=TRUE)))
    }
    for(dyn_group in levels(plot_df$dyn_group)){
      sub_order = get_sub_order(dyn_group)
      dyn_group_idx = which(cluster_order == dyn_group)
      if(length(sub_order) > 0){
        if(dyn_group_idx < length(pos)){
          pos[(dyn_group_idx+1):length(pos)] <- pos[(dyn_group_idx+1):length(pos)] + length(sub_order) - 1
        }
        pos <- append(pos, seq(pos[dyn_group_idx], pos[dyn_group_idx]+length(sub_order)-1, 1), after=dyn_group_idx) 
      }
      else {
        pos[(dyn_group_idx+1):length(pos)] <- pos[(dyn_group_idx+1):length(pos)] -2
      }
      pos <- pos[-dyn_group_idx]
      cluster_order <- append(cluster_order, sub_order, after=dyn_group_idx)
      cluster_order <- cluster_order[-dyn_group_idx]
    }
    plot_df$dyn <- factor(pos[as.numeric(factor(plot_df$dyn, levels = cluster_order))], levels = seq(1, max(pos)))
    cols = setNames(cols[cluster_order], pos)
  }
  else {
    plot_df$dyn <- factor(plot_df$dyn, levels = levels(reorder(plot_df$dyn, plot_df[,paste(type)], median, decreasing=TRUE)))
  }
  
  violin_plot <- ggviolin(plot_df, x = "dyn", y = paste(type), fill = "dyn",
           add = "boxplot", add.params = list(fill = "white")) + scale_fill_manual(values=cols) + rotate_x_text(45) + guides(fill="none") + ylab(paste(axis_titles[[type]])) + theme(axis.title.x=element_blank())
  if(grouped){violin_plot <- violin_plot + scale_x_discrete(breaks = pos, labels=cluster_order, drop=FALSE)}
  if(only_violin){return(violin_plot)}
  ecdf_plot <- ggecdf(plot_df, x = paste(type), color = "dyn") + scale_colour_manual(values=cols) + guides(color = "none") + xlab(paste(axis_titles[[type]]))
  if(legend){
    violin_plot <- violin_plot + guides(fill = guide_legend(nrow = 2))
    ecdf_plot <- ecdf_plot + guides(color = guide_legend(nrow = 2))
  }
  if(type=="cg.density"){title<-""}else{title <- paste(timepoint)}
  arrangeGrob(violin_plot, ecdf_plot, top=title, ncol=2) #grid.arrange for direct plotting
}

plot_balloonplot <- function(x_dyns, y_dyns, limit=50, sort_similar=FALSE, remove_empty=FALSE, remove_dyns=c(), x_caption="", y_caption=""){
  
  genesAbsolute <- lapply(y_dyns, function(y){as.data.frame(lapply(x_dyns, function(x){length(intersect(y, x))}))})
  genesAbsolute <- do.call(rbind, genesAbsolute)
  
  fisher_p <- function(row, column){
    mat = genesAbsolute
    fisher_mat <- matrix(c(mat[row, column], sum(mat[-row, column]), sum(mat[row, -column]), sum(mat[-row, -column])), nrow=2)
    p_val <- fisher.test(fisher_mat)$p.value
    return(max(p_val, 10^(-limit)))
  }
  
  p_mat = data.frame(matrix(NA, nrow = nrow(genesAbsolute), ncol = ncol(genesAbsolute)), row.names = rownames(genesAbsolute))
  colnames(p_mat) = colnames(genesAbsolute)
  
  for (row in 1:nrow(p_mat)){
    for (col in 1:ncol(p_mat)){
      p_mat[row, col] <- fisher_p(row, col)
    }
  }
  
  p_mat <- matrix(p.adjust(as.matrix(p_mat), method = "fdr"), ncol=ncol(p_mat), dimnames=list(rownames(p_mat), colnames(p_mat)))
  p_mat = -log10(p_mat)
  p_mat[genesAbsolute == 0] <- NA
  
  if(sort_similar){
    p_mat[p_mat < (-log10(0.05))] <- 0
    y_dist<-dist(p_mat, method="euclidean")
    y_clust<-hclust(y_dist, method='complete')
  
    x_dist<-dist(t(p_mat), method="euclidean")
    x_clust<-hclust(x_dist, method='complete')

    p_mat <- p_mat[y_clust$order, x_clust$order]
  }
  
  p_mat[p_mat < (-log10(1e-2))] <- NA

  if(remove_empty){
    p_mat <- p_mat[rowSums(is.na(p_mat)) != ncol(p_mat), ] 
  }
  p_mat <- p_mat[!rownames(p_mat) %in% remove_dyns, !colnames(p_mat) %in% remove_dyns]
  
  my_cols <- rev(c("#DC267F", "#FE6100", "#FFB000", "#FFFFFF"))
  ggballoonplot(p_mat, fill="value") +
                        theme_light() + labs(x="", y="") + labs(x=x_caption, y=y_caption) + rotate_x_text(45) +
                        scale_fill_gradientn(colors = my_cols, limits=c(0,limit), values=c(0, 0.6, 0.7, 1), breaks=c(2, 10, 25, 35, limit), labels = c(1e-2, 1e-10, 1e-25, 1e-35, 10^(-limit))) + guides(size=guide_legend(title="adjusted p value"), fill=guide_legend(title="adjusted p value")) + scale_size(range = c(0, 10), limits=c(0, limit), breaks =c(2, 10, 25, 35, limit), labels = c(1e-2, 1e-10, 1e-25, 1e-35, 10^(-limit)))
}

assignChromosomeRegion <- function(peaks.RD, exon, TSS, utr5, utr3, proximal.promoter.cutoff = c(
                                     upstream = 2000,
                                     downstream = 100
                                   ), immediate.downstream.cutoff = c(
                                     upstream = 0,
                                     downstream = 1000
                                   ), nucleotideLevel = FALSE, precedence = NULL,
                                   TxDb = NULL) {
  if (!is.null(TxDb)) {
    if (!inherits(TxDb, "TxDb")) {
      stop("TxDb must be an object of TxDb, \n                     try\n?TxDb\tto see more info.")
    }
    if (!inherits(peaks.RD, c("GRanges"))) {
      stop("peaks.RD must be a GRanges object.")
    }
    if (!all(c("upstream", "downstream") %in% names(proximal.promoter.cutoff))) {
      stop("proximal.promoter.cutoff must contain elements upstream and downstream")
    }
    if (!all(c("upstream", "downstream") %in% names(immediate.downstream.cutoff))) {
      stop("immediate.downstream.cutoff must contain elements upstream and downstream")
    }
    if (!is.null(precedence)) {
      if (!all(precedence %in% c(
        "Exons", "Introns", "fiveUTRs",
        "threeUTRs", "Promoters", "immediateDownstream"
      ))) {
        stop("precedence must be a combination of \n                         Exons, Introns, fiveUTRs, threeUTRs, \n                         Promoters, immediateDownstream")
      }
    }
    ignore.strand <- all(as.character(strand(peaks.RD)) ==
      "*")
    exons <- exons(TxDb, columns = NULL)
    introns <- unique(unlist(intronsByTranscript(TxDb)))
    fiveUTRs <- unique(unlist(fiveUTRsByTranscript(TxDb)))
    threeUTRs <- unique(unlist(threeUTRsByTranscript(TxDb)))
    transcripts <- unique(transcripts(TxDb, columns = NULL))
    options(warn = -1)
    try({
      promoters <- unique(promoters(TxDb,
        upstream = proximal.promoter.cutoff["upstream"],
        downstream = proximal.promoter.cutoff["downstream"]
      ))
      immediateDownstream <- unique(downstreams(transcripts,
        upstream = immediate.downstream.cutoff["upstream"],
        downstream = immediate.downstream.cutoff["downstream"]
      ))
      promoters <- GenomicRanges::trim(promoters)
      immediateDownstream <- GenomicRanges::trim(immediateDownstream)
    })
    microRNAs <- tryCatch(microRNAs(TxDb), error = function(e) {
      return(NULL)
    })
    tRNAs <- tryCatch(tRNAs(TxDb), error = function(e) {
      return(NULL)
    })
    options(warn = 0)
    annotation <- list(
      exons, introns, fiveUTRs, threeUTRs,
      promoters, immediateDownstream
    )
    if (!is.null(microRNAs)) {
      annotation <- c(annotation, microRNAs = microRNAs)
    }
    if (!is.null(tRNAs)) {
      annotation <- c(annotation, tRNAs = tRNAs)
    }
    annotation <- lapply(annotation, function(.anno) {
      mcols(.anno) <- NULL
      .anno
    })
    names(annotation)[1:6] <- c(
      "Exons", "Introns", "fiveUTRs",
      "threeUTRs", "Promoters", "immediateDownstream"
    )
    # peaks.RD <- formatSeqnames(peaks.RD, exons)
    peaks.RD <- unique(peaks.RD)
    annotation <- GRangesList(annotation)
    newAnno <- c(unlist(annotation))
    if (ignore.strand) {
      newAnno.rd <- newAnno
      strand(newAnno.rd) <- "*"
      newAnno.rd <- reduce(trim(newAnno.rd))
      Intergenic.Region <- gaps(newAnno.rd, end = seqlengths(TxDb))
      Intergenic.Region <- Intergenic.Region[strand(Intergenic.Region) ==
        "*"]
    } else {
      newAnno.rd <- reduce(trim(newAnno))
      Intergenic.Region <- gaps(newAnno.rd, end = seqlengths(TxDb))
      Intergenic.Region <- Intergenic.Region[strand(Intergenic.Region) !=
        "*"]
    }
    if (!all(seqlevels(peaks.RD) %in% seqlevels(newAnno))) {
      warning("peaks.RD has sequence levels not in TxDb.")
      sharedlevels <- intersect(seqlevels(newAnno), seqlevels(peaks.RD))
      peaks.RD <- keepSeqlevels(peaks.RD, sharedlevels,
        pruning.mode = "coarse"
      )
    }
    mcols(peaks.RD) <- NULL
    if (!is.null(precedence)) {
      annotation <- annotation[unique(c(precedence, names(annotation)))]
    }
    names(Intergenic.Region) <- NULL
    annotation$Intergenic.Region <- Intergenic.Region
    anno.names <- names(annotation)
    ol.anno <- findOverlaps(peaks.RD, annotation, ignore.strand = ignore.strand)
    if (nucleotideLevel) {
      jaccardIndex <- unlist(lapply(annotation, function(.ele) {
        intersection <- intersect(.ele, peaks.RD, ignore.strand = ignore.strand)
        union <- union(.ele, peaks.RD, ignore.strand = ignore.strand)
        sum(as.numeric(width(intersection))) / sum(as.numeric(width(union)))
      }))
      jaccardIndex <- jaccardIndex[anno.names]
      names(jaccardIndex) <- anno.names
      jaccardIndex[is.na(jaccardIndex)] <- 0
      newAnno <- unlist(annotation)
      newAnno$source <- rep(names(annotation), lengths(annotation))
      newAnno.disjoin <- disjoin(newAnno,
        with.revmap = TRUE,
        ignore.strand = ignore.strand
      )
      if (!is.null(precedence)) {
        revmap <- cbind(
          from = unlist(newAnno.disjoin$revmap),
          to = rep(seq_along(newAnno.disjoin), lengths(newAnno.disjoin$revmap))
        )
        revmap <- revmap[order(revmap[, "to"], revmap[
          ,
          "from"
        ]), , drop = FALSE]
        revmap <- revmap[!duplicated(revmap[, "to"]), ,
          drop = FALSE
        ]
        newAnno.disjoin$source <- newAnno[revmap[, "from"]]$source
      } else {
        revmap <- unlist(newAnno.disjoin$revmap)
        newAnno.disjoin <- rep(newAnno.disjoin, lengths(newAnno.disjoin$revmap))
        newAnno.disjoin$source <- newAnno[revmap]$source
      }
      ol.anno <- findOverlaps(peaks.RD, newAnno.disjoin,
        ignore.strand = ignore.strand
      )
      queryHits <- peaks.RD[queryHits(ol.anno)]
      subjectHits <- newAnno.disjoin[subjectHits(ol.anno)]
      totalLen <- sum(as.numeric(width(peaks.RD)))
      queryHits.list <- split(queryHits, subjectHits$source)
      lens <- unlist(lapply(queryHits.list, function(.ele) sum(as.numeric(width(unique(.ele))))))
      percentage <- 100 * lens / totalLen
    } else {
      ol.anno.splited <- split(queryHits(ol.anno), anno.names[subjectHits(ol.anno)])
      jaccardIndex <- unlist(lapply(anno.names, function(.name) {
        union <- length(annotation[[.name]]) + length(peaks.RD) -
          length(unique(subjectHits(findOverlaps(peaks.RD,
            annotation[[.name]],
            ignore.strand = ignore.strand
          ))))
        intersection <- length(ol.anno.splited[[.name]])
        intersection / union
      }))
      names(jaccardIndex) <- anno.names
      ol.anno <- as.data.frame(ol.anno)
      ol.anno.splited <- split(ol.anno, ol.anno[, 2])
      hasAnnoHits <- do.call(rbind, ol.anno.splited[names(ol.anno.splited) !=
        as.character(length(annotation))])
      hasAnnoHits <- unique(hasAnnoHits[, 1])
      ol.anno <- ol.anno[!(ol.anno[, 2] == length(annotation) &
        (ol.anno[, 1] %in% hasAnnoHits)), ]
      if (!is.null(precedence)) {
        ol.anno <- ol.anno[!duplicated(ol.anno[, 1]), ]
      }
      subjectHits <- anno.names[ol.anno[, 2]]
      counts <- table(subjectHits)
      percentage <- 100 * counts / length(peaks.RD)
    }
    len <- length(anno.names) - length(percentage)
    if (len > 0) {
      tobeadd <- rep(0, len)
      names(tobeadd) <- anno.names[!anno.names %in% names(percentage)]
      percentage <- c(percentage, tobeadd)
    }
    percentage <- percentage[anno.names]
    return(list(percentage = percentage, jaccard = jaccardIndex))
  } else {
    message("Please try to use TxDb next time. Try\n\n                    ?TxDb\tto see more info.")
    annotationList <- list(exon, TSS, utr5, utr3)
    names(annotationList) <- c("Exon", "TSS", "UTR5", "UTR3")
    status <- lapply(annotationList, function(.ele) {
      if (!inherits(.ele, "GRanges")) {
        stop("Annotation of exon, TSS, utr5, utr3 must \n                         be objects of GRanges.")
      }
    })
    if (!inherits(peaks.RD, "GRanges")) {
      stop("peaks.RD must be a GRanges object.")
    }
    ann.peaks <- annotatePeakInBatch(peaks.RD, AnnotationData = TSS)
    ann.peaks <- ann.peaks[!is.na(ann.peaks$distancetoFeature)]
    upstream <- ann.peaks[ann.peaks$insideFeature == "upstream" |
      (ann.peaks$distancetoFeature < 0 & ann.peaks$insideFeature ==
        "overlapStart" & abs(ann.peaks$distancetoFeature) >
        ann.peaks$shortestDistance) | ann.peaks$insideFeature ==
      "includeFeature" | (ann.peaks$distancetoFeature >=
      0 & ann.peaks$insideFeature == "overlapStart" & ann.peaks$distancetoFeature ==
      ann.peaks$shortestDistance)]
    proximal.promoter.n <- length(upstream[upstream$distancetoFeature >=
      -proximal.promoter.cutoff | upstream$shortestDistance <=
      proximal.promoter.cutoff])
    enhancer.n <- length(upstream) - proximal.promoter.n
    downstream <- ann.peaks[ann.peaks$insideFeature == "downstream"]
    immediateDownstream.n <- length(downstream[downstream$distancetoFeature <=
      immediate.downstream.cutoff, ])
    enhancer.n <- enhancer.n + dim(downstream[downstream$distancetoFeature >
      immediate.downstream.cutoff, ])
    inside.peaks <- ann.peaks[ann.peaks$insideFeature ==
      "inside" | ann.peaks$insideFeature == "overlapEnd" |
      (ann.peaks$insideFeature == "overlapStart" & ann.peaks$distancetoFeature >=
        0 & ann.peaks$distancetoFeature != ann.peaks$shortestDistance) |
      (ann.peaks$insideFeature == "overlapStart" & ann.peaks$distancetoFeature <
        0 & abs(ann.peaks$distancetoFeature) == ann.peaks$shortestDistance)]
    ann.utr5.peaks <- annotatePeakInBatch(inside.peaks, AnnotationData = utr5)
    proximal.promoter.n <- proximal.promoter.n + length(ann.utr5.peaks[ann.utr5.peaks$insideFeature ==
      "upstream"])
    utr5.n <- length(ann.utr5.peaks[ann.utr5.peaks$insideFeature %in%
      c("includeFeature", "inside") | (ann.utr5.peaks$insideFeature ==
      "overlapStart" & ann.utr5.peaks$distancetoFeature >=
      0 & ann.utr5.peaks$distancetoFeature != ann.utr5.peaks$shortestDistance) |
      (ann.utr5.peaks$insideFeature == "overlapStart" &
        ann.utr5.peaks$distancetoFeature < 0 & abs(ann.utr5.peaks$distancetoFeature) ==
        ann.utr5.peaks$shortestDistance) | (ann.utr5.peaks$insideFeature ==
      "overlapEnd" & ann.utr5.peaks$strand == "+" & abs(start(ann.utr5.peaks) -
      ann.utr5.peaks$end_position) >= (end(ann.utr5.peaks) -
      ann.utr5.peaks$end_position)) | (ann.utr5.peaks$insideFeature ==
      "overlapEnd" & ann.utr5.peaks$strand == "-" & abs(end(ann.utr5.peaks) -
      ann.utr5.peaks$start_position) >= abs(start(ann.utr5.peaks) -
      ann.utr5.peaks$start_position))])
    proximal.promoter.n <- proximal.promoter.n + length(ann.utr5.peaks[(ann.utr5.peaks$insideFeature ==
      "overlapStart" & ann.utr5.peaks$distancetoFeature >=
      0 & ann.utr5.peaks$distancetoFeature == ann.utr5.peaks$shortestDistance) |
      (ann.utr5.peaks$insideFeature == "overlapStart" &
        ann.utr5.peaks$distancetoFeature < 0 & abs(ann.utr5.peaks$distancetoFeature) !=
        ann.utr5.peaks$shortestDistance)])
    downstream.utr5 <- ann.utr5.peaks[ann.utr5.peaks$insideFeature ==
      "downstream" | (ann.utr5.peaks$insideFeature == "overlapEnd" &
      ann.utr5.peaks$strand == "+" & abs(start(ann.utr5.peaks) -
      ann.utr5.peaks$end_position) < (end(ann.utr5.peaks) -
      ann.utr5.peaks$end_position)) | (ann.utr5.peaks$insideFeature ==
      "overlapEnd" & ann.utr5.peaks$strand == "-" & abs(end(ann.utr5.peaks) -
      ann.utr5.peaks$start_position) < abs(start(ann.utr5.peaks) -
      ann.utr5.peaks$start_position))]
    ann.utr3.peaks <- annotatePeakInBatch(downstream.utr5,
      AnnotationData = utr3
    )
    utr3.n <- length(ann.utr3.peaks[ann.utr3.peaks$insideFeature %in%
      c(
        "includeFeature", "overlapStart", "overlapEnd",
        "inside"
      )])
    rest.peaks <- ann.utr3.peaks[ann.utr3.peaks$insideFeature %in%
      c("downstream", "upstream")]
    ann.rest.peaks <- annotatePeakInBatch(rest.peaks, AnnotationData = exon)
    intron.n <- length(ann.rest.peaks[ann.rest.peaks$insideFeature %in%
      c("downstream", "upstream")])
    exon.n <- length(ann.rest.peaks) - intron.n
    total <- length(peaks.RD) / 100
    list(
      Exons = exon.n / total, Introns = intron.n / total,
      fiveUTRs = utr5.n / total, threeUTRs = utr3.n / total,
      Promoters = proximal.promoter.n / total, immediate.Downstream = immediateDownstream.n / total,
      Intergenic.Region = enhancer.n / total
    )
  }
}