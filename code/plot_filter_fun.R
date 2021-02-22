################################################################################
#### Project: TP Network
#### Title:   Function | Wrangle | Format ASV data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    9 February 2021
#### ---------------------------------------------------------------------------

plot_filter_fun <- function(seqData){
  ## Rank abundances ----
  # Melt by filtering method and aggregate at phylum level
  allMelted <- lapply(
    seqData,
    function(x){
      x %>%
        # make phyloseq object long
        psmelt %>%
        # sum abundances at class level
        aggregate(Abundance ~ Class, ., sum)
    }
  )
  # Join data
  abundances <- plyr::join_all(
    allMelted, 
    by = "Class", 
    type = "left"
  )
  # rename columns
  colnames(abundances) <- c("Class", names(seqData))
  # format
  forPlot <- abundances %>%
    # replace NAs
    replace(is.na(.), 0) %>%
    # change missing phyla to "unknown" and reorder by abundance
    mutate(Class = ifelse(Class == "", "Unknown", Class)) %>%
    mutate(Class = fct_reorder(Class, rawRA, .desc = T))
  # Build base plot
  basePlot <- ggplot(forPlot) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.y = element_text(size = 6)) +
    coord_flip() +
    scale_y_log10() +
    geom_segment() +
    geom_point(shape = 21, fill = "white") +
    ylab(expression(paste("Lo", g[10], " abundance")))
  # Add data aesthetics
  rankAbundances <- list(
    raw = basePlot %+% aes(x = Class, y = rawCounts, xend = Class, yend = 0) %+% xlab(""),
    rab = basePlot %+% aes(x = Class, y = rawRA, xend = Class, yend = 0) %+% xlab(""),
    rag = basePlot %+% aes(x = Class, y = rawGMPR, xend = Class, yend = 0) %+% xlab(""),
    pcr = basePlot %+% aes(x = Class, y = percRA, xend = Class, yend = 0) %+% xlab(""),
    pcg = basePlot %+% aes(x = Class, y = percGMPR, xend = Class, yend = 0) %+% xlab(""),
    pvr = basePlot %+% aes(x = Class, y = prevRA, xend = Class, yend = 0) %+% xlab(""),
    pvg = basePlot %+% aes(x = Class, y = prevGMPR, xend = Class, yend = 0) %+% xlab("")
  )
  ## Do NMDS ----
  nmdsPlots <- lapply(
    seqData,
    function(x){
      # remove NAs (only extreme filtering - should not occur)
      complete <- colSums(is.na(x@otu_table)) == 0
      noNAs <- x@otu_table[, complete]
      # calculate distance and generate NMDS
      dist <- vegan::vegdist(t(noNAs))
      nmds <- vegan::metaMDS(dist, k = 2, trymax = 200)
      # make bi plot data (my function)
      plotter <- make_biplot_data(
        PC1 = nmds$points[, 1], 
        PC2 = nmds$points[, 2],
        groupVar = x@sam_data$site[complete])
      # plot to object and return
      out <- ggplot(plotter) +
        theme_bw() +
        theme(panel.grid = element_blank()) +
        guides(col = "none") +
        aes(col = groupVar, x = PC1, y = PC2, xend = meanPC1, yend = meanPC2) +
        geom_point() +
        geom_segment() +
        labs(x = "MDS1", y = "MDS2")
      return(out)
    }
  )
  ## Collate plots ----
  # generate filenames
  files <- paste0("./plots/qc_filtering/", "fun_", names(seqData), ".pdf")
  # plot
  for(i in 1:length(files)){
    pdf(file = files[i], width = mm2in(240), height = mm2in(100))
    print(cowplot::plot_grid(
      rankAbundances[[i]], nmdsPlots[[i]],
      align = "h", rel_widths = c(2, 1)
    ))
    dev.off()
  }
}
