library(knitr)
library(kableExtra)
library(dplyr)
library(alakazam)
library(shazam)
library(stringr)

library(data.table)
library(igraph)
library(gplots)
library(circlize)
library(UpSetR)
library(gtools)

theme_set(theme_bw(base_family = "ArialMT") +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(family="ArialMT")))

outdir <- "repertoire_per_patient/repertoire_comparison"
# setwd to results folder (containing alakazam, shazam, etc. folders)
### Read all the tables as produced by the pipeline in the current folder and joins them together in the df_all dataframe
#all_files <- system(paste0("find '",datadir,"' -name '*germ-pass.tsv'"), intern=T)
dir.create(outdir)

df_all <- read.csv("../pipeline_results_per_patient/repertoire_analysis/repertoire_comparison/all_data.tsv", sep="\t")

# Fix swap of samples
#df_clad4 <- df_all[which(df_all$source == "CLAD4"),]

df_all <- df_all[!(df_all$source == "CLAD4"),]

#df_clad4$correct_extract_time <- sapply(df_clad4$extract_time, function(x) str_replace(as.character(x), "6months", "fix"))
#df_clad4$correct_extract_time <- sapply(df_clad4$correct_extract_time, function(x) str_replace(as.character(x), "baseline", "6months"))
#df_clad4$correct_extract_time <- sapply(df_clad4$correct_extract_time, function(x) str_replace(as.character(x), "fix", "baseline"))
#df_clad4$extract_time <- df_clad4$correct_extract_time
#df_clad4$correct_extract_time <- NULL

#df_all <- rbind(df_all, df_clad4)

# Remove underscores in these columns
df_all$treatment <- sapply(df_all$treatment, function(x) str_replace(as.character(x), "_", ""))
df_all$source <- sapply(df_all$source, function(x) str_replace(as.character(x), "_", ""))
df_all$extract_time <- sapply(df_all$extract_time, function(x) str_replace(as.character(x), "_", ""))
df_all$population <- sapply(df_all$population, function(x) str_replace(as.character(x), "_", ""))
# Annotate sample and samplepop (sample + population) by add ing all the conditions
df_all$sample <- as.factor(paste(df_all$treatment, df_all$extract_time, df_all$source, sep="_"))
df_all$sample_pop <- as.factor(paste(df_all$treatment, df_all$extract_time, df_all$source, df_all$population, sep="_"))

df_subset <- df_all[!(df_all$source == "CLAD4"),]

# Splitting data in a per patient basis
df_subset <- df_subset[,c("treatment", "extract_time", "source", "population",
                          "clone_id",
                          "v_call", "d_call", "j_call", "junction_length",
                          "sample", "sample_pop")]


df_pat <- split(df_subset, df_subset$source)
#df_pat$CLAD4<-NULL
dir.create(paste(outdir,"Clone_overlap", sep = "/"))

all_clonedf_tp <- data.frame()
all_seqdf_tp <- data.frame()
all_clonedf <- data.frame()
all_seqdf <- data.frame()
all_clonedf_pop_tp <- data.frame()
all_seqdf_pop_tp <- data.frame()

for (i in c(1:length(df_pat))) {
  patdir_overlap <- paste(outdir,"Clone_overlap",names(df_pat)[i], sep="/")
  dir.create(patdir_overlap)
  
  ###################################################################
  ## Plot chordplot comparison between time points
  ###################################################################
  df_time <- split(df_pat[[i]], df_pat[[i]]$extract_time)
  
  # Need to calculate vennplot to be able to count the clones without overlaps to other time points
  vennplot <- venn(list(unique(df_time[[1]]$clone_id), unique(df_time[[2]]$clone_id), unique(df_time[[3]]$clone_id)), names = names(df_time))
  count_clones <- countClones(df_pat[[i]])
  
  combin_tp <- expand.grid(c("baseline"), c("6months","12months"))
  clonedf_tp <- rbind(combin_tp, 
                      data.frame(Var1=c("6months"), Var2=c("12months")))
  names(clonedf_tp) <- c("from","to")
  
  seqdf_tp <- clonedf_tp
  
  lenintersects = numeric(0)
  seqsintersects = numeric(0)
  for (j in c(1:nrow(clonedf_tp))){
    inter <- intersect(df_time[[which(grepl(clonedf_tp[j,1], names(df_time)))]]$clone_id,
                       df_time[[which(grepl(clonedf_tp[j,2], names(df_time)))]]$clone_id)
    
    clones_subset <- count_clones[which(count_clones$clone_id %in% as.character(inter)),]
    
    lenintersects <- c(lenintersects, length(inter))
    seqsintersects <- c(seqsintersects, sum(clones_subset$seq_count))
  }
  
  clonedf_tp$value <- lenintersects
  seqdf_tp$value <- seqsintersects
  
  self_comb <- data.frame(from = names(df_time), to = names(df_time))
  self_clonedf_tp <- self_comb
  self_seqdf_tp <- self_comb
  
  lenintersects <- numeric(0)
  seqsintersects <- numeric(0)
  for (tp in self_comb$from){
    inter <- attributes(vennplot)[["intersections"]][[tp]]
    clones_subset <- count_clones[which(count_clones$clone_id %in% as.character(inter)),]
    
    lenintersects <- c(lenintersects, length(inter))
    seqsintersects <- c(seqsintersects, sum(clones_subset$seq_count))
  }
  self_clonedf_tp$value <- lenintersects
  self_seqdf_tp$value <- seqsintersects
  
  clonedf_tp <- rbind(clonedf_tp, self_clonedf_tp)
  seqdf_tp <- rbind(seqdf_tp, self_seqdf_tp)
  

  # Plot chordplots between time points
  grid.col = c("#a6cee3", "#b2df8a", "#fb9a99")
  names(grid.col) <- names(df_time)
  
  
  # Clone overlap plot
  svg(filename = paste(patdir_overlap,"/Clone_overlap_comparison_time_points_",
                       df_time[[1]]$treatment[1], "_",
                       df_time[[1]]$source[1], ".svg", sep=""))
  chordDiagram(clonedf_tp,
               grid.col = grid.col,
               self.link = 1,
               transparency = 0.3,
               annotationTrack="grid",
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(clonedf_tp))))))
  circos.track(track.index = 1,
               panel.fun = function(x, y) {
                 circos.text(CELL_META$xcenter,
                             CELL_META$ylim[2],
                             CELL_META$sector.index,
                             adj = c(0, 0.5))
               }, bg.border = NA)
  title(paste("clone overlap", df_time[[1]]$treatment[1], df_time[[1]]$source[1]), cex = 0.8)
  circos.clear()
  dev.off()
  
  png(filename = paste(patdir_overlap, "/Clone_overlap_comparison_time_points_", 
                       df_time[[1]]$treatment[1], "_", df_time[[1]]$source[1], ".png", sep=""), 
      width=15, height=15, units = "cm", res = 300)
  chordDiagram(clonedf_tp,
               grid.col = grid.col,
               self.link = 1,
               transparency = 0.3,
               annotationTrack="grid",
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(clonedf_tp))))))
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[2], CELL_META$sector.index,
                adj = c(0, 0.5))
  }, bg.border = NA)
  title(paste("clone_id OVERLAP", df_time[[1]]$treatment[1], df_time[[1]]$source[1]), cex = 0.8)
  circos.clear()
  dev.off()
  
  # Sequences overlap plot
  svg(filename = paste(patdir_overlap,"/Clone_seqN_overlap_comparison_time_points_", 
                       df_time[[1]]$treatment[1], "_", 
                       df_time[[1]]$source[1], ".svg", sep=""))
  chordDiagram(seqdf_tp,
               grid.col = grid.col,
               self.link = 1,
               transparency = 0.3,
               annotationTrack="grid",
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(seqdf_tp))))))
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[2], CELL_META$sector.index,
                adj = c(0, 0.5))
  }, bg.border = NA)
  title(paste("CLONE SEQ NUM OVERLAP", df_time[[1]]$treatment[1], df_time[[1]]$source[1]), cex = 0.8)
  circos.clear()
  dev.off()
  
  png(filename = paste(patdir_overlap, "/Clone_seqN_overlap_comparison_time_points_", 
                       df_time[[1]]$treatment[1], "_", df_time[[1]]$source[1], ".png", sep=""), 
      width=15, height=15, units = "cm", res = 300)
  chordDiagram(seqdf_tp,
               grid.col = grid.col,
               self.link = 1,
               transparency = 0.3,
               annotationTrack="grid",
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(seqdf_tp))))))
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[2], CELL_META$sector.index,
                adj = c(0, 0.5))
  }, bg.border = NA)
  title(paste("CLONE SEQ NUM OVERLAP", df_time[[1]]$treatment[1], df_time[[1]]$source[1]), cex = 0.8)
  circos.clear()
  dev.off()
  
  # Upset plots per time point
  
  listinput <- list()
  for (k in c(1:length(df_time))){
    listinput <- append(listinput, list(df_time[[k]]$clone_id))
  }    
  names(listinput) <- names(df_time)
  
  
  svg(filename = paste(patdir_overlap,"/Set_plot_timepoint", df_time[[1]]$treatment[1], "_",df_time[[1]]$source[1], ".svg", sep="" ),
      width = 20*0.4, height=10*0.4)
  print(upset( fromList(listinput), 
               nsets = length(listinput), 
               nintersects = NA,
               keep.order=TRUE,
               order.by = "freq",
               #group.by = "sets",
               point.size = 3.5, 
               line.size=2, 
               mainbar.y.label = "Clone intersections", 
               sets.x.label = "Clones per time point"))
  
  dev.off()
  
  png(filename = paste(patdir_overlap,"/Set_plot_timepoint", df_time[[1]]$treatment[1], "_",df_time[[1]]$source[1], ".png", sep="" ),
      res=600,width = 15, height=10, units="cm")
  print(upset( fromList(listinput), 
               nsets = length(listinput), 
               nintersects = NA,
               keep.order=TRUE,
               order.by = "freq",
               #group.by = "sets",
               point.size = 3.5, 
               line.size=2, 
               mainbar.y.label = "Clone intersections", 
               sets.x.label = "Clones per time point"))
  
  dev.off()
  
  
   # Saving both tables
  clonedf_tp$source <- rep(df_time[[1]]$source[1], nrow(clonedf_tp))
  seqdf_tp$source <- rep(df_time[[1]]$source[1], nrow(seqdf_tp))
  
  # Saving both tables
  write.table(clonedf_tp, file = paste(patdir_overlap,"/Clone_overlap_comparison_time_points_", df_time[[1]]$treatment[1], "_", df_time[[1]]$source[1], ".tsv", sep=""), sep = "\t", quote = F, row.names = F)
  write.table(seqdf_tp, file = paste(patdir_overlap,"/Clone_seqN_overlap_comparison_time_points_", df_time[[1]]$treatment[1], "_", df_time[[1]]$source[1], ".tsv", sep=""), sep = "\t", quote = F, row.names = F)
  
  all_clonedf_tp <- rbind(all_clonedf_tp, clonedf_tp)
  all_seqdf_tp <- rbind(all_seqdf_tp, seqdf_tp)
  
  ###################################################################
  ## Plot chordplot comparison time points per patient and population
  ###################################################################
  
  df_pat[[i]]$time_pop <- as.factor(paste(df_pat[[i]]$extract_time, df_pat[[i]]$population, sep="_"))
  df_pop_time <- split(df_pat[[i]], df_pat[[i]]$time_pop)
  
  count_clones <- countClones(df_pat[[i]])
  
  ## Calculating overlaps between time points per patient and population
  combin <- expand.grid(unique(df_pat[[i]]$extract_time), unique(df_pat[[i]]$population))
  combin$names <- apply(combin[,c("Var1", "Var2")], 1, paste, collapse = "_")
  
  baselines <- subset(combin, combin$Var1 == "baseline")
  six_months <- subset(combin, combin$Var1 == "6months")
  other <- subset(combin, combin$Var1 != "baseline")
  twelve_months <- subset(combin, combin$Var1 == "12months")
  
  clonedfb <- expand.grid(levels(as.factor(baselines$names)), levels(as.factor(other$names)))
  clonedf612 <- expand.grid(levels(as.factor(six_months$names)), levels(as.factor(twelve_months$names)))
  clonedf <- rbind(clonedfb, clonedf612)
  colnames(clonedf) <- c("from","to")
  seqdf <- clonedf
  
  lenintersects = numeric(0)
  seqsintersects = numeric(0)
  for (j in c(1:nrow(clonedf))){
    inter <- intersect(df_pop_time[[which(grepl(clonedf[j,1], names(df_pop_time)))]]$clone_id,
                       df_pop_time[[which(grepl(clonedf[j,2], names(df_pop_time)))]]$clone_id)
    
    clones_subset <- count_clones[which(count_clones$clone_id %in% as.character(inter)),]
    
    lenintersects <- c(lenintersects, length(inter))
    seqsintersects <- c(seqsintersects, sum(clones_subset$seq_count))
  }
  
  clonedf$value <- lenintersects
  seqdf$value <- seqsintersects
  

  grid.col = c("#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#1f78b4", "#33a02c", "#e31a1c", "#ff7f00","#1f78b4", "#33a02c", "#e31a1c", "#ff7f00")
  names(grid.col) <- append(levels(clonedf$from), levels(as.factor(twelve_months$names)))
  
  
  # Clone overlap plot
  svg(filename = paste(patdir_overlap,"/Clone_overlap_comparison_time_points_population",
                       df_pop_time[[1]]$treatment[1], "_",
                       df_pop_time[[1]]$source[1], ".svg", sep=""))
  chordDiagram(clonedf,
               grid.col = grid.col,
               self.link = 1,
               transparency = 0.3,
               annotationTrack="grid",
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(clonedf))))))
  circos.track(track.index = 1,
               panel.fun = function(x, y) {
                 circos.text(CELL_META$xcenter,
                             CELL_META$ylim[2],
                             CELL_META$sector.index,
                             adj = c(0, 0.5))
               }, bg.border = NA)
  title(paste("clone overlap", df_pop_time[[1]]$treatment[1], df_pop_time[[1]]$source[1]), cex = 0.8)
  circos.clear()
  dev.off()
  
  png(filename = paste(patdir_overlap, "/Clone_overlap_comparison_time_points_population", df_pop_time[[1]]$treatment[1], "_", df_pop_time[[1]]$source[1], ".png", sep=""), width=15, height=15, units = "cm", res = 300)
  chordDiagram(clonedf,
               grid.col = grid.col,
               self.link = 1,
               transparency = 0.3,
               annotationTrack="grid",
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(clonedf))))))
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[2], CELL_META$sector.index,
                adj = c(0, 0.5))
  }, bg.border = NA)
  title(paste("clone_id OVERLAP", df_pop_time[[1]]$treatment[1], df_pop_time[[1]]$source[1]), cex = 0.8)
  circos.clear()
  dev.off()
  
  # Sequences overlap plot
  svg(filename = paste(patdir_overlap,"/Clone_seqN_overlap_comparison_time_points_population", df_pop_time[[1]]$treatment[1], "_", df_pop_time[[1]]$source[1], ".svg", sep=""))
  chordDiagram(seqdf,
               grid.col = grid.col,
               self.link = 1,
               transparency = 0.3,
               annotationTrack="grid",
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(seqdf))))))
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[2], CELL_META$sector.index,
                adj = c(0, 0.5))
  }, bg.border = NA)
  title(paste("CLONE SEQ NUM OVERLAP", df_pop_time[[1]]$treatment[1], df_pop_time[[1]]$source[1]), cex = 0.8)
  circos.clear()
  dev.off()
  
  png(filename = paste(patdir_overlap, "/Clone_seqN_overlap_comparison_time_points_population", df_pop_time[[1]]$treatment[1], "_", df_pop_time[[1]]$source[1], ".png", sep=""), width=15, height=15, units = "cm", res = 300)
  chordDiagram(seqdf,
               grid.col = grid.col,
               self.link = 1,
               transparency = 0.3,
               annotationTrack="grid",
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(seqdf))))))
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[2], CELL_META$sector.index,
                adj = c(0, 0.5))
    }, bg.border = NA)
  title(paste("CLONE SEQ NUM OVERLAP", df_pop_time[[1]]$treatment[1], df_pop_time[[1]]$source[1]), cex = 0.8)
  circos.clear()
  dev.off()
  
  # Upset plots
  
  listinput <- list()
  for (k in c(1:length(df_pop_time))){
    listinput <- append(listinput, list(df_pop_time[[k]]$clone_id))
  }    
  names(listinput) <- names(df_pop_time)
  
  
  
  # Filtering out only intersects between time points (didnt work)
  interlist <- list()
  for (j in c(1:nrow(clonedf))){
    inter <- list(as.character(clonedf[j,1]), as.character(clonedf[j,2]))
    interlist <- append(interlist,list(inter))
  }
  
  svg(filename = paste(patdir_overlap,"/Set_plot_population", df_pop_time[[1]]$treatment[1], "_",df_pop_time[[1]]$source[1], ".svg", sep="" ),
      width = 60*0.4, height=18*0.4)
  print(upset( fromList(listinput), 
               nsets = length(listinput), 
               nintersects = NA,
               keep.order=TRUE,
               order.by = "freq",
               #group.by = "sets",
               point.size = 3.5, 
               line.size=2, 
               mainbar.y.label = "Clone intersections", 
               sets.x.label = "Clones per population"))
  
  dev.off()
  
  pdf(file = paste(patdir_overlap,"/Set_plot_population", df_pop_time[[1]]$treatment[1], "_",df_pop_time[[1]]$source[1], ".pdf", sep="" ),
      width=60*0.4, height=18*0.4)
  print(upset( fromList(listinput), 
               nsets = length(listinput), 
               nintersects = NA,
               keep.order=TRUE,
               order.by = "freq",
               #group.by = "sets",
               point.size = 3.5, 
               line.size=2, 
               mainbar.y.label = "Clone intersections", 
               sets.x.label = "Clones per population"))
  dev.off()
  
  png(filename = paste(patdir_overlap,"/Set_plot_population", df_pop_time[[1]]$treatment[1],"_",df_pop_time[[1]]$source[1], ".png", sep=""), 
      res = 600, width = 60, height=18, units="cm")
  print(upset( fromList(listinput), 
               nsets = length(listinput),
               nintersects = NA,
               keep.order=TRUE,
               order.by = "freq",
               #group.by = "sets",
               point.size = 3.5,
               line.size=2,
               mainbar.y.label = "Clone intersections",
               sets.x.label = "Clones per population"))
  dev.off()
  
  png(filename = paste(patdir_overlap,"/Set_plot_population", df_pop_time[[1]]$treatment[1],"_",df_pop_time[[1]]$source[1], "_60intersects.png", sep=""), 
      res = 600, width = 30, height=18, units="cm")
  print(upset( fromList(listinput), 
               nsets = length(listinput),
               nintersects = 60,
               keep.order=TRUE,
               order.by = "freq",
               #group.by = "sets",
               point.size = 3.5,
               line.size=2,
               mainbar.y.label = "Clone intersections",
               sets.x.label = "Clones per population"))
  dev.off()
  
  
  # Saving both tables
  clonedf$source <- rep(df_pop_time[[1]]$source[1], nrow(clonedf))
  seqdf$source <- rep(df_pop_time[[1]]$source[1], nrow(seqdf))
  
  write.table(clonedf, file = paste(patdir_overlap,"/Clone_overlap_comparison_time_points_population", df_pop_time[[1]]$treatment[1], "_", df_pop_time[[1]]$source[1], ".tsv", sep=""), sep = "\t", quote = F, row.names = F)
  write.table(seqdf, file = paste(patdir_overlap,"/Clone_seqN_overlap_comparison_time_points_population", df_pop_time[[1]]$treatment[1], "_", df_pop_time[[1]]$source[1], ".tsv", sep=""), sep = "\t", quote = F, row.names = F)
  
  all_clonedf <- rbind(all_clonedf, clonedf)
  all_seqdf <- rbind(all_seqdf, seqdf)
  
  ######################################
  # Chordplots for Memory population
  ######################################
  
  combin <- expand.grid(unique(df_pat[[i]]$extract_time), unique(df_pat[[i]]$population))
  combin$names <- apply(combin[,c("Var1", "Var2")], 1, paste, collapse = "_")
  
  other <- subset(combin, combin$Var1 != "baseline")
  others = c("baseline_N","baseline_P","baseline_DN",levels(as.factor(other$names)))
  
  clonedf <- expand.grid(c("baseline_M","6months_M"), others)
  clonedf <- subset(clonedf, !(clonedf$Var1 == "6months_M" & clonedf$Var2 == "6months_M"))
  colnames(clonedf) <- c("from","to")
  seqdf <- clonedf
  
  lenintersects = numeric(0)
  seqsintersects = numeric(0)
  for (j in c(1:nrow(clonedf))){
    inter <- intersect(df_pop_time[[which(grepl(clonedf[j,1], names(df_pop_time)))]]$clone_id,
                       df_pop_time[[which(grepl(clonedf[j,2], names(df_pop_time)))]]$clone_id)
    
    clones_subset <- count_clones[which(count_clones$clone_id %in% as.character(inter)),]
    
    lenintersects <- c(lenintersects, length(inter))
    seqsintersects <- c(seqsintersects, sum(clones_subset$seq_count))
  }
  
  
  clonedf$value <- lenintersects
  seqdf$value <- seqsintersects
  
  
  #grid.col = c("#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#1f78b4", "#33a02c", "#e31a1c", "#ff7f00","#1f78b4", "#33a02c", "#e31a1c", "#ff7f00")
  #names(grid.col) <- append(levels(clonedf$from), levels(as.factor(twelve_months$names)))
  
  
  # Clone overlap plot
  svg(filename = paste(patdir_overlap,"/Clone_overlap_comparison_time_points_memory_population",
                       df_pop_time[[1]]$treatment[1], "_",
                       df_pop_time[[1]]$source[1], ".svg", sep=""))
  chordDiagram(clonedf,
               grid.col = grid.col,
               self.link = 1,
               transparency = 0.3,
               annotationTrack="grid",
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(clonedf))))))
  circos.track(track.index = 1,
               panel.fun = function(x, y) {
                 circos.text(CELL_META$xcenter,
                             CELL_META$ylim[2],
                             CELL_META$sector.index,
                             adj = c(0, 0.5))
               }, bg.border = NA)
  title(paste("clone overlap", df_pop_time[[1]]$treatment[1], df_pop_time[[1]]$source[1]), cex = 0.8)
  circos.clear()
  dev.off()
  
  png(filename = paste(patdir_overlap, "/Clone_overlap_comparison_time_points_memory_population", df_pop_time[[1]]$treatment[1], "_", df_pop_time[[1]]$source[1], ".png", sep=""), width=15, height=15, units = "cm", res = 300)
  chordDiagram(clonedf,
               grid.col = grid.col,
               self.link = 1,
               transparency = 0.3,
               annotationTrack="grid",
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(clonedf))))))
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[2], CELL_META$sector.index,
                adj = c(0, 0.5))
  }, bg.border = NA)
  title(paste("clone_id OVERLAP", df_pop_time[[1]]$treatment[1], df_pop_time[[1]]$source[1]), cex = 0.8)
  circos.clear()
  dev.off()
  
  # Sequences overlap plot
  svg(filename = paste(patdir_overlap,"/Clone_seqN_overlap_comparison_time_points_memory_population", df_pop_time[[1]]$treatment[1], "_", df_pop_time[[1]]$source[1], ".svg", sep=""))
  chordDiagram(seqdf,
               grid.col = grid.col,
               self.link = 1,
               transparency = 0.3,
               annotationTrack="grid",
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(seqdf))))))
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[2], CELL_META$sector.index,
                adj = c(0, 0.5))
  }, bg.border = NA)
  title(paste("CLONE SEQ NUM OVERLAP", df_pop_time[[1]]$treatment[1], df_pop_time[[1]]$source[1]), cex = 0.8)
  circos.clear()
  dev.off()
  
  png(filename = paste(patdir_overlap, "/Clone_seqN_overlap_comparison_time_points_memory_population", df_pop_time[[1]]$treatment[1], "_", df_pop_time[[1]]$source[1], ".png", sep=""), width=15, height=15, units = "cm", res = 300)
  chordDiagram(seqdf,
               grid.col = grid.col,
               self.link = 1,
               transparency = 0.3,
               annotationTrack="grid",
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(seqdf))))))
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[2], CELL_META$sector.index,
                adj = c(0, 0.5))
  }, bg.border = NA)
  title(paste("CLONE SEQ NUM OVERLAP", df_pop_time[[1]]$treatment[1], df_pop_time[[1]]$source[1]), cex = 0.8)
  circos.clear()
  dev.off()
  

  ###############################################################
  # Chordplots overlap populations for each time point separately
  ################################################################

  df_TP <- split(df_pat[[i]], df_pat[[i]]$extract_time)

  # Plots per patient and time point - overlap populations
  for (n in c(1:length(df_TP))) {
      df_pop <- split(df_TP[[n]], df_TP[[n]]$population)
      vennplot <- venn(list(unique(df_pop[[1]]$clone_id), unique(df_pop[[2]]$clone_id), unique(df_pop[[3]]$clone_id), unique(df_pop[[4]]$clone_id)), names = names(df_pop))

      listInput <- list(df_pop[[1]]$clone_id, df_pop[[2]]$clone_id, df_pop[[3]]$clone_id, df_pop[[4]]$clone_id)
      names(listInput) <- names(df_pop)
    

      # Upset plots
      pdf(file = paste(patdir_overlap,"/Set_plot_time_points", df_pop[[1]]$treatment[1], "_",df_pop[[1]]$extract_time[1], "_",df_pop[[1]]$source[1], ".pdf", sep="" ))
      print(upset(fromList(listInput), 
                  group.by = "sets", 
                  order.by="freq", 
                  point.size = 3.5, 
                  line.size=2, 
                  mainbar.y.label = "Clone intersections", 
                  sets.x.label = "Clones per population"))
      dev.off()

      png(filename = paste(patdir_overlap,"/Set_plot_time_points", df_pop[[1]]$treatment[1], "_",df_pop[[1]]$extract_time[1], "_",df_pop[[1]]$source[1], ".png", sep=""), res = 600, width = 15, height=10, units = "cm")
      print(upset(fromList(listInput), 
                  order.by="freq", 
                  group.by = "sets", 
                  point.size = 3.5, 
                  line.size=2, 
                  mainbar.y.label = "Clone intersections", 
                  sets.x.label = "Clones per population"))
      dev.off()

      pops_vec = c("N", "DN", "M", "P")
      combin <- data.frame(from=combinations(4,2,pops_vec,repeats.allowed=F)[,1], to=combinations(4,2,pops_vec,repeats.allowed=F)[,2])
      
      clonedf_pop_tp <- combin
      seqdf_pop_tp <- combin

      lenintersects = numeric(0)
      seqsintersects = numeric(0)
      for (j in c(1:nrow(clonedf_pop_tp))){
        inter <- intersect(df_pop[[which(grepl(paste0("^",clonedf_pop_tp[j,1]), names(df_pop)))]]$clone_id,
                           df_pop[[which(grepl(paste0("^",clonedf_pop_tp[j,2]), names(df_pop)))]]$clone_id)

        clones_subset <- count_clones[which(count_clones$clone_id %in% as.character(inter)),]

        lenintersects <- c(lenintersects, length(inter))
        seqsintersects <- c(seqsintersects, sum(clones_subset$seq_count))
      }

      clonedf_pop_tp$value <- lenintersects
      seqdf_pop_tp$value <- seqsintersects


      self_comb_pop_tp <- data.frame(from = names(df_pop), to = names(df_pop))
      self_clonedf_pop_tp <- self_comb_pop_tp
      self_seqdf_pop_tp <- self_comb_pop_tp

      lenintersects <- numeric(0)
      seqsintersects <- numeric(0)
      for (pop in self_comb_pop_tp$from){
        inter <- attributes(vennplot)[["intersections"]][[pop]]
        clones_subset <- count_clones[which(count_clones$clone_id %in% as.character(inter)),]

        lenintersects <- c(lenintersects, length(inter))
        seqsintersects <- c(seqsintersects, sum(clones_subset$seq_count))
      }
      self_clonedf_pop_tp$value <- lenintersects
      self_seqdf_pop_tp$value <- seqsintersects

      clonedf_pop_tp <- rbind(clonedf_pop_tp, self_clonedf_pop_tp)
      seqdf_pop_tp <- rbind(seqdf_pop_tp, self_seqdf_pop_tp)
      
      grid.col = c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00")
      names(grid.col) = c("DN","M","N","P")


      # Plots clone overlap
      svg(filename = paste(patdir_overlap,"/Clone_overlap_comparison_population_", df_pop[[1]]$treatment[1], "_", df_pop[[1]]$extract_time[1], "_", df_pop[[1]]$source[1], ".svg", sep=""))
      chordDiagram(clonedf_pop_tp,
                   grid.col = grid.col,
                   self.link = 1,
                   transparency = 0.3,
                   annotationTrack="grid",
                   preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(clonedf_pop_tp))))))
      circos.track(track.index = 1, panel.fun = function(x, y) {
                    circos.text(CELL_META$xcenter, CELL_META$ylim[2], CELL_META$sector.index,
                    adj = c(0, 0.5))
                    }, bg.border = NA)
      title(paste("clone_id OVERLAP", df_pop[[1]]$treatment[1], df_pop[[1]]$source[1], df_pop[[1]]$extract_time[1]), cex = 0.8)
      circos.clear()
      dev.off()

      png(filename = paste(patdir_overlap,"/Clone_overlap_comparison_population_", df_pop[[1]]$treatment[1], "_", df_pop[[1]]$extract_time[1], "_", df_pop[[1]]$source[1], ".png", sep=""), res = 600, width = 15, height=10, units = "cm")
      chordDiagram(clonedf_pop_tp,
                   grid.col = grid.col,
                   self.link = 1,
                   transparency = 0.3,
                   annotationTrack="grid",
                   preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(clonedf_pop_tp))))))
      circos.track(track.index = 1, panel.fun = function(x, y) {
        circos.text(CELL_META$xcenter, CELL_META$ylim[2], CELL_META$sector.index,
                    adj = c(0, 0.5))
      }, bg.border = NA)
      title(paste("clone_id OVERLAP", df_pop[[1]]$treatment[1], df_pop[[1]]$source[1], df_pop[[1]]$extract_time[1]), cex = 0.8)
      circos.clear()
      dev.off()

      # Plots clone sequence numbers overlap
      svg(filename = paste(patdir_overlap,"/Clone_seqN_overlap_comparison_population_", df_pop[[1]]$treatment[1], "_", df_pop[[1]]$extract_time[1], "_", df_pop[[1]]$source[1], ".svg", sep=""))
      chordDiagram(seqdf_pop_tp,
                   grid.col = grid.col,
                   self.link = 1,
                   transparency = 0.3,
                   annotationTrack="grid",
                   preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(seqdf_pop_tp))))))
      circos.track(track.index = 1, panel.fun = function(x, y) {
        circos.text(CELL_META$xcenter, CELL_META$ylim[2], CELL_META$sector.index,
                    adj = c(0, 0.5))
      }, bg.border = NA)
      title(paste("CLONE SEQ NUM OVERLAP", df_pop[[1]]$treatment[1], df_pop[[1]]$source[1], df_pop[[1]]$extract_time[1]), cex = 0.8)
      circos.clear()
      dev.off()

      png(filename = paste(patdir_overlap,"/Clone_seqN_overlap_comparison_population_", df_pop[[1]]$treatment[1], "_", df_pop[[1]]$extract_time[1], "_", df_pop[[1]]$source[1], ".png", sep=""), res = 600, width = 15, height=10, units = "cm")
      chordDiagram(seqdf_pop_tp,
                   grid.col = grid.col,
                   self.link = 1,
                   transparency = 0.3,
                   annotationTrack="grid",
                   preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(seqdf_pop_tp))))))
      circos.track(track.index = 1, panel.fun = function(x, y) {
        circos.text(CELL_META$xcenter, CELL_META$ylim[2], CELL_META$sector.index,
                    adj = c(0, 0.5))
      }, bg.border = NA)
      title(paste("CLONE SEQ NUM OVERLAP", df_pop[[1]]$treatment[1], df_pop[[1]]$source[1], df_pop[[1]]$extract_time[1]), cex = 0.8)
      circos.clear()
      dev.off()
      
      seqdf_pop_tp$time_point <- rep(df_pop[[1]]$extract_time[1], nrow(seqdf_pop_tp))
      seqdf_pop_tp$patient <- rep(df_pop[[1]]$source[1], nrow(seqdf_pop_tp))
      clonedf_pop_tp$time_point <- rep(df_pop[[1]]$extract_time[1], nrow(clonedf_pop_tp))
      clonedf_pop_tp$patient <- rep(df_pop[[1]]$source[1], nrow(clonedf_pop_tp))
      
      all_clonedf_pop_tp <- rbind(all_clonedf_pop_tp, clonedf_pop_tp)
      all_seqdf_pop_tp <- rbind(all_seqdf_pop_tp, seqdf_pop_tp)
  }
  
}

write.table(all_clonedf_pop_tp, file = paste(outdir,"/Clone_overlap","/Clone_overlap_comparison_population_individual_tp.tsv", sep=""), sep = "\t", quote = F, row.names = F)
write.table(all_seqdf_pop_tp, file = paste(outdir,"/Clone_overlap","/Clone_seqN_overlap_comparison_population_individual_tp.tsv", sep=""), sep = "\t", quote = F, row.names = F)
write.table(all_clonedf, file = paste(outdir,"/Clone_overlap", "/Clone_overlap_comparison_time_points_population.tsv", sep=""), sep = "\t", quote = F, row.names = F)
write.table(all_seqdf, file = paste(outdir,"/Clone_overlap","/Clone_seqN_overlap_comparison_time_points_population.tsv", sep=""), sep = "\t", quote = F, row.names = F)
write.table(all_clonedf_tp, file = paste(outdir,"/Clone_overlap", "/Clone_overlap_comparison_time_points.tsv", sep=""), sep = "\t", quote = F, row.names = F)
write.table(all_seqdf_tp, file = paste(outdir,"/Clone_overlap","/Clone_seqN_overlap_comparison_time_points.tsv", sep=""), sep = "\t", quote = F, row.names = F)


all_clonedf_tp$fromto <- paste(all_clonedf_tp$from,"_",all_clonedf_tp$to, sep = "")
all_clonedf_tp <- all_clonedf_tp[!(all_clonedf_tp$fromto %in% c("baseline_baseline","6months_6months","12months_12months")),]
all_clonedf_tp$fromto <- factor(all_clonedf_tp$fromto, levels=c("baseline_6months","6months_12months","baseline_12months"))
ggplot(all_clonedf_tp, aes(x=fromto, y=value)) +
  geom_boxplot() +
  geom_jitter(aes(color=source), width = 0.2) +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
ggsave(paste(outdir,"/Clone_overlap/Clone_overlap_boxplot_allvalues.png", sep = ""), device="png",
       width = 12, height = 8, units = "cm")

all_clonedf_tp_subset <- all_clonedf_tp[,c("fromto", "source", "value")]
all_clonedf_tp_transformed <- dcast(data=as.data.table(all_clonedf_tp_subset), formula=source~fromto, fun.aggregate=sum, value.var="value")
all_clonedf_tp_transformed$delta_6 <- all_clonedf_tp_transformed$`6months_12months` / all_clonedf_tp_transformed$baseline_6months

ggplot(all_clonedf_tp_transformed, aes(x=source,y=delta_6)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = 1, color="blue") +
  ggtitle("overlaps 6-12mo / overlaps B-6mo") +
  ylab("Overlap number difference") + xlab("") 
ggsave(paste(outdir,"/Clone_overlap/Clone_overlap_ratio_between_6months12months_and_6monthsbaseline.png",sep=""), device="png",
       width = 12, height = 8, units = "cm")

write.table(all_clonedf_tp_transformed, file=paste(outdir,"/Clone_overlap/Clone_overlap_comparison_timepoints_unmelted.tsv", sep = ""), sep = "\t", quote=F, row.names = F)

