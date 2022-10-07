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
library(cowplot)
library(ggpubr)

theme_set(theme_bw(base_family = "ArialMT") +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(family="ArialMT")))

# Run inside improved_reports/repertoire_per_patient
outdir <- "repertoire_comparison_bootstrapped"
dir.create(outdir)

## Reading clonal abundance table, which has been bootstrapped.
df_ab_pat <- read.csv("repertoire_comparison/Abundance/Clonal_abundance_data_cell_population.tsv", sep = "\t")

# Remove patient CLAD4 as it is missing a time point
df_subset <- df_ab_pat[!(df_ab_pat$patient == "CLAD4"),]


df_pat <- split(df_subset, df_subset$patient)
df_pat$CLAD4<-NULL
dir.create(paste(outdir,"Clone_overlap", sep = "/"))

all_clonedf_tp <- data.frame()

for (i in c(1:length(df_pat))) {
  patdir_overlap <- paste(outdir,"Clone_overlap",names(df_pat)[i], sep="/")
  dir.create(patdir_overlap)

  ###################################################################
  ## Plot chordplot comparison between time points
  ###################################################################
  df_time <- split(df_pat[[i]], df_pat[[i]]$time_point)

  # Need to calculate vennplot to be able to count the clones without overlaps to other time points
  vennplot <- venn(list(unique(df_time[[1]]$clone_id), unique(df_time[[2]]$clone_id), unique(df_time[[3]]$clone_id)), names = names(df_time))
  
  # Create dataframe with desired combinations of overlaps from and to
  combin_tp <- expand.grid(c("baseline"), c("6months","12months"))
  clonedf_tp <- rbind(combin_tp,
                      data.frame(Var1=c("6months"), Var2=c("12months")))
  names(clonedf_tp) <- c("from","to")


  lenintersects = numeric(0)
  for (j in c(1:nrow(clonedf_tp))){
    # intersect (inclusive) of clone IDs between the two time points.
    inter <- base::intersect(df_time[[as.character(clonedf_tp[j,1])]]$clone_id,
                        df_time[[as.character(clonedf_tp[j,2])]]$clone_id)
    
    # Getting those clones
    clones_subset <- df_pat[[i]][which(df_pat[[i]]$clone_id %in% as.character(inter)),]
  
    # Number of clones on those 2 time points
    lenintersects <- c(lenintersects, length(inter))

  }

  clonedf_tp$value <- lenintersects

  # Count the clones without overlaps to other time points
  self_comb <- data.frame(from = names(df_time), to = names(df_time))
  self_clonedf_tp <- self_comb

  lenintersects <- numeric(0)
  for (tp in self_comb$from){
    inter <- attributes(vennplot)[["intersections"]][[tp]]
    clones_subset <- df_pat[[i]][which(df_pat[[i]]$clone_id %in% as.character(inter)),]
    
    lenintersects <- c(lenintersects, length(inter))

  }
  self_clonedf_tp$value <- lenintersects

  clonedf_tp <- rbind(clonedf_tp, self_clonedf_tp)


  # Plot chordplots between time points
  grid.col = c("#a6cee3", "#b2df8a", "#fb9a99")
  names(grid.col) <- names(df_time)


  # Clone overlap plot
  svg(filename = paste(patdir_overlap,"/Clone_overlap_comparison_time_points_",
                        df_time[[1]]$treatment[1], "_",
                        df_time[[1]]$patient[1], ".svg", sep=""))
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
  title(paste("clone overlap", df_time[[1]]$treatment[1], df_time[[1]]$patient[1]), cex = 0.8)
  circos.clear()
  dev.off()

  png(filename = paste(patdir_overlap, "/Clone_overlap_comparison_time_points_",
                        df_time[[1]]$treatment[1], "_", df_time[[1]]$patient[1], ".png", sep=""),
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
  title(paste("clone_id OVERLAP", df_time[[1]]$treatment[1], df_time[[1]]$patient[1]), cex = 0.8)
  circos.clear()
  dev.off()

  # Upset plots per time point

  listinput <- list()
  for (k in c(1:length(df_time))){
    listinput <- append(listinput, list(df_time[[k]]$clone_id))
  }
  names(listinput) <- names(df_time)


  svg(filename = paste(patdir_overlap,"/Set_plot_timepoint", df_time[[1]]$treatment[1], "_",df_time[[1]]$patient[1], ".svg", sep="" ),
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

  png(filename = paste(patdir_overlap,"/Set_plot_timepoint", df_time[[1]]$treatment[1], "_",df_time[[1]]$patient[1], ".png", sep="" ),
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
  clonedf_tp$patient <- rep(df_time[[1]]$patient[1], nrow(clonedf_tp))

  # Saving both tables
  write.table(clonedf_tp, file = paste(patdir_overlap,"/Clone_overlap_comparison_time_points_", df_time[[1]]$treatment[1], "_", df_time[[1]]$patient[1], ".tsv", sep=""), sep = "\t", quote = F, row.names = F)
  all_clonedf_tp <- rbind(all_clonedf_tp, clonedf_tp)
}

write.table(all_clonedf_tp, file = paste(outdir,"/Clone_overlap", "/Clone_overlap_comparison_time_points.tsv", sep=""), sep = "\t", quote = F, row.names = F)


df_ab_pop <- read.csv("./repertoire_comparison/Abundance/Clonal_abundance_data_cell_population.tsv", sep = "\t")

# Remove patient CLAD4 as it is missing a time point
df_subset <- df_ab_pop[!(df_ab_pop$patient == "CLAD4"),]


df_pat <- split(df_subset, df_subset$patient)
df_pat$CLAD4<-NULL
dir.create(paste(outdir,"Clone_overlap", sep = "/"))


all_clonedf <- data.frame()
all_clonedf_m <- data.frame()

for (i in c(1:length(df_pat))) {
  patdir_overlap <- paste(outdir,"Clone_overlap",names(df_pat)[i], sep="/")
  ## Plot chordplot comparison time points per patient and population
  df_pat[[i]]$time_pop <- as.factor(paste(df_pat[[i]]$time_point, df_pat[[i]]$population, sep="_"))
  df_pop_time <- split(df_pat[[i]], df_pat[[i]]$time_pop)


  ## Calculating overlaps between time points per patient and population
  combin <- expand.grid(unique(df_pat[[i]]$time_point), unique(df_pat[[i]]$population))
  combin$names <- apply(combin[,c("Var1", "Var2")], 1, paste, collapse = "_")

  baselines <- subset(combin, combin$Var1 == "baseline")
  six_months <- subset(combin, combin$Var1 == "6months")
  other <- subset(combin, combin$Var1 != "baseline")
  twelve_months <- subset(combin, combin$Var1 == "12months")

  clonedfb <- expand.grid(levels(as.factor(baselines$names)), levels(as.factor(other$names)))
  clonedf612 <- expand.grid(levels(as.factor(six_months$names)), levels(as.factor(twelve_months$names)))
  clonedf <- rbind(clonedfb, clonedf612)
  colnames(clonedf) <- c("from","to")


  lenintersects = numeric(0)
  for (j in c(1:nrow(clonedf))){
    inter <- intersect(df_pop_time[[which(grepl(clonedf[j,1], names(df_pop_time)))]]$clone_id,
                        df_pop_time[[which(grepl(clonedf[j,2], names(df_pop_time)))]]$clone_id)

    clones_subset <- clones_subset <- df_pat[[i]][which(df_pat[[i]]$clone_id %in% as.character(inter)),]

    lenintersects <- c(lenintersects, length(inter))
  }

  clonedf$value <- lenintersects


  grid.col = c("#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#1f78b4", "#33a02c", "#e31a1c", "#ff7f00","#1f78b4", "#33a02c", "#e31a1c", "#ff7f00")
  names(grid.col) <- append(levels(clonedf$from), levels(as.factor(twelve_months$names)))


  # Clone overlap plot
  svg(filename = paste(patdir_overlap,"/Clone_overlap_comparison_time_points_population",
                        df_pop_time[[1]]$treatment[1], "_",
                        df_pop_time[[1]]$patient[1], ".svg", sep=""))
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
  title(paste("clone overlap", df_pop_time[[1]]$treatment[1], df_pop_time[[1]]$patient[1]), cex = 0.8)
  circos.clear()
  dev.off()

  png(filename = paste(patdir_overlap, "/Clone_overlap_comparison_time_points_population", df_pop_time[[1]]$treatment[1], "_", df_pop_time[[1]]$patient[1], ".png", sep=""), width=15, height=15, units = "cm", res = 300)
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
  title(paste("clone_id OVERLAP", df_pop_time[[1]]$treatment[1], df_pop_time[[1]]$patient[1]), cex = 0.8)
  circos.clear()
  dev.off()


  # Chordplots for Memory population
  
  clonedf_m <- clonedf
  clonedf_m <- clonedf_m[which(clonedf_m$from %in% c("baseline_M","6months_M")),]
  clonedf_m <- clonedf_m[which(clonedf_m$to %in% c("6months_M","12months_M")),]
  
  self_comb_m <- data.frame(from = c("baseline_M","6months_M","12months_M"), to = c("baseline_M","6months_M","12months_M"))
  self_clonedf_m <- self_comb_m
  
  baseline_M <- unique(df_pop_time$baseline_M$clone_id)
  six_months_M <- unique(df_pop_time$`6months_M`$clone_id)
  twelve_months_M <- unique(df_pop_time$`12months_M`$clone_id)
  
  self_clonedf_m$value <- c(
  length(setdiff(baseline_M, union(six_months_M, twelve_months_M))),
  length(setdiff(six_months_M, union(baseline_M, twelve_months_M))),
  length(setdiff(twelve_months_M, union(six_months_M, baseline_M)))
)
  
  clonedf_m <- rbind(clonedf_m, self_clonedf_m)
  
  grid.col = c("#a6cee3", "#b2df8a", "#fb9a99")
  names(grid.col) <- c("12months_M","6months_M","baseline_M")


  # Clone overlap plot
  svg(filename = paste(patdir_overlap,"/Clone_overlap_comparison_time_points_memory_population",
                        df_pop_time[[1]]$treatment[1], "_",
                        df_pop_time[[1]]$patient[1], ".svg", sep=""))
  chordDiagram(clonedf_m,
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
  title(paste("clone overlap", df_pop_time[[1]]$treatment[1], df_pop_time[[1]]$patient[1]), cex = 0.8)
  circos.clear()
  dev.off()

  png(filename = paste(patdir_overlap, "/Clone_overlap_comparison_time_points_memory_population", df_pop_time[[1]]$treatment[1], "_", df_pop_time[[1]]$patient[1], ".png", sep=""), width=15, height=15, units = "cm", res = 300)
  chordDiagram(clonedf_m,
                grid.col = grid.col,
                self.link = 1,
                transparency = 0.3,
                annotationTrack="grid",
                preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(clonedf))))))
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[2], CELL_META$sector.index,
                adj = c(0, 0.5))
  }, bg.border = NA)
  title(paste("clone_id OVERLAP", df_pop_time[[1]]$treatment[1], df_pop_time[[1]]$patient[1]), cex = 0.8)
  circos.clear()
  dev.off()


  # Saving both tables
  clonedf$patient <- rep(df_pop_time[[1]]$patient[1], nrow(clonedf))
  clonedf_m$patient  <- rep(df_pop_time[[1]]$patient[1], nrow(clonedf_m))

  write.table(clonedf, file = paste(patdir_overlap,"/Clone_overlap_comparison_time_points_population", df_pop_time[[1]]$treatment[1], "_", df_pop_time[[1]]$patient[1], ".tsv", sep=""), sep = "\t", quote = F, row.names = F)


  # Upset plots

  listinput_m <- list()
  for (k in c("baseline_M","6months_M","12months_M")){
    listinput_m <- append(listinput_m, list(df_pop_time[[k]]$clone_id))
  }
  names(listinput_m) <- c("baseline_M","6months_M","12months_M")


  svg(filename = paste(patdir_overlap,"/Set_plot_population", df_pop_time[[1]]$treatment[1], "_",df_pop_time[[1]]$patient[1], ".svg", sep="" ),
      width = 20*0.4, height=10*0.4)
  print(upset( fromList(listinput_m),
                nsets = length(listinput_m),
                nintersects = NA,
                keep.order=TRUE,
                order.by = "freq",
                #group.by = "sets",
                point.size = 3.5,
                line.size=2,
                mainbar.y.label = "Clone intersections",
                sets.x.label = "Clones per population"))

  dev.off()

  pdf(file = paste(patdir_overlap,"/Set_plot_population", df_pop_time[[1]]$treatment[1], "_",df_pop_time[[1]]$patient[1], ".pdf", sep="" ),
      width=20*0.4, height=10*0.4)
  print(upset( fromList(listinput_m),
                nsets = length(listinput_m),
                nintersects = NA,
                keep.order=TRUE,
                order.by = "freq",
                #group.by = "sets",
                point.size = 3.5,
                line.size=2,
                mainbar.y.label = "Clone intersections",
                sets.x.label = "Clones per population"))
  dev.off()

  png(filename = paste(patdir_overlap,"/Set_plot_population", df_pop_time[[1]]$treatment[1],"_",df_pop_time[[1]]$patient[1], ".png", sep=""),
      res=600,width = 15, height=10, units="cm")
  print(upset( fromList(listinput_m),
                nsets = length(listinput_m),
                nintersects = NA,
                keep.order=TRUE,
                order.by = "freq",
                #group.by = "sets",
                point.size = 3.5,
                line.size=2,
                mainbar.y.label = "Clone intersections",
                sets.x.label = "Clones per population"))
  dev.off()

  all_clonedf <- rbind(all_clonedf, clonedf)
  all_clonedf_m <- rbind(all_clonedf_m, clonedf_m)
}

write.table(all_clonedf, file = paste(outdir,"/Clone_overlap", "/Clone_overlap_comparison_time_points_population.tsv", sep=""), sep = "\t", quote = F, row.names = F)
write.table(all_clonedf_m, file = paste(outdir,"/Clone_overlap", "/Clone_overlap_comparison_time_points_memory_population.tsv", sep=""), sep = "\t", quote = F, row.names = F)

# Boxplot overlaps comparison

all_clonedf_tp$fromto <- paste(all_clonedf_tp$from,"_",all_clonedf_tp$to, sep = "")
all_clonedf_tp <- all_clonedf_tp[!(all_clonedf_tp$fromto %in% c("baseline_baseline","6months_6months","12months_12months")),]
all_clonedf_tp$fromto <- factor(all_clonedf_tp$fromto, levels=c("baseline_6months","6months_12months","baseline_12months"))

my_comparisons <- list( c("baseline_6months", "6months_12months"), c("6months_12months", "baseline_12months"), c("baseline_6months", "baseline_12months") )
boxplot_overlaps <- ggpubr::ggboxplot(all_clonedf_tp, x="fromto", y="value",
                                           add = "jitter", legend="none", 
                                          line.color = "gray", line.size = 0.4,
                                           short.panel.labs = T) +
  scale_x_discrete(breaks=c("baseline_6months","6months_12months","baseline_12months"),
                   labels=c("B-6M", "6M-12M", "B-12M")) +
  xlab("") + ylab("Number overlaps") +
  stat_compare_means(comparisons = my_comparisons, paired=T, method = "wilcox", label="p.signif") +
  stat_compare_means(paired=T, label="p.format", label.y = 35000)
boxplot_overlaps

ggsave(plot = boxplot_overlaps, paste(outdir,"/Clone_overlap/Clone_overlap_boxplot_allvalues.png", sep = ""), device="png", width = 10, height = 7, units = "cm")
ggsave(plot = boxplot_overlaps, paste(outdir,"/Clone_overlap/Clone_overlap_boxplot_allvalues.pdf", sep = ""), device="pdf", width = 10, height = 7, units = "cm")

all_clonedf_tp_nob12 <- all_clonedf_tp[!(all_clonedf_tp$fromto == "baseline_12months"),]

my_comparisons <- list( c("baseline_6months", "6months_12months"))
boxplot_overlaps_paired <- ggpubr::ggpaired(all_clonedf_tp_nob12, x="fromto", y="value",
                                      add = "jitter", legend="none", 
                                      line.color = "gray", line.size = 0.4,
                                      short.panel.labs = T) +
  scale_x_discrete(breaks=c("baseline_6months","6months_12months"),
                   labels=c("B-6M", "6M-12M")) +
  xlab("") + ylab("Number overlaps") +
  stat_compare_means(comparisons = my_comparisons, paired=T, method = "wilcox", label="p.signif")
boxplot_overlaps_paired

ggsave(plot = boxplot_overlaps_paired, paste(outdir,"/Clone_overlap/Clone_overlap_boxplot_paired.png", sep = ""), device="png", width = 8, height = 6.5, units = "cm")
ggsave(plot = boxplot_overlaps_paired, paste(outdir,"/Clone_overlap/Clone_overlap_boxplot_paired.pdf", sep = ""), device="pdf", width = 8, height = 6.5, units = "cm")


all_clonedf_tp_subset <- all_clonedf_tp[,c("fromto", "patient", "value")]
all_clonedf_tp_transformed <- dcast(data=as.data.table(all_clonedf_tp_subset), formula=patient~fromto, fun.aggregate=sum, value.var="value")
all_clonedf_tp_transformed$delta_6 <- all_clonedf_tp_transformed$`6months_12months` / all_clonedf_tp_transformed$baseline_6months

ggplot(all_clonedf_tp_transformed, aes(x=patient,y=delta_6)) +
  geom_bar(stat="identity") +
  ggtitle("Increase in number of clone overlaps (# overlaps 6-12mo / # overlaps B-6mo)") +
  ylab("Overlap number difference") + xlab("")
ggsave(paste(outdir,"/Clone_overlap/Clone_overlap_difference_between_6months12months_and_6monthsbaseline.png",sep=""), device="png")

write.table(all_clonedf_tp_transformed, file=paste(outdir,"/Clone_overlap/Clone_overlap_comparison_timepoints_unmelted.tsv", sep = ""), sep = "\t", quote=F, row.names = F)


## Boxplots overlap comparison for memory cells

all_clonedf_m$fromto <- paste(all_clonedf_m$from,"_",all_clonedf_m$to, sep = "")
all_clonedf_m$fromto <- factor(all_clonedf_m$fromto, levels=c("baseline_M_6months_M","6months_M_12months_M","baseline_M_12months_M"))

ggplot(all_clonedf_m, aes(x=fromto, y=value)) +
  geom_boxplot()
  geom_jitter(aes(color=patient), width = 0.2)
ggsave(paste(outdir,"/Clone_overlap/Clone_overlap_boxplot_memory_cells_allvalues.png", sep = ""), device="png")

my_comparisons <- list( c("baseline_M_6months_M", "6months_M_12months_M"), c("6months_M_12months_M", "baseline_M_12months_M"), c("baseline_M_6months_M", "baseline_M_12months_M") )
boxplot_overlaps_m <- ggpubr::ggboxplot(all_clonedf_m, x="fromto", y="value",
                                      add = "jitter", legend="none", 
                                      line.color = "gray", line.size = 0.4,
                                      short.panel.labs = T) +
  scale_x_discrete(breaks=c("baseline_M_6months_M","6months_M_12months_M","baseline_M_12months_M"),
                   labels=c("B-6M", "6M-12M", "B-12M")) +
  xlab("") + ylab("Number overlaps") +
  stat_compare_means(comparisons = my_comparisons, paired=T, method = "wilcox", label="p.signif") +
  stat_compare_means(paired=T, label="p.format", label.y = 4500)
boxplot_overlaps_m

ggsave(plot = boxplot_overlaps_m, paste(outdir,"/Clone_overlap/Clone_overlap_boxplot_memory_allvalues.png", sep = ""), device="png", width = 10, height = 7, units = "cm")
ggsave(plot = boxplot_overlaps_m, paste(outdir,"/Clone_overlap/Clone_overlap_boxplot_memory_allvalues.pdf", sep = ""), device="pdf", width = 10, height = 7, units = "cm")

all_clonedf_m_nob12 <- all_clonedf_m[!(all_clonedf_m$fromto == "baseline_M_12months_M"),]


all_clonedf_m_subset <- all_clonedf_m[,c("fromto", "patient", "value")]
all_clonedf_m_transformed <- dcast(data=as.data.table(all_clonedf_m_subset), formula=patient~fromto, fun.aggregate=sum, value.var="value")
all_clonedf_m_transformed$delta_6 <- all_clonedf_m_transformed$`6months_M_12months_M` / all_clonedf_m_transformed$baseline_M_6months_M

ggplot(all_clonedf_m_transformed, aes(x=patient,y=delta_6)) +
  geom_bar(stat="identity") +
  ggtitle("Increase in number of clone overlaps (# overlaps 6-12mo / # overlaps B-6mo)") +
  ylab("Overlap number difference") + xlab("")
ggsave(paste(outdir,"/Clone_overlap/Clone_overlap_difference_between_6months12months_and_6monthsbaseline_memory.png",sep=""), device="png")

write.table(all_clonedf_m_transformed, file=paste(outdir,"/Clone_overlap/Clone_overlap_comparison_timepoints_unmelted_memory.tsv", sep = ""), sep = "\t", quote=F, row.names = F)

