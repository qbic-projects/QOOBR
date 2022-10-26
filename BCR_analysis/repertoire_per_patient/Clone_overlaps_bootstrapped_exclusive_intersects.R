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
library(ggupset)
library(ggpubr)

traceback()

theme_set(theme_bw(base_family = "ArialMT") +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(family="ArialMT")))

# Run inside repertoire_per_patient
outdir <- "Clone_overlaps_exclusive_intersects"
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

  #
  # # Need to calculate vennplot to be able to count the clones without overlaps to other time points
  # vennplot <- venn(list(unique(df_time[[1]]$clone_id), unique(df_time[[2]]$clone_id), unique(df_time[[3]]$clone_id)), names = names(df_time))
  #
  # # Create dataframe with desired combinations of overlaps from and to
  combin_tp <- expand.grid(c("baseline"), c("6months","12months"))
  clonedf_tp <- rbind(combin_tp,
                      data.frame(Var1=c("6months"), Var2=c("12months")))
  names(clonedf_tp) <- c("from","to")

  self_comb <- data.frame(from = c("baseline","6months","12months"), to = c("baseline","6months","12months"))
  self_clonedf_tp <- self_comb

  clonedf_tp <- rbind(clonedf_tp, self_clonedf_tp)

  baseline <- unique(df_time[["baseline"]]$clone_id)
  sixmonths <- unique(df_time[["6months"]]$clone_id)
  twemonths <- unique(df_time[["12months"]]$clone_id)

  exclusive_intersects = c(
    length(setdiff(intersect(baseline, sixmonths), twemonths)),
    length(setdiff(intersect(baseline, twemonths), sixmonths)),
    length(setdiff(intersect(sixmonths, twemonths), baseline)),
    length(setdiff(baseline, union(sixmonths, twemonths))),
    length(setdiff(sixmonths, union(baseline, twemonths))),
    length(setdiff(twemonths, union(sixmonths, baseline)))
  )

  clonedf_tp$value <- exclusive_intersects

  # Upset plots per time point

  listinput <- list()
  for (k in c(1:length(df_time))){
    listinput <- append(listinput, list(df_time[[k]]$clone_id))
  }
  names(listinput) <- names(df_time)


  ## Adding overlap of 3 time points
  clonedf_tp = clonedf_tp %>% mutate(across(c(from,to), as.character)) %>%
                              mutate(across(c(value), as.numeric))

  overlap_3tp <- length(intersect(baseline, intersect(sixmonths, twemonths)))
  clonedf_tp[nrow(clonedf_tp)+1,] <- c("baseline","6months_12months",overlap_3tp)


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
all_clonedf_n <- data.frame()
all_clonedf_p <- data.frame()
all_clonedf_dn <- data.frame()

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


   # Calculate overlaps for Memory population
  # -------------------------------------------

  clonedf_m <- clonedf
  clonedf_m <- clonedf_m[which(clonedf_m$from %in% c("baseline_M","6months_M")),]
  clonedf_m <- clonedf_m[which(clonedf_m$to %in% c("6months_M","12months_M")),]

  self_comb_m <- data.frame(from = c("baseline_M","6months_M","12months_M"), to = c("baseline_M","6months_M","12months_M"))
  self_clonedf_m <- self_comb_m

  clonedf_m <- rbind(clonedf_m, self_clonedf_m)

  baseline_M <- unique(df_pop_time$baseline_M$clone_id)
  six_months_M <- unique(df_pop_time$`6months_M`$clone_id)
  twelve_months_M <- unique(df_pop_time$`12months_M`$clone_id)

  exclusive_intersects_m = c(
    length(setdiff(intersect(baseline_M, twelve_months_M), six_months_M)),
    length(setdiff(intersect(baseline_M, six_months_M), twelve_months_M)),
    length(setdiff(intersect(six_months_M, twelve_months_M), baseline_M)),
    length(setdiff(baseline_M, union(six_months_M, twelve_months_M))),
    length(setdiff(six_months_M, union(baseline_M, twelve_months_M))),
    length(setdiff(twelve_months_M, union(six_months_M, baseline_M)))
  )

  clonedf_m$value <- exclusive_intersects_m

  grid.col = c("#a6cee3", "#b2df8a", "#fb9a99")
  names(grid.col) <- c("12months_M","6months_M","baseline_M")



  ## Adding overlap of 3 time points
  clonedf_m = clonedf_m %>% mutate(across(c(from,to), as.character)) %>%
                            mutate(across(c(value), as.numeric))

  overlap_3_m <- length(intersect(intersect(baseline_M, twelve_months_M), six_months_M))
  clonedf_m[nrow(clonedf_m)+1,] <- c("baseline_M","6months_M_12months_M",overlap_3_m)

  # Calculate overlaps for Plasmablast population
  # -------------------------------------------

  clonedf_n <- clonedf
  clonedf_n <- clonedf_n[which(clonedf_n$from %in% c("baseline_N","6months_N")),]
  clonedf_n <- clonedf_n[which(clonedf_n$to %in% c("6months_N","12months_N")),]

  self_comb_n <- data.frame(from = c("baseline_N","6months_N","12months_N"), to = c("baseline_N","6months_N","12months_N"))
  self_clonedf_n <- self_comb_n

  clonedf_n <- rbind(clonedf_n, self_clonedf_n)

  baseline_N <- unique(df_pop_time$baseline_N$clone_id)
  six_months_N <- unique(df_pop_time$`6months_N`$clone_id)
  twelve_months_N <- unique(df_pop_time$`12months_N`$clone_id)

  exclusive_intersects_n = c(
    length(setdiff(intersect(baseline_N, twelve_months_N), six_months_N)),
    length(setdiff(intersect(baseline_N, six_months_N), twelve_months_N)),
    length(setdiff(intersect(six_months_N, twelve_months_N), baseline_N)),
    length(setdiff(baseline_N, union(six_months_N, twelve_months_N))),
    length(setdiff(six_months_N, union(baseline_N, twelve_months_N))),
    length(setdiff(twelve_months_N, union(six_months_N, baseline_N)))
  )

  clonedf_n$value <- exclusive_intersects_n

  grid.col = c("#a6cee3", "#b2df8a", "#fb9a99")
  names(grid.col) <- c("12months_N","6months_N","baseline_N")



  ## Adding overlap of 3 time points
  clonedf_n = clonedf_n %>% mutate(across(c(from,to), as.character)) %>%
    mutate(across(c(value), as.numeric))

  overlap_3_n <- length(intersect(intersect(baseline_N, twelve_months_N), six_months_N))
  clonedf_n[nrow(clonedf_n)+1,] <- c("baseline_N","6months_N_12months_N",overlap_3_n)

  # Calculate overlaps for Plasmablast population
  # -------------------------------------------

  clonedf_p <- clonedf
  clonedf_p <- clonedf_p[which(clonedf_p$from %in% c("baseline_P","6months_P")),]
  clonedf_p <- clonedf_p[which(clonedf_p$to %in% c("6months_P","12months_P")),]

  self_comb_p <- data.frame(from = c("baseline_P","6months_P","12months_P"), to = c("baseline_P","6months_P","12months_P"))
  self_clonedf_p <- self_comb_p

  clonedf_p <- rbind(clonedf_p, self_clonedf_p)

  baseline_P <- unique(df_pop_time$baseline_P$clone_id)
  six_months_P <- unique(df_pop_time$`6months_P`$clone_id)
  twelve_months_P <- unique(df_pop_time$`12months_P`$clone_id)

  exclusive_intersects_p = c(
    length(setdiff(intersect(baseline_P, twelve_months_P), six_months_P)),
    length(setdiff(intersect(baseline_P, six_months_P), twelve_months_P)),
    length(setdiff(intersect(six_months_P, twelve_months_P), baseline_P)),
    length(setdiff(baseline_P, union(six_months_P, twelve_months_P))),
    length(setdiff(six_months_P, union(baseline_P, twelve_months_P))),
    length(setdiff(twelve_months_P, union(six_months_P, baseline_P)))
  )

  clonedf_p$value <- exclusive_intersects_p

  grid.col = c("#a6cee3", "#b2df8a", "#fb9a99")
  names(grid.col) <- c("12months_P","6months_P","baseline_P")



  ## Adding overlap of 3 time points
  clonedf_p = clonedf_p %>% mutate(across(c(from,to), as.character)) %>%
    mutate(across(c(value), as.numeric))

  overlap_3_p <- length(intersect(intersect(baseline_P, twelve_months_P), six_months_P))
  clonedf_p[nrow(clonedf_p)+1,] <- c("baseline_P","6months_P_12months_P",overlap_3_p)

  # Calculate overlaps for double negative population
  # -------------------------------------------

  clonedf_dn <- clonedf
  clonedf_dn <- clonedf_dn[which(clonedf_dn$from %in% c("baseline_DN","6months_DN")),]
  clonedf_dn <- clonedf_dn[which(clonedf_dn$to %in% c("6months_DN","12months_DN")),]

  self_comb_dn <- data.frame(from = c("baseline_DN","6months_DN","12months_DN"), to = c("baseline_DN","6months_DN","12months_DN"))
  self_clonedf_dn <- self_comb_dn

  clonedf_dn <- rbind(clonedf_dn, self_clonedf_dn)

  baseline_DN <- unique(df_pop_time$baseline_DN$clone_id)
  six_months_DN <- unique(df_pop_time$`6months_DN`$clone_id)
  twelve_months_DN <- unique(df_pop_time$`12months_DN`$clone_id)

  exclusive_intersects_dn = c(
    length(setdiff(intersect(baseline_DN, twelve_months_DN), six_months_DN)),
    length(setdiff(intersect(baseline_DN, six_months_DN), twelve_months_DN)),
    length(setdiff(intersect(six_months_DN, twelve_months_DN), baseline_DN)),
    length(setdiff(baseline_DN, union(six_months_DN, twelve_months_DN))),
    length(setdiff(six_months_DN, union(baseline_DN, twelve_months_DN))),
    length(setdiff(twelve_months_DN, union(six_months_DN, baseline_DN)))
  )

  clonedf_dn$value <- exclusive_intersects_dn

  grid.col = c("#a6cee3", "#b2df8a", "#fb9a99")
  names(grid.col) <- c("12months_DN","6months_DN","baseline_DN")



  ## Adding overlap of 3 time points
  clonedf_dn = clonedf_dn %>% mutate(across(c(from,to), as.character)) %>%
    mutate(across(c(value), as.numeric))

  overlap_3_dn <- length(intersect(intersect(baseline_DN, twelve_months_DN), six_months_DN))
  clonedf_dn[nrow(clonedf_dn)+1,] <- c("baseline_DN","6months_DN_12months_DN",overlap_3_dn)

  # ---------------------------------------------------------
  # Saving both tables
  clonedf$patient <- rep(df_pop_time[[1]]$patient[1], nrow(clonedf))
  clonedf_m$patient  <- rep(df_pop_time[[1]]$patient[1], nrow(clonedf_m))
  clonedf_n$patient  <- rep(df_pop_time[[1]]$patient[1], nrow(clonedf_n))
  clonedf_p$patient  <- rep(df_pop_time[[1]]$patient[1], nrow(clonedf_p))
  clonedf_dn$patient  <- rep(df_pop_time[[1]]$patient[1], nrow(clonedf_dn))



  write.table(clonedf, file = paste(patdir_overlap,"/Clone_overlap_comparison_time_points_population", df_pop_time[[1]]$treatment[1], "_", df_pop_time[[1]]$patient[1], ".tsv", sep=""), sep = "\t", quote = F, row.names = F)


  #all_clonedf <- rbind(all_clonedf, clonedf)
  all_clonedf_m <- rbind(all_clonedf_m, clonedf_m)
  all_clonedf_n <- rbind(all_clonedf_n, clonedf_n)
  all_clonedf_p <- rbind(all_clonedf_p, clonedf_p)
  all_clonedf_dn <- rbind(all_clonedf_dn, clonedf_dn)


}

#write.table(all_clonedf, file = paste(outdir,"/Clone_overlap", "/Clone_overlap_comparison_time_points_population.tsv", sep=""), sep = "\t", quote = F, row.names = F)
write.table(all_clonedf_m, file = paste(outdir,"/Clone_overlap", "/Clone_overlap_comparison_time_points_memory_population.tsv", sep=""), sep = "\t", quote = F, row.names = F)
write.table(all_clonedf_n, file = paste(outdir,"/Clone_overlap", "/Clone_overlap_comparison_time_points_naive_population.tsv", sep=""), sep = "\t", quote = F, row.names = F)
write.table(all_clonedf_p, file = paste(outdir,"/Clone_overlap", "/Clone_overlap_comparison_time_points_plasma_population.tsv", sep=""), sep = "\t", quote = F, row.names = F)
write.table(all_clonedf_dn, file = paste(outdir,"/Clone_overlap", "/Clone_overlap_comparison_time_points_doubleneg_population.tsv", sep=""), sep = "\t", quote = F, row.names = F)

# Boxplot overlaps comparison


all_clonedf_tp$fromto <- paste(all_clonedf_tp$from,"_",all_clonedf_tp$to, sep = "")
all_clonedf_tp <- all_clonedf_tp[!(all_clonedf_tp$fromto %in% c("baseline_baseline","6months_6months","12months_12months")),]
all_clonedf_tp$fromto <- factor(all_clonedf_tp$fromto, levels=c("baseline_6months","6months_12months","baseline_12months", "baseline_6months_12months"))
all_clonedf_tp$value <- as.numeric(all_clonedf_tp$value)
# Add new row with overlaps of 3 TPs

my_comparisons <- list( c("baseline_6months", "6months_12months"),
                        c("baseline_6months_12months","baseline_6months")
                        )
boxplot_overlaps <- ggpubr::ggboxplot(all_clonedf_tp, x="fromto", y="value",
                                           add = "jitter", legend="none",
                                          line.color = "gray", line.size = 0.4,
                                           short.panel.labs = T) +
  scale_x_discrete(breaks=c("baseline_6months","6months_12months","baseline_12months", "baseline_6months_12months"),
                   labels=c("B-6M", "6M-12M", "B-12M", "B-6M-12M")) +
  xlab("") + ylab("Number overlaps") +
  stat_compare_means(comparisons = my_comparisons, paired=T, method = "wilcox", label="p.signif", hide.ns = T) +
 stat_compare_means(paired=T, label="p.format", label.y = 35000)
boxplot_overlaps

ggsave(plot = boxplot_overlaps, paste(outdir,"/Clone_overlap/Clone_overlap_boxplot_allvalues.png", sep = ""), device="png", width = 10, height = 7, units = "cm")
ggsave(plot = boxplot_overlaps, paste(outdir,"/Clone_overlap/Clone_overlap_boxplot_allvalues.pdf", sep = ""), device="pdf", width = 10, height = 7, units = "cm")



## Boxplots overlap comparison for memory cells

all_clonedf_m$fromto <- paste(all_clonedf_m$from,"_",all_clonedf_m$to, sep = "")
all_clonedf_m <- all_clonedf_m[!(all_clonedf_m$fromto %in% c("baseline_M_baseline_M","6months_M_6months_M","12months_M_12months_M")),]
all_clonedf_m$fromto <- factor(all_clonedf_m$fromto, levels=c("baseline_M_6months_M","6months_M_12months_M","baseline_M_12months_M", "baseline_M_6months_M_12months_M"))

all_clonedf_m$value <- as.numeric(all_clonedf_m$value)

my_comparisons <- list( c("baseline_M_6months_M", "6months_M_12months_M"),
                        c("baseline_M_6months_M_12months_M", "baseline_M_6months_M"),
                        c("baseline_M_6months_M_12months_M", "6months_M_12months_M")
                  )
boxplot_overlaps_m <- ggpubr::ggboxplot(all_clonedf_m, x="fromto", y="value",
                                      add = "jitter", legend="none",
                                      line.color = "gray", line.size = 0.4,
                                      short.panel.labs = T) +
  scale_x_discrete(breaks=c("baseline_M_6months_M","6months_M_12months_M","baseline_M_12months_M", "baseline_M_6months_M_12months_M"),
                   labels=c("B-6M", "6M-12M", "B-12M", "B-6M-12M")) +
  xlab("") + ylab("Number overlaps") +
  stat_compare_means(comparisons = my_comparisons, paired=T, method = "wilcox", label="p.signif", hide.ns = T) +
  stat_compare_means(paired=T, label="p.format", label.y = 3500)
boxplot_overlaps_m

ggsave(plot = boxplot_overlaps_m, paste(outdir,"/Clone_overlap/Clone_overlap_boxplot_memory_allvalues.png", sep = ""), device="png", width = 10, height = 7, units = "cm")
ggsave(plot = boxplot_overlaps_m, paste(outdir,"/Clone_overlap/Clone_overlap_boxplot_memory_allvalues.pdf", sep = ""), device="pdf", width = 10, height = 7, units = "cm")


## Boxplots overlap comparison for naive cells

all_clonedf_n$fromto <- paste(all_clonedf_n$from,"_",all_clonedf_n$to, sep = "")
all_clonedf_n <- all_clonedf_n[!(all_clonedf_n$fromto %in% c("baseline_N_baseline_N","6months_N_6months_N","12months_N_12months_N")),]
all_clonedf_n$fromto <- factor(all_clonedf_n$fromto, levels=c("baseline_N_6months_N","6months_N_12months_N","baseline_N_12months_N", "baseline_N_6months_N_12months_N"))

all_clonedf_n$value <- as.numeric(all_clonedf_n$value)

my_comparisons <- list( c("baseline_N_6months_N", "6months_N_12months_N"),
                        c("baseline_N_6months_N", "baseline_N_12months_N"),
                        c("baseline_N_6months_N_12months_N", "baseline_N_6months_N"),
                        c("baseline_N_6months_N_12months_N", "6months_N_12months_N"),
                        c("baseline_N_6months_N_12months_N", "baseline_N_12months_N")
)
boxplot_overlaps_n <- ggpubr::ggboxplot(all_clonedf_n, x="fromto", y="value",
                                        add = "jitter", legend="none",
                                        line.color = "gray", line.size = 0.4,
                                        short.panel.labs = T) +
  scale_x_discrete(breaks=c("baseline_N_6months_N","6months_N_12months_N","baseline_N_12months_N", "baseline_N_6months_N_12months_N"),
                   labels=c("B-6M", "6M-12M", "B-12M", "B-6M-12M")) +
  xlab("") + ylab("Number overlaps") +
  stat_compare_means(comparisons = my_comparisons, paired=T, method = "wilcox", label="p.signif", hide.ns = T) +
  stat_compare_means(paired=T, label="p.format", label.y = 3500)
boxplot_overlaps_n

ggsave(plot = boxplot_overlaps_n, paste(outdir,"/Clone_overlap/Clone_overlap_boxplot_naive_allvalues.png", sep = ""), device="png", width = 10, height = 7, units = "cm")
ggsave(plot = boxplot_overlaps_n, paste(outdir,"/Clone_overlap/Clone_overlap_boxplot_naive_allvalues.pdf", sep = ""), device="pdf", width = 10, height = 7, units = "cm")

## Boxplots overlap comparison for plasmablasts cells

all_clonedf_p$fromto <- paste(all_clonedf_p$from,"_",all_clonedf_p$to, sep = "")
all_clonedf_p <- all_clonedf_p[!(all_clonedf_p$fromto %in% c("baseline_P_baseline_P","6months_P_6months_P","12months_P_12months_P")),]
all_clonedf_p$fromto <- factor(all_clonedf_p$fromto, levels=c("baseline_P_6months_P","6months_P_12months_P","baseline_P_12months_P", "baseline_P_6months_P_12months_P"))

all_clonedf_p$value <- as.numeric(all_clonedf_p$value)

my_comparisons <- list( c("baseline_P_6months_P", "6months_P_12months_P"),
                        c("baseline_P_6months_P", "baseline_P_12months_P"),
                        c("baseline_P_6months_P_12months_P", "baseline_P_6months_P"),
                        c("baseline_P_6months_P_12months_P", "6months_P_12months_P"),
                        c("baseline_P_6months_P_12months_P", "baseline_P_12months_P")
)
boxplot_overlaps_p <- ggpubr::ggboxplot(all_clonedf_p, x="fromto", y="value",
                                        add = "jitter", legend="none",
                                        line.color = "gray", line.size = 0.4,
                                        short.panel.labs = T) +
  scale_x_discrete(breaks=c("baseline_P_6months_P","6months_P_12months_P","baseline_P_12months_P", "baseline_P_6months_P_12months_P"),
                   labels=c("B-6M", "6M-12M", "B-12M", "B-6M-12M")) +
  xlab("") + ylab("Number overlaps") +
  stat_compare_means(comparisons = my_comparisons, paired=T, method = "wilcox", label="p.signif", hide.ns = T) +
  stat_compare_means(paired=T, label="p.format", label.y = 3500)
boxplot_overlaps_p

ggsave(plot = boxplot_overlaps_p, paste(outdir,"/Clone_overlap/Clone_overlap_boxplot_plasmablasts_allvalues.png", sep = ""), device="png", width = 10, height = 7, units = "cm")
ggsave(plot = boxplot_overlaps_p, paste(outdir,"/Clone_overlap/Clone_overlap_boxplot_plasmablasts_allvalues.pdf", sep = ""), device="pdf", width = 10, height = 7, units = "cm")

## Boxplots overlap comparison for double_negative cells

all_clonedf_dn$fromto <- paste(all_clonedf_dn$from,"_",all_clonedf_dn$to, sep = "")
all_clonedf_dn <- all_clonedf_dn[!(all_clonedf_dn$fromto %in% c("baseline_DN_baseline_DN","6months_DN_6months_DN","12months_DN_12months_DN")),]
all_clonedf_dn$fromto <- factor(all_clonedf_dn$fromto, levels=c("baseline_DN_6months_DN","6months_DN_12months_DN","baseline_DN_12months_DN", "baseline_DN_6months_DN_12months_DN"))

all_clonedf_dn$value <- as.numeric(all_clonedf_dn$value)

my_comparisons <- list( c("baseline_DN_6months_DN", "6months_DN_12months_DN"),
                        c("baseline_DN_6months_DN", "baseline_DN_12months_DN"),
                        c("baseline_DN_6months_DN_12months_DN", "baseline_DN_6months_DN"),
                        c("baseline_DN_6months_DN_12months_DN", "6months_DN_12months_DN"),
                        c("baseline_DN_6months_DN_12months_DN", "baseline_DN_12months_DN")
)
boxplot_overlaps_dn <- ggpubr::ggboxplot(all_clonedf_dn, x="fromto", y="value",
                                        add = "jitter", legend="none",
                                        line.color = "gray", line.size = 0.4,
                                        short.panel.labs = T) +
  scale_x_discrete(breaks=c("baseline_DN_6months_DN","6months_DN_12months_DN","baseline_DN_12months_DN", "baseline_DN_6months_DN_12months_DN"),
                   labels=c("B-6M", "6M-12M", "B-12M", "B-6M-12M")) +
  xlab("") + ylab("Number overlaps") +
  stat_compare_means(comparisons = my_comparisons, paired=T, method = "wilcox", label="p.signif", hide.ns = T) +
  stat_compare_means(paired=T, label="p.format", label.y = 2000)
boxplot_overlaps_dn

ggsave(plot = boxplot_overlaps_dn, paste(outdir,"/Clone_overlap/Clone_overlap_boxplot_double_negative_allvalues.png", sep = ""), device="png", width = 10, height = 7, units = "cm")
ggsave(plot = boxplot_overlaps_dn, paste(outdir,"/Clone_overlap/Clone_overlap_boxplot_double_negative_allvalues.pdf", sep = ""), device="pdf", width = 10, height = 7, units = "cm")
