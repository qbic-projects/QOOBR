# Peptides table
library(dplyr)
library(stringr)
library(ggplot2)
library(data.table)
theme_set(theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

options(error = function() traceback())

all_patients_tab <- data.frame()
all_patients_pep_tab_subset <- data.frame()
# Run inside proteomics folder

for (pat in c("CLAD1","CLAD2","CLAD3","CLAD5","CLAD6","CLAD7","CLAD8")) {
  pep_tab <- read.table(paste0("results/",pat,"/combined/txt/peptides.txt"), sep="\t", header = T)
  seq_tab <- read.table(paste0("translation/",pat,"_igblast_results.tab"), sep = "\t", header = T)
  seq_tab$sequence_id_vdj <- paste0(seq_tab$sequence_id, seq_tab$v_call, seq_tab$d_call, seq_tab$j_call)

  # Remove a duplicate sample
  if (pat == "CLAD2"){
    pep_tab$Identification.type.CLAD_P2_12M_1 <- NULL
    pep_tab$Intensity.CLAD_P2_12M_1 <- NULL
    pep_tab$Experiment.CLAD_P2_12M_1 <- NULL
  }

  # Annotate where in the sequence (cdr3, cdr1, cdr2, frw1...) does the peptide hit.
  peptide_regions <- c()
  pep_seq_col <- c()
  prot_id_col <- c()
  prot_id_novdj_col <- c()
  pep_region_1 <- c()
  pep_region_2 <- c()
  pep_region_1_aa <- c()
  pep_region_2_aa <- c()
  pep_region_naa <- c()
  cdr3_col <- c()
  cdr3_aa_col <- c()
  cdr2_col <- c()
  cdr2_aa_col <- c()
  cdr1_col <- c()
  cdr1_aa_col <- c()
  number_matches <- c()

  for (p in c(1:nrow(pep_tab))) {

      pep_seq = pep_tab$Sequence[p]
      pep_prots <- pep_tab$Proteins[p]

      pep_prots <- str_split(pep_prots,";",simplify=T)
      number_matches <- append(number_matches, length(pep_prots))

      # if peptides matches more than 10 sequences then skip peptide
      if (length(pep_prots) > 100){
        next
      }

      # For each of the protein matches per peptide
      for (m in c(1:length(pep_prots))) {
          pep_prot = pep_prots[m]

          # Extract sequence ID without VDJ ID.
          prot_id_novdj = str_extract(pep_prot, pattern = "[0-9]{1,3}[A-Z]{8}")
          prot_id = pep_prot

          # Search for sequence ID in the AIRR table
          seq_prot = seq_tab[which(seq_tab$sequence_id_vdj == prot_id),] %>%
                      distinct()

          # Identify the aligned aa to reference
          seq_prot_aa = seq_prot$sequence_alignment_aa

          # Locate where the peptide is located in the original sequence
          pep_location = str_locate(seq_prot_aa, pep_seq)
          pep_start_nt = (pep_location[1]-1)*3+1

          # if peptide not found, go to next iteration
          if (is.na(pep_start_nt)) {
            next
          }

          pep_end_nt = (pep_location[2])*3

          # Identify in which Ig region the peptide matches

          regions = c("fwr1_start","fwr1_end","cdr1_start","cdr1_end","fwr2_start","fwr2_end","cdr2_start","cdr2_end","fwr3_start","fwr3_end","cdr3_start","cdr3_end")
          regions = regions[!is.na(seq_prot[regions])]


          if (seq_prot$v_sequence_start > 1) { alignment_start = seq_prot$v_sequence_start } else { alignment_start = 0 }

          # If peptide position is higher than available regions positions, skip match
          if ( pep_start_nt >= seq_prot[regions[length(regions)]]  ) {
            next
          }

          # Search for region where the peptide starts
          r=1
          while ( pep_start_nt >= (seq_prot[regions[r]] -1 - alignment_start) ) {
            r = r + 1
            if (r > length(regions)){
              break
            }
          }

          if (r > length(regions)){
            next
          }

          # search for region where the peptide ends
          pep_region_1 <- append(pep_region_1, regions[r])
          if (pep_end_nt >= seq_prot[regions[r]] - alignment_start) { # peptide goes over to next region
            naa_region1 = as.numeric(((seq_prot[regions[r]] +1 - alignment_start)-pep_start_nt)/3)
            pep_region_1_aa <- append(pep_region_1_aa, naa_region1)
          } else { # peptide stays within this region
            naa_region1 = nchar(pep_seq)
            pep_region_1_aa <- append(pep_region_1_aa, naa_region1)
          }

          if ( r < length(regions) ) {
            if ( !(pep_end_nt <= (seq_prot[regions[r+1]] - alignment_start) ) ) {
              pep_region_2 <- append(pep_region_2, regions[r+1])
              naa_region2 = as.numeric((pep_end_nt - (seq_prot[regions[r+1]] - alignment_start -1))/3)
              pep_region_2_aa <- append(pep_region_2_aa, naa_region2)
            } else {
              pep_region_2 <- append(pep_region_2, "")
              pep_region_2_aa <- append(pep_region_2_aa, "")
              naa_region2 = 0
            }
          } else {
            pep_region_2 <- append(pep_region_2, "")
            pep_region_2_aa <- append(pep_region_2_aa, "")
            naa_region2 = 0
          }
          pep_seq_col <- append(pep_seq_col, pep_seq)
          prot_id_col <- append(prot_id_col, prot_id)
          prot_id_novdj_col <- append(prot_id_novdj_col, prot_id_novdj)


      }

  }

  pep_tab$match_number <- number_matches

  pep_tab$baseline <- select(pep_tab, contains("BL"))[,3]
  pep_tab$`6months` <- select(pep_tab, contains("6M"))[,3]
  pep_tab$`12months` <- select(pep_tab, contains("12M"))[,3]
  pep_tab_subset <- pep_tab[,c("Sequence","Mass", "Proteins", "Start.position", "End.position", "PEP", "Score", "match_number","baseline","6months","12months")]

  pep_tab_subset <- melt(data.table(pep_tab_subset), id.vars= c("Sequence","Mass", "Proteins", "Start.position", "End.position", "PEP", "Score", "match_number"), measure.vars=c("baseline","6months","12months"), value.name="Intensity", variable.name = "time_point")


  ggplot(pep_tab, aes(match_number)) +
    geom_histogram(binwidth = 1) +
    scale_x_continuous(limits = c(0,100), breaks = c(1,10,20,30,40,50,60,70,80,90,100)) +
    xlab("Number of protein matches") + ylab("Peptide number")
  ggsave(file = paste0("results/",pat,"/peptide_match_distribution.png"), width = 10, height = 6, units = "cm")

  write.table(pep_tab_subset, file = paste0("results/",pat,"/original_peptide_protein_matches.tsv"), sep = "\t", quote = F, row.names = F)
  pep_tab_subset$subject <- rep(pat,nrow(pep_tab_subset))
  all_patients_pep_tab_subset <- rbind(all_patients_pep_tab_subset,pep_tab_subset)

  df <- data.frame(peptide_id = pep_seq_col,
                    protein_id = prot_id_novdj_col,
                    protein_id_vdj = prot_id_col,
                    peptide_region_1 = pep_region_1,
                    peptide_region_2 = pep_region_2,
                    peptide_region_1_aa = pep_region_1_aa,
                    peptide_region_2_aa = pep_region_2_aa)

  write.table(df, file = paste0("results/",pat,"/peptide_protein_region_matches.tsv"), sep = "\t", quote = F, row.names = F)

  merged_df_peptab <- merge(x = df, y=pep_tab_subset, by.x = "peptide_id", by.y="Sequence", all.x=T)


  # Annotation of overlapping to cdrs and how many aas are overlapping
  merged_df_peptab$cdr3 <- c(grepl("cdr3", merged_df_peptab$peptide_region_1 ) | grepl("cdr3", merged_df_peptab$peptide_region_2))
  merged_df_peptab$cdr2 <- c(grepl("cdr2", merged_df_peptab$peptide_region_1)  | grepl("cdr2", merged_df_peptab$peptide_region_2))
  merged_df_peptab$cdr1 <- c(grepl("cdr1", merged_df_peptab$peptide_region_1)  | grepl("cdr1", merged_df_peptab$peptide_region_2))

  merged_df_peptab$cdr3_aa <- c("")
  merged_df_peptab$cdr2_aa <- c("")
  merged_df_peptab$cdr1_aa <- c("")

  merged_df_peptab[which(grepl("cdr3", merged_df_peptab$peptide_region_1 )),]$cdr3_aa <- merged_df_peptab[which(grepl("cdr3", merged_df_peptab$peptide_region_1 )),c("peptide_region_1_aa")]
  merged_df_peptab[which(grepl("cdr3", merged_df_peptab$peptide_region_2 )),]$cdr3_aa <- merged_df_peptab[which(grepl("cdr3", merged_df_peptab$peptide_region_2 )),c("peptide_region_2_aa")]
  merged_df_peptab[which(grepl("cdr2", merged_df_peptab$peptide_region_1 )),]$cdr2_aa <- merged_df_peptab[which(grepl("cdr2", merged_df_peptab$peptide_region_1 )),c("peptide_region_1_aa")]
  merged_df_peptab[which(grepl("cdr2", merged_df_peptab$peptide_region_2 )),]$cdr2_aa <- merged_df_peptab[which(grepl("cdr2", merged_df_peptab$peptide_region_2 )),c("peptide_region_2_aa")]
  merged_df_peptab[which(grepl("cdr1", merged_df_peptab$peptide_region_1 )),]$cdr1_aa <- merged_df_peptab[which(grepl("cdr1", merged_df_peptab$peptide_region_1 )),c("peptide_region_1_aa")]
  merged_df_peptab[which(grepl("cdr1", merged_df_peptab$peptide_region_2 )),]$cdr1_aa <- merged_df_peptab[which(grepl("cdr1", merged_df_peptab$peptide_region_2 )),c("peptide_region_2_aa")]

  merged_df_peptab$region <- c("other")
  merged_df_peptab[which(merged_df_peptab$cdr3 ),]$region <- "cdr3"
  merged_df_peptab[which(merged_df_peptab$cdr2 ),]$region <- "cdr2"
  merged_df_peptab[which(merged_df_peptab$cdr1 ),]$region <- "cdr1"

  merged_df_peptab$cdr1_aa <- as.numeric(merged_df_peptab$cdr1_aa)
  merged_df_peptab$cdr2_aa <- as.numeric(merged_df_peptab$cdr2_aa)
  merged_df_peptab$cdr3_aa <- as.numeric(merged_df_peptab$cdr3_aa)

  merged_df_peptab$cdr1_aa[is.na(merged_df_peptab$cdr1_aa)] <- 0
  merged_df_peptab$cdr2_aa[is.na(merged_df_peptab$cdr2_aa)] <- 0
  merged_df_peptab$cdr3_aa[is.na(merged_df_peptab$cdr3_aa)] <- 0

  df_naa_plot <- melt(data.table(merged_df_peptab), id.vars= c("peptide_id", "protein_id_vdj"), measure.vars=c("cdr1_aa","cdr2_aa","cdr3_aa"), value.name="number_aa")
  df_naa_plot$number_aa <- as.numeric(df_naa_plot$number_aa)

  pr <- ggplot(merged_df_peptab, aes(x = region)) +
          geom_bar() +
          xlab("")
  ggsave(filename = paste0("results/",pat,"/peptide_regions_distribution.pdf"), plot=pr, width = 6, height = 6, units = "cm")
  ggsave(filename = paste0("results/",pat,"/peptide_regions_distribution.png"), plot=pr, width = 6, height = 6, units = "cm")


  pr_aa <- ggplot(df_naa_plot, aes(x=variable, y=number_aa, group=variable)) +
            geom_boxplot() +
            xlab("")
  ggsave(filename = paste0("results/",pat,"/peptide_regions_number_aa_distribution.pdf"), plot=pr_aa, width = 6, height = 5, units = "cm")
  ggsave(filename = paste0("results/",pat,"/peptide_regions_number_aa_distribution.png"), plot=pr_aa, width = 6, height = 5, units = "cm")

  merged_df_peptab$cdr_aa <- merged_df_peptab$cdr1_aa + merged_df_peptab$cdr2_aa + merged_df_peptab$cdr3_aa
  merged_df_peptab <- merged_df_peptab[,c("peptide_id","protein_id","protein_id_vdj","peptide_region_1","peptide_region_2","peptide_region_1_aa", "peptide_region_2_aa", "cdr3", "cdr2", "cdr1", "cdr3_aa", "cdr2_aa", "cdr1_aa", "cdr_aa", "region", "Mass", "Proteins", "Start.position", "End.position", "PEP", "Score", "match_number","time_point","Intensity")]
  merged_df_peptab <-  merged_df_peptab %>% filter(Intensity > 0)
  write.table(merged_df_peptab, file = paste0("results/",pat,"/peptide_protein_match_region_merged_tables.tsv"), sep = "\t", quote = F, row.names = F)

  merged_df_peptab$subject_id <- pat
  all_patients_tab <- rbind(all_patients_tab, merged_df_peptab)


}

write.table(all_patients_tab, file = paste0("results/","all_patients_peptide_protein_match_region_merged_tables.tsv"), sep = "\t", quote = F, row.names = F)
write.table(all_patients_pep_tab_subset, file = paste0("results/","all_patients_peptide_protein_match_unfiltered.tsv"), sep = "\t", quote = F, row.names = F)

# Plot Number of protein matches per peptide for all patients together
ggplot(all_patients_tab, aes(match_number)) +
  geom_histogram(binwidth = 1) +
  scale_x_continuous(limits = c(0,100), breaks = c(1,10,20,30,40,50,60,70,80,90,100)) +
  xlab("Number of protein matches") + ylab("Peptide number")
ggsave(file = paste0("results/peptide_match_all_patients_distribution_after_filtering.png"), width = 10, height = 6, units = "cm")
ggsave(file = paste0("results/peptide_match_all_patients_distribution_after_filtering.pdf"), width = 10, height = 6, units = "cm")


ggplot(all_patients_pep_tab_subset, aes(match_number)) +
  geom_histogram(binwidth = 1) +
  scale_x_continuous(limits = c(0,100), breaks = c(1,10,20,30,40,50,60,70,80,90,100)) +
  xlab("Number of protein matches") + ylab("Peptide number")
ggsave(file = paste0("results/peptide_match_all_patients_distribution_before_filtering.png"), width = 10, height = 6, units = "cm")
ggsave(file = paste0("results/peptide_match_all_patients_distribution_before_filtering.pdf"), width = 10, height = 6, units = "cm")

