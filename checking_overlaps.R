df_all <- read.csv("../../pipeline_results_per_patient/repertoire_analysis/repertoire_comparison/all_data.tsv", sep="\t")

df_CLAD2 <- df_all[which(df_all$source=="CLAD2"),]
df_CLAD2_baseline <- df_CLAD2[(df_CLAD2$extract_time=="baseline"),]
df_CLAD2_12months <- df_CLAD2[(df_CLAD2$extract_time=="12months"),]
df_CLAD2_6months <- df_CLAD2[(df_CLAD2$extract_time=="6months"),]

length(unique(df_CLAD2_baseline$clone_id))
length(unique(df_CLAD2_12months$clone_id))
length(unique(df_CLAD2_6months$clone_id))
length(intersect(df_CLAD2_baseline,df_CLAD2_6months))

df_ab <- read.csv("./repertoire_comparison/Abundance/Clonal_abundance_data_subject.tsv", sep = "\t")
df_ab_CLAD2 <- df_ab[(df_ab$patient=="CLAD2"),]
df_ab_CLAD2_baseline <- df_ab_CLAD2[(df_ab_CLAD2$time_point == "baseline"),]
df_ab_CLAD2_6months <- df_ab_CLAD2[(df_ab_CLAD2$time_point == "6months"),]
df_ab_CLAD2_12months <- df_ab_CLAD2[(df_ab_CLAD2$time_point == "12months"),]

length(unique(df_ab_CLAD2_baseline$clone_id))
length(unique(df_ab_CLAD2_6months$clone_id))
length(unique(df_ab_CLAD2_12months$clone_id))
length(intersect(df_ab_CLAD2_baseline,df_ab_CLAD2_6months))
