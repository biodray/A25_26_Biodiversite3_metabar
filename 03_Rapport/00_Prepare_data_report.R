
library(tidyverse)
library(here)
library(Hmisc)

source(file.path(here::here(), "01_Code", "Functions", "get.value.R"))

# Déclare ici tes objets nécessaires: 'LOCUS', etc.
# Suppose que 'data.info' existe déjà (sinon charge-le)

data.info <- read_csv(here("00_Data", "00_FileInfos", "SeqInfo.csv"))
LOCUS <- data.info$Loci %>% unique()
#projets <- data.info$ID_subproject %>% unique() %>% str_subset(pattern = "ALL", negate = T)


projets <- stringr::str_split(get.value("group.metabaR", file = file.path(here::here(), "Options.txt")), pattern = ";")[[1]]
# Stats

# Étape 1: charger toutes les tables indépendamment du sous-projet
cutadapt.res <- read_csv(here("00_Data", "02a_Cutadapt", "log", "Cutadapt_Stats.csv"))
esv.reads.res <- read_csv(here("02_Results", "02_Filtrations", "ESVtab_Stats_Nreads.csv"))
names(esv.reads.res)[4] <- "ESVdada2"

# Fusionner pour tous les locus
dada2.summary <- map_dfr(LOCUS, function(l) {
  read_csv(here("00_Data", "02b_Filtered_dada2", "log", paste0(l, "_dada2_summary.csv"))) %>%
    mutate(Loci = l)
})
dada2.summary <- dada2.summary %>%
  mutate(Trim = reads.out,
         Adapt = reads.in
  ) %>%
  select(-c(reads.in, reads.out))

tag.summary <- map_dfr(LOCUS, function(l) {
  read_csv(here("00_Data", "04_ESVcorrected", paste0("ESVtab_Stats_postTagjump_", l, ".csv"))) %>% mutate(Loci = l)
})


esv.motus.res <- read_csv(here("02_Results", "02_Filtrations", "ESVtab_Stats_NESV.csv"))
names(esv.motus.res)[5] <- "ESVdada2"



for (proj in projets) {
  
  metabar.summary <- data.frame()
  
  for(l in LOCUS){
    
    
    path <- here::here("00_Data", "04_ESVcorrected", paste0("ESVtab_Stats_postMetabaR_", l, "_", proj, ".csv"))
    
    if(file.exists(path)){
      
      metabar.int <- readr::read_csv(path) %>% mutate(Loci = ifelse(Loci == "1.2e+161", "12S160", Loci))
      metabar.summary <- bind_rows(metabar.summary, metabar.int)    
      
    }
    
  }
  
  
  nreads.res <- cutadapt.res %>%
    left_join(
      dada2.summary %>%
        mutate(ID_sample = ID %>% str_remove(paste(paste0("_", LOCUS), collapse = "|"))) %>%
        select(ID_sample, Loci, Trim)
    ) %>%
    left_join(esv.reads.res) %>%
    left_join(tag.summary %>% select(ID_sample = sample_id, ESVtagjump = nb_reads.tagjump, Loci)) %>%
    left_join(metabar.summary %>% select(ID_sample = sample_id, ESVfinal = nb_reads_postmetabaR, Loci)) %>%
    pivot_longer(c(Raw, Adapt, Trim, Merge,  ESVlength, ESVdada2, ESVtagjump, ESVfinal),
                 names_to = "Step", values_to = "Nreads") %>%
    left_join(
      data.info %>%
        select(ID_sample, Loci, ID_project, Sample_type, Inhibition),
      by = c("ID_sample", "Loci")
    ) %>%
    mutate(Nreads = ifelse(is.na(Nreads), 0, Nreads),
           Nreads = ifelse((Sample_type != "ECH" & Step == "ESVfinal"), NA, Nreads),
           Nreads1 = Nreads + 1,
           Step = factor(Step, levels = c("Raw", "Adapt", "Trim", "Merge", "ESVlength", "ESVdada2", "ESVtagjump", "ESVfinal"))
    )
  
  saveRDS(nreads.res, here("03_Rapport/99_Prepared", paste0("nreads_res_", proj, ".rds")))
  
  
  nmotus.res <- esv.motus.res %>%
    left_join(tag.summary %>% select(ID_sample = sample_id, ESVtagjump = nb_motus.tagjump, Loci)) %>%
    left_join(metabar.summary %>% select(ID_sample = sample_id, ESVfinal = nb_motus_postmetabaR, Loci)) %>%
    pivot_longer(c(Merge,ESVlength, ESVdada2, ESVtagjump, ESVfinal), names_to = "Step", values_to = "Nmotus") %>%
    left_join(data.info %>% select(ID_sample, Loci, ID_project, Sample_type),
              by = c("ID_sample", "Loci")) %>%
    mutate(
      Nmotus = ifelse(is.na(Nmotus), 0, Nmotus),
      Nmotus = ifelse((Sample_type != "ECH" & Step == "ESVfinal"), NA, Nmotus),
      Nmotus1 = Nmotus + 1,
      Step = factor(Step, levels = c("Raw", "Adapt", "Trim", "Merge", "ESVlength", "ESVdada2", "ESVtagjump", "ESVfinal"))
    )
  
  saveRDS(nmotus.res, here("03_Rapport/99_Prepared", paste0("nmotus_res_", proj, ".rds")))
  
}



# Control negative and positive

ESV.table.control.long <- vector("list", length(LOCUS))
names(ESV.table.control.long) <- LOCUS

for (l in LOCUS) {
  MOTUs.control.table <- data.table::fread(here("00_Data", "04_ESVcorrected", paste0("MOTUs.Metabarinfo.postTagjump_", l, "_ALL.csv"))) %>% 
    dplyr::select(ESV, Taxon, phylum, not_a_max_conta, not_an_exclude_taxa, not_detected_in_control)
  ESV.control.table <- data.table::fread(here("00_Data", "04_ESVcorrected", paste0("ESVtab.postTagjump_", l, "_ALL.csv")))
  names(ESV.control.table)[1] <- "ID_sample"
  ESV.control.table <- ESV.control.table %>% dplyr::left_join(data.info %>% select(ID_sample, ID_subproject, Sample_type) %>% 
                                                                distinct(.keep_all = T),
                                                              relationship = "many-to-one") %>% 
    dplyr::filter( Sample_type %nin% c("ECH", "POOL")) 
  ESV.control.table.long <- ESV.control.table %>%
    pivot_longer(-c(ID_sample, ID_subproject, Sample_type), names_to = "ID", values_to = "Nreads") %>%
    dplyr::filter(Nreads > 0) %>% 
    left_join(MOTUs.control.table, by = c("ID" = "ESV")) %>%
    mutate(Loci = l,
           Taxon= ifelse(is.na(Taxon), "Unassigned", Taxon),
           Category = ifelse(not_a_max_conta == FALSE & not_an_exclude_taxa == TRUE, "Contaminant",
                             ifelse(not_a_max_conta == FALSE & not_an_exclude_taxa == FALSE, "Contaminant + exclusion taxa list",
                                    ifelse(not_a_max_conta == TRUE & not_an_exclude_taxa == FALSE, "Exclusion taxa list",       
                                           ifelse(not_a_max_conta == TRUE & not_an_exclude_taxa== TRUE, "Good MOTU", "Undefined??")))))
  
  
  ESV.table.control.long[[l]] <- ESV.control.table.long
  
}  

# Fusionner tous les locus si tu veux un seul fichier par projet

ESV.table.control.long.df <- bind_rows(ESV.table.control.long)

# Selectionner les projets

for (proj in projets) {
  
  saveRDS(  ESV.table.control.long.df %>% dplyr::filter(ID_subproject %in% c(proj,  "ALL", "All", "all")), here("03_Rapport/99_Prepared", paste0("ESVtable_control_long_", proj, ".rds")))
  
}






