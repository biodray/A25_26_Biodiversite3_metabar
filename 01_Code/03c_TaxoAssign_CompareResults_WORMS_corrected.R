# Info --------------------------------------------------------------------

# Compare all taxonomic assignment that were performed
# and create a final taxonomic assignment result
#
# Audrey Bourret
# 2024
#

# Library -----------------------------------------------------------------

library(stringr)
library(readr)
library(magrittr)
library(dplyr)
library(ggplot2)
library(tidyverse)

source(file.path(here::here(), "01_Code", "Functions", "get.value.R"))
`%nin%` = Negate(`%in%`)

# Load assignment results -------------------------------------------------

# This script will check across all folders available
res.dirs <- list.dirs("02_Results/03_TaxoAssign", full.names = F)[-1] 
res.dirs

# Compile the results
res <- tibble()

for(i in res.dirs){

res.file <- list.files(file.path("02_Results/03_TaxoAssign",i), pattern = "RES.all", full.names = T) %>% 
              str_subset("wSamples", negate = T)  

res.int <- read_csv(res.file) %>% mutate(Folder = i)

if("QueryAccVer" %in% names(res.int)){
  res.int <- res.int %>% dplyr::rename(ESV = QueryAccVer )
}

res <- bind_rows(res, res.int)

}

write_csv(res , "02_Results/03_TaxoAssign/Assignements.ALL.csv")

# Graphics ----------------------------------------------------------------

graph.ESV <- res %>% dplyr::filter(Levels %in% c("species")) %>% 
  dplyr::group_by(Taxon, phylum, Levels, Method, Threshold, Folder, RefSeq ) %>% 
  dplyr::summarize(N_ESV = n()) %>%
  mutate(CAT = paste(Method, Threshold)) %>% 
  ggplot(aes(y = Taxon, x = CAT))+
  geom_tile(aes(fill = N_ESV)) +
  scale_fill_viridis_c() +
  facet_grid(phylum~RefSeq + Method, scale = "free", space = "free")+
theme_bw() +
  labs(title = "Visual comparison detected species",
       x = "Assignement method")+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 8),
        axis.text.y = element_text(face = "italic", size = 8, vjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.y = element_text(angle = 0, size = 8), 
        panel.spacing.y=unit(0.1, "lines"),
        legend.position = "right",
        legend.text = element_text(size = 8),
        panel.grid = element_blank()
  )

graph.ESV

graph.loci <- res %>% dplyr::filter(Script %in%  c("Blast.local",   "Blast.local3"),
                                    phylum == "Chordata") %>% 
  dplyr::group_by(Taxon, order, Levels, Method, Threshold, Folder, RefSeq ) %>% 
  summarise(N_loci = length(unique(Loci))) %>% 
  mutate(CAT = paste(Method, Threshold)) %>% 
  ggplot(aes(y = Taxon, x = CAT))+
  geom_tile(aes(fill = N_loci)) +
  scale_fill_viridis_c() +
  facet_grid(order~RefSeq + Method, scale = "free", space = "free")+
  theme_bw(base_size = 8) +
  labs(title = "Visual comparison detected species",
       x = "Assignement method")+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 8),
        axis.text.y = element_text(face = "italic", size = 8, vjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.y = element_text(angle = 0, size = 8), 
        panel.spacing.y=unit(0.1, "lines"),
        legend.position = "right",
        legend.text = element_text(size = 8),
        panel.grid = element_blank()
  )

graph.loci


ggsave(filename = file.path("02_Results/03_TaxoAssign", "Comparison_ALL_nESV.png"), plot =  graph.ESV, height = 12, width = 10)

ggsave(filename = file.path("02_Results/03_TaxoAssign", "Comparison_ALL_nLoci.png"), plot =  graph.loci, height = 12, width = 10)


# Final dataset -----------------------------------------------------------

res %>% group_by(Script, Method, Threshold) %>% summarise(N = n())


res.comp <- res %>% mutate(Assign = paste(RefSeq, Method, Threshold, sep = "_")) %>% 
  dplyr::filter(Threshold %in% c(95,99), Method == "TOP") %>% 
  dplyr::select(ESV, Taxon, Assign) %>% 
  distinct(.keep_all = T) %>% 
  pivot_wider(names_from = Assign, values_from = Taxon) 


# Marine

# The first method that will be used (usually a local DB)
method.1 <- get.value("assign.1")
method.1

# The second method that will be used to "fill" the hole
method.2 <- get.value("assign.2")
method.2

# Subset the dataset with the preferred assignment method
# minor change for the 18S
res.1 <- res %>% dplyr::filter(Loci != "18SV9",
                               Method ==  stringr::str_split(method.1, pattern = ";")[[1]][2] ,
                               Threshold == stringr::str_split(method.1, pattern = ";")[[1]][3],
                               Script  == stringr::str_split(method.1, pattern = ";")[[1]][1]) %>% 
  mutate(Assignment.method = method.1)

res.1 %>% group_by(Loci) %>% summarise(N = n())

# Subset the dataset with the second assignment method
res.2 <- res %>% dplyr::filter(Loci != "18SV9",
                               Method ==  stringr::str_split(method.2, pattern = ";")[[1]][2] ,
                               Threshold == stringr::str_split(method.2, pattern = ";")[[1]][3],
                               Script  == stringr::str_split(method.2, pattern = ";")[[1]][1]) %>% 
  mutate(Assignment.method = method.2)

res.2 %>% group_by(Loci) %>% summarise(N = n())

res.2 %>% dplyr::filter(ESV %nin% res.1$ESV) %>%  group_by(Loci) %>% summarise(N = n())


res.18S <- res %>% dplyr::filter(Loci == "18SV9",
                                 Method ==  "LCA" ,
                                 Threshold == 95,
                                 Script  == stringr::str_split(method.2, pattern = ";")[[1]][1]) %>% 
  mutate(Assignment.method = "Blast;LCA;95")

res.18S %>%  group_by(Loci) %>% summarise(N = n())

res.marine.final <- bind_rows(res.1, 
                       res.2 %>% dplyr::filter(ESV %nin% res.1$ESV),
                       res.18S)
res.marine.final


# Check that no duplicated values persisted
table(duplicated(res.marine.final$ESV))
res.marine.final <- res.marine.final[!duplicated(res.marine.final$ESV),]
# Save the results

write_csv(res.marine.final,  "02_Results/03_TaxoAssign/Assignements.marine.not.corrected.csv")


# NT
res.NT.final <- bind_rows(res.2,
                              res.18S)
res.NT.final


# Check that no duplicated values persisted
table(duplicated(res.NT.final$ESV))
res.NT.final <- res.NT.final[!duplicated(res.NT.final$ESV),]
# Save the results

write_csv(res.NT.final,  "02_Results/03_TaxoAssign/Assignements.NT.not.corrected.csv")



# Estuary

# The first method that will be used (usually a local DB)
method.1 <- get.value("assign.1")
method.1 <- "Blast.local3;TOP;95"

# The second method that will be used to "fill" the hole
method.2 <- get.value("assign.2")
method.2

# Subset the dataset with the preferred assignment method
# minor change for the 18S
res$Script %>% table()

res.1 <- res %>% dplyr::filter(Loci != "18SV9",
  Method ==  stringr::str_split(method.1, pattern = ";")[[1]][2] ,
  Threshold == stringr::str_split(method.1, pattern = ";")[[1]][3],
  Script  == stringr::str_split(method.1, pattern = ";")[[1]][1]) %>% 
  mutate(Assignment.method = method.1)

res.1 %>% group_by(Loci) %>% summarise(N = n())

# Subset the dataset with the second assignment method
res.2 <- res %>% dplyr::filter(Loci != "18SV9",
                               Method ==  stringr::str_split(method.2, pattern = ";")[[1]][2] ,
                               Threshold == stringr::str_split(method.2, pattern = ";")[[1]][3],
                               Script  == stringr::str_split(method.2, pattern = ";")[[1]][1]) %>% 
  mutate(Assignment.method = method.2)

res.2 %>% group_by(Loci) %>% summarise(N = n())

res.2 %>% dplyr::filter(ESV %nin% res.1$ESV) %>%  group_by(Loci) %>% summarise(N = n())


res.estuary.final <- bind_rows(res.1, 
                                  res.2 %>% dplyr::filter(ESV %nin% res.1$ESV),
                                  res.18S)
res.estuary.final


# Check that no duplicated values persisted
table(duplicated(res.estuary.final$ESV))
res.estuary.final <- res.estuary.final[!duplicated(res.estuary.final$ESV),]
# Save the results

write_csv(res.estuary.final,  "02_Results/03_TaxoAssign/Assignements.estuary.not.corrected.csv")




graph.ESV.final <- res.NT.final %>% dplyr::filter(Levels %in% c("species", "genus")) %>% 
  dplyr::group_by(Taxon, phylum, Levels, Loci, Assignment.method ) %>% 
  summarise(N_ESV = n()) %>% 
  dplyr::mutate(Assignment.method = factor(Assignment.method)) %>% 
   ggplot(aes(y = Taxon, x = Assignment.method))+
  geom_tile(aes(fill = N_ESV)) +
  scale_fill_viridis_c() +
  facet_grid(phylum~Loci, scale = "free", space = "free")+
  theme_bw() +
  labs(title = "Detected species combining 2 methods",
       x = "Assignement method")+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 8),
        axis.text.y = element_text(face = "italic", size = 8, vjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.y = element_text(angle = 0, size = 8), 
        panel.spacing.y=unit(0.1, "lines"),
        legend.position = "right",
        legend.text = element_text(size = 8),
        panel.grid = element_blank()
  )

graph.ESV.final

ggsave(filename = file.path("02_Results/03_TaxoAssign", "Comparison_freshwater_NT_nESV.png"), plot =  graph.ESV.final, height = 12, width = 10)



graph.ESV.final <- res.marine.final %>% dplyr::filter(Levels %in% c("species", "genus")) %>% 
  dplyr::group_by(Taxon, phylum, Levels, Loci, Assignment.method ) %>% 
  summarise(N_ESV = n()) %>% 
  dplyr::mutate(Assignment.method = factor(Assignment.method)) %>% 
  ggplot(aes(y = Taxon, x = Assignment.method))+
  geom_tile(aes(fill = N_ESV)) +
  scale_fill_viridis_c() +
  facet_grid(phylum~Loci, scale = "free", space = "free")+
  theme_bw() +
  labs(title = "Detected species combining 2 methods",
       x = "Assignement method")+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 8),
        axis.text.y = element_text(face = "italic", size = 8, vjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.y = element_text(angle = 0, size = 8), 
        panel.spacing.y=unit(0.1, "lines"),
        legend.position = "right",
        legend.text = element_text(size = 8),
        panel.grid = element_blank()
  )

graph.ESV.final

ggsave(filename = file.path("02_Results/03_TaxoAssign", "Comparison_Marine_Final_nESV.png"), plot =  graph.ESV.final, height = 12, width = 10)


# Check and correct taxonomy for WoRMS ------------------------------
#import function validateWORMS
source(file.path(here::here(), "01_Code", "Functions", "fct_WoRMS.R"))

#load Final assignation file
assignations <- bind_rows(readr::read_csv(file.path("02_Results/03_TaxoAssign/Assignements.marine.not.corrected.csv")),
                          readr::read_csv(file.path("02_Results/03_TaxoAssign/Assignements.NT.not.corrected.csv")),
                          readr::read_csv(file.path("02_Results/03_TaxoAssign/Assignements.estuary.not.corrected.csv")))

#remove duplicate Taxon and NA
assignations.unique <- assignations %>% distinct(Taxon, .keep_all = T) %>% 
  select(species,Taxon, Levels, genus, family, order, class, phylum, kingdom, ESV)
assignations.unique <- assignations.unique %>% drop_na(Taxon)

#check if a taxon or a vector of taxon complies with WORMS ID
list.taxon <- assignations.unique$Taxon

valid.taxo <- data.frame()

sep <- seq(from = 1, to = length(list.taxon), by = 10) #separate list to avoid error with taxa without entry in WORMS
sep

for(x in sep){
  cat("\nWe now process row number", x)
  y <- x + 9
  
  if(y > length(list.taxon)){
    y <- length(list.taxon)
  }
  
  valid.taxo.int <- validateWORMS(list.taxon[c(x:y)])
  
  valid.taxo <- bind_rows(valid.taxo, valid.taxo.int)
  
}

readr::write_csv(valid.taxo, file = file.path(here::here(), "02_Results/03_TaxoAssign/Assignements_WoRMS_valid.csv"))

#select and check only the taxa with the "unaccepted" status
unacc.taxo <- valid.taxo %>% filter(status %in% c("unaccepted", "accepted | unaccepted"))

valid.taxo.f <- valid.taxo %>% separate(valid_name, c("valid_name", "B"),sep = "\\|", remove = TRUE) %>% #sep is a regular expression so the | needs to be escaped with \\
  separate(kingdom,  c("kingdom", "C"),sep = "\\|",remove=TRUE) %>%
  separate(phylum,  c("phylum", "D"),sep = "\\|", remove =TRUE) %>%
  separate(class,  c("class", "E"),sep = "\\|", remove=TRUE) %>%
  separate(order,  c("order", "F"),sep = "\\|", remove=TRUE) %>%
  separate(family, c("family", "G"),sep = "\\|", remove =TRUE) %>%
  select(-c(B,C,D,E,F,G)) %>%
  mutate(species = case_when(rank == "Species" ~ valid_name, rank != "Species" ~ NA) ) %>%
  dplyr::rename(Taxon = scientificname)


#add WORMS into the original Assignement file
assignations.marine.valid <- res.marine.final %>%   left_join(
  valid.taxo.f %>%
    select(Taxon, AphiaID, authority, status, valid_name,
           kingdom, phylum, class, order, family, genus, species),
  by = "Taxon",
  suffix = c(".ori", ".worms")
) %>%
  # valid_name peut contenir "Genus species Author, year" ; on récupère le premier mot (Genus)
  tidyr::separate(valid_name, into = c("A", NA), remove = FALSE) %>%
  mutate(
    # Recalcule du genus à partir du valid_name uniquement quand c'est pertinent
    genus_from_validname = case_when(
      Levels %in% c("species", "genus") ~ A,
      TRUE ~ NA_character_
    ),
    # Pour chaque rang, si le validé existe -> on le prend, sinon on garde l'original
    kingdom = coalesce(kingdom.worms, kingdom.ori),
    phylum  = coalesce(phylum.worms,  phylum.ori),
    class   = coalesce(class.worms,   class.ori),
    order   = coalesce(order.worms,   order.ori),
    family  = coalesce(family.worms,  family.ori),
    # genus : priorité au genus recalculé depuis valid_name (quand applicable), 
    # puis au genus validé (si présent), sinon à l'original
    genus   = coalesce(genus_from_validname, genus.worms, genus.ori),
    species = coalesce(species.worms, species.ori),
    # Le libellé principal : si validé connu, on prend valid_name ; sinon on garde Taxon d'origine
    Taxon   = coalesce(valid_name, Taxon)
  ) %>%
  # On garde les métadonnées utiles du côté validé
  mutate(
    AphiaID  = AphiaID,   # si tu veux prioriser un AphiaID validé
    authority = authority,
    status    = status
  ) #%>%
# Nettoyage : on enlève les colonnes intermédiaires
#  select(
#    # champs finaux
#    
#    Taxon, AphiaID, authority, status,
#    kingdom, phylum, class, order, family, genus, species,
# si tu souhaites garder une trace des originaux pour audit, dé-commente ceci :
# kingdom_ori, phylum_ori, class_ori, order_ori, family_ori, genus_ori, species_ori
#  ) %>%
# Si tu veux renommer pour conserver ton schéma initial côté "ori" en mémoire :
# (sinon, saute cette étape)
#  identity()
dup.ESV <- assignations.marine.valid %>% dplyr::filter(duplicated(ESV)) %>% pull(ESV)

assignations.marine.valid %>% dplyr::filter(ESV %in% dup.ESV) %>% View()

assignations.marine.valid %>% distinct(.keep_all = T) %>% nrow()
assignations %>% nrow()

assignations.marine.valid <- assignations.marine.valid %>% distinct(.keep_all = T)

readr::write_csv(assignations.marine.valid, file = file.path(here::here(), "02_Results/03_TaxoAssign/Assignements.marine.Final.csv"))


assignations.NT.valid <- res.NT.final %>%
  left_join(
    valid.taxo.f %>%
      select(Taxon, AphiaID, authority, status, valid_name,
             kingdom, phylum, class, order, family, genus, species),
    by = "Taxon",
    suffix = c(".ori", ".worms")
  ) %>%
  # valid_name peut contenir "Genus species Author, year" ; on récupère le premier mot (Genus)
  tidyr::separate(valid_name, into = c("A", NA), remove = FALSE) %>%
  mutate(
    # Recalcule du genus à partir du valid_name uniquement quand c'est pertinent
    genus_from_validname = case_when(
      Levels %in% c("species", "genus") ~ A,
      TRUE ~ NA_character_
    ),
    # Pour chaque rang, si le validé existe -> on le prend, sinon on garde l'original
    kingdom = coalesce(kingdom.worms, kingdom.ori),
    phylum  = coalesce(phylum.worms,  phylum.ori),
    class   = coalesce(class.worms,   class.ori),
    order   = coalesce(order.worms,   order.ori),
    family  = coalesce(family.worms,  family.ori),
    # genus : priorité au genus recalculé depuis valid_name (quand applicable), 
    # puis au genus validé (si présent), sinon à l'original
    genus   = coalesce(genus_from_validname, genus.worms, genus.ori),
    species = coalesce(species.worms, species.ori),
    # Le libellé principal : si validé connu, on prend valid_name ; sinon on garde Taxon d'origine
    Taxon   = coalesce(valid_name, Taxon)
  ) %>%
  # On garde les métadonnées utiles du côté validé
  mutate(
    AphiaID  = AphiaID,   # si tu veux prioriser un AphiaID validé
    authority = authority,
    status    = status
  ) #%>%
# Nettoyage : on enlève les colonnes intermédiaires
#  select(
#    # champs finaux
#    
#    Taxon, AphiaID, authority, status,
#    kingdom, phylum, class, order, family, genus, species,
# si tu souhaites garder une trace des originaux pour audit, dé-commente ceci :
# kingdom_ori, phylum_ori, class_ori, order_ori, family_ori, genus_ori, species_ori
#  ) %>%
# Si tu veux renommer pour conserver ton schéma initial côté "ori" en mémoire :
# (sinon, saute cette étape)
#  identity()
dup.ESV <- assignations.NT.valid %>% dplyr::filter(duplicated(ESV)) %>% pull(ESV)

assignations.NT.valid %>% dplyr::filter(ESV %in% dup.ESV) %>% View()

assignations.NT.valid %>% distinct(.keep_all = T) %>% nrow()
assignations %>% nrow()

assignations.NT.valid <- assignations.NT.valid %>% distinct(.keep_all = T)


readr::write_csv(assignations.NT.valid, file = file.path(here::here(), "02_Results/03_TaxoAssign/Assignements.NT.Final.csv"))



# https://www.marinespecies.org/aphia.php?p=manual to check the significance of WoRMS status. 




assignations.estuary.valid <- res.estuary.final %>%
  left_join(
    valid.taxo.f %>%
      select(Taxon, AphiaID, authority, status, valid_name,
             kingdom, phylum, class, order, family, genus, species),
    by = "Taxon",
    suffix = c(".ori", ".worms")
  ) %>%
  # valid_name peut contenir "Genus species Author, year" ; on récupère le premier mot (Genus)
  tidyr::separate(valid_name, into = c("A", NA), remove = FALSE) %>%
  mutate(
    # Recalcule du genus à partir du valid_name uniquement quand c'est pertinent
    genus_from_validname = case_when(
      Levels %in% c("species", "genus") ~ A,
      TRUE ~ NA_character_
    ),
    # Pour chaque rang, si le validé existe -> on le prend, sinon on garde l'original
    kingdom = coalesce(kingdom.worms, kingdom.ori),
    phylum  = coalesce(phylum.worms,  phylum.ori),
    class   = coalesce(class.worms,   class.ori),
    order   = coalesce(order.worms,   order.ori),
    family  = coalesce(family.worms,  family.ori),
    # genus : priorité au genus recalculé depuis valid_name (quand applicable), 
    # puis au genus validé (si présent), sinon à l'original
    genus   = coalesce(genus_from_validname, genus.worms, genus.ori),
    species = coalesce(species.worms, species.ori),
    # Le libellé principal : si validé connu, on prend valid_name ; sinon on garde Taxon d'origine
    Taxon   = coalesce(valid_name, Taxon)
  ) %>%
  # On garde les métadonnées utiles du côté validé
  mutate(
    AphiaID  = AphiaID,   # si tu veux prioriser un AphiaID validé
    authority = authority,
    status    = status
  ) #%>%
# Nettoyage : on enlève les colonnes intermédiaires
#  select(
#    # champs finaux
#    
#    Taxon, AphiaID, authority, status,
#    kingdom, phylum, class, order, family, genus, species,
# si tu souhaites garder une trace des originaux pour audit, dé-commente ceci :
# kingdom_ori, phylum_ori, class_ori, order_ori, family_ori, genus_ori, species_ori
#  ) %>%
# Si tu veux renommer pour conserver ton schéma initial côté "ori" en mémoire :
# (sinon, saute cette étape)
#  identity()
dup.ESV <- assignations.estuary.valid %>% dplyr::filter(duplicated(ESV)) %>% pull(ESV)

assignations.estuary.valid %>% dplyr::filter(ESV %in% dup.ESV) %>% View()

assignations.estuary.valid %>% distinct(.keep_all = T) %>% nrow()
assignations %>% nrow()

assignations.estuary.valid <- assignations.estuary.valid %>% distinct(.keep_all = T)


readr::write_csv(assignations.estuary.valid, file = file.path(here::here(), "02_Results/03_TaxoAssign/Assignements.estuary.Final.csv"))



