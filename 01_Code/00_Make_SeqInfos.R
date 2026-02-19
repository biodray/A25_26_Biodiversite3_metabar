# Info --------------------------------------------------------------------

# code to extract information from excel/ACCESS files and GQ file to generate the SeqInfo.csv file

# M Chevrinais & A Bourret
# 2024-04-18
# Last Modified : 2024-08-14 (AB)

# Library -----------------------------------------------------------------
library(tidyverse)
library(dplyr)

# Upload excel/ACCESS file####
#packages####
#list.of.packages <- c("RODBC", "tidyverse", "remotes", "xlsx")
#new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
#if(length(new.packages)) install.packages(new.packages)

library(RODBC)
library(remotes)
library(xlsx)

# Package maison pour faciliter l'utilisation de R pour accéder aux données
#remotes::install_github("GenomicsMLI-DFO/BD-LabGeno-IML_gestion",
#                       upgrade = c("never")) # Can be changed if you want to dependencies to be updated

library(BDLG.gestion)

# Aller chercher les métadonnées dans la db ACCESS####
# Étape 1 : Vérifier la connection ----------------------------------------

# Cette fonction vérifies seulement qu'on est capable d'établit une connection avec LabGeno
# Si ça ne fonctionne pas, vérifier qu'on a accès au S, et que le ODBC est bien configuré

test_DB()

# Étape 2 : Lister les données disponibles --------------------------------

# Cette fonction permet de listes les tables et les requêtes disponibles pour le téléchargement
list_DB()

# Étape 3 : Télécharger les données qui nous intéressent --------------------

# Ensuite, une fois qu'on connait le nom des tables et/ou des requêtes qui nous intéressent,
# on peut les télécharger et les assigner à des objets qui vont se retrouver dans notre
# environnement R (mémoire). Ici on va comparer la table Analyse_externe à la requête Metadata_Analyses_Externes

metadata <- load_DB("R_Metadata_Metabarcoding_ALL")

#data2 <- load_DB("Metadata_Analyses_Externes")

#Formatter les métadonnées pour avoir toutes les colonnes sur un seul df####
#select

metadata |> names()
metadata |> group_by(Projet_GQ, No_soumission_GQ) |> summarise(N = n())

metadata <- metadata %>% filter(Projet_GQ %in% "A25_26_Biodiversite_3") |> distinct(.keep_all = T)

metadata |> group_by(Nom_projet_librairie_ADNe, Nom_sous_projet_librairies_ADNe, Loci) |>
  summarise(N = n()) # |> View()

metadata  |> dplyr::filter(Nom_sous_projet_librairies_ADNe == "Twells_Innu_Nation", Type_echantillon_librairie_ADNe == "FNC" ) |> View()

metadata$Loci |> table()

metadata <- metadata |> dplyr::mutate(Loci = ifelse(Loci == "COI_Elbrecht_B1", "COIB1",
                                                    ifelse(Loci == "18SV9_Stoeck", "18SV9", Loci)),
                                      Loci = str_remove_all(Loci, "_"))


#metadata$Loci[metadata$Loci=="COI_Leray"] <- "COI"
#metadata$Loci[metadata$Loci=="16SVert"] <- "Vert16S"
metadata$ID_Nanuq |> unique() |> length()
metadata$Numero_unique_extrait_ADNe |> unique() |> length()
metadata$Numero_unique_extrait |> unique() |> length()

c(metadata$Numero_unique_extrait_ADNe, metadata$Numero_unique_extrait) |> unique() |> length()

metadata |> group_by(Nom_projet_librairie_ADNe, Loci) |> summarise(N = n())
metadata |> group_by(Nom_projet_librairie_ADNe, Nom_sous_projet_librairies_ADNe, Loci) |> summarise(N = n())

metadata |> group_by(Nom_projet_librairie_ADNe, ID_Nanuq, Numero_unique_extrait_ADNe, Numero_unique_extrait) |> summarise(N = n()) |>
  arrange(N) |> View()

metadata |> names()

lib.ADN <- metadata |> dplyr::filter(Nom_projet_librairie_ADNe == "Contenus_Pborealis") |>
  dplyr::select(ID_sample = ID_Nanuq,
                ID_extrait = Numero_unique_extrait,
                Sample_type = Type_echantillon_librairie_ADNe,
                ID_project = Nom_projet_librairie_ADNe,
                ID_subproject = Nom_commun,
                Loci,
                Inhibition = Inhibition_ADNe.1,
                Cq_PCR1,
                Nombre_cycles_PCR1,
                Kit_extraction,
                ID_plate = Set_index,
                ID_well = Puits_index,
                Index_i5,
                Index_i7,
                Region = Region_echantillonnage,
                Site = Lieu_echantillonnage,
                Latitude = Latitude_echantillonnage_DD,
                Longitude = Longitude_echantillonnage_DD,
                Nom_commun,
                ID_sample_ori = Numero_reception_specimen
  )

lib.ADNe <- metadata |> #dplyr::filter(Nom_projet_librairie_ADNe !=  "Contenus_Pborealis") |>
  dplyr::select(ID_sample = Numero_unique_extrait_ADNe, #ID_Nanuq,
                ID_extrait = Numero_unique_extrait_ADNe,
                Sample_type = Type_echantillon_librairie_ADNe,
                ID_project = Nom_projet_librairie_ADNe,
                ID_subproject = Nom_sous_projet_librairies_ADNe,
                Loci,
                Inhibition = Inhibition_ADNe,
                Cq_PCR1,
                Nombre_cycles_PCR1,
                Kit_extraction = Kit_extraction_ADNe,
                ID_plate = Set_index,
                ID_well = Puits_index,
                Index_i5,
                Index_i7,

                Region = Region_ADNe,
                Site = Site_echantillonnage_ADNe,
                Sampling_site = Nom_site_reception_ADNe,
                Latitude = Latitude_ADNe_DD,
                Longitude = Longitude_ADNe_DD,
                Trait_echantillonnage = Trait_echantillonnage_ADNe,
                ID_sample_ori = Nom_reception_echantillon_ADNe,
                Annee_echantillonnage = Annee_echantillonnage_ADNe,
                Mois_echantillonnage = Mois_echantillonnage_ADNe,
                Jour_echantillonnage = Jour_echantillonnage_ADNe,
                Heure_prelevement = Heure_prelevement_echantillon_ADNe,
                Volume_filtre_L,
                Systeme_filtration,
                Type_filtre,
                Pore_filtre_um,
                Taille_filtre_mm

  ) |> mutate(ID_subproject = ifelse(ID_subproject == "Twells_Innu_Nation", paste(ID_subproject, Annee_echantillonnage, sep = "_"), ID_subproject))

lib.total <- bind_rows(lib.ADNe )



# Donnees GQ --------------------------------------------------------------


GQ.data <- readxl::read_excel("00_Data/00_FileInfos/A25_26_Biodiversite3_NovaSeqReadSet_2026-02-12.xlsx") |>
  mutate(Run2 = paste(Run, Region, sep = "."))
GQ.data |> summary()


##save metadata file
#write.csv(metadata, "00_Data/00_FileInfos/user_input/metadata.csv", row.names=TRUE)

data.info  <- lib.total|> left_join(GQ.data |> dplyr::select(ID_sample = `Library Name`, Run = Run2, File = `Filename Prefix`)) |>
  mutate(ID_subproject = ifelse(is.na(ID_subproject), "ALL", ID_subproject))

View(data.info)

data.info  <- data.info |> mutate(#tag_fwd = Index_i7, tag_rev = Index_i5,
  #Station = Nom_reception_station_ADNe,
  #Site_echantillonnage = Site_echantillonnage_ADNe,
  #ID_puit = Puits_index,
  ID_plate = ifelse(ID_plate == "A", 1,
                    ifelse(ID_plate == "B", 2,
                           ifelse(ID_plate == "C", 3,
                                  ifelse(ID_plate == "D", 4, 99
                                  )))))

data.info |>  nrow()
lib.total |> nrow()

data.info |> group_by(ID_project, ID_subproject, ) |> summarise(N = n())

data.info <- data.info |> distinct(.keep_all = T)

readr::write_csv(data.info, "00_Data/00_FileInfos/SeqInfo.csv")
