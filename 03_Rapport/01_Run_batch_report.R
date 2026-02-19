# Script pour faire les rapports automatisés

rmarkdown::render(
  input = here::here("03_Rapport/Automatic_Report_ALL_BDA.Rmd"),
  output_file = here::here("03_Rapport/R25-26_016_Metabarcoding_ADNe_Projet_BDA_Biodiversity_COI_12SMiFish_12S160_Geneviève_Faille.pdf")
)

rmarkdown::render(
  input = here::here("03_Rapport/Automatic_Report_ALL_Biodiversite_Parc_Marin.Rmd"),
  output_file = here::here("03_Rapport/R25-26_015_Metabarcoding_ADNe_Projet_Parc_marin_Biodiversity_COI_12SMiFish_12S160_Samuel_Turgeon.pdf")
)

rmarkdown::render(
  input = here::here("03_Rapport/Automatic_Report_ALL_Twells_Innu_Nation_2024.Rmd"),
  output_file = here::here("03_Rapport/R25-26_014a_Metabarcoding_ADNe_Projet_AMP_Terre-Neuve_Biodiversity_COI_12SMiFish_12S160_18SV9_rbcL_vertes+rouges_Terri_Wells_2024.pdf")
)

rmarkdown::render(
  input = here::here("03_Rapport/Automatic_Report_ALL_Twells_Innu_Nation_2025.Rmd"),
  output_file = here::here("03_Rapport/R25-26_014b_Metabarcoding_ADNe_Projet_AMP_Terre-Neuve_Biodiversity_COI_12SMiFish_12S160_18SV9_rbcL_vertes+rouges_Terri_Wells_2025.pdf")
)

