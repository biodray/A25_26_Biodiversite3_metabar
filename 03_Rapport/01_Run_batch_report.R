# Script pour faire les rapports automatisés


rmarkdown::render(
  input = here::here("03_Rapport/Automatic_Report_ALL_Anse_a_Beaufils.Rmd"),
  output_file = here::here("03_Rapport/R25-26_005_Metabarcoding_ADNe_Projet_JHill_Anse_a_beaufils_Biodiversity_COI_12SMiFish_12S160_Jaclyn_Hill.pdf")
)


rmarkdown::render(
  input = here::here("03_Rapport/Automatic_Report_ALL_NGSL.Rmd"),
  output_file = here::here("03_Rapport/R25-26_008_Metabarcoding_ADNe_Projet_Cabot_Biodiversity_COI_12SMiFish_12S160_Genevieve_Parent.pdf")
)

rmarkdown::render(
  input = here::here("03_Rapport/Automatic_Report_ALL_Nmelanostomus_Parc_Canada.Rmd"),
  output_file = here::here("03_Rapport/R25-26_003_Metabarcoding_ADNe_Projet_Gobie_Canal_Chambly_Biodiversity_COI_12SMiFish_12S160_Unio_Adele_Hoarau.pdf")
)

rmarkdown::render(
  input = here::here("03_Rapport/Automatic_Report_ALL_PPO_Charcot.Rmd"),
  output_file = here::here("03_Rapport/R25-26_007_Metabarcoding_ADNe_Projet_PPO_StPierre_Miquelon_Charcot_Biodiversity_COI_12SMiFish_12S160_Yanick_Gendreau.pdf")
)

rmarkdown::render(
  input = here::here("03_Rapport/Automatic_Report_ALL_Restauration_Restigouche.Rmd"),
  output_file = here::here("03_Rapport/R25-26_006_Metabarcoding_ADNe_Projet_JHill_Restigouche_Biodiversity_COI_12SMiFish_12S160_Jaclyn_Hill.pdf")
)


