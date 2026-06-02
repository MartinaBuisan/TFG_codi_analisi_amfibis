############################################################

# SEGUIMENT D'AMFIBIS ALS PARCS DE LA DIBA
# Anàlisi de patrons reproductors (postes, larves i metamòrfics)

############################################################

# Autora: Martina Buisán Rodríguez
# Treball Final de Grau: 
#   Efectes de la sequera en la comunitat d'amfibis de la XPN de la DIBA

# Descripció:
#   Aquest script genera heatmaps per a totes les espècies monitoritzades,
#   mostrant la intensitat reproductora (postes, larves, metamòrfics)
#   per parc, any i temporada.

# Variables representades:
#   - Postes (N_Postes_num)
#   - Larves (N_Larves_num)
#   - Metamòrfics (N_Metam_num)

############################################################

# NETEJA DE L'ENTORN I CÀRREGA DE PAQUETS

rm(list = ls())

library(tidyverse)
library(stringr)
library(patchwork)

############################################################
# 1. CARREGA I PREPARACIÓ DE DADES
############################################################

inds <- read.csv("taulaubis_postes_script2.csv", stringsAsFactors = FALSE)

df <- inds %>%
  mutate(
    park   = str_sub(CodiMostreig, 1, 3),
    year   = as.numeric(str_sub(CodiMostreig, 7, 10)),
    season = as.numeric(str_sub(CodiMostreig, 12, 13)),
    species = CodiEspecie
  ) %>%
  group_by(park, year, season, species) %>%
  summarise(
    postes = max(N_Postes_num, na.rm = TRUE),
    larves = max(N_Larves_num, na.rm = TRUE),
    metam  = max(N_Metam_num,  na.rm = TRUE),
    .groups = "drop"
  )

# Completar combinacions inexistents amb zeros
df_full <- df %>%
  complete(
    park, year, season, species,
    fill = list(postes = 0, larves = 0, metam = 0)
  )

############################################################
# 2. LLISTA D’ESPÈCIES I NOM CIENTÍFIC
############################################################

especies <- c("Aalm","Bspi","Ecal","Hmer","Pper","Ppun","Ssal",
              "Tmar","Ialp","Lhel","Pcul","Rtem","Dpic")

gspecie <- c("A. almogavarii","B. spinosus","E. calamita","H. meridionalis",
             "P. perezi","P. punctatus","S. salamandra",
             "T. marmoratus","I. alpestris","L. helveticus",
             "P. cultripes","R. temporaria","D. pictus")

nomcomplet <- data.frame(sp = especies, gspecie = gspecie)

############################################################
# 3. FUNCIÓ PER GENERAR HEATMAPS
############################################################

plot_heatmap <- function(df_sp, variable, titulo) {
  
  df_sp <- df_sp %>%
    mutate(year_season = paste(year, sprintf("%02d", season), sep = "_")) %>%
    arrange(year, season)
  
  orden_x <- df_sp %>%
    distinct(year, season, year_season) %>%
    arrange(year, season) %>%
    pull(year_season)
  
  df_sp$year_season <- factor(df_sp$year_season, levels = orden_x)
  
  ggplot(df_sp, aes(x = year_season, y = park, fill = .data[[variable]])) +
    geom_tile() +
    scale_fill_gradient2(
      low = "pink",
      mid = "white",
      high = "red",
      midpoint = median(df_sp[[variable]], na.rm = TRUE),
      name = variable,
      na.value = "white"
    ) +
    labs(
      title = titulo,
      x = "Any - Temporada",
      y = "Parc"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5),
      panel.grid = element_blank()
    )
}

############################################################
# 4. HEATMAPS PER A TOTES LES ESPÈCIES 
############################################################

for(i in seq_along(especies)) {
  
  sp <- especies[i]
  df_sp <- df_full %>% filter(species == sp)
  
  print(plot_heatmap(df_sp, "postes",
                     paste("Heatmap postes -", nomcomplet$gspecie[i])))
  
  print(plot_heatmap(df_sp, "larves",
                     paste("Heatmap larves -", nomcomplet$gspecie[i])))
  
  print(plot_heatmap(df_sp, "metam",
                     paste("Heatmap metamòrfics -", nomcomplet$gspecie[i])))
}

############################################################
# 5. AGRUPACIÓ D’ESPÈCIES PER HIDROPERÍODE
############################################################

grup_permanents <- c("Pper","Tmar","Rtem")
grup_intermedis <- c("Ssal","Aalm","Bspi","Hmer","Pcul")
grup_efimers <- c("Ppun","Ecal")
grup_poc_mostrejats <- c("Ialp","Lhel","Dpic")

############################################################
# 6. CREACIÓ DE PLOTS PER ESPÈCIE 
############################################################

plots_por_especie <- list()

for(i in seq_along(especies)) {
  
  sp <- especies[i]
  nom <- nomcomplet$gspecie[i]
  
  df_sp <- df_full %>% filter(species == sp)
  
  plots_por_especie[[sp]] <- list(
    postes = plot_heatmap(df_sp, "postes", paste("Postes -", nom)),
    larves = plot_heatmap(df_sp, "larves", paste("Larves -", nom)),
    metam  = plot_heatmap(df_sp, "metam",  paste("Metamòrfics -", nom))
  )
}

############################################################
# 7. QUADRÍCULES PER GRUPS D’ESPÈCIES 
############################################################

# PERMANENTS
cuadricula_permaments <-
  wrap_plots(lapply(grup_permanents, function(sp) plots_por_especie[[sp]]$postes), nrow = 1) /
  wrap_plots(lapply(grup_permanents, function(sp) plots_por_especie[[sp]]$larves), nrow = 1) /
  wrap_plots(lapply(grup_permanents, function(sp) plots_por_especie[[sp]]$metam),  nrow = 1) +
  plot_annotation(title = "Espècies de hidroperíode permanent")

cuadricula_permaments

# INTERMEDIS
cuadricula_intermedis <-
  wrap_plots(lapply(grup_intermedis, function(sp) plots_por_especie[[sp]]$postes), nrow = 1) /
  wrap_plots(lapply(grup_intermedis, function(sp) plots_por_especie[[sp]]$larves), nrow = 1) /
  wrap_plots(lapply(grup_intermedis, function(sp) plots_por_especie[[sp]]$metam),  nrow = 1) +
  plot_annotation(title = "Espècies de hidroperíode intermedi")

cuadricula_intermedis

# EFÍMERS
cuadricula_efimers <-
  wrap_plots(lapply(grup_efimers, function(sp) plots_por_especie[[sp]]$postes), nrow = 1) /
  wrap_plots(lapply(grup_efimers, function(sp) plots_por_especie[[sp]]$larves), nrow = 1) /
  wrap_plots(lapply(grup_efimers, function(sp) plots_por_especie[[sp]]$metam),  nrow = 1) +
  plot_annotation(title = "Espècies de hidroperíode efímer")

cuadricula_efimers

# POC MOSTREJATS
cuadricula_poc_mostrejats <-
  wrap_plots(lapply(grup_poc_mostrejats, function(sp) plots_por_especie[[sp]]$postes), nrow = 1) /
  wrap_plots(lapply(grup_poc_mostrejats, function(sp) plots_por_especie[[sp]]$larves), nrow = 1) /
  wrap_plots(lapply(grup_poc_mostrejats, function(sp) plots_por_especie[[sp]]$metam),  nrow = 1) +
  plot_annotation(title = "Espècies poc mostrejades")

cuadricula_poc_mostrejats
