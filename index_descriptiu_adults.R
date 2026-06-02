############################################################

# SEGUIMENT D'AMFIBIS ALS PARCS DE LA DIBA
# Anàlisi de tendències poblacionals (adults)

############################################################

# Autora: Martina Buisán Rodríguez
# Treball Final de Grau: 
#   Efectes de la sequera en la comunitat d'amfibis de la XPN de la DIBA

# Descripció:
  # Aquest script calcula tendències poblacionals 
  # per a espècies amb dades insuficients per aplicar models TRIM

# Espècies analitzades:
# - Lissotriton helveticus
# - Pelobates cultripes
# - Rana temporaria
# - Discoglossus pictus

# Motiu exclusió TRIM:
#   1) Massa zeros per lloc
#   2) Pocs positius per any o lloc
#   3) TRIM elimina gairebé tots els llocs → models no convergeixen

# Decisió final:
#   → Índex simple = (nº individus any) / (nº llocs visitats)

############################################################

# NETEJA DE L'ENTORN I CÀRREGA DE PAQUETS

rm(list=ls())
library(tidyverse)
library(rtrim)

############################################################
# 1. CARREGA I PREPARACIÓ DE DADES
############################################################

# Dades d’individus
inds <- read.csv("totmerged_DIBA_2025.csv", stringsAsFactors = TRUE) %>%
  mutate(
    adults  = N_FemAdult + N_MasAdult + N_IndetAdult,
    reclut  = ifelse(N_Metam == 0, "no", "yes"),
    reprod  = ifelse(N_Postes != 0 | N_Larves != 0, "yes", "no"),
    juveni  = N_Juvenils
  ) %>%
  select(Parc, CodiBassa, CodiMostreig, CodiEspecie, TipusMostreig,
         adults, juveni, reclut, reprod, Hidroperiode)

# Dades d’observacions (visites)
obsr <- read.csv("DadesObserv_DIBA_2025.csv", stringsAsFactors = TRUE) %>%
  mutate(
    Ocasio = str_sub(CodiMostreig, 7, 23),
    site   = str_sub(CodiMostreig, 1, 5)
  ) %>%
  select(CodiMostreig, Ocasio, site)

# Filtrar llocs amb ≥1 visita
nvisites <- table(obsr$site)
obsr <- obsr %>% filter(site %in% names(nvisites[nvisites > 0]))

############################################################
# 2. OBJECTES D’INTERÈS
############################################################

especies  <- c("Lhel", "Pcul", "Rtem", "Dpic")
gspecie   <- c("L. helveticus", "P. cultripes", "R. temporaria", "D. pictus")

nomcomplet <- tibble(sp = especies, gspecie = gspecie)

############################################################
# 3. PREPARACIÓ PER ESPÈCIE: omplir zeros i estructurar dades
############################################################

# Comptatge màxim d’adults per mostreig (visual + acústic)
trobats <- inds %>%
  group_by(CodiMostreig, CodiEspecie, TipusMostreig) %>%
  summarise(adults = sum(adults), .groups = "drop") %>%
  group_by(CodiMostreig, CodiEspecie) %>%
  summarise(adults = max(adults), .groups = "drop")

# Crear objectes per espècie
for (s in especies) {
  df <- trobats %>% filter(CodiEspecie == s)
  
  sp_df <- obsr %>%
    left_join(df, by = "CodiMostreig") %>%
    mutate(
      CodiEspecie = ifelse(is.na(CodiEspecie), s, CodiEspecie),
      adults      = ifelse(is.na(adults), 0, adults),
      site        = factor(str_sub(CodiMostreig, 1, 5)),
      year        = as.numeric(str_sub(CodiMostreig, 7, 10)),
      season      = factor(str_sub(CodiMostreig, 12, 13)),
      park        = factor(str_sub(CodiMostreig, 1, 3)),
      count       = adults,
      species     = CodiEspecie
    ) %>%
    select(site, year, count, park, season, species)
  
  assign(s, sp_df)
}

datalist <- mget(especies)

############################################################
# 4. MILLORA I AGREGA PER ANY I LLOC
############################################################

datalist2 <- lapply(datalist, function(x) {
  x %>%
    mutate(
      habitat = park,
      site_year_season = interaction(site, year, season, sep = "_")
    )
})

datalistyear <- lapply(datalist2, function(x) {
  anual <- x %>%
    group_by(site, year, habitat) %>%
    summarise(count = max(count), .groups = "drop")
  
  basses_sense <- anual %>%
    group_by(site) %>%
    summarise(max_count = max(count), .groups = "drop") %>%
    filter(max_count == 0) %>%
    pull(site)
  
  anual %>% filter(!site %in% basses_sense)
})

############################################################
# 5. RESUM DE CENSOS PER ESPÈCIE
############################################################

resumcensos_sp <- lapply(datalist2, count_summary)

resumcensos_df <- tibble(
  species         = especies,
  sites           = sapply(resumcensos_sp, `[[`, "sites"),
  positive_counts = sapply(resumcensos_sp, `[[`, "positive_counts"),
  total_observed  = sapply(resumcensos_sp, `[[`, "total_observed")
)

print(resumcensos_df)

############################################################
# 6. CÀLCUL I GRÀFICS D’ÍNDEX SIMPLE
############################################################

triats <- c("Lhel","Pcul","Rtem","Dpic")

for(nom in triats) {
  
  comptatges <- datalist2[[nom]]
  
  idx_simple <- comptatges %>%
    group_by(year) %>%
    summarise(
      total = sum(count, na.rm = TRUE),
      n_sites = n_distinct(site),
      index = total / n_sites,
      .groups = "drop"
    ) %>%
    arrange(year)
  
  assign(paste0("index_simple_", nom), idx_simple)
  
  # --- MOSTRAR GRÀFIC ---
  print(
    ggplot(idx_simple, aes(year, index)) +
      geom_line(colour = "darkred", linewidth = 1) +
      geom_point(colour = "darkred", size = 2) +
      labs(
        x = "Any",
        y = "Índex simple (individus / lloc)",
        title = paste("Índex simple -", nom)
      ) +
      theme_bw()
  )
  
  cat("Índex simple calculat i gràfic mostrat per:", nom, "\n")
}

############################################################
# 7. GRÀFIC CONJUNT DE LES 4 ESPÈCIES
############################################################

index_simple_all <- lapply(triats, function(sp) {
  df <- get(paste0("index_simple_", sp))
  df$species <- sp
  df
}) %>% bind_rows() %>%
  left_join(nomcomplet, by = c("species" = "sp"))

# --- MOSTRAR GRÀFIC ---
print(
  ggplot(index_simple_all,
         aes(x = year, y = index, colour = gspecie, group = gspecie)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    labs(
      x = "Any",
      y = "Índex simple (individus / lloc)",
      colour = "Espècie",
      title = "Índex simple anual de les 4 espècies"
    ) +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
)

############################################################
# 8. GRÀFICS PER PARC I BASSA (només llocs amb deteccions)
############################################################

for(nom in triats) {
  
  dades_sp <- datalistyear[[nom]]
  
  basses_positives <- dades_sp %>%
    group_by(site) %>%
    summarise(max_count = max(count, na.rm = TRUE)) %>%
    filter(max_count > 0) %>%
    pull(site)
  
  if(length(basses_positives) == 0) {
    message("No hi ha basses positives per a ", nom)
    next
  }
  
  dades_pos <- dades_sp %>% filter(site %in% basses_positives)
  parcs_positius <- sort(unique(dades_pos$habitat))
  
  for(parc in parcs_positius) {
    
    dades_parc <- dades_pos %>% filter(habitat == parc)
    if(nrow(dades_parc) == 0) next
    
    # --- MOSTRAR GRÀFIC ---
    print(
      ggplot(dades_parc,
             aes(x = year, y = count, group = site, colour = site)) +
        geom_line(linewidth = 1) +
        geom_point(size = 2) +
        labs(
          title = paste0("Sèrie temporal - ", nom, " - Parc ", parc),
          subtitle = "Només basses on s'ha detectat l'espècie",
          x = "Any",
          y = "Nombre d'individus (màxim anual per bassa)",
          colour = "Bassa"
        ) +
        scale_x_continuous(breaks = sort(unique(dades_parc$year))) +
        theme_bw() +
        theme(panel.grid.minor = element_blank())
    )
  }
}
