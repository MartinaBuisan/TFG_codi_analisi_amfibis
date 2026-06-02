############################################################

# SEGUIMENT D'AMFIBIS ALS PARCS DE LA DIBA
# Tendències poblacionals amb models TRIM (adults)

############################################################

# Autora: Martina Buisán Rodríguez
# Treball Final de Grau:
#   Efectes de la sequera en la comunitat d'amfibis de la XPN de la DIBA

# Descripció:
#   Aquest script aplica models TRIM a les espècies que
#   SÍ compleixen els requisits mínims de dades
# Espècies analitzades:
# -  A. almogavarii, B. spinosus, E. calamita 
# - H. meridionalis, P. perezi, P. punctatus
# - S. salamandra, T. marmoratus, I. alpestris

############################################################

# NETEJA DE L'ENTORN I CÀRREGA DE PAQUETS

rm(list = ls())

library(tidyverse)
library(dplyr)
library(rtrim)

############################################################
# 1. CARREGA I PREPARACIÓ DE DADES
############################################################

inds <- read.csv("totmerged_DIBA_2025.csv", stringsAsFactors = TRUE) %>%
  mutate(
    adults = N_FemAdult + N_MasAdult + N_IndetAdult,
    reclut = ifelse(N_Metam == 0, "no", "yes"),
    reprod = ifelse(N_Postes != 0 | N_Larves != 0, "yes", "no"),
    juveni = N_Juvenils
  ) %>%
  select(Parc, CodiBassa, CodiMostreig, CodiEspecie, TipusMostreig,
         adults, juveni, reclut, reprod, Hidroperiode)

obsr <- read.csv("DadesObserv_DIBA_2025.csv", stringsAsFactors = TRUE) %>%
  mutate(
    Ocasio = str_sub(CodiMostreig, 7, 23),
    site   = str_sub(CodiMostreig, 1, 5)
  ) %>%
  select(CodiMostreig, Ocasio, site)

# Filtrar llocs amb mínim 3 visites
nvisites <- table(obsr$site)
obsr <- obsr %>% filter(site %in% names(nvisites[nvisites > 2]))

############################################################
# 2. OBJECTES D’INTERÈS
############################################################

especies <- c("Aalm","Bspi","Ecal","Hmer","Pper","Ppun","Ssal","Tmar","Ialp")
gspecie  <- c("A. almogavarii","B. spinosus","E. calamita","H. meridionalis",
              "P. perezi","P. punctatus","S. salamandra",
              "T. marmoratus","I. alpestris")

nomcomplet <- data.frame(sp = especies, gspecie = gspecie)
rownames(nomcomplet) <- especies

############################################################
# 3. PREPARACIÓ PER ESPÈCIE
############################################################

trobatsadults_tmo <- aggregate(adults ~ CodiMostreig + CodiEspecie + TipusMostreig,
                               data = inds, sum)

trobatsadults_max <- aggregate(adults ~ CodiMostreig + CodiEspecie,
                               data = trobatsadults_tmo, max)

for(s in especies) {
  ss1 <- filter(trobatsadults_max, CodiEspecie == s)
  
  ss2 <- left_join(obsr, ss1, "CodiMostreig") %>%
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
  
  assign(s, ss2)
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
  anual <- aggregate(count ~ site + year + habitat, data = x, max)
  anual <- arrange(anual, site)
  
  maxobs <- tapply(anual$count, anual$site, max)
  sense <- names(maxobs[maxobs == 0])
  
  anual %>% filter(!site %in% sense)
})

############################################################
# 5. RESUM DE CENSOS
############################################################

resumcensos_sp <- lapply(datalist2, count_summary)

resumcensos_df <- data.frame(
  sites           = sapply(resumcensos_sp, `[[`, "sites"),
  positive_counts = sapply(resumcensos_sp, `[[`, "positive_counts"),
  total_observed  = sapply(resumcensos_sp, `[[`, "total_observed")
)

print(resumcensos_df)

triats <- especies

############################################################
# 6. EXEMPLE TRIM (canviar a espècie d'interès)
############################################################

comptatges_Ialp <- datalist2[["Ialp"]]

print(
  trim(data = comptatges_Ialp,
       count ~ site + (year + season),
       overdisp = TRUE, serialcor = FALSE)
)

############################################################
# 7. TRIMS PER A TOTES LES ESPÈCIES 
############################################################

for(nom in triats) {
  
  comptatges <- datalist2[[nom]]
  
  # Model TRIM
  trim1 <- trim(
    data = comptatges,
    count ~ site + (year + season),
    overdisp = TRUE,
    serialcor = FALSE
  )

  # Resums
  resum_comptatge <- count_summary(comptatges)
  resum_trim      <- summary(trim1)
  conclusio_stat  <- overall(trim1)
  
  # 1. Heatmap imputats
  print(
    heatmap(trim1, what = "imputed",
            main = nomcomplet[nom, "gspecie"], las = 1)
  )
  
  # 2. Observats vs imputats
  tiobs <- totals(trim1, "imputed", obs = TRUE)
  print(
    plot(tiobs,
         main = nomcomplet[nom, "gspecie"],
         leg.pos = "bottomright",
         ylab = "Observats vs. Observats + Imputats",
         xlab = "Any")
  )

  # 3. Modelitzats vs imputats
  tf <- totals(trim1, "fitted")
  ti <- totals(trim1, "imputed")
  print(
    plot(tf, ti,
         names = c("Modelitzats", "Imputats + Observats"),
         main = nomcomplet[nom, "gspecie"],
         leg.pos = "bottomright",
         ylab = "Imputats + Observats vs. Modelitzats",
         xlab = "Any")
  )
  
  # 4. Index formal
  idx1 <- index(trim1, method = "formal")
  print(
    plot(idx1,
         main = nomcomplet[nom, "gspecie"])
  )
  
  # 5. Index scaled
  idx2 <- index(trim1, method = "scaled")
  print(
    plot(idx2,
         main = nomcomplet[nom, "gspecie"])
  )

  # Guardar resultats a l’entorn (no a fitxers)
  assign(
    paste0("trim_", nom),
    list(
      nom = nom,
      resum_comptatge = resum_comptatge,
      modeltrim = trim1,
      resum_trim = resum_trim,
      conclusio_stat = conclusio_stat,
      imputats_observats = tiobs,
      estimats = tf,
      index1 = idx1,
      index2 = idx2
    )
  )
}
