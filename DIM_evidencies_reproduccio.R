############################################################

# SEGUIMENT D'AMFIBIS ALS PARCS DE LA DIBA
# Anàlisi de les evidències de reproducció (postes, larves i metamòrfics)
# mitjançant el model DIM proposat per MacKenzie (2003)

############################################################

# Autora: Martina Buisán Rodríguez
# Treball Final de Grau: 
#   Efectes de la sequera en la comunitat d'amfibis de la XPN de la DIBA

# Descripció:
#   Aquest script implementa un model d'ocupació multiseason (colext)
#   per a les espècies que SÍ compleixen els requisits mínims de dades.

#   Inclou:
#     - Preparació de dades (detecció i covariables)
#     - Construcció de matrius y, obsCovs, siteCovs i yearlySiteCovs
#     - Model base i model complet
#     - Bucle automàtic de selecció de models (fins a 32.768 combinacions)
#     - Selecció dels millors models per espècie

############################################################

# NETEJA DE L'ENTORN I CÀRREGA DE PAQUETS

rm(list = ls())

library(unmarked)
library(MuMIn)
library(texreg)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(zoo)
library(ggplot2)
library(patchwork)
library(scales)
library(reshape2)
library(fastDummies)

############################################################
# 1. CARREGA I PREPARACIÓ DE DADES
############################################################

ap   <- read.csv("TOTAL_binari.csv", stringsAsFactors = FALSE)
covs <- read.csv("DadesObserv_DIBA_2025.csv", stringsAsFactors = FALSE)

############################################################
# 2. CORRECCIONS AP
############################################################

ap <- subset(ap, select = -NA.)
ap <- ap %>% mutate(across(everything(), ~replace_na(.x, 0)))

colSums(ap[, 4:ncol(ap)])

############################################################
# 3. CORRECCIÓ DE COVARIABLES
############################################################

### 3.1 ECELS — imputació temporal
covs$Data <- as.Date(covs$Data)
covs <- covs[order(covs$Data), ]

down <- na.locf(covs$ECELS, na.rm = FALSE)
up   <- na.locf(covs$ECELS, fromLast = TRUE, na.rm = FALSE)

ECELS_final <- covs$ECELS

for (i in which(is.na(covs$ECELS))) {
  prev_idx <- max(which(!is.na(covs$ECELS[1:i])))
  next_candidates <- which(!is.na(covs$ECELS[i:nrow(covs)]))
  if (length(next_candidates) == 0) {
    ECELS_final[i] <- down[i]
    next
  }
  next_idx <- i - 1 + min(next_candidates)
  dist_prev <- abs(as.numeric(covs$Data[i] - covs$Data[prev_idx]))
  dist_next <- abs(as.numeric(covs$Data[next_idx] - covs$Data[i]))
  ECELS_final[i] <- ifelse(dist_prev <= dist_next, down[i], up[i])
}

covs$ECELS <- ECELS_final

### 3.2 HIDROPERIODE — mode per bassa
safe_mode <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

covs <- covs %>%
  group_by(CodiBassa) %>%
  mutate(
    group_mode = safe_mode(Hidroperiode),
    fill_value = ifelse(is.na(group_mode), "B", group_mode),
    Hidroperiode = coalesce(Hidroperiode, fill_value)
  ) %>%
  select(-group_mode, -fill_value)

### 3.3 SUPERFÍCIE
covs <- covs %>%
  mutate(
    Superf_temp = Superfmax_m2 %>%
      as.character() %>%
      map_dbl(~{
        x <- .
        if (is.na(x)) x <- ""
        if (str_detect(tolower(x), "tram")) return(2)
        nums <- str_extract_all(x, "\\d+\\.?\\d*")[[1]]
        if (length(nums) == 0) return(NA_real_)
        nums <- as.numeric(nums)
        if (length(nums) == 1) return(nums)
        return(mean(nums))
      })
  )
covs <- covs %>%
  mutate(
    Superfmax_m2 = case_when(
      is.na(Superf_temp) & CodiBassa == "GUI09" ~ 90,
      is.na(Superf_temp) & CodiBassa == "GUI12" ~ 130,
      is.na(Superf_temp) & CodiBassa == "STL28" ~ 160,
      is.na(Superf_temp) & CodiBassa == "STL31" ~ 220,
      is.na(Superf_temp) & CodiBassa == "STL32" ~ 420,
      TRUE ~ Superf_temp
    )
  ) %>%
  select(-Superf_temp)

### 3.4 PEIXOS / CRANC / FAUNA DOMÈSTICA
covs$PresenciaPeixosExotics[is.na(covs$PresenciaPeixosExotics)] <- "B"
covs$PresenciaCranc[is.na(covs$PresenciaCranc)] <- "B"
covs$Presencia_fauna_dom[is.na(covs$Presencia_fauna_dom)] <- "B"

### 3.5 CUBETA — interpolació temporal
covs <- covs %>%
  arrange(CodiBassa, Data) %>%
  group_by(CodiBassa) %>%
  mutate(
    Data_num = as.numeric(Data) + seq_along(Data) * 1e-6,
    CubetaAmbAigua_percent = na.approx(CubetaAmbAigua_percent,
                                       x = Data_num,
                                       na.rm = FALSE,
                                       rule = 2),
    CubetaAmbAigua_percent = na.locf(CubetaAmbAigua_percent, na.rm = FALSE),
    CubetaAmbAigua_percent = na.locf(CubetaAmbAigua_percent, fromLast = TRUE, na.rm = FALSE)
  ) %>%
  select(-Data_num) %>%
  ungroup()

### 3.6 TERBOLESA
covs$Terbolesa[is.na(covs$Terbolesa)] <- "A"

### 3.7 COLUMNA OBSERVABLE
covs <- covs %>%
  mutate(
    percen_Columnaobservable = ifelse(
      is.na(percen_Columnaobservable) | percen_Columnaobservable == "-",
      100,
      percen_Columnaobservable
    ),
    percen_Columnaobservable = str_replace(percen_Columnaobservable, "<", ""),
    percen_Columnaobservable = as.numeric(as.character(percen_Columnaobservable))
  )

############################################################
# 4. RECOMPTES I FILTRATGE DE VISITES
############################################################

covs <- covs %>%
  mutate(
    CodiBassa = str_extract(CodiMostreig, "^[A-Z]+\\d+"),
    Any = str_extract(CodiMostreig, "(?<=_)\\d{4}(?=_)") %>% as.numeric()
  )

visites_any <- covs %>%
  group_by(CodiBassa, Any) %>%
  summarise(N_visites = n(), .groups = "drop")

tots_anys <- expand.grid(
  CodiBassa = unique(covs$CodiBassa),
  Any = 2021:2025
)

visites_any_complet <- tots_anys %>%
  left_join(visites_any, by = c("CodiBassa", "Any")) %>%
  mutate(N_visites = ifelse(is.na(N_visites), 0, N_visites))

visites_totals <- visites_any_complet %>%
  group_by(CodiBassa) %>%
  summarise(Visites_5anys = sum(N_visites))

umbral <- 6
basses_filtrades <- visites_totals %>%
  filter(Visites_5anys > umbral) %>%
  pull(CodiBassa)

covs_filtrat <- covs %>% filter(CodiBassa %in% basses_filtrades)
ap_filtrat   <- ap   %>% filter(CodiBassa %in% basses_filtrades)

############################################################
# 5. DUPLICACIÓ DE CODIS (4 visites per any)
############################################################

DUPLICATCODIS <- function(df) {
  df_sep <- df %>%
    separate(CodiMostreig, into = c("bassa", "year", "visit"), sep = "_") %>%
    mutate(
      year = as.integer(year),
      visit = sprintf("%02d", as.integer(visit))
    )
  
  basses <- sort(unique(df_sep$bassa))
  
  disseny <- expand_grid(
    bassa = basses,
    year = 2021:2025,
    visit = sprintf("%02d", 1:4)
  ) %>%
    unite("CodiMostreig", bassa, year, visit, sep = "_")
  
  df_sep <- df_sep %>%
    unite("CodiMostreig", bassa, year, visit, sep = "_")
  
  df_complet <- disseny %>%
    left_join(df_sep, by = "CodiMostreig")
  
  return(df_complet)
}

ap   <- DUPLICATCODIS(ap_filtrat)
covs <- DUPLICATCODIS(covs_filtrat)  

############################################################
# 6. AFEGIR PARC I CODIBASSA
############################################################

rellenar_codis <- function(df) {
  df %>%
    mutate(
      Parc = if_else(is.na(Parc), substr(CodiMostreig, 1, 3), Parc),
      CodiBassa = if_else(is.na(CodiBassa), substr(CodiMostreig, 1, 5), CodiBassa)
    )
}

ap   <- rellenar_codis(ap)
covs <- rellenar_codis(covs)

############################################################
# 7. AFEGIR SALABRE I PRESÈNCIA DE L’ESPÈCIE
############################################################

especie <- "Pper" #Canviar per Ssal, Aalm, Ppun, Bspi, Pper si escau

taulasalabre <- read.csv("totmerged_DIBA_2025.csv", stringsAsFactors = FALSE)

taulasalabre <- taulasalabre %>%
  mutate(
    TipusMostreig = ifelse(is.na(TipusMostreig), "A", TipusMostreig)
  )

taula_sp <- taulasalabre %>%
  group_by(CodiMostreig) %>%
  summarise(
    TipusMostreig = case_when(
      any(TipusMostreig == "S") ~ "B",
      any(TipusMostreig %in% c("A", "V")) ~ "A",
      TRUE ~ "A"
    ),
    .groups = "drop"
  )

covs <- covs %>%
  left_join(taula_sp, by = "CodiMostreig")

ap_sp <- ap %>%
  select(CodiMostreig, Parc, all_of(especie))

covs <- covs %>% distinct(CodiMostreig, .keep_all = TRUE)

dat <- covs %>%
  left_join(ap_sp, by = "CodiMostreig") %>%
  mutate(
    site = CodiBassa
  )

############################################################
# 8. CREAR ÍNDEXS DE MOSTREIG
############################################################

dat <- dat %>%
  mutate(
    any = as.numeric(str_extract(CodiMostreig, "\\d{4}")),
    primary = as.numeric(factor(any, ordered = TRUE)),
    secondary = as.numeric(str_extract(CodiMostreig, "\\d{2}$"))
  )

Primary  <- length(unique(dat$primary))
Secondary <- max(dat$secondary)
Sites    <- length(unique(dat$site))

############################################################
# 9. MATRIU DE DETECCIÓ
############################################################

dat <- dat %>%
  mutate(occasion = (primary - 1) * Secondary + secondary)

y_wide <- dat %>%
  group_by(site, occasion) %>%
  slice(1) %>%
  ungroup() %>%
  select(site, occasion, all_of(especie)) %>%
  pivot_wider(id_cols = site, names_from = occasion, values_from = all_of(especie)) %>%
  arrange(site)

sites <- y_wide$site
y_mat <- y_wide %>% select(-site) %>% as.matrix()

############################################################
# 10. AGRUPACIÓ DE NIVELLS
############################################################

dat <- dat %>%
  mutate(
    ZONA = case_when(
      Parc.x %in% c("FOX", "OLE", "GRF") ~ "SUD",
      Parc.x %in% c("GUI", "MTQ") ~ "NORD",
      Parc.x %in% c("PSL", "STL") ~ "CENTRE",
      TRUE ~ NA_character_
    ),
    DepredadorsInvasors = case_when(
      PresenciaPeixosExotics == "A" | PresenciaCranc == "A" ~ "A",
      TRUE ~ "B"
    )
  )

############################################################
# 11. OBSERVATION COVARIATES
############################################################

make_obs_matrix <- function(var){
  dat %>%
    group_by(site, occasion) %>%
    slice(1) %>%
    ungroup() %>%
    select(site, occasion, all_of(var)) %>%
    pivot_wider(id_cols = site, names_from = occasion, values_from = all_of(var)) %>%
    arrange(site) %>%
    select(-site) %>%
    as.matrix()
}

obsCovs_list <- list(
  terbol = make_obs_matrix("Terbolesa"),
  columnaob = make_obs_matrix("percen_Columnaobservable"),
  salabre = make_obs_matrix("TipusMostreig")
)

############################################################
# 12. SITE COVARIATES
############################################################

siteCovs_df <- dat %>%
  group_by(site) %>%
  summarise(
    ZONA = first(na.omit(ZONA)),
    Hidroperiode = first(na.omit(Hidroperiode)),
    Superfmax_m2 = first(na.omit(Superfmax_m2))
  )

############################################################
# 13. YEARLY SITE COVARIATES
############################################################

dat <- dat %>%
  group_by(CodiBassa, any) %>%
  mutate(
    ECELS_MITJANA = {
      m <- mean(ECELS, na.rm = TRUE)
      if (is.nan(m)) NA_real_ else m
    },
    CUBETA_MITJANA = {
      m <- mean(CubetaAmbAigua_percent, na.rm = TRUE)
      if (is.nan(m)) NA_real_ else m
    }
  ) %>%
  ungroup()

dat <- dat %>%
  group_by(CodiBassa, any) %>%
  mutate(
    depredadors_invasors = if_else(any(DepredadorsInvasors == "A", na.rm = TRUE), "A", "B")
  ) %>%
  ungroup()

covs_yearly <- dat %>%
  distinct(CodiBassa, any, ECELS_MITJANA, CUBETA_MITJANA, depredadors_invasors) %>%
  arrange(CodiBassa, any)

ECELS_mat <- covs_yearly %>%
  select(CodiBassa, any, ECELS_MITJANA) %>%
  pivot_wider(names_from = any, values_from = ECELS_MITJANA) %>%
  select(-CodiBassa) %>%
  as.matrix()

CUBETA_mat <- covs_yearly %>%
  select(CodiBassa, any, CUBETA_MITJANA) %>%
  pivot_wider(names_from = any, values_from = CUBETA_MITJANA) %>%
  select(-CodiBassa) %>%
  as.matrix()

imputa_temporal <- function(vec) {
  for (i in seq_along(vec)) {
    if (is.na(vec[i])) {
      dist <- abs(seq_along(vec) - i)
      dist[is.na(vec)] <- Inf
      nearest <- which.min(dist)
      vec[i] <- vec[nearest]
    }
  }
  return(vec)
}

CUBETA_mat <- t(apply(CUBETA_mat, 1, imputa_temporal))
ECELS_mat  <- t(apply(ECELS_mat, 1, imputa_temporal))

depredadors_mat <- covs_yearly %>%
  select(CodiBassa, any, depredadors_invasors) %>%
  pivot_wider(names_from = any, values_from = depredadors_invasors) %>%
  select(-CodiBassa) %>%
  as.matrix()

season_mat <- matrix(
  rep(1:Primary, each = Sites),
  nrow = Sites,
  ncol = Primary,
  byrow = FALSE
)

yearlySiteCovs_list <- list(
  seas = season_mat,
  ECELS_MITJANA = ECELS_mat,
  CUBETA_MITJANA = CUBETA_mat,
  depredadors_invasors = depredadors_mat
)
  
############################################################
# 14. CREACIÓ DE L’UNMARKEDMULTFRAME
############################################################

vv_unmarkedMultFrame <- unmarkedMultFrame(
  y = y_mat,
  numPrimary = Primary,
  obsCovs = obsCovs_list,
  siteCovs = siteCovs_df,
  yearlySiteCovs = yearlySiteCovs_list
)

summary(vv_unmarkedMultFrame)

############################################################
# 15. MODEL BASE I MODEL COMPLET
############################################################

model_base <- colext(
  psiformula = ~ 1,
  gammaformula = ~ 1,
  epsilonformula = ~ 1,
  pformula = ~ 1,
  data = vv_unmarkedMultFrame
)

summary(model_base)

model_complet <- colext(
  psiformula = ~ ZONA + Hidroperiode + Superfmax_m2,
  gammaformula = ~ ECELS_MITJANA + CUBETA_MITJANA + depredadors_invasors,
  epsilonformula = ~ ECELS_MITJANA + CUBETA_MITJANA + depredadors_invasors,
  pformula = ~ terbol + columnaob + salabre,
  data = vv_unmarkedMultFrame
)

summary(model_complet)

############################################################
# 16. BUCLE DE SELECCIÓ AUTOMÀTICA DE MODELS
############################################################

psi_covs     <- c("ZONA", "Hidroperiode", "Superfmax_m2")
gamma_covs   <- c("ECELS_MITJANA", "CUBETA_MITJANA", "depredadors_invasors")
epsilon_covs <- c("ECELS_MITJANA", "CUBETA_MITJANA", "depredadors_invasors")
p_covs       <- c("terbol", "columnaob", "salabre")

all_formulas_from_covs <- function(covs) {
  if (length(covs) == 0) return("~ 1")
  mat <- expand.grid(rep(list(c(FALSE, TRUE)), length(covs)))
  colnames(mat) <- covs
  forms <- apply(mat, 1, function(row) {
    sel <- covs[as.logical(row)]
    if (length(sel) == 0) "~ 1"
    else paste("~", paste(sel, collapse = " + "))
  })
  unique(forms)
}

psi_forms     <- all_formulas_from_covs(psi_covs)
gamma_forms   <- all_formulas_from_covs(gamma_covs)
epsilon_forms <- all_formulas_from_covs(epsilon_covs)
p_forms       <- all_formulas_from_covs(p_covs)

### Funció de filtratge robust
is_good_unmarked_model <- function(fit) {
  if (inherits(fit, "try-error")) return(FALSE)
  if (!inherits(fit, "unmarkedFit")) return(FALSE)
  
  opt <- fit@opt
  if (is.null(opt$conv) || opt$conv != 0) return(FALSE)
  
  coefs <- coef(fit)
  if (any(!is.finite(coefs))) return(FALSE)
  if (any(abs(coefs) > 20)) return(FALSE)
  
  H <- opt$hessian
  if (is.null(H)) return(FALSE)
  if (any(!is.finite(H))) return(FALSE)
  if (abs(det(H)) < 1e-6) return(FALSE)
  
  V <- try(solve(H), silent = TRUE)
  if (inherits(V, "try-error")) return(FALSE)
  if (any(!is.finite(diag(V)))) return(FALSE)
  
  TRUE
}

### Bucle
results <- list()
k <- 1

for (psi_f in psi_forms) {
  for (g_f in gamma_forms) {
    for (e_f in epsilon_forms) {
      for (p_f in p_forms) {
        
        fit <- try(
          colext(
            psiformula = as.formula(psi_f),
            gammaformula = as.formula(g_f),
            epsilonformula = as.formula(e_f),
            pformula = as.formula(p_f),
            data = vv_unmarkedMultFrame
          ),
          silent = TRUE
        )
        if (!is_good_unmarked_model(fit)) next
        
        aic_val <- AIC(fit)
        
        results[[k]] <- list(
          fit = fit,
          spec = list(
            psi = psi_f,
            gamma = g_f,
            epsilon = e_f,
            p = p_f
          ),
          AIC = aic_val
        )
        k <- k + 1
      }
    }
  }
}
if (length(results) == 0) {
  warning("Cap model vàlid trobat.")
} else {
  aics <- sapply(results, function(x) x$AIC)
  ord <- order(aics)
  top_n <- min(10, length(results))
  best10 <- results[ord[1:top_n]]
  
  best10_summary <- lapply(best10, function(x) {
    list(
      AIC = x$AIC,
      psi = x$spec$psi,
      gamma = x$spec$gamma,
      epsilon = x$spec$epsilon,
      p = x$spec$p
    )
  })
  
  print(best10_summary)
}

############################################################
# 17. 10 MILLORS MODELS PER SSAL SEGONS L'ESCALA AIC
############################################################

m1 <- colext(
  psiformula      = ~ 1,
  gammaformula    = ~ ECELS_MITJANA + CUBETA_MITJANA,
  epsilonformula  = ~ ECELS_MITJANA,
  pformula        = ~ columnaob + salabre,
  data = vv_unmarkedMultFrame
)
m1

m2 <- colext(
  psiformula      = ~ Superfmax_m2,
  gammaformula    = ~ ECELS_MITJANA + CUBETA_MITJANA,
  epsilonformula  = ~ ECELS_MITJANA,
  pformula        = ~ columnaob + salabre,
  data = vv_unmarkedMultFrame
)
m2

m3 <- colext(
  psiformula      = ~ 1,
  gammaformula    = ~ ECELS_MITJANA + CUBETA_MITJANA,
  epsilonformula  = ~ depredadors_invasors,
  pformula        = ~ columnaob + salabre,
  data = vv_unmarkedMultFrame
)
m3

m4 <- colext(
  psiformula      = ~ Superfmax_m2,
  gammaformula    = ~ ECELS_MITJANA + CUBETA_MITJANA,
  epsilonformula  = ~ depredadors_invasors,
  pformula        = ~ columnaob + salabre,
  data = vv_unmarkedMultFrame
)
m4

m5 <- colext(
  psiformula      = ~ 1,
  gammaformula    = ~ ECELS_MITJANA + CUBETA_MITJANA,
  epsilonformula  = ~ ECELS_MITJANA + depredadors_invasors,
  pformula        = ~ columnaob + salabre,
  data = vv_unmarkedMultFrame
)
m5

m6 <- colext(
  psiformula      = ~ 1,
  gammaformula    = ~ ECELS_MITJANA + CUBETA_MITJANA,
  epsilonformula  = ~ 1,
  pformula        = ~ columnaob + salabre,
  data = vv_unmarkedMultFrame
)
m6 # millora significativament el model amb 0.09955092 Pr > Chisq

m7 <- colext(
  psiformula      = ~ Superfmax_m2,
  gammaformula    = ~ ECELS_MITJANA + CUBETA_MITJANA,
  epsilonformula  = ~ ECELS_MITJANA + depredadors_invasors,
  pformula        = ~ columnaob + salabre,
  data = vv_unmarkedMultFrame
)
m7

m8 <- colext(
  psiformula      = ~ Superfmax_m2,
  gammaformula    = ~ ECELS_MITJANA + CUBETA_MITJANA,
  epsilonformula  = ~ 1,
  pformula        = ~ columnaob + salabre,
  data = vv_unmarkedMultFrame
)
m8

m9 <- colext(
  psiformula      = ~1,
  gammaformula    = ~ ECELS_MITJANA + CUBETA_MITJANA,
  epsilonformula  = ~ depredadors_invasors,
  pformula        = ~ terbol + columnaob + salabre,
  data = vv_unmarkedMultFrame
)
m9

m10 <- colext(
  psiformula      = ~1,
  gammaformula    = ~ ECELS_MITJANA,
  epsilonformula  = ~ ECELS_MITJANA + depredadors_invasors,
  pformula        = ~ columnaob + salabre,
  data = vv_unmarkedMultFrame
)
m10

############################################################
# 18. 10 MILLORS MODELS PER AALM SEGONS L'ESCALA AIC
############################################################

m1 <- colext(
  psiformula      = ~ ZONA,
  gammaformula    = ~ ECELS_MITJANA + depredadors_invasors,
  epsilonformula  = ~ CUBETA_MITJANA,
  pformula        = ~ terbol + columnaob + salabre,
  data = vv_unmarkedMultFrame
)
m1

m2 <- colext(
  psiformula      = ~ 1,
  gammaformula    = ~ ECELS_MITJANA + depredadors_invasors,
  epsilonformula  = ~ CUBETA_MITJANA,
  pformula        = ~ terbol + columnaob + salabre,
  data = vv_unmarkedMultFrame
)
m2

m3 <- colext(
  psiformula      = ~ ZONA,
  gammaformula    = ~ ECELS_MITJANA + depredadors_invasors,
  epsilonformula  = ~ ECELS_MITJANA + CUBETA_MITJANA,
  pformula        = ~ terbol + columnaob + salabre,
  data = vv_unmarkedMultFrame
)
m3

m4 <- colext(
  psiformula      = ~ 1,
  gammaformula    = ~ ECELS_MITJANA + depredadors_invasors,
  epsilonformula  = ~ ECELS_MITJANA + CUBETA_MITJANA,
  pformula        = ~ terbol + columnaob + salabre,
  data = vv_unmarkedMultFrame
)
m4

m5 <- colext(
  psiformula      = ~ ZONA + Superfmax_m2,
  gammaformula    = ~ ECELS_MITJANA + depredadors_invasors,
  epsilonformula  = ~ CUBETA_MITJANA,
  pformula        = ~ terbol + columnaob + salabre,
  data = vv_unmarkedMultFrame
)
m5

m6 <- colext(
  psiformula      = ~ ZONA + Superfmax_m2,
  gammaformula    = ~ ECELS_MITJANA + depredadors_invasors,
  epsilonformula  = ~ ECELS_MITJANA + CUBETA_MITJANA,
  pformula        = ~ terbol + columnaob + salabre,
  data = vv_unmarkedMultFrame
)
m6

m7 <- colext(
  psiformula      = ~ ZONA,
  gammaformula    = ~ ECELS_MITJANA,
  epsilonformula  = ~ CUBETA_MITJANA,
  pformula        = ~ terbol + columnaob + salabre,
  data = vv_unmarkedMultFrame
)
m7

m8 <- colext(
  psiformula      = ~ 1,
  gammaformula    = ~ ECELS_MITJANA + depredadors_invasors,
  epsilonformula  = ~ CUBETA_MITJANA + depredadors_invasors,
  pformula        = ~ terbol + columnaob + salabre,
  data = vv_unmarkedMultFrame
)
m8

m9 <- colext(
  psiformula      = ~ Superfmax_m2,
  gammaformula    = ~ ECELS_MITJANA + depredadors_invasors,
  epsilonformula  = ~ CUBETA_MITJANA,
  pformula        = ~ terbol + columnaob + salabre,
  data = vv_unmarkedMultFrame
)
m9

m10 <- colext(
  psiformula      = ~ 1,
  gammaformula    = ~ ECELS_MITJANA,
  epsilonformula  = ~ CUBETA_MITJANA,
  pformula        = ~ terbol + columnaob + salabre,
  data = vv_unmarkedMultFrame
)
m10

############################################################
# 19. 10 MILLORS MODELS PER PPUN SEGONS L'ESCALA AIC
############################################################

m1 <- colext(
  psiformula      = ~ ZONA + Hidroperiode,
  gammaformula    = ~ 1,
  epsilonformula  = ~ ECELS_MITJANA + CUBETA_MITJANA,
  pformula        = ~ columnaob + salabre,
  data = vv_unmarkedMultFrame
)
m1

m2 <- colext(
  psiformula      = ~ ZONA,
  gammaformula    = ~ 1,
  epsilonformula  = ~ ECELS_MITJANA + CUBETA_MITJANA,
  pformula        = ~ columnaob + salabre,
  data = vv_unmarkedMultFrame
)
m2

m3 <- colext(
  psiformula      = ~ ZONA + Hidroperiode,
  gammaformula    = ~ ECELS_MITJANA,
  epsilonformula  = ~ ECELS_MITJANA + CUBETA_MITJANA,
  pformula        = ~ columnaob + salabre,
  data = vv_unmarkedMultFrame
)
m3

m4 <- colext(
  psiformula      = ~ ZONA + Hidroperiode,
  gammaformula    = ~ depredadors_invasors,
  epsilonformula  = ~ ECELS_MITJANA + CUBETA_MITJANA,
  pformula        = ~ columnaob + salabre,
  data = vv_unmarkedMultFrame
)
m4

m5 <- colext(
  psiformula      = ~ ZONA + Hidroperiode,
  gammaformula    = ~ depredadors_invasors,
  epsilonformula  = ~ ECELS_MITJANA + CUBETA_MITJANA + depredadors_invasors,
  pformula        = ~ columnaob + salabre,
  data = vv_unmarkedMultFrame
)
m5

m6 <- colext(
  psiformula      = ~ ZONA + Hidroperiode,
  gammaformula    = ~ CUBETA_MITJANA,
  epsilonformula  = ~ ECELS_MITJANA + CUBETA_MITJANA,
  pformula        = ~ columnaob + salabre,
  data = vv_unmarkedMultFrame
)
m6

m7 <- colext(
  psiformula      = ~ ZONA + Superfmax_m2,
  gammaformula    = ~ 1,
  epsilonformula  = ~ ECELS_MITJANA + CUBETA_MITJANA,
  pformula        = ~ columnaob + salabre,
  data = vv_unmarkedMultFrame
)
m7

m8 <- colext(
  psiformula      = ~ ZONA + Hidroperiode,
  gammaformula    = ~ ECELS_MITJANA + depredadors_invasors,
  epsilonformula  = ~ CUBETA_MITJANA,
  pformula        = ~ columnaob + salabre,
  data = vv_unmarkedMultFrame
)
m8

m9 <- colext(
  psiformula      = ~ ZONA + Hidroperiode + Superfmax_m2,
  gammaformula    = ~ 1,
  epsilonformula  = ~ ECELS_MITJANA + CUBETA_MITJANA,
  pformula        = ~ columnaob + salabre,
  data = vv_unmarkedMultFrame
)
m9

m10 <- colext(
  psiformula      = ~ ZONA,
  gammaformula    = ~ ECELS_MITJANA,
  epsilonformula  = ~ ECELS_MITJANA + CUBETA_MITJANA,
  pformula        = ~ columnaob + salabre,
  data = vv_unmarkedMultFrame
)
m10

############################################################
# 20. 10 MILLORS MODELS PER BSPI SEGONS L'ESCALA AIC
############################################################

m1 <- colext(
  psiformula      = ~ Superfmax_m2,
  gammaformula    = ~ CUBETA_MITJANA,
  epsilonformula  = ~ CUBETA_MITJANA,
  pformula        = ~ terbol,
  data = vv_unmarkedMultFrame
)
m1

m2 <- colext(
  psiformula      = ~ Superfmax_m2,
  gammaformula    = ~ CUBETA_MITJANA,
  epsilonformula  = ~ CUBETA_MITJANA,
  pformula        = ~ terbol + salabre,
  data = vv_unmarkedMultFrame
)
m2

m3 <- colext(
  psiformula      = ~ Superfmax_m2,
  gammaformula    = ~ CUBETA_MITJANA + depredadors_invasors,
  epsilonformula  = ~ CUBETA_MITJANA,
  pformula        = ~ terbol,
  data = vv_unmarkedMultFrame
)
m3

m4 <- colext(
  psiformula      = ~ Superfmax_m2,
  gammaformula    = ~ CUBETA_MITJANA,
  epsilonformula  = ~ ECELS_MITJANA + CUBETA_MITJANA,
  pformula        = ~ terbol,
  data = vv_unmarkedMultFrame
)
m4

m5 <- colext(
  psiformula      = ~ Superfmax_m2,
  gammaformula    = ~ CUBETA_MITJANA,
  epsilonformula  = ~ CUBETA_MITJANA + depredadors_invasors,
  pformula        = ~ terbol,
  data = vv_unmarkedMultFrame
)
m5

m6 <- colext(
  psiformula      = ~ Superfmax_m2,
  gammaformula    = ~ ECELS_MITJANA + CUBETA_MITJANA,
  epsilonformula  = ~ CUBETA_MITJANA,
  pformula        = ~ terbol,
  data = vv_unmarkedMultFrame
)
m6

m7 <- colext(
  psiformula      = ~ Hidroperiode + Superfmax_m2,
  gammaformula    = ~ CUBETA_MITJANA,
  epsilonformula  = ~ CUBETA_MITJANA,
  pformula        = ~ terbol,
  data = vv_unmarkedMultFrame
)
m7

m8 <- colext(
  psiformula      = ~ Superfmax_m2,
  gammaformula    = ~ CUBETA_MITJANA,
  epsilonformula  = ~ CUBETA_MITJANA,
  pformula        = ~ terbol + columnaob + salabre,
  data = vv_unmarkedMultFrame
)
m8

m9 <- colext(
  psiformula      = ~ Superfmax_m2,
  gammaformula    = ~ CUBETA_MITJANA + depredadors_invasors,
  epsilonformula  = ~ CUBETA_MITJANA,
  pformula        = ~ terbol + salabre,
  data = vv_unmarkedMultFrame
)
m9

m10 <- colext(
  psiformula      = ~ Superfmax_m2,
  gammaformula    = ~ CUBETA_MITJANA,
  epsilonformula  = ~ ECELS_MITJANA + CUBETA_MITJANA,
  pformula        = ~ terbol + salabre,
  data = vv_unmarkedMultFrame
)
m10

############################################################
# 21. 10 MILLORS MODELS PER PPER SEGONS L'ESCALA AIC
############################################################

m1 <- colext(
  psiformula      = ~ 1,
  gammaformula    = ~ CUBETA_MITJANA,
  epsilonformula  = ~ 1,
  pformula        = ~ columnaob + salabre,
  data = vv_unmarkedMultFrame
)
m1

m2 <- colext(
  psiformula      = ~ 1,
  gammaformula    = ~ CUBETA_MITJANA,
  epsilonformula  = ~ 1,
  pformula        = ~ terbol + columnaob + salabre,
  data = vv_unmarkedMultFrame
)
m2

m3 <- colext(
  psiformula      = ~ Superfmax_m2,
  gammaformula    = ~ CUBETA_MITJANA,
  epsilonformula  = ~ ECELS_MITJANA,
  pformula        = ~ terbol + salabre,
  data = vv_unmarkedMultFrame
)
m3

m4 <- colext(
  psiformula      = ~ 1,
  gammaformula    = ~ CUBETA_MITJANA,
  epsilonformula  = ~ ECELS_MITJANA,
  pformula        = ~ terbol + columnaob + salabre,
  data = vv_unmarkedMultFrame
)
m4

m5 <- colext(
  psiformula      = ~ Superfmax_m2,
  gammaformula    = ~ CUBETA_MITJANA,
  epsilonformula  = ~ CUBETA_MITJANA,
  pformula        = ~ salabre,
  data = vv_unmarkedMultFrame
)
m5

m6 <- colext(
  psiformula      = ~ Superfmax_m2,
  gammaformula    = ~ CUBETA_MITJANA,
  epsilonformula  = ~ CUBETA_MITJANA,
  pformula        = ~ terbol + salabre,
  data = vv_unmarkedMultFrame
)
m6

m7 <- colext(
  psiformula      = ~ 1,
  gammaformula    = ~ 1,
  epsilonformula  = ~ CUBETA_MITJANA,
  pformula        = ~ salabre,
  data = vv_unmarkedMultFrame
)
m7

m8 <- colext(
  psiformula      = ~ 1,
  gammaformula    = ~ 1,
  epsilonformula  = ~ 1,
  pformula        = ~ salabre,
  data = vv_unmarkedMultFrame
)
m8

m9 <- colext(
  psiformula      = ~ 1,
  gammaformula    = ~ 1,
  epsilonformula  = ~ CUBETA_MITJANA,
  pformula        = ~ terbol + salabre,
  data = vv_unmarkedMultFrame
)
m9

m10 <- colext(
  psiformula      = ~ 1,
  gammaformula    = ~ CUBETA_MITJANA,
  epsilonformula  = ~ CUBETA_MITJANA,
  pformula        = ~ columnaob + salabre,
  data = vv_unmarkedMultFrame
)
m10

############################################################
# 22. COMPARACIÓ DE MODELS MITJANÇANT LRT
############################################################

fl_valid <- fitList(
  model_base = model_base,
  model_complet = model_complet,
  m1 = m1,
  m2 = m2,
  m3 = m3,
  m4 = m4,
  m5 = m5,
  m6 = m6,
  m7 = m7,
  m8 = m8,
  m9 = m9,
  m10 = m10
)
modSel(fl_valid)

#Likelihood Ratio Test
LRT(m1,m2)

############################################################
# 23. MILLORS MODELS PER ESPÈCIE (SSAL, AALM, PPUN, BSPI, PPER)
############################################################

bestSSAL <- colext(
  psiformula = ~ 1,
  gammaformula = ~ ECELS_MITJANA + CUBETA_MITJANA,
  epsilonformula = ~ ECELS_MITJANA,
  pformula = ~ columnaob + salabre,
  data = vv_unmarkedMultFrame
)

bestAALM <- colext(
  psiformula = ~ ZONA,
  gammaformula = ~ ECELS_MITJANA + depredadors_invasors,
  epsilonformula = ~ CUBETA_MITJANA,
  pformula = ~ terbol + columnaob + salabre,
  data = vv_unmarkedMultFrame
)

bestPPUN <- colext(
  psiformula = ~ ZONA + Hidroperiode,
  gammaformula = ~ 1,
  epsilonformula = ~ ECELS_MITJANA + CUBETA_MITJANA,
  pformula = ~ columnaob + salabre,
  data = vv_unmarkedMultFrame
)

bestBSPI <- colext(
  psiformula = ~ Superfmax_m2,
  gammaformula = ~ CUBETA_MITJANA,
  epsilonformula = ~ CUBETA_MITJANA,
  pformula = ~ terbol,
  data = vv_unmarkedMultFrame
)

bestPPER <- colext(
  psiformula = ~ 1,
  gammaformula = ~ CUBETA_MITJANA,
  epsilonformula = ~ 1,
  pformula = ~ columnaob + salabre,
  data = vv_unmarkedMultFrame
)

############################################################
# 18. GRÀFICS: SALAMANDRA SALAMANDRA 
############################################################

mod <- bestSSAL

# Seqüències
ecels_seq <- seq(
  min(vv_unmarkedMultFrame@yearlySiteCovs$ECELS_MITJANA, na.rm = TRUE),
  max(vv_unmarkedMultFrame@yearlySiteCovs$ECELS_MITJANA, na.rm = TRUE),
  length.out = 100
)

cubeta_seq <- seq(
  min(vv_unmarkedMultFrame@yearlySiteCovs$CUBETA_MITJANA, na.rm = TRUE),
  max(vv_unmarkedMultFrame@yearlySiteCovs$CUBETA_MITJANA, na.rm = TRUE),
  length.out = 100
)

cubeta_vals <- quantile(cubeta_seq, probs = c(0.1, 0.5, 0.9))

# Gamma
gamma_df <- bind_rows(lapply(cubeta_vals, function(cu) {
  newdat <- data.frame(ECELS_MITJANA = ecels_seq, CUBETA_MITJANA = cu)
  pred <- predict(mod, type = "col", newdata = newdat)
  data.frame(
    ECELS_MITJANA = ecels_seq,
    estimate = pred$Predicted,
    lower = pred$lower,
    upper = pred$upper,
    CUBETA = round(cu, 2)
  )
}))

GRAFICgammaSSAL <- ggplot(gamma_df,
                          aes(x = ECELS_MITJANA, y = estimate,
                              color = factor(CUBETA), fill = factor(CUBETA))) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25) +
  geom_line(linewidth = 1.2) +
  labs(x = "ECELS", y = "Colonització") +
  theme_classic()

# Extinció
ext_df <- {
  newdat <- data.frame(ECELS_MITJANA = ecels_seq)
  pred <- predict(mod, type = "ext", newdata = newdat)
  data.frame(
    ECELS_MITJANA = ecels_seq,
    estimate = pred$Predicted,
    lower = pred$lower,
    upper = pred$upper
  )
}

GRAFICextSSAL <- ggplot(ext_df, aes(ECELS_MITJANA, estimate)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25) +
  geom_line(linewidth = 1.2, color = "darkred") +
  labs(x = "ECELS", y = "Extinció") +
  theme_classic()

# Ocupació dinàmica
re <- ranef(mod)
state_mat <- bup(re, stat = "mean")

psi_df <- data.frame(
  Year = seq_len(ncol(state_mat)),
  Occupancy = colMeans(state_mat)
)

GRAFOcupacióDinSSAL <- ggplot(psi_df, aes(Year, Occupancy)) +
  geom_line(linewidth = 1.2, color = "darkblue") +
  geom_point(size = 2) +
  labs(x = "Any", y = "Ocupació") +
  theme_classic()

# Plot final
final_plotSSAL <- (GRAFICgammaSSAL / (GRAFICextSSAL | GRAFOcupacióDinSSAL)) +
  plot_annotation(tag_levels = "A")

print(final_plotSSAL)

############################################################
# 19. GRÀFICS: ALYTES ALMOGAVARII
############################################################

mod <- bestAALM

# Seqüències
ecels_seq <- seq(
  min(vv_unmarkedMultFrame@yearlySiteCovs$ECELS_MITJANA, na.rm = TRUE),
  max(vv_unmarkedMultFrame@yearlySiteCovs$ECELS_MITJANA, na.rm = TRUE),
  length.out = 100
)

cubeta_seq <- seq(
  min(vv_unmarkedMultFrame@yearlySiteCovs$CUBETA_MITJANA, na.rm = TRUE),
  max(vv_unmarkedMultFrame@yearlySiteCovs$CUBETA_MITJANA, na.rm = TRUE),
  length.out = 100
)

cubeta_vals <- quantile(cubeta_seq, probs = c(0.1, 0.5, 0.9))

# Gamma
gamma_df <- bind_rows(lapply(cubeta_vals, function(cu) {
  newdat <- data.frame(ECELS_MITJANA = ecels_seq, CUBETA_MITJANA = cu)
  pred <- predict(mod, type = "col", newdata = newdat)
  data.frame(
    ECELS_MITJANA = ecels_seq,
    estimate = pred$Predicted,
    lower = pred$lower,
    upper = pred$upper,
    CUBETA = round(cu, 2)
  )
}))

GRAFICgammaAALM <- ggplot(gamma_df,
                          aes(x = ECELS_MITJANA, y = estimate,
                              color = factor(CUBETA), fill = factor(CUBETA))) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25) +
  geom_line(linewidth = 1.2) +
  labs(x = "ECELS", y = "Colonització") +
  theme_classic()

# Extinció
ext_df <- {
  newdat <- data.frame(ECELS_MITJANA = ecels_seq)
  pred <- predict(mod, type = "ext", newdata = newdat)
  data.frame(
    ECELS_MITJANA = ecels_seq,
    estimate = pred$Predicted,
    lower = pred$lower,
    upper = pred$upper
  )
}

GRAFICextAALM <- ggplot(ext_df, aes(ECELS_MITJANA, estimate)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25) +
  geom_line(linewidth = 1.2, color = "darkred") +
  labs(x = "ECELS", y = "Extinció") +
  theme_classic()

# Ocupació dinàmica
re <- ranef(mod)
state_mat <- bup(re, stat = "mean")

psi_df <- data.frame(
  Year = seq_len(ncol(state_mat)),
  Occupancy = colMeans(state_mat)
)

GRAFOcupacióDinAALM <- ggplot(psi_df, aes(Year, Occupancy)) +
  geom_line(linewidth = 1.2, color = "darkblue") +
  geom_point(size = 2) +
  labs(x = "Any", y = "Ocupació") +
  theme_classic()

# Plot final
final_plotAALM <- (GRAFICgammaAALM / (GRAFICextAALM | GRAFOcupacióDinAALM)) +
  plot_annotation(tag_levels = "A")

print(final_plotAALM)

############################################################
# 20. GRÀFICS: PELODYTES PUNCTATUS
############################################################

mod <- bestPPUN

# Seqüències
ecels_seq <- seq(
  min(vv_unmarkedMultFrame@yearlySiteCovs$ECELS_MITJANA, na.rm = TRUE),
  max(vv_unmarkedMultFrame@yearlySiteCovs$ECELS_MITJANA, na.rm = TRUE),
  length.out = 100
)

cubeta_seq <- seq(
  min(vv_unmarkedMultFrame@yearlySiteCovs$CUBETA_MITJANA, na.rm = TRUE),
  max(vv_unmarkedMultFrame@yearlySiteCovs$CUBETA_MITJANA, na.rm = TRUE),
  length.out = 100
)

cubeta_vals <- quantile(cubeta_seq, probs = c(0.1, 0.5, 0.9))

# Gamma (PPUN no té gamma → només extinció)
# Extinció
ext_df <- bind_rows(lapply(cubeta_vals, function(cu) {
  newdat <- data.frame(ECELS_MITJANA = ecels_seq, CUBETA_MITJANA = cu)
  pred <- predict(mod, type = "ext", newdata = newdat)
  data.frame(
    ECELS_MITJANA = ecels_seq,
    estimate = pred$Predicted,
    lower = pred$lower,
    upper = pred$upper,
    CUBETA = round(cu, 2)
  )
}))

GRAFICextPPUN <- ggplot(ext_df,
                        aes(x = ECELS_MITJANA, y = estimate,
                            color = factor(CUBETA), fill = factor(CUBETA))) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25) +
  geom_line(linewidth = 1.2) +
  labs(x = "ECELS", y = "Extinció") +
  theme_classic()

# Ocupació inicial
site_covs <- vv_unmarkedMultFrame@siteCovs
zones <- levels(site_covs$ZONA)
hidro_levels <- levels(site_covs$Hidroperiode)

psi_df <- bind_rows(lapply(zones, function(z) {
  bind_rows(lapply(hidro_levels, function(h) {
    newdat <- data.frame(
      ZONA = factor(z, levels = zones),
      Hidroperiode = factor(h, levels = hidro_levels)
    )
    pred <- predict(mod, type = "psi", newdata = newdat)
    data.frame(
      ZONA = z,
      Hidroperiode = h,
      estimate = pred$Predicted,
      lower = pred$lower,
      upper = pred$upper
    )
  }))
}))

GRAFICpsiPPUN <- ggplot(psi_df,
                        aes(x = ZONA, y = estimate, color = Hidroperiode)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.2,
                position = position_dodge(width = 0.5)) +
  labs(x = "Zona", y = "Ocupació inicial") +
  theme_classic()

# Ocupació dinàmica
re <- ranef(mod)
state_mat <- bup(re, stat = "mean")

psi_dyn <- data.frame(
  Year = seq_len(ncol(state_mat)),
  Occupancy = colMeans(state_mat)
)

GRAFOcupacióDinPPUN <- ggplot(psi_dyn, aes(Year, Occupancy)) +
  geom_line(linewidth = 1.2, color = "darkblue") +
  geom_point(size = 2) +
  labs(x = "Any", y = "Ocupació") +
  theme_classic()

# Plot final
final_plotPPUN <- (GRAFICextPPUN / (GRAFICpsiPPUN | GRAFOcupacióDinPPUN)) +
  plot_annotation(tag_levels = "A")

print(final_plotPPUN)

############################################################
# 21. GRÀFICS: BUFO SPINOSUS
############################################################

mod <- bestBSPI

# Seqüències
cubeta_seq <- seq(
  min(vv_unmarkedMultFrame@yearlySiteCovs$CUBETA_MITJANA, na.rm = TRUE),
  max(vv_unmarkedMultFrame@yearlySiteCovs$CUBETA_MITJANA, na.rm = TRUE),
  length.out = 100
)

# Gamma
gamma_df <- {
  newdat <- data.frame(CUBETA_MITJANA = cubeta_seq)
  pred <- predict(mod, type = "col", newdata = newdat)
  data.frame(
    CUBETA_MITJANA = cubeta_seq,
    estimate = pred$Predicted,
    lower = pred$lower,
    upper = pred$upper
  )
}

GRAFICgammaBSPI <- ggplot(gamma_df, aes(CUBETA_MITJANA, estimate)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25, fill = "lightblue") +
  geom_line(linewidth = 1.2, color = "blue4") +
  labs(x = "% Cubeta plena", y = "Colonització") +
  theme_classic()

# Extinció
ext_df <- {
  newdat <- data.frame(CUBETA_MITJANA = cubeta_seq)
  pred <- predict(mod, type = "ext", newdata = newdat)
  data.frame(
    CUBETA_MITJANA = cubeta_seq,
    estimate = pred$Predicted,
    lower = pred$lower,
    upper = pred$upper
  )
}

GRAFICextBSPI <- ggplot(ext_df, aes(CUBETA_MITJANA, estimate)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25, fill = "grey70") +
  geom_line(linewidth = 1.2, color = "darkred") +
  labs(x = "% Cubeta plena", y = "Extinció") +
  theme_classic()

# Ocupació inicial
superf_seq <- seq(
  min(vv_unmarkedMultFrame@siteCovs$Superfmax_m2, na.rm = TRUE),
  max(vv_unmarkedMultFrame@siteCovs$Superfmax_m2, na.rm = TRUE),
  length.out = 100
)

psi_df <- {
  newdat <- data.frame(Superfmax_m2 = superf_seq)
  pred <- predict(mod, type = "psi", newdata = newdat)
  data.frame(
    Superfmax_m2 = superf_seq,
    estimate = pred$Predicted,
    lower = pred$lower,
    upper = pred$upper
  )
}

GRAFICpsiBSPI <- ggplot(psi_df, aes(Superfmax_m2, estimate)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25, fill = "lightgreen") +
  geom_line(linewidth = 1.2, color = "darkgreen") +
  labs(x = "Superfície (m²)", y = "Ocupació inicial") +
  theme_classic()

# Ocupació dinàmica
re <- ranef(mod)
state_mat <- bup(re, stat = "mean")

psi_dyn <- data.frame(
  Year = seq_len(ncol(state_mat)),
  Occupancy = colMeans(state_mat)
)

GRAFOcupacióDinBSPI <- ggplot(psi_dyn, aes(Year, Occupancy)) +
  geom_line(linewidth = 1.2, color = "darkblue") +
  geom_point(size = 2) +
  labs(x = "Any", y = "Ocupació") +
  theme_classic()

# Plot final
final_plotBSPI <- (GRAFICgammaBSPI / (GRAFICextBSPI | GRAFICpsiBSPI | GRAFOcupacióDinBSPI)) +
  plot_annotation(tag_levels = "A")

print(final_plotBSPI)

############################################################
# 22. GRÀFICS: PELOPHYLAX PEREZI
############################################################

mod <- bestPPER

# Seqüències
cubeta_seq <- seq(
  min(vv_unmarkedMultFrame@yearlySiteCovs$CUBETA_MITJANA, na.rm = TRUE),
  max(vv_unmarkedMultFrame@yearlySiteCovs$CUBETA_MITJANA, na.rm = TRUE),
  length.out = 100
)

# Gamma
gamma_df <- {
  newdat <- data.frame(CUBETA_MITJANA = cubeta_seq)
  pred <- predict(mod, type = "col", newdata = newdat)
  data.frame(
    CUBETA_MITJANA = cubeta_seq,
    estimate = pred$Predicted,
    lower = pred$lower,
    upper = pred$upper
  )
}

GRAFICgammaPPER <- ggplot(gamma_df, aes(CUBETA_MITJANA, estimate)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25, fill = "lightblue") +
  geom_line(linewidth = 1.2, color = "blue4") +
  labs(x = "% Cubeta plena", y = "Colonització") +
  theme_classic()

# Extinció
pred_ext <- predict(mod, type = "ext")

ext_df <- data.frame(
  parameter = "Extinció",
  estimate = pred_ext$Predicted,
  lower = pred_ext$lower,
  upper = pred_ext$upper
)

GRAFICextPPER <- ggplot(ext_df, aes(parameter, estimate)) +
  geom_point(size = 3, color = "darkred") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  labs(x = "", y = "Extinció") +
  theme_classic()

# Ocupació inicial
pred_psi <- predict(mod, type = "psi")

psi_df <- data.frame(
  parameter = "Ocupació inicial",
  estimate = pred_psi$Predicted,
  lower = pred_psi$lower,
  upper = pred_psi$upper
)

GRAFICpsiPPER <- ggplot(psi_df, aes(parameter, estimate)) +
  geom_point(size = 3, color = "darkgreen") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  labs(x = "", y = "Ocupació inicial") +
  theme_classic()

# Ocupació dinàmica
re <- ranef(mod)
state_mat <- bup(re, stat = "mean")

psi_dyn <- data.frame(
  Year = seq_len(ncol(state_mat)),
  Occupancy = colMeans(state_mat)
)

GRAFOcupacióDinPPER <- ggplot(psi_dyn, aes(Year, Occupancy)) +
  geom_line(linewidth = 1.2, color = "darkblue") +
  geom_point(size = 2) +
  labs(x = "Any", y = "Ocupació") +
  theme_classic()

# Plot final
final_plotPPER <- (GRAFICgammaPPER | GRAFOcupacióDinPPER) +
  plot_annotation(tag_levels = "A")

print(final_plotPPER)

############################################################
# 23. GRÀFICS: MATRIU DE CORRELACIONS
############################################################

# Correlació ECELS vs CUBETA
print(cor(dat$ECELS_MITJANA, dat$CUBETA_MITJANA, use = "complete.obs"))
print(cor.test(dat$ECELS_MITJANA, dat$CUBETA_MITJANA))

# Gràfic ECELS vs CUBETA
print(
  ggplot(dat, aes(ECELS_MITJANA, CUBETA_MITJANA)) +
    geom_point() +
    geom_smooth(method = "lm") +
    labs(title = "Relació ECELS vs % Cubeta plena")
)

# Matriu de correlacions
vars <- c(
  "ZONA", "Hidroperiode", "Superfmax_m2",
  "ECELS_MITJANA", "CUBETA_MITJANA", "depredadors_invasors",
  "Terbolesa", "percen_Columnaobservable", "TipusMostreig"
)

dat_sense_NA <- na.omit(dat)
sub <- dat_sense_NA[, vars]

sub_dummy <- fastDummies::dummy_cols(
  sub,
  select_columns = c("ZONA", "Hidroperiode", "depredadors_invasors"),
  remove_first_dummy = TRUE
)

nums <- sub_dummy[sapply(sub_dummy, is.numeric)]
nums <- nums[, sapply(nums, function(x) sd(x, na.rm = TRUE) != 0)]

corr <- cor(nums, use = "complete.obs")
corr_melt <- reshape2::melt(corr)

print(
  ggplot(corr_melt, aes(Var1, Var2, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
      low = "blue", high = "red", mid = "white",
      midpoint = 0, limit = c(-1, 1), name = "Correlació"
    ) +
    geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_fixed() +
    labs(title = "Matriu de correlacions entre covariables")
)
