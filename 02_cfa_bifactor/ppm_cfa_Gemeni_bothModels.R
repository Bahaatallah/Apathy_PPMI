#!/usr/bin/env Rscript
# ==========================================================================
#  CFA PIPELINE: Correlated domains, second-order g, unidim g, & Bifactor
# ==========================================================================
#  Models fitted: correlated factors (cfa_att_exec, cfa_memory, cfa_visuo), second-order g
#  (so_g and its domains), unidimensional g (cfa_g) and the bifactor model (bifactor_g and
#  its specific factors). The manuscript's global cognition score is bifactor_g, the general
#  factor of the bifactor model; the other scores are exported for reference.
# ==========================================================================

# --- 0. Packages -----------------------------------------------------------
required_pkgs <- c("lavaan", "semTools", "semPlot", "dplyr", "readr",
                   "tidyr", "readxl", "ggplot2", "patchwork", "corrplot")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

cat("=== Comprehensive CFA for PPMI Cognitive Data ===\n\n")

# ==========================================================================
#  1. FILE PATHS
# ==========================================================================
basePath     <- file.path(Sys.getenv("PPMI_CLINICAL_ROOT", "."), "ppmi_data")
curatedFile  <- file.path(basePath, "PPMI_Curated_Data_Cut_Public_20260223.xlsx")
sheetName    <- "20260105"

downloadPath <- file.path(basePath, "download")
gdsFile      <- file.path(downloadPath, "Geriatric_Depression_Scale__Short_Version__16Mar2026.csv")
mocaFile     <- file.path(downloadPath, "Montreal_Cognitive_Assessment__MoCA__16Mar2026.csv")

outputDir    <- file.path(basePath, "MATLAB_Results")
if (!dir.exists(outputDir)) dir.create(outputDir, recursive = TRUE)

# ==========================================================================
#  2. LOAD DATA
# ==========================================================================
cat("[1] Loading curated dataset + raw instrument CSVs...\n")

df_main <- read_excel(curatedFile, sheet = sheetName)
df_main$PATNO <- as.character(df_main$PATNO)

if (file.exists(gdsFile)) {
  gdsRaw <- read_csv(gdsFile, show_col_types = FALSE) %>% mutate(PATNO = as.character(PATNO))
} else { gdsRaw <- NULL }

if (file.exists(mocaFile)) {
  mocaRaw <- read_csv(mocaFile, show_col_types = FALSE) %>% mutate(PATNO = as.character(PATNO))
} else { mocaRaw <- NULL }

# ==========================================================================
#  3. COMPUTE COMPOSITE SCORES (With Purified Anatomical Mapping)
# ==========================================================================
cat("[2] Computing composite scores from raw items...\n")

safe_numeric <- function(x) suppressWarnings(as.numeric(x))
nansum_if_complete <- function(mat, min_valid) {
  n_valid <- rowSums(!is.na(mat))
  s <- rowSums(mat, na.rm = TRUE)
  s[n_valid < min_valid] <- NA
  s
}

# --- 3A. GDS ---
if (!is.null(gdsRaw)) {
  gdsItemCols <- c("GDSSATIS", "GDSDROPD", "GDSEMPTY", "GDSBORED", "GDSGSPIR",
                   "GDSAFRAD", "GDSHAPPY", "GDSHLPLS", "GDSHOME",  "GDSMEMRY",
                   "GDSALIVE", "GDSWRTLS", "GDSENRGY", "GDSHOPLS", "GDSBETER")
  reverseIdx  <- c(1, 5, 7, 11, 13)
  
  gdsScored <- matrix(NA_real_, nrow = nrow(gdsRaw), ncol = 15)
  for (it in seq_along(gdsItemCols)) {
    col_name <- gdsItemCols[it]
    if (col_name %in% names(gdsRaw)) {
      rawV <- safe_numeric(gdsRaw[[col_name]])
      gdsScored[, it] <- ifelse(it %in% reverseIdx, 1 - rawV, rawV)
    }
  }
  
  gdsRaw$gds_apathy_3item <- nansum_if_complete(gdsScored[, c(2, 9, 13)], 2)
  gdsRaw$depression_factor <- nansum_if_complete(gdsScored[, c(1, 3, 4, 5, 6, 7, 8, 11, 14, 15)], 7)
  
  gds_merge <- gdsRaw %>% select(PATNO, EVENT_ID, gds_apathy_3item, depression_factor) %>% distinct(PATNO, EVENT_ID, .keep_all = TRUE)
  df_main <- df_main %>% left_join(gds_merge, by = c("PATNO", "EVENT_ID"), suffix = c("", ".gds"))
  
  if ("depression_factor.gds" %in% names(df_main)) { df_main$depression_factor <- coalesce(df_main$depression_factor.gds, df_main$depression_factor); df_main$depression_factor.gds <- NULL }
  if ("gds_apathy_3item.gds" %in% names(df_main)) { df_main$gds_apathy_3item <- coalesce(df_main$gds_apathy_3item.gds, df_main$gds_apathy_3item); df_main$gds_apathy_3item.gds <- NULL }
}

# --- 3B. MoCA (Purified) ---
if (!is.null(mocaRaw)) {
  mc <- mocaRaw %>% mutate(
    PATNO = as.character(PATNO),
    moca_visuo_draw = safe_numeric(MCACUBE) + safe_numeric(MCACLCKC) + safe_numeric(MCACLCKN) + safe_numeric(MCACLCKH),
    moca_trail      = safe_numeric(MCAALTTM),
    moca_fluency    = safe_numeric(MCAVF),
    moca_abstr      = safe_numeric(MCAABSTR),
    moca_naming     = safe_numeric(MCALION) + safe_numeric(MCARHINO) + safe_numeric(MCACAMEL),
    moca_repetition = safe_numeric(MCASNTNC),
    moca_attn       = safe_numeric(MCAFDS) + safe_numeric(MCABDS) + safe_numeric(MCAVIGIL) + safe_numeric(MCASER7),
    moca_memory     = safe_numeric(MCAREC1) + safe_numeric(MCAREC2) + safe_numeric(MCAREC3) + safe_numeric(MCAREC4) + safe_numeric(MCAREC5),
    moca_orientation= safe_numeric(MCADATE) + safe_numeric(MCAMONTH) + safe_numeric(MCAYR) + safe_numeric(MCADAY) + safe_numeric(MCAPLACE) + safe_numeric(MCACITY)
  )
  
  existing_moca_cols <- grep("^moca_", names(df_main), value = TRUE)
  df_main <- df_main %>% select(-any_of(existing_moca_cols))
  moca_clean <- mc %>% select(PATNO, EVENT_ID, starts_with("moca_")) %>% distinct(PATNO, EVENT_ID, .keep_all = TRUE)
  df_main <- df_main %>% left_join(moca_clean, by = c("PATNO", "EVENT_ID"))
}

# ==========================================================================
#  4. FILTER TO PD COHORT & PREPARE CFA INDICATORS
# ==========================================================================
cat("[3] Filtering to PD cohort and standardising indicators...\n")

T_pd <- df_main %>% 
  filter(COHORT == "1") %>%
  mutate(
    np_lns            = safe_numeric(lns),
    np_sdmt           = safe_numeric(SDMTOTAL),
    np_semflu         = safe_numeric(VLTANIM),
    np_tmta           = safe_numeric(TMT_A),
    np_tmtb           = safe_numeric(TMT_B),
    np_hvlt_recall    = safe_numeric(DVT_TOTAL_RECALL),
    np_hvlt_retention = safe_numeric(hvlt_retention),
    np_hvlt_recog     = safe_numeric(DVT_RECOG_DISC_INDEX),
    np_jlo            = safe_numeric(bjlot),
    YEAR_num          = safe_numeric(YEAR)
  ) %>%
  mutate(np_tmta_inv = -1 * np_tmta, np_tmtb_inv = -1 * np_tmtb)

cfa_vars <- c("moca_visuo_draw", "moca_trail", "moca_fluency", "moca_abstr", 
              "moca_naming", "moca_repetition", "moca_attn", "moca_memory", "moca_orientation",
              "np_lns", "np_sdmt", "np_semflu", "np_tmta_inv", "np_tmtb_inv",
              "np_hvlt_recall", "np_hvlt_retention", "np_hvlt_recog", "np_jlo")

T_pd_bl <- T_pd %>% filter(YEAR_num == 0)

for (v in cfa_vars) {
  if (v %in% names(T_pd)) {
    mu  <- mean(T_pd_bl[[v]], na.rm = TRUE)
    sig <- sd(T_pd_bl[[v]], na.rm = TRUE)
    T_pd[[v]] <- (T_pd[[v]] - mu) / sig
  }
}
T_pd_bl <- T_pd %>% filter(YEAR_num == 0)

# ==========================================================================
#  5. CFA MODEL SPECIFICATIONS
# ==========================================================================
cat("[4] Specifying CFA models (correlated + second-order are primary)...\n\n")

# --- PRIMARY Model 1: Correlated Domains (Clinical Subcomponents) ---
correlated_model <- '
  att_exec =~ np_lns + np_sdmt + np_tmta_inv + np_tmtb_inv + np_semflu + moca_attn + moca_trail + moca_fluency + moca_abstr
  memory   =~ np_hvlt_recall + np_hvlt_retention + np_hvlt_recog + moca_memory
  visuo    =~ np_jlo + moca_visuo_draw
'

# --- PRIMARY Model 2: Second-order general factor on top of correlated domains ---
# Higher-order g extracts what is shared among the three domains; gives a
# well-identified global cognition score alongside the three domain factors.
# Note: with exactly three first-order factors, this model is statistically
# equivalent in fit to the correlated-factors model (df = 0 at 2nd order).
second_order_model <- '
  att_exec =~ np_lns + np_sdmt + np_tmta_inv + np_tmtb_inv + np_semflu + moca_attn + moca_trail + moca_fluency + moca_abstr
  memory   =~ np_hvlt_recall + np_hvlt_retention + np_hvlt_recog + moca_memory
  visuo    =~ np_jlo + moca_visuo_draw

  g =~ att_exec + memory + visuo
'

# --- REFERENCE Model A: Unidimensional (single g across all indicators, incl. floating MoCA items) ---
unidim_model <- '
  g =~ np_lns + np_sdmt + np_tmta_inv + np_tmtb_inv + np_semflu + moca_attn + moca_trail + moca_fluency + moca_abstr +
       np_hvlt_recall + np_hvlt_retention + np_hvlt_recog + moca_memory +
       np_jlo + moca_visuo_draw +
       moca_naming + moca_repetition + moca_orientation
'

# --- REFERENCE Model B: ORIGINAL Bifactor (naive, preserved for comparison) ---
bifactor_model <- '
  g =~ np_lns + np_sdmt + np_tmta_inv + np_tmtb_inv + np_semflu + moca_attn + moca_trail + moca_fluency + moca_abstr +
       np_hvlt_recall + np_hvlt_retention + np_hvlt_recog + moca_memory +
       np_jlo + moca_visuo_draw +
       moca_naming + moca_repetition + moca_orientation

  att_exec =~ np_lns + np_sdmt + np_tmta_inv + np_tmtb_inv + np_semflu + moca_attn + moca_trail + moca_fluency + moca_abstr
  memory   =~ np_hvlt_recall + np_hvlt_retention + np_hvlt_recog + moca_memory
  visuo    =~ np_jlo + moca_visuo_draw

  g ~~ 0*att_exec
  g ~~ 0*memory
  g ~~ 0*visuo
  att_exec ~~ 0*memory
  att_exec ~~ 0*visuo
  memory   ~~ 0*visuo
'

# --- PRIMARY Model 3: RESPECIFIED Bifactor ---
# Symmetric bifactor with three targeted respecifications to obtain a
# well-identified, fit-stable solution:
#   (1) TMT-A <-> TMT-B residual covariance is freed (absorbs the shared
#       pencil-and-paper / timing method variance that otherwise inflates
#       the att/exec specific factor).
#   (2) HVLT recall <-> retention residual covariance is freed (retention is
#       algebraically derived from recall, so they share method variance).
#   (3) The two visuospatial indicators (Benton JLO, MoCA visuoconstruction)
#       are constrained to load equally on the visuo specific factor. With
#       only two indicators this is required for identification. In practice
#       JLO carries the bulk of the signal; the specific factor primarily
#       captures visuospatial variance above and beyond g.
# All three specific factors are orthogonal to g and to each other
# (symmetric bifactor constraints retained).
bifactor_respec_model <- '
  g =~ np_lns + np_sdmt + np_tmta_inv + np_tmtb_inv + np_semflu + moca_attn + moca_trail + moca_fluency + moca_abstr +
       np_hvlt_recall + np_hvlt_retention + np_hvlt_recog + moca_memory +
       np_jlo + moca_visuo_draw +
       moca_naming + moca_repetition + moca_orientation

  att_exec =~ np_lns + np_sdmt + np_tmta_inv + np_tmtb_inv + np_semflu + moca_attn + moca_trail + moca_fluency + moca_abstr
  memory   =~ np_hvlt_recall + np_hvlt_retention + np_hvlt_recog + moca_memory
  visuo    =~ v*np_jlo + v*moca_visuo_draw

  np_tmta_inv    ~~ np_tmtb_inv
  np_hvlt_recall ~~ np_hvlt_retention

  g ~~ 0*att_exec
  g ~~ 0*memory
  g ~~ 0*visuo
  att_exec ~~ 0*memory
  att_exec ~~ 0*visuo
  memory   ~~ 0*visuo
'

# ==========================================================================
#  6. FIT MODELS (STRICTLY ON BASELINE DATA)
# ==========================================================================
cat("[5] Fitting CFA models on BASELINE data (Year 0)...\n")

fit_opts <- list(estimator = "MLR", missing = "fiml", std.lv = TRUE)

# PRIMARY
fit_corr        <- tryCatch({ do.call(cfa, c(list(model = correlated_model,   data = T_pd_bl), fit_opts)) }, error = function(e) { NULL })
fit_secondorder <- tryCatch({ do.call(cfa, c(list(model = second_order_model, data = T_pd_bl), fit_opts)) }, error = function(e) { NULL })

# PRIMARY bifactor (respecified)
fit_bifactor_respec <- tryCatch({
  suppressWarnings(do.call(cfa, c(list(model = bifactor_respec_model, data = T_pd_bl), fit_opts)))
}, error = function(e) { NULL })

# REFERENCE / archival
fit_uni  <- tryCatch({ do.call(cfa, c(list(model = unidim_model, data = T_pd_bl), fit_opts)) }, error = function(e) { NULL })
fit_bifactor <- tryCatch({
  suppressWarnings(do.call(cfa, c(list(model = bifactor_model, data = T_pd_bl), fit_opts)))
}, error = function(e) { NULL })

# Quick fit-summary printout so the user sees primary vs reference at a glance
safe_fit <- function(fit, idx = c("cfi.robust","tli.robust","rmsea.robust","srmr","aic","bic")) {
  if (is.null(fit) || !lavInspect(fit, "converged")) return(setNames(rep(NA, length(idx)), idx))
  fm <- tryCatch(fitMeasures(fit, idx), error = function(e) setNames(rep(NA, length(idx)), idx))
  fm
}
fit_tbl <- rbind(
  Correlated        = safe_fit(fit_corr),
  SecondOrder       = safe_fit(fit_secondorder),
  Bifactor_respec   = safe_fit(fit_bifactor_respec),
  Bifactor_original = safe_fit(fit_bifactor),
  Unidim            = safe_fit(fit_uni)
)
cat("\n--- Baseline model fit (robust indices) ---\n")
print(round(fit_tbl, 3))
write_csv(cbind(Model = rownames(fit_tbl), as.data.frame(round(fit_tbl, 3))),
          file.path(outputDir, "CFA_model_fit_comparison.csv"))

# ==========================================================================
#  7. EXTRACT FACTOR SCORES (PROJECT ONTO LONGITUDINAL)
# ==========================================================================
cat("\n[6] Extracting factor scores for full longitudinal cohort...\n")

# --- PRIMARY: Correlated-factors domain scores ---
if (is.null(fit_corr) || !lavInspect(fit_corr, "converged")) {
  stop("Correlated-factors model failed to converge. Check data distribution.")
}
fs_domains <- as.data.frame(lavPredict(fit_corr, newdata = T_pd, method = "regression", type = "lv"))
names(fs_domains) <- c("cfa_att_exec", "cfa_memory", "cfa_visuo")

# --- PRIMARY: Second-order g (higher-order global score) ---
# Rationale: with 3 first-order factors, the only unique output of the
# second-order model (vs. the correlated model) is the higher-order g score.
# The first-order domain scores from the second-order model are either
# identical to the correlated-model domain scores (well-identified case) or
# pathological (if any domain->g loading is at the 1.0 boundary). We therefore
# export only so_g here and rely on cfa_att_exec / cfa_memory / cfa_visuo for
# the domain scores.
if (!is.null(fit_secondorder) && lavInspect(fit_secondorder, "converged")) {
  so_loadings <- standardizedSolution(fit_secondorder) %>%
    dplyr::filter(op == "=~", lhs == "g") %>%
    dplyr::select(domain = rhs, est.std)
  heywood <- any(abs(so_loadings$est.std) > 0.99)
  if (heywood) {
    cat("    Second-order g converged but with Heywood-boundary loading(s); ",
        "only so_g is exported (domain residuals are not identified here).\n", sep = "")
  } else {
    cat("    Second-order model converged cleanly. so_g exported.\n")
  }
  fs_so_all <- as.data.frame(lavPredict(fit_secondorder, newdata = T_pd,
                                        method = "regression", type = "lv"))
  fs_so <- data.frame(so_g = fs_so_all[["g"]])
} else {
  cat("    Second-order model failed to converge. NAs written for so_g.\n")
  fs_so <- data.frame(so_g = rep(NA_real_, nrow(T_pd)))
}

# --- REFERENCE: Unidimensional g ---
if (!is.null(fit_uni) && lavInspect(fit_uni, "converged")) {
  fs_global <- as.data.frame(lavPredict(fit_uni, newdata = T_pd,
                                        method = "regression", type = "lv"))
  names(fs_global) <- "cfa_g"
} else {
  fs_global <- data.frame(cfa_g = NA_real_)[rep(1, nrow(T_pd)), , drop = FALSE]
}

# --- PRIMARY: Original bifactor (g + 3 orthogonal specific factors) ---
# These scores preserve the original bifactor_* column values used in
# downstream work. The bifactor_g score itself is well-defined despite the
# model's vcov singularity (the singularity affects parameter standard
# errors, not the point-estimate factor score predictions).
if (!is.null(fit_bifactor) && lavInspect(fit_bifactor, "converged")) {
  cat("    Original bifactor converged. g + 3 specifics exported as bifactor_*.\n")
  fs_bf <- as.data.frame(lavPredict(fit_bifactor, newdata = T_pd,
                                    method = "regression", type = "lv"))
  fs_bf <- fs_bf[, c("g", "att_exec", "memory", "visuo"), drop = FALSE]
  names(fs_bf) <- c("bifactor_g", "bifactor_att_exec",
                    "bifactor_memory", "bifactor_visuo")
} else {
  cat("    Original bifactor failed. NAs written for bifactor_* columns.\n")
  fs_bf <- data.frame(bifactor_g = NA_real_, bifactor_att_exec = NA_real_,
                      bifactor_memory = NA_real_, bifactor_visuo = NA_real_)[rep(1, nrow(T_pd)), ]
}

# --- ALTERNATIVE: Respecified bifactor (kept available under bifactor_respec_*) ---
# TMT residual free + HVLT residual free + visuo loadings equality-constrained.
# Fits cleanly (CFI 0.998) but specific factors carry little reliable variance.
# Do not swap into bifactor_* without informing downstream analyses.
if (!is.null(fit_bifactor_respec) && lavInspect(fit_bifactor_respec, "converged")) {
  fs_bf_respec <- as.data.frame(lavPredict(fit_bifactor_respec, newdata = T_pd,
                                           method = "regression", type = "lv"))
  fs_bf_respec <- fs_bf_respec[, c("g", "att_exec", "memory", "visuo"), drop = FALSE]
  names(fs_bf_respec) <- c("bifactor_respec_g", "bifactor_respec_att_exec",
                           "bifactor_respec_memory", "bifactor_respec_visuo")
} else {
  fs_bf_respec <- data.frame(bifactor_respec_g = NA_real_, bifactor_respec_att_exec = NA_real_,
                             bifactor_respec_memory = NA_real_, bifactor_respec_visuo = NA_real_)[rep(1, nrow(T_pd)), ]
}

# --- PRIMARY: Observed composite scores (mean of z-scored indicators) -----
# Rationale: these are the simplest, most transparent domain scores. They
# capture ALL of the domain's variance (both g-related and specific), and
# are intended to be used in downstream regressions with an explicit global
# cognition covariate (e.g., bifactor_g) so that reported effects reflect
# domain contributions above and beyond general cognition. Minimum valid
# items per composite: att_exec >= 5/9, memory >= 3/4, visuo >= 1/2,
# whole-battery >= 10/18.
obs_items <- list(
  att_exec   = c("np_lns","np_sdmt","np_tmta_inv","np_tmtb_inv","np_semflu",
                 "moca_attn","moca_trail","moca_fluency","moca_abstr"),
  memory     = c("np_hvlt_recall","np_hvlt_retention","np_hvlt_recog","moca_memory"),
  visuo      = c("np_jlo","moca_visuo_draw"),
  whole      = c("np_lns","np_sdmt","np_tmta_inv","np_tmtb_inv","np_semflu",
                 "moca_attn","moca_trail","moca_fluency","moca_abstr",
                 "np_hvlt_recall","np_hvlt_retention","np_hvlt_recog","moca_memory",
                 "np_jlo","moca_visuo_draw",
                 "moca_naming","moca_repetition","moca_orientation")
)
min_valid <- list(att_exec = 5, memory = 3, visuo = 1, whole = 10)

mean_composite <- function(df, items, min_valid) {
  M <- as.matrix(df[, items, drop = FALSE])
  n_valid <- rowSums(!is.na(M))
  m <- rowMeans(M, na.rm = TRUE)
  m[n_valid < min_valid] <- NA_real_
  m
}

fs_obs <- data.frame(
  obs_g        = mean_composite(T_pd, obs_items$whole,    min_valid$whole),
  obs_att_exec = mean_composite(T_pd, obs_items$att_exec, min_valid$att_exec),
  obs_memory   = mean_composite(T_pd, obs_items$memory,   min_valid$memory),
  obs_visuo    = mean_composite(T_pd, obs_items$visuo,    min_valid$visuo)
)

# --- Residualised observed domain composites (regressed on bifactor_g) ----
# Each domain composite is residualised against bifactor_g (the primary global
# cognition score: omega-H = 0.76, all 18 indicators, FIML handles missing
# data so coverage is ~97% vs ~79% for obs_g). Resulting scores are orthogonal
# to bifactor_g by construction -- they can be used directly as predictors
# without adding g as a covariate. Between themselves the residualised domains
# retain some correlation (real construct overlap above what g captures);
# include the other domains as covariates if you need full mutual adjustment.
resid_vs_g <- function(y, x) {
  out <- rep(NA_real_, length(y))
  ok  <- is.finite(y) & is.finite(x)
  if (sum(ok) < 30) return(out)
  fit <- lm(y[ok] ~ x[ok])
  out[ok] <- residuals(fit)
  out
}
fs_obs$obs_att_exec_resid <- resid_vs_g(fs_obs$obs_att_exec, fs_bf$bifactor_g)
fs_obs$obs_memory_resid   <- resid_vs_g(fs_obs$obs_memory,   fs_bf$bifactor_g)
fs_obs$obs_visuo_resid    <- resid_vs_g(fs_obs$obs_visuo,    fs_bf$bifactor_g)

# Cronbach's alpha on the baseline sample for each composite (reliability
# check for the observed domain scores).
cronbach_alpha <- function(df_bl, items) {
  M <- as.matrix(df_bl[, items, drop = FALSE])
  M <- M[complete.cases(M), , drop = FALSE]
  if (nrow(M) < 50 || ncol(M) < 2) return(NA_real_)
  k <- ncol(M); var_items <- apply(M, 2, var); var_total <- var(rowSums(M))
  (k / (k - 1)) * (1 - sum(var_items) / var_total)
}
alpha_tbl <- data.frame(
  composite = names(obs_items),
  n_items   = sapply(obs_items, length),
  alpha     = sapply(obs_items, function(it) cronbach_alpha(T_pd_bl, it))
)
alpha_tbl$alpha <- round(alpha_tbl$alpha, 3)
cat("\n--- Observed composite Cronbach's alpha (baseline) ---\n")
print(alpha_tbl)
write_csv(alpha_tbl, file.path(outputDir, "Observed_composite_alpha.csv"))

# Combine: PRIMARY (bifactor + observed composites + correlated domains),
# then ALTERNATIVE (respecified bifactor + second-order g) and REFERENCE (unidim g).
fscores_all <- bind_cols(fs_bf, fs_obs, fs_so, fs_domains, fs_global, fs_bf_respec)


# ==========================================================================
#  7b. RELIABILITY (omega) AND FACTOR CORRELATIONS FOR PRIMARY MODEL
# ==========================================================================
cat("\n--- Reliability (coefficient omega) for correlated-factors model ---\n")
# Compute omega manually from standardized loadings and residuals -- robust
# against compRelSEM failures when SE/vcov cannot be inverted.
omega_manual <- function(fit) {
  std <- standardizedSolution(fit)
  load_df <- std %>% filter(op == "=~")
  resid_df <- std %>% filter(op == "~~", lhs == rhs,
                             !lhs %in% unique(load_df$lhs))
  out <- sapply(unique(load_df$lhs), function(fac) {
    lam   <- load_df$est.std[load_df$lhs == fac]
    items <- load_df$rhs[load_df$lhs == fac]
    theta <- resid_df$est.std[resid_df$lhs %in% items]
    (sum(lam))^2 / ((sum(lam))^2 + sum(theta))
  })
  out
}
omega_tbl <- tryCatch(omega_manual(fit_corr), error = function(e) NULL)
if (!is.null(omega_tbl)) {
  omega_df <- data.frame(Factor = names(omega_tbl),
                         omega  = round(as.numeric(omega_tbl), 3))
  print(omega_df)
  write_csv(omega_df, file.path(outputDir, "CFA_correlated_reliability.csv"))
}

cat("\n--- First-order factor intercorrelations (phi) ---\n")
phi <- lavInspect(fit_corr, "cor.lv")
print(round(phi, 3))
write_csv(cbind(Factor = rownames(phi), as.data.frame(round(phi, 3))),
          file.path(outputDir, "CFA_correlated_phi.csv"))

if (!is.null(fit_secondorder) && lavInspect(fit_secondorder, "converged")) {
  cat("\n--- Second-order loadings (domain -> g) ---\n")
  so_load <- standardizedSolution(fit_secondorder) %>%
    dplyr::filter(op == "=~", lhs == "g") %>%
    dplyr::select(from = lhs, to = rhs, std_loading = est.std, se, pvalue) %>%
    dplyr::mutate(std_loading = round(std_loading, 3))
  print(so_load)
  write_csv(so_load, file.path(outputDir, "CFA_secondorder_loadings.csv"))
}

# ==========================================================================
#  7c. ORIGINAL BIFACTOR (the manuscript's bifactor_g):
#      STANDARDISED LOADINGS AND RELIABILITY INDICES
# --------------------------------------------------------------------------
#  Model: fit_bifactor (bifactor_model, "REFERENCE Model B: ORIGINAL
#  Bifactor"), whose g is exported as bifactor_g. NOT the respecified
#  bifactor (fit_bifactor_respec -> bifactor_respec_*).
#  Explicit formulas on the completely standardised solution
#  (standardizedSolution(fit_bifactor)). For indicator i (18 indicators):
#  lg_i = loading on g; lk_i = loading on specific factor k in {att_exec,
#  memory, visuo} (0 if i does not load on k; each indicator loads on at
#  most one specific factor); th_i = standardised residual variance
#  (= 1 - lg_i^2 - lk_i^2: all factors orthogonal, no residual covariances).
#  S = sum over the indicators named; I_k = indicators of specific factor k.
#    omega (total) = [ (S lg)^2 + sum_k (S_{I_k} lk)^2 ] / D
#                    D = (S lg)^2 + sum_k (S_{I_k} lk)^2 + S th     (all 18)
#    omega-H       = (S lg)^2 / D
#    ECV           = S lg^2 / ( S lg^2 + sum_k S_{I_k} lk^2 )
#    PUC           = 1 - sum_k C(n_k, 2) / C(18, 2), n_k = size of I_k
#                    (share of indicator pairs that share no specific
#                    factor; the g-only indicators share none)
#    omega-S_k     = [ (S_{I_k} lg)^2 + (S_{I_k} lk)^2 ] / D_k
#                    D_k = (S_{I_k} lg)^2 + (S_{I_k} lk)^2 + S_{I_k} th
#    omega-HS_k    = (S_{I_k} lk)^2 / D_k
#  Baseline (YEAR 0) convergent validity of the bifactor_g factor score
#  (Pearson r, pairwise complete):
#    r(g, MoCA total) with the curated `moca` (education-adjusted MoCA total)
#    r(g, mean of the 18 standardised indicators) with obs_g (mean of the
#      18 z-scored indicators, >= 10 valid, as defined in section 7)
#  Bifactor_std_loadings.csv: indicator, g, att_exec, memory, visuo,
#    resid_std (layout of the earlier file).
#  Bifactor_reliability.csv: layout of the earlier file (Factor, omega,
#    omega_hierarchical, ECV); g row: omega = omega total, omega_hierarchical
#    = omega-H; specific-factor rows: omega = omega-S_k, omega_hierarchical
#    = omega-HS_k. Added columns: PUC, the two baseline correlations with
#    their n, the fitted n, a note on the indices and the model name.
# ==========================================================================
cat("\n--- [MOCA-FIX] Original bifactor (bifactor_g): standardised loadings and reliability ---\n")
if (!is.null(fit_bifactor) && lavInspect(fit_bifactor, "converged")) {
  bf_std  <- standardizedSolution(fit_bifactor, se = FALSE, zstat = FALSE,
                                  pvalue = FALSE, ci = FALSE)
  bf_ind  <- lavNames(fit_bifactor, "ov")
  bf_fac  <- c("g", "att_exec", "memory", "visuo")
  bf_lam  <- bf_std[bf_std$op == "=~", ]
  bf_load <- data.frame(indicator = bf_ind, stringsAsFactors = FALSE)
  for (f in bf_fac) {
    lf <- bf_lam[bf_lam$lhs == f, ]
    bf_load[[f]] <- ifelse(bf_ind %in% lf$rhs, lf$est.std[match(bf_ind, lf$rhs)], 0)
  }
  bf_th <- bf_std[bf_std$op == "~~" & bf_std$lhs == bf_std$rhs & bf_std$lhs %in% bf_ind, ]
  bf_load$resid_std <- bf_th$est.std[match(bf_ind, bf_th$lhs)]
  stopifnot(length(bf_ind) == 18, all(bf_load$g != 0),
            all(rowSums(bf_load[, bf_fac[-1]] != 0) <= 1))
  if (any(bf_load$resid_std < 0) || any(abs(as.matrix(bf_load[, bf_fac])) >= 1))
    cat("    WARNING: Heywood case in the original bifactor standardised solution.\n")
  write_csv(bf_load, file.path(outputDir, "Bifactor_std_loadings.csv"))

  spec_lab <- c(att_exec = "Att/Exec (specific)", memory = "Memory (specific)",
                visuo = "Visuospatial (specific)")
  S_g   <- sum(bf_load$g)
  S_k   <- sapply(names(spec_lab), function(k) sum(bf_load[[k]]))
  D_all <- S_g^2 + sum(S_k^2) + sum(bf_load$resid_std)
  bf_omega_total <- (S_g^2 + sum(S_k^2)) / D_all
  bf_omega_H     <- S_g^2 / D_all
  bf_ECV <- sum(bf_load$g^2) /
            (sum(bf_load$g^2) + sum(sapply(names(spec_lab), function(k) sum(bf_load[[k]]^2))))
  n_k    <- sapply(names(spec_lab), function(k) sum(bf_load[[k]] != 0))
  bf_PUC <- 1 - sum(choose(n_k, 2)) / choose(nrow(bf_load), 2)
  bf_omega_S <- bf_omega_HS <- setNames(numeric(length(spec_lab)), names(spec_lab))
  for (k in names(spec_lab)) {
    ik  <- bf_load[[k]] != 0
    a   <- sum(bf_load$g[ik]); b <- sum(bf_load[[k]][ik])
    D_k <- a^2 + b^2 + sum(bf_load$resid_std[ik])
    bf_omega_S[k]  <- (a^2 + b^2) / D_k
    bf_omega_HS[k] <- b^2 / D_k
  }
  bl_rows <- which(T_pd$YEAR_num == 0)
  g_bl    <- fs_bf$bifactor_g[bl_rows]
  moca_bl <- safe_numeric(T_pd$moca[bl_rows])
  obsg_bl <- fs_obs$obs_g[bl_rows]
  ok_m <- is.finite(g_bl) & is.finite(moca_bl)
  ok_o <- is.finite(g_bl) & is.finite(obsg_bl)
  r_g_moca <- if (sum(ok_m) > 2) cor(g_bl[ok_m], moca_bl[ok_m]) else NA_real_
  r_g_obsg <- if (sum(ok_o) > 2) cor(g_bl[ok_o], obsg_bl[ok_o]) else NA_real_
  bf_rel <- data.frame(
    Factor             = c("General (g)", unname(spec_lab)),
    omega              = round(c(bf_omega_total, bf_omega_S), 3),
    omega_hierarchical = round(c(bf_omega_H, bf_omega_HS), 3),
    ECV                = c(round(bf_ECV, 3), NA, NA, NA),
    PUC                = c(round(bf_PUC, 3), NA, NA, NA),
    r_g_MoCA_total_BL  = c(round(r_g_moca, 3), NA, NA, NA),
    n_g_MoCA_total_BL  = c(sum(ok_m), NA, NA, NA),
    r_g_mean18z_BL     = c(round(r_g_obsg, 3), NA, NA, NA),
    n_g_mean18z_BL     = c(sum(ok_o), NA, NA, NA),
    n_fit_BL           = c(lavInspect(fit_bifactor, "nobs"), NA, NA, NA),
    index_note         = c("omega = omega total; omega_hierarchical = omega-H",
                           rep("omega = omega-S (subscale); omega_hierarchical = omega-HS (subscale)", 3)),
    model              = "Bifactor_original (fit_bifactor; exported as bifactor_g)",
    stringsAsFactors   = FALSE
  )
  print(bf_load %>% mutate(across(where(is.numeric), ~ round(.x, 3))), row.names = FALSE)
  print(bf_rel[, 1:10], row.names = FALSE)
  write_csv(bf_rel, file.path(outputDir, "Bifactor_reliability.csv"))
} else {
  cat("    Original bifactor did not converge: Bifactor_std_loadings.csv and Bifactor_reliability.csv NOT written.\n")
}

# ==========================================================================
#  X. LONGITUDINAL MEASUREMENT INVARIANCE
#     Tests on the bifactor scoring model (primary) and the unidimensional
#     model on the same 18 indicators (sensitivity).
# ==========================================================================
cat("\n========================================\n")
cat("  LONGITUDINAL MEASUREMENT INVARIANCE\n")
cat("========================================\n\n")

time_col <- "YEAR_num"
time_points <- sort(unique(T_pd[[time_col]]))
t0 <- time_points[1]
t1_idx <- which(time_points >= 2)[1]
if (is.na(t1_idx)) t1_idx <- length(time_points)
t1 <- time_points[t1_idx]

df_mi <- T_pd %>% filter(.data[[time_col]] %in% c(t0, t1)) %>% mutate(time_grp = factor(.data[[time_col]]))

get_lmi_fit <- function(fit, name) {
  if (is.null(fit) || !lavInspect(fit, "converged")) return(data.frame(Model = name, cfi = NA, tli = NA, rmsea = NA))
  fm <- fitMeasures(fit, c("cfi.robust", "tli.robust", "rmsea.robust"))
  data.frame(Model = name, cfi = round(fm[1], 3), tli = round(fm[2], 3), rmsea = round(fm[3], 3))
}

get_lmi_fit_flex <- function(fit, name) {
  # Prefer robust indices; fall back to standard CFI/TLI/RMSEA if robust
  # versions can't be computed (occurs with bifactor multi-group + MLR
  # when the baseline-model information matrix is singular).
  if (is.null(fit) || !lavInspect(fit, "converged")) {
    return(data.frame(Model = name, cfi = NA, tli = NA, rmsea = NA))
  }
  fm_robust <- tryCatch(fitMeasures(fit, c("cfi.robust","tli.robust","rmsea.robust")),
                        error = function(e) c(NA, NA, NA))
  if (any(is.na(fm_robust))) {
    fm_std <- tryCatch(fitMeasures(fit, c("cfi","tli","rmsea")),
                       error = function(e) c(NA, NA, NA))
    return(data.frame(Model = name, cfi = round(fm_std[1], 3),
                      tli = round(fm_std[2], 3), rmsea = round(fm_std[3], 3)))
  }
  data.frame(Model = name, cfi = round(fm_robust[1], 3),
             tli = round(fm_robust[2], 3), rmsea = round(fm_robust[3], 3))
}

run_mi <- function(model_syntax, model_name, df_mi) {
  cat(sprintf("\n  >> Invariance test: %s\n", model_name))
  out <- tryCatch({
    mi_config <- suppressWarnings(cfa(model = model_syntax, data = df_mi, group = "time_grp",
                                      estimator = "MLR", missing = "fiml", std.lv = TRUE))
    mi_metric <- suppressWarnings(cfa(model = model_syntax, data = df_mi, group = "time_grp",
                                      group.equal = "loadings",
                                      estimator = "MLR", missing = "fiml", std.lv = TRUE))
    mi_scalar <- suppressWarnings(cfa(model = model_syntax, data = df_mi, group = "time_grp",
                                      group.equal = c("loadings","intercepts"),
                                      estimator = "MLR", missing = "fiml", std.lv = TRUE))
    rbind(
      cbind(Model = model_name, Step = "Configural", get_lmi_fit_flex(mi_config, "Configural")[, -1]),
      cbind(Model = model_name, Step = "Metric",     get_lmi_fit_flex(mi_metric, "Metric")[, -1]),
      cbind(Model = model_name, Step = "Scalar",     get_lmi_fit_flex(mi_scalar, "Scalar")[, -1])
    )
  }, error = function(e) {
    cat(sprintf("    %s invariance failed: %s\n", model_name, conditionMessage(e)))
    data.frame(Model = model_name, Step = c("Configural","Metric","Scalar"),
               cfi = NA_real_, tli = NA_real_, rmsea = NA_real_)
  })
  print(out, row.names = FALSE)
  if (all(!is.na(out$cfi))) {
    d1 <- out$cfi[1] - out$cfi[2]; d2 <- out$cfi[2] - out$cfi[3]
    cat(sprintf("    Delta-CFI configural->metric: %+.4f  %s\n", d1, ifelse(abs(d1) < 0.010, "PASS", "FAIL")))
    cat(sprintf("    Delta-CFI metric->scalar:     %+.4f  %s\n", d2, ifelse(abs(d2) < 0.010, "PASS", "FAIL")))
  }
  out
}

if (min(sum(df_mi$time_grp == t0), sum(df_mi$time_grp == t1)) >= 100) {
  # Invariance tests for the actual scoring models:
  #   * Bifactor (PRIMARY): the model from which `bifactor_g` is extracted.
  #   * Unidimensional (SENSITIVITY): same 18-indicator system, single-factor
  #     parameterisation -- corroborates that the indicator measurement
  #     structure underlying `bifactor_g` is stable across time.
  mi_bf       <- run_mi(bifactor_model,        "Bifactor",         df_mi)
  mi_bf_respec <- run_mi(bifactor_respec_model, "Bifactor_respec",  df_mi)
  mi_uni      <- run_mi(unidim_model,          "Unidim",           df_mi)

  # ---- Partial measurement invariance (unidim model on all 18 indicators) --
  # Strategy: fit metric (all loadings equal), then iteratively free the
  # single loading with the largest modification index (largest evidence of
  # cross-time non-invariance) until Delta-CFI(configural -> partial-metric)
  # falls below 0.010. Report the freed indicators. Repeat for scalar
  # (intercepts), conditional on the final partial-metric solution.
  cat("\n\n  >> PARTIAL invariance test: Unidim (free the most non-invariant loadings/intercepts)\n")

  fit_uni_config <- suppressWarnings(cfa(unidim_model, data = df_mi, group = "time_grp",
                                         estimator = "MLR", missing = "fiml", std.lv = TRUE))
  cfi_config <- get_lmi_fit_flex(fit_uni_config, "Configural")$cfi

  fit_uni_metric <- suppressWarnings(cfa(unidim_model, data = df_mi, group = "time_grp",
                                         group.equal = "loadings",
                                         estimator = "MLR", missing = "fiml", std.lv = TRUE))
  cfi_metric_full <- get_lmi_fit_flex(fit_uni_metric, "Metric")$cfi

  partial_loadings <- character(0)
  cur_metric <- fit_uni_metric
  cur_cfi    <- cfi_metric_full
  max_releases <- 6  # cap to avoid pathological loops
  while (!is.na(cur_cfi) && (cfi_config - cur_cfi) > 0.010 && length(partial_loadings) < max_releases) {
    mi_df <- modificationIndices(cur_metric, op = "=~", free.remove = FALSE) %>%
      filter(group == 1) %>%   # constraints applied across both groups
      arrange(desc(mi))
    if (nrow(mi_df) == 0) break
    target  <- mi_df$rhs[1]
    if (target %in% partial_loadings) break
    partial_loadings <- c(partial_loadings, target)
    cat(sprintf("    Freeing loading: %s  (MI = %.1f)\n", target, mi_df$mi[1]))

    cur_metric <- suppressWarnings(cfa(unidim_model, data = df_mi, group = "time_grp",
                                       group.equal     = "loadings",
                                       group.partial   = paste0("g=~", partial_loadings),
                                       estimator = "MLR", missing = "fiml", std.lv = TRUE))
    cur_cfi <- get_lmi_fit_flex(cur_metric, "PartialMetric")$cfi
    cat(sprintf("    -> Partial Metric CFI = %.3f, Delta = %.4f\n",
                cur_cfi, cfi_config - cur_cfi))
  }
  fit_uni_partial_metric <- cur_metric
  cfi_partial_metric     <- cur_cfi

  partial_intercepts <- character(0)
  fit_uni_scalar <- suppressWarnings(cfa(unidim_model, data = df_mi, group = "time_grp",
                                         group.equal   = c("loadings","intercepts"),
                                         group.partial = paste0("g=~", partial_loadings),
                                         estimator = "MLR", missing = "fiml", std.lv = TRUE))
  cur_scalar <- fit_uni_scalar
  cur_cfi_s  <- get_lmi_fit_flex(cur_scalar, "Scalar")$cfi
  while (!is.na(cur_cfi_s) && (cfi_partial_metric - cur_cfi_s) > 0.010 &&
         length(partial_intercepts) < max_releases) {
    mi_df <- modificationIndices(cur_scalar, op = "~1", free.remove = FALSE) %>%
      filter(group == 1) %>% arrange(desc(mi))
    if (nrow(mi_df) == 0) break
    target <- mi_df$lhs[1]
    if (target %in% partial_intercepts) break
    partial_intercepts <- c(partial_intercepts, target)
    cat(sprintf("    Freeing intercept: %s  (MI = %.1f)\n", target, mi_df$mi[1]))

    free_intercepts <- paste0(partial_intercepts, "~1")
    cur_scalar <- suppressWarnings(cfa(unidim_model, data = df_mi, group = "time_grp",
                                       group.equal   = c("loadings","intercepts"),
                                       group.partial = c(paste0("g=~", partial_loadings), free_intercepts),
                                       estimator = "MLR", missing = "fiml", std.lv = TRUE))
    cur_cfi_s <- get_lmi_fit_flex(cur_scalar, "PartialScalar")$cfi
    cat(sprintf("    -> Partial Scalar CFI = %.3f, Delta = %.4f\n",
                cur_cfi_s, cfi_partial_metric - cur_cfi_s))
  }

  cat("\n  Partial invariance summary:\n")
  cat(sprintf("    Configural CFI       = %.3f\n", cfi_config))
  cat(sprintf("    Partial Metric CFI   = %.3f   (loadings freed: %s)\n",
              cfi_partial_metric,
              if (length(partial_loadings)) paste(partial_loadings, collapse = ", ") else "none"))
  cat(sprintf("    Partial Scalar CFI   = %.3f   (intercepts freed: %s)\n",
              cur_cfi_s,
              if (length(partial_intercepts)) paste(partial_intercepts, collapse = ", ") else "none"))

  partial_row <- data.frame(
    Model = "Unidim_partial",
    Step  = c("Configural", "PartialMetric", "PartialScalar"),
    cfi   = c(cfi_config, cfi_partial_metric, cur_cfi_s),
    tli   = c(get_lmi_fit_flex(fit_uni_config, "")$tli,
              get_lmi_fit_flex(fit_uni_partial_metric, "")$tli,
              get_lmi_fit_flex(cur_scalar, "")$tli),
    rmsea = c(get_lmi_fit_flex(fit_uni_config, "")$rmsea,
              get_lmi_fit_flex(fit_uni_partial_metric, "")$rmsea,
              get_lmi_fit_flex(cur_scalar, "")$rmsea)
  )

  partial_meta <- data.frame(
    Model = c("Unidim_partial_freed_loadings","Unidim_partial_freed_intercepts"),
    Step  = c(paste(partial_loadings, collapse = "|"),
              paste(partial_intercepts, collapse = "|")),
    cfi = NA_real_, tli = NA_real_, rmsea = NA_real_
  )

  mi_all <- bind_rows(mi_bf, mi_bf_respec, mi_uni, partial_row, partial_meta)
  write_csv(mi_all, file.path(outputDir, "CFA_measurement_invariance.csv"))
}

# ==========================================================================
#  8. EXPORT TO MATLAB
# ==========================================================================
cat("\n[7] Exporting factor scores...\n")

export_df <- T_pd %>%
  select(PATNO, any_of(c("EVENT_ID", "YEAR_num"))) %>%
  bind_cols(fscores_all) %>%
  rename(YEAR = YEAR_num)

export_path <- file.path(outputDir, "cfa_factor_scores.csv")
write_csv(export_df, export_path)
cat(sprintf("\n    -> %s successfully written for MATLAB.\n", export_path))
cat("\nColumn guide:\n",
    "  PRIMARY -- latent bifactor scores (original, preserved):\n",
    "    bifactor_g                                -> latent global cognition (all 18 indicators)\n",
    "    bifactor_att_exec / bifactor_memory /     -> specific factors (orthogonal to g\n",
    "    bifactor_visuo                              and each other by design).\n",
    "  PRIMARY -- observed composites (RECOMMENDED for domain-level analyses):\n",
    "    obs_g                                     -> mean of all 18 z-scored indicators\n",
    "    obs_att_exec / obs_memory / obs_visuo     -> mean of z-scored domain indicators\n",
    "    obs_att_exec_resid / obs_memory_resid /   -> domain composites residualised on\n",
    "    obs_visuo_resid                              bifactor_g (orthogonal to g;\n",
    "                                                 ready-to-use as predictors).\n",
    "  ALTERNATIVE:\n",
    "    bifactor_respec_*                         -> respecified bifactor (CFI 0.998,\n",
    "                                                 TMT + HVLT residuals freed, visuo\n",
    "                                                 loadings equality-constrained).\n",
    "    so_g                                      -> second-order g (15-indicator subset)\n",
    "    cfa_att_exec / cfa_memory / cfa_visuo     -> correlated-factors domains\n",
    "  REFERENCE:\n",
    "    cfa_g                                     -> unidimensional g\n",
    sep = "")