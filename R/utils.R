# =============================================================================
# utils.R — Shared constants, palettes, theme, and functions
# Source this file at the top of every analysis Rmd:
#   source(here::here("R/utils.R"))
# =============================================================================

library(tidyverse)
library(vegan)
library(phyloseq)
library(purrr)

# -----------------------------------------------------------------------------
# Constants
# -----------------------------------------------------------------------------

SEED <- 711

# Rarefaction depths — one named constant per analysis level
RAREFY_DEPTH_ALL    <- 435   # whole-community rarefaction
RAREFY_DEPTH_CURTO  <- 165   # Curtobacterium-only
RAREFY_DEPTH_SUBCL  <- 161   # subclade-level avgdist / NMDS
RAREFY_DEPTH_SUBCL2 <- 162   # centroid pipeline
RAREFY_DEPTH_TOP4   <- 143   # top-4-subclades euclidean
N_ITER              <- 100   # rarefaction iterations

# Data paths — update if the project moves
DATA_DIR    <- here::here("data")
CA_CURTO_ROOT <- here::here("data", "external")
RESULTS_DIR <- here::here("results")
RDS_DIR     <- file.path(RESULTS_DIR, "rds")

# -----------------------------------------------------------------------------
# Color palettes
# -----------------------------------------------------------------------------

reserve_colors <- c(
  "EMR"   = "darkgoldenrod4",
  "MOJ"   = "#FA8072",
  "COP"   = "#104E8B",
  "SGWK"  = "#D1A33D",
  "BLOK"  = "#AD6F3B",
  "BOMR"  = "royalblue2",
  "McLL"  = "#7CCD7C",
  "SNARL" = "#C84248",
  "LOMA"  = "#548B54"
)

inoc_colors <- c(
  "EMR"   = "darkgoldenrod4",
  "MOJ"   = "#FA8072",
  "COP"   = "#104E8B",
  "SGWK"  = "#D1A33D",
  "BLOK"  = "#AD6F3B",
  "BOMR"  = "royalblue2",
  "McLL"  = "#7CCD7C",
  "SNARL" = "#C84248",
  "LOMA"  = "#548B54",
  "none"  = "gray51"
)

ecosystem_colors <- c(
  "Oak Woodland"   = "#8B5A2B",
  "Desert"         = "#B7410E",
  "Coastal Marine" = "#1E6091",
  "Chaparral"      = "#6B8E23",
  "Conifer Forest" = "#2E8B57",
  "Grassland"      = "#DAA520",
  "Montane"        = "#778899"
)

colors <- c(
  "#CBD588", "#5F7FC7", "orange", "#DA5724", "#508578",
  "#AD6F3B", "#673770", "#CD9BCD", "#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F", "#D1A33D", "#AD6F3B", "#8A7C64", "#599861", "red",
  "blue", "black", "darkgoldenrod4", "darkcyan", "gray", "#7ECFE8", "yellow",
  "#8B0000", "#006400", "#ADD8E6", "#A020F0", "#FFC0CB", "#483D8B", "#8B4500",
  "#FA8072", "#00A400", "#7171FF", "#B0B038", "#00FF36", "#FF69B4",
  "#CAFF70", "#CD6600", "#B4EEB4", "#104E8B", "#E6E6FA"
)

path_colors <- c(
  LOMA_clim_eff = "#023903",
  LOMA_litt_eff = "#38913A",
  LOMA_int_eff  = "#98C698",
  away_clim_eff = "#AC64B9",
  away_litt_eff = "#D5ADDD"
)

# -----------------------------------------------------------------------------
# ggplot2 theme
# -----------------------------------------------------------------------------

apatheme <- theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_blank(),
    axis.line        = element_line(),
    text             = element_text(family = "Arial"),
    axis.title       = element_text(size = 12),
    axis.text        = element_text(size = 10),
    legend.text      = element_text(size = 10),
    legend.title     = element_text(size = 11),
    strip.text       = element_text(size = 12)
  )

# -----------------------------------------------------------------------------
# Convergence path definitions
# Each path = c(S0_label, S1_label, R_label) in i_r_l grammar.
# i_r_l = "<i_code>|<reserve_cat> & <s_cat>"
# -----------------------------------------------------------------------------

paths <- list(
  # Tests: does LOMA climate signature shift community?
  LOMA_clim_eff = c("LOMA|LOMA & LOMA", "LOMA|away & LOMA",     "away|away & LOMA"),
  # Tests: does non-LOMA litter shift community away from LOMA?
  LOMA_litt_eff = c("LOMA|LOMA & LOMA", "LOMA|LOMA & non-LOMA", "away|LOMA & non-LOMA"),
  # Tests: combined climate x litter interaction
  LOMA_int_eff  = c("LOMA|LOMA & LOMA", "LOMA|away & non-LOMA", "away|away & non-LOMA"),
  # Tests: does LOMA climate shift the away community toward LOMA?
  away_clim_eff = c("away|away & non-LOMA", "away|LOMA & non-LOMA",  "LOMA|LOMA & non-LOMA"),
  # Tests: does LOMA litter shift the away community toward LOMA?
  away_litt_eff = c("away|away & non-LOMA", "away|away & LOMA",      "LOMA|away & LOMA")
)

# Path table for NMDS trajectory plotting (includes source column)
plot_paths <- tribble(
  ~source, ~path,          ~step, ~i_r_l,
  "LOMA", "Climate",        0L,   "LOMA|LOMA & LOMA",
  "LOMA", "Climate",        1L,   "LOMA|away & LOMA",
  "LOMA", "Climate",        2L,   "away|away & LOMA",
  "LOMA", "Litter",         0L,   "LOMA|LOMA & LOMA",
  "LOMA", "Litter",         1L,   "LOMA|LOMA & non-LOMA",
  "LOMA", "Litter",         2L,   "away|LOMA & non-LOMA",
  "LOMA", "Interaction",    0L,   "LOMA|LOMA & LOMA",
  "LOMA", "Interaction",    1L,   "LOMA|away & non-LOMA",
  "LOMA", "Interaction",    2L,   "away|away & non-LOMA",
  "away", "Climate",        0L,   "away|away & non-LOMA",
  "away", "Climate",        1L,   "away|LOMA & non-LOMA",
  "away", "Climate",        2L,   "LOMA|LOMA & non-LOMA",
  "away", "Litter",         0L,   "away|away & non-LOMA",
  "away", "Litter",         1L,   "away|away & LOMA",
  "away", "Litter",         2L,   "LOMA|away & LOMA"
)

# Path table for relative-abundance trajectory plots
plot_paths_RA <- tribble(
  ~source, ~path,          ~step, ~i_r_l,
  "LOMA", "Climate",        0L,   "LOMA|LOMA & LOMA",
  "LOMA", "Climate",        1L,   "LOMA|away & LOMA",
  "LOMA", "Climate",        2L,   "away|away & non-LOMA",
  "LOMA", "Litter",         0L,   "LOMA|LOMA & LOMA",
  "LOMA", "Litter",         1L,   "LOMA|LOMA & non-LOMA",
  "LOMA", "Litter",         2L,   "away|away & non-LOMA",
  "LOMA", "Interaction",    0L,   "LOMA|LOMA & LOMA",
  "LOMA", "Interaction",    1L,   "LOMA|away & non-LOMA",
  "LOMA", "Interaction",    2L,   "away|away & non-LOMA",
  "away", "Climate",        0L,   "away|away & non-LOMA",
  "away", "Climate",        1L,   "away|LOMA & non-LOMA",
  "away", "Climate",        2L,   "LOMA|LOMA & LOMA",
  "away", "Litter",         0L,   "away|away & non-LOMA",
  "away", "Litter",         1L,   "away|away & LOMA",
  "away", "Litter",         2L,   "LOMA|LOMA & LOMA"
)

# -----------------------------------------------------------------------------
# Metadata helpers
# -----------------------------------------------------------------------------

#' Add pair label column ("X vs LOMA") to a metadata data frame.
add_pair_label <- function(df) {
  df |>
    mutate(pair = case_when(
      reserve_code != "LOMA" & i_code == "LOMA" ~ paste(reserve_code, "vs", "LOMA"),
      reserve_code != "LOMA" & i_code != "LOMA" ~ paste(reserve_code, "vs", "LOMA"),
      reserve_code == "LOMA" & i_code != "LOMA" ~ paste(i_code,       "vs", "LOMA"),
      TRUE                                       ~ paste(s_code,       "vs", "LOMA")
    ))
}

#' Add reserve_litter and i_r_l columns (inoculum x reserve x litter context).
add_irl_label <- function(df) {
  df |>
    mutate(
      reserve_litter = case_when(
        reserve_code == "LOMA" & s_code == "LOMA" ~ "LOMA & LOMA",
        reserve_code == "LOMA" & s_code != "LOMA" ~ "LOMA & non-LOMA",
        reserve_code != "LOMA" & s_code == "LOMA" ~ "away & LOMA",
        reserve_code != "LOMA" & s_code != "LOMA" ~ "away & non-LOMA"
      ),
      i_r_l = case_when(
        i_code == "LOMA" ~ paste(i_code, reserve_litter, sep = "|"),
        i_code != "LOMA" ~ paste("away",  reserve_litter, sep = "|")
      )
    )
}

#' Duplicate LOMA-vs-LOMA rows so they appear in every away-vs-LOMA pair context.
#' @param df     Data frame with a `pair` column from add_pair_label().
#' @param others Away reserve codes (character vector, not including "LOMA").
expand_loma_rows <- function(df, others) {
  loma_rows <- df |> filter(i_code == "LOMA", pair == "LOMA vs LOMA")
  new_rows  <- loma_rows |>
    crossing(other = others) |>
    mutate(pair = paste(other, "vs LOMA")) |>
    select(-other)
  bind_rows(df, new_rows)
}

# -----------------------------------------------------------------------------
# Distance matrix helpers
# -----------------------------------------------------------------------------

#' Mean pairwise distance within a set of samples (upper triangle only).
mean_within <- function(D, ids) {
  ids <- intersect(ids, rownames(D))
  if (length(ids) < 2) return(NA_real_)
  M <- D[ids, ids, drop = FALSE]
  mean(M[upper.tri(M)], na.rm = TRUE)
}

#' Mean pairwise distance between two sets of samples.
mean_cross <- function(D, I, J) {
  I <- intersect(I, rownames(D)); J <- intersect(J, rownames(D))
  if (length(I) == 0 || length(J) == 0) return(NA_real_)
  mean(D[I, J, drop = FALSE], na.rm = TRUE)
}

#' Within- and between-group average distances for a labelled distance matrix.
#' @return List: between (GxG), within (G vec), n, pairs_between, pairs_within.
avg_block_dists_plus <- function(D, groups) {
  stopifnot(is.matrix(D),
            isTRUE(all(rownames(D) == colnames(D))),
            length(groups) == nrow(D))
  lev <- levels(factor(groups)); G <- length(lev)

  between       <- matrix(NA_real_, G, G, dimnames = list(lev, lev))
  pairs_between <- matrix(0L,       G, G, dimnames = list(lev, lev))
  within        <- rep(NA_real_, G); names(within)       <- lev
  pairs_within  <- rep(0L,       G); names(pairs_within) <- lev
  n_in          <- setNames(integer(G), lev)

  for (i in seq_len(G)) {
    Ii <- which(groups == lev[i]); n_in[i] <- length(Ii)
    if (n_in[i] >= 2) {
      m               <- D[Ii, Ii, drop = FALSE]
      within[i]       <- mean(m[upper.tri(m)], na.rm = TRUE)
      pairs_within[i] <- choose(n_in[i], 2)
    }
    for (j in i:G) {
      Jj <- which(groups == lev[j])
      if (i == j) {
        between[i, i]       <- within[i]
        pairs_between[i, i] <- pairs_within[i]
      } else {
        between[i, j]       <- between[j, i]       <- mean(D[Ii, Jj], na.rm = TRUE)
        pairs_between[i, j] <- pairs_between[j, i] <- length(Ii) * length(Jj)
      }
    }
  }
  list(between = between, within = within, n = n_in,
       pairs_between = pairs_between, pairs_within = pairs_within)
}

# -----------------------------------------------------------------------------
# Alpha diversity
# -----------------------------------------------------------------------------

#' Average Shannon diversity and richness via repeated rarefaction.
#' @param otu_table      Counts matrix (rows = samples, cols = OTUs).
#' @param n_iterations   Number of rarefaction iterations.
#' @param sequence_depth Rarefaction depth.
#' @param set_seed       Integer seed.
#' @return data.frame with Richness and Shannon columns, rownames = sample IDs.
calc_alpha_diversity <- function(otu_table,
                                 n_iterations   = 100,
                                 sequence_depth,
                                 set_seed       = NULL) {
  if (!requireNamespace("vegan", quietly = TRUE))
    stop("Package 'vegan' is required.")
  if (!is.null(set_seed)) set.seed(set_seed)

  otu_mat    <- as.matrix(otu_table)
  alpha_list <- vector("list", n_iterations)
  S_list     <- vector("list", n_iterations)

  for (i in seq_len(n_iterations)) {
    r             <- vegan::rrarefy(otu_mat, sample = sequence_depth)
    alpha_list[[i]] <- vegan::diversity(r, index = "shannon")
    S_list[[i]]     <- vegan::specnumber(r)
  }

  data.frame(
    Richness = round(rowMeans(do.call(cbind, S_list))),
    Shannon  = rowMeans(do.call(cbind, alpha_list)),
    row.names = rownames(otu_mat)
  )
}

# -----------------------------------------------------------------------------
# Rarefied Euclidean distances
# -----------------------------------------------------------------------------

#' Average Euclidean distance matrix over repeated rarefactions.
#' @param otu_table      Counts matrix (rows = samples, cols = OTUs).
#' @param n_iterations   Number of rarefaction iterations.
#' @param sequence_depth Rarefaction depth.
#' @param set_seed       Integer seed.
#' @return A `dist` object.
calc_rarefaction_euclidean <- function(otu_table,
                                       n_iterations   = 100,
                                       sequence_depth,
                                       set_seed       = NULL) {
  if (!requireNamespace("vegan", quietly = TRUE))
    stop("Package 'vegan' is required.")
  if (!is.null(set_seed)) set.seed(set_seed)

  otu_mat        <- as.matrix(otu_table)
  euclidean_list <- vector("list", n_iterations)

  for (i in seq_len(n_iterations)) {
    r                  <- suppressWarnings(vegan::rrarefy(otu_mat, sample = sequence_depth))
    euclidean_list[[i]] <- as.matrix(dist(r, method = "euclidean"))
  }

  euc_array  <- array(unlist(euclidean_list),
                      dim = c(nrow(otu_mat), nrow(otu_mat), n_iterations))
  avg_matrix <- apply(euc_array, c(1, 2), mean)
  rownames(avg_matrix) <- rownames(otu_mat)
  colnames(avg_matrix) <- rownames(otu_mat)
  as.dist(avg_matrix)
}

# -----------------------------------------------------------------------------
# Centroid pipeline
# -----------------------------------------------------------------------------

#' Average Bray-Curtis centroid distance matrix via repeated rarefaction.
#' @param comm       samples x taxa count matrix (rownames = sample IDs).
#' @param meta       data.frame with id_col and group_cols.
#' @param id_col     Column name for sample IDs in meta.
#' @param group_cols Columns defining a treatment cell.
#' @param depth      Rarefaction depth (NULL = min library size).
#' @param iterations Number of rarefaction iterations.
#' @param method     Dissimilarity for vegdist.
#' @param seed       Random seed.
#' @param sep        Separator for collapsing group labels.
#' @param return_centroids If TRUE, also return averaged centroid profiles.
centroid_dm_rarefy <- function(comm, meta,
                               id_col           = "original_label",
                               group_cols       = c("reserve_code", "s_code", "i_code", "Bag_type"),
                               depth            = NULL,
                               iterations       = 100,
                               method           = "bray",
                               seed             = NULL,
                               sep              = " | ",
                               return_centroids = TRUE) {
  if (!is.null(seed)) set.seed(seed)

  meta2 <- meta[, c(id_col, group_cols), drop = FALSE]
  stopifnot(all(meta2[[id_col]] %in% rownames(comm)))
  comm2 <- comm[meta2[[id_col]], , drop = FALSE]

  libs  <- rowSums(comm2)
  if (is.null(depth)) depth <- min(libs)
  stopifnot(depth > 0)

  grp         <- interaction(meta2[, group_cols], drop = TRUE, sep = sep)
  n_per_group <- setNames(as.integer(table(grp)), names(table(grp)))

  one_iter <- function() {
    comm_r <- vegan::rrarefy(comm2, sample = depth)
    C <- rowsum(comm_r, grp) / as.vector(table(grp))
    if (ncol(C) > 1) C <- C[, colSums(C) > 0, drop = FALSE]
    list(C = C, D = as.matrix(vegan::vegdist(C, method = method)))
  }

  iters     <- purrr::map(seq_len(iterations), ~ one_iter())
  all_names <- Reduce(union, lapply(iters, function(x) rownames(x$D)))

  align_D <- function(M) {
    M2 <- matrix(NA_real_, length(all_names), length(all_names),
                 dimnames = list(all_names, all_names))
    M2[rownames(M), colnames(M)] <- M; M2
  }
  Ds    <- lapply(iters, function(x) align_D(x$D))
  num   <- Reduce(`+`, lapply(Ds, function(x) { x[is.na(x)] <- 0; x }))
  denom <- Reduce(`+`, lapply(Ds, function(x) !is.na(x)))
  Davg  <- num / denom

  out <- list(dist_mat = Davg, depth = depth, labels = all_names,
              sep = sep, n_per_group = n_per_group)

  if (return_centroids) {
    all_taxa <- Reduce(union, lapply(iters, function(x) colnames(x$C)))
    align_C  <- function(C) {
      C2 <- matrix(NA_real_, length(all_names), length(all_taxa),
                   dimnames = list(all_names, all_taxa))
      C2[rownames(C), colnames(C)] <- C; C2
    }
    Cs     <- lapply(iters, function(x) align_C(x$C))
    Cnum   <- Reduce(`+`, lapply(Cs, function(x) { x[is.na(x)] <- 0; x }))
    Cdenom <- Reduce(`+`, lapply(Cs, function(x) !is.na(x)))
    out$centroids_avg <- Cnum / Cdenom
    out$taxa          <- colnames(out$centroids_avg)
  }
  out
}

#' Extract paired centroid distances relative to LOMA from a centroid_dm_rarefy result.
extract_loma_pairs <- function(res,
                               group_cols = c("reserve_code", "s_code", "i_code", "Bag_type"),
                               loma       = "LOMA",
                               sep        = " | ") {
  M    <- res$dist_mat
  labs <- rownames(M)
  stopifnot(!is.null(labs), identical(labs, colnames(M)))

  parts         <- stringr::str_split_fixed(labs, stringr::fixed(sep), n = length(group_cols))
  colnames(parts) <- group_cols
  idx           <- tibble::tibble(row = labs, !!!as.data.frame(parts))

  reserve_col   <- group_cols[1]
  scenario_cols <- group_cols[-1]
  away          <- dplyr::filter(idx, .data[[reserve_col]] != loma)
  home          <- dplyr::filter(idx, .data[[reserve_col]] == loma)
  paired        <- dplyr::inner_join(away, home, by = scenario_cols, suffix = c("_R", "_LOMA"))

  tibble::tibble(
    reserve   = paired[[paste0(reserve_col, "_R")]],
    !!!rlang::set_names(paired[scenario_cols], scenario_cols),
    comp_dist = M[cbind(paired$row_R, paired$row_LOMA)]
  )
}

# -----------------------------------------------------------------------------
# Convergence path machinery
# -----------------------------------------------------------------------------

#' Build lookup: for each (pair, i_r_l), store a list of sample IDs.
build_pair_lookup <- function(meta_ext) {
  meta_ext |>
    dplyr::group_by(pair, i_r_l) |>
    dplyr::summarise(ids = list(original_label), .groups = "drop")
}

#' Safely retrieve sample IDs for (pair, i_r_l) from the lookup table.
get_ids <- function(lookup, p, lab) {
  idx <- which(lookup$pair == p & lookup$i_r_l == lab)
  if (length(idx) == 0) character(0) else unlist(lookup$ids[idx], use.names = FALSE)
}

#' Compute convergence metrics for one named path across all pairs.
#' @param D_full        Full square distance matrix.
#' @param lookup        Output of build_pair_lookup().
#' @param path_name     Label string for this path.
#' @param path_steps    c(S0_label, S1_label, R_label).
#' @param R0_label      i_r_l label for the native-away reference.
#' @param zero_R_within Treat within-R distance as 0?
#' @param zero_R0_within Treat within-R0 distance as 0?
compute_path <- function(D_full, lookup, path_name, path_steps,
                         R0_label       = "away|away & non-LOMA",
                         zero_R_within  = FALSE,
                         zero_R0_within = FALSE) {
  stopifnot(is.matrix(D_full), isTRUE(all(rownames(D_full) == colnames(D_full))))
  stopifnot(length(path_steps) == 3)

  S0 <- path_steps[1]; S1 <- path_steps[2]; R <- path_steps[3]

  purrr::map_dfr(sort(unique(lookup$pair)), function(p) {
    ids_S0 <- get_ids(lookup, p, S0)
    ids_S1 <- get_ids(lookup, p, S1)
    ids_R  <- get_ids(lookup, p, R)
    ids_R0 <- get_ids(lookup, p, R0_label)
    if (any(c(length(ids_S0), length(ids_S1), length(ids_R), length(ids_R0)) == 0))
      return(tibble::tibble())

    dS0_R  <- mean_cross(D_full, ids_S0, ids_R)
    dS1_R  <- mean_cross(D_full, ids_S1, ids_R)
    dR_R   <- if (zero_R_within) 0 else mean_within(D_full, ids_R)
    D_S0S1 <- mean_cross(D_full, ids_S0, ids_S1)
    D_S1R  <- mean_cross(D_full, ids_S1, ids_R)
    delta1 <- dS0_R - dS1_R
    delta2 <- dS1_R - dR_R
    eta1   <- if (is.na(D_S0S1) || D_S0S1 == 0) NA_real_ else delta1 / D_S0S1
    eta2   <- if (is.na(D_S1R)  || D_S1R  == 0) NA_real_ else delta2 / D_S1R
    F_s1   <- if (is.na(dS0_R)  || dS0_R  == 0) NA_real_ else 1 - (dS1_R / dS0_R)

    dS0_R0 <- mean_cross(D_full, ids_S0, ids_R0)
    dS1_R0 <- mean_cross(D_full, ids_S1, ids_R0)
    dR_R0  <- mean_cross(D_full, ids_R,  ids_R0)
    dR0_R0 <- if (zero_R0_within) 0 else mean_within(D_full, ids_R0)
    D_S1R0 <- mean_cross(D_full, ids_S1, ids_R0)
    delta1_R0 <- dS0_R0 - dS1_R0
    delta2_R0 <- dS1_R0 - dR_R0
    total_delta_R0 <- dS0_R0 - dR_R0
    eta1_R0 <- if (is.na(D_S0S1) || D_S0S1 == 0) NA_real_ else delta1_R0 / D_S0S1
    eta2_R0 <- if (is.na(D_S1R0) || D_S1R0 == 0) NA_real_ else delta2_R0 / D_S1R0
    F_s1_R0 <- if (is.na(dS0_R0) || dS0_R0 == 0) NA_real_ else 1 - (dS1_R0 / dS0_R0)
    native_shift_frac <- if (is.na(dS0_R0) || dS0_R0 == 0) NA_real_ else 1 - (dR_R0 / dS0_R0)

    tibble::tibble(
      pair = p, path = path_name,
      S0 = S0, S1 = S1, R = R, R0 = R0_label,
      n_S0 = length(ids_S0), n_S1 = length(ids_S1),
      n_R  = length(ids_R),  n_R0 = length(ids_R0),
      d_to_R_S0 = dS0_R, d_to_R_S1 = dS1_R, d_to_R_R = dR_R,
      delta1 = delta1, delta2 = delta2, total_delta = dS0_R - dR_R,
      D_S0_S1 = D_S0S1, D_S1_R = D_S1R,
      eta1 = eta1, eta2 = eta2, F_s1 = F_s1,
      d_to_R0_S0 = dS0_R0, d_to_R0_S1 = dS1_R0, d_to_R0_R = dR_R0, d_to_R0_R0 = dR0_R0,
      delta1_R0 = delta1_R0, delta2_R0 = delta2_R0, total_delta_R0 = total_delta_R0,
      D_S1_R0 = D_S1R0, eta1_R0 = eta1_R0, eta2_R0 = eta2_R0, F_s1_R0 = F_s1_R0,
      dR_R0 = dR_R0, native_shift_frac = native_shift_frac
    )
  }) |> dplyr::arrange(pair, path)
}

#' Run all 5 convergence paths and bind results.
compute_all_paths <- function(D_full, meta_ext,
                              R0_label       = "away|away & non-LOMA",
                              zero_R_within  = FALSE,
                              zero_R0_within = FALSE) {
  lookup <- build_pair_lookup(meta_ext)
  dplyr::bind_rows(
    compute_path(D_full, lookup, "LOMA_clim_eff", paths$LOMA_clim_eff,
                 R0_label = R0_label,         zero_R_within = zero_R_within, zero_R0_within = zero_R0_within),
    compute_path(D_full, lookup, "LOMA_litt_eff", paths$LOMA_litt_eff,
                 R0_label = R0_label,         zero_R_within = zero_R_within, zero_R0_within = zero_R0_within),
    compute_path(D_full, lookup, "LOMA_int_eff",  paths$LOMA_int_eff,
                 R0_label = R0_label,         zero_R_within = zero_R_within, zero_R0_within = zero_R0_within),
    compute_path(D_full, lookup, "away_clim_eff", paths$away_clim_eff,
                 R0_label = "LOMA|LOMA & LOMA", zero_R_within = zero_R_within, zero_R0_within = zero_R0_within),
    compute_path(D_full, lookup, "away_litt_eff", paths$away_litt_eff,
                 R0_label = "LOMA|LOMA & LOMA", zero_R_within = zero_R_within, zero_R0_within = zero_R0_within)
  ) |> dplyr::arrange(pair)
}

#' Report which (pair, label) combinations are missing from the lookup table.
check_paths <- function(lookup, paths) {
  tibble::tibble(pair = unique(lookup$pair)) |>
    tidyr::crossing(label = unlist(paths, use.names = FALSE)) |>
    dplyr::mutate(has = purrr::pmap_lgl(list(pair, label),
                                        ~ any(lookup$pair == ..1 & lookup$i_r_l == ..2))) |>
    dplyr::filter(!has) |>
    dplyr::arrange(pair, label)
}

#' Compute pairwise distances for all treatment comparisons within each pair context.
calculate_comp_dists <- function(distance_matrix, membership,
                                 pairs    = NULL,
                                 comp_ids = NULL,
                                 min_n    = 1) {
  D_full <- as.matrix(distance_matrix)
  stopifnot(is.matrix(D_full), is.character(rownames(D_full)))

  mem <- membership |>
    filter(grp %in% c("A", "B")) |>
    mutate(grp = factor(grp, levels = c("A", "B")))
  if (!is.null(pairs))    mem <- mem |> filter(pair    %in% pairs)
  if (!is.null(comp_ids)) mem <- mem |> filter(comp_id %in% comp_ids)

  comp_lookup <- mem |>
    distinct(pair, comp_id, grp, i_r_l) |>
    pivot_wider(names_from = grp, values_from = i_r_l, names_prefix = "t") |>
    rename(t1 = tA, t2 = tB)

  keys <- mem |> distinct(pair, comp_id) |> arrange(pair, comp_id)

  pmap_dfr(keys, function(pair, comp_id) {
    subm <- mem |> filter(pair == !!pair, comp_id == !!comp_id)
    ids  <- intersect(subm$original_label, rownames(D_full))
    if (length(ids) < 2) return(tibble())

    Dp  <- D_full[ids, ids, drop = FALSE]
    grp <- factor(subm$grp[match(ids, subm$original_label)], levels = c("A", "B"))
    nA  <- sum(grp == "A"); nB <- sum(grp == "B")
    if (nA < min_n || nB < min_n) return(tibble())

    res   <- avg_block_dists_plus(Dp, grp)
    t_row <- comp_lookup |> filter(pair == !!pair, comp_id == !!comp_id)
    tibble(
      pair, comp_id,
      t1 = if (nrow(t_row) > 0) t_row$t1[1] else NA_character_,
      t2 = if (nrow(t_row) > 0) t_row$t2[1] else NA_character_,
      n_A = as.integer(res$n["A"]), n_B = as.integer(res$n["B"]),
      mean_between   = res$between["A", "B"],
      mean_within_A  = res$between["A", "A"],
      mean_within_B  = res$between["B", "B"],
      pairs_between  = as.integer(res$pairs_between["A", "B"]),
      pairs_within_A = as.integer(res$pairs_within["A"]),
      pairs_within_B = as.integer(res$pairs_within["B"])
    )
  }) |> arrange(pair, comp_id)
}
