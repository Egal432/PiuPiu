# =============================================================================
# Flight initiation distance (FID) of dune birds - simplified analysis
# -----------------------------------------------------------------------------
# Data      : Birds_of_the_dune.xlsx (sheets "Behaviour" and "Code")
# Response  : FID = Distance at which birds took flight (m)
# Predictors: number of people, flock size (number of birds), species
# Excluded  : all non-flight responses (W, L, S) and rows without a response;
#             Time_settle is not used.
# Terns     : Brandseeschwalbe (Ts), Flussseeschwalbe (Sh/"Fluss") and the
#             unresolved "Schwalbe" are pooled into one group "Schwalbe".
#
# THE KEY PROBLEM: species and flock size are confounded (terns and plovers
# occur in big flocks, pipits and gulls mostly alone). A pooled model therefore
# mixes (i) the effect of flock size WITHIN a species with (ii) plain
# differences BETWEEN species. We separate the two:
#
#   MODEL A  pooled, all flights incl. mixed flocks, no species term   (n = 72)
#            -> overall ASSOCIATION across all birds; may be confounded
#   MODEL B  species as FIXED effect, flights with known species        (n = 70)
#            -> effect of flock size / people WITHIN species (main test)
#   MODEL C  within- vs between-species decomposition of flock size     (n = 70)
#            -> are the within- and between-species slopes similar?
# Mixed flocks ("Flock", species unknown) can only be used in model A.
# =============================================================================


# 0. Packages & settings ------------------------------------------------------
pkgs <- c("tidyverse", "readxl", "lme4", "lmerTest", "performance",
          "patchwork", "broom", "broom.mixed", "lmtest", "sandwich", "ggtext")
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

data_file     <- "Birds_of_the_dune.xlsx"   # path to your Excel file
out_dir       <- "fid_simple_output"        # figures + tables are saved here
min_n_species <- 5                          # species with >= this many flights get
# species-level tests
dir.create(out_dir, showWarnings = FALSE)
theme_set(theme_minimal(base_size = 12))
save_fig <- function(p, name, w = 8, h = 5)
  ggsave(file.path(out_dir, name), p, width = w, height = h, dpi = 300, bg = "white")


# 1. Import & data preparation ------------------------------------------------
raw <- read_excel(data_file, sheet = "Behaviour") |>
  mutate(row_id = row_number())

codes <- read_excel(data_file, sheet = "Code")
names(codes) <- trimws(names(codes))              # header is "Code " (trailing space)
# Latin names come from the Code sheet. Four spellings in that sheet are wrong and
# are corrected here (please fix them in the Excel file too).
latin_fix <- c("Coleus monedula"     = "Coloeus monedula",
               "Soateria mollissima" = "Somateria mollissima",
               "Triga totanus"       = "Tringa totanus",
               "Anas platyrhychos"   = "Anas platyrhynchos")
codes <- transmute(codes, code = trimws(Code), name_ger = Name_ger,
                   name_sc = dplyr::recode(trimws(Name_sc), !!!latin_fix))
latin <- deframe(select(codes, code, name_sc))

# "Fluss" = Flussseeschwalbe (code Sh). Ts + Sh + "Schwalbe" -> one group.
# label    = plain text (tables); label_md = with italics for the plots (ggtext)
raw <- raw |>
  mutate(code    = dplyr::recode(str_trim(Code), Fluss = "Sh"),
         species = if_else(code %in% c("Ts", "Sh", "Schwalbe"), "Schwalbe", code),
         label   = case_when(species == "Schwalbe" ~ "Terns (S. hirundo + T. sandvicensis)",
                             species == "Flock"    ~ "Mixed flock",
                             TRUE ~ coalesce(unname(latin[species]), species)),
         label_md = case_when(species == "Schwalbe" ~ "Terns (*S. hirundo* + *T. sandvicensis*)",
                              species == "Flock"    ~ "Mixed flock",
                              TRUE ~ paste0("*", label, "*")))

# Rough "approach event" ID (only used for one robustness check): consecutive
# rows with identical distance and number of people are assumed to belong to the
# same approach. The data contain no event ID, so this is a heuristic.
raw <- raw |>
  mutate(event = cumsum(c(TRUE, diff(Distance) != 0 | diff(Num_people) != 0)))

cat("\n--- Responses in the raw data (NA = not recorded) ---\n")
print(table(raw$Response, useNA = "ifany"))

# All flights, INCLUDING mixed flocks  -> model A, descriptives, figures
fid_all <- raw |>
  filter(Response == "F") |>
  transmute(row_id, event, species, label, label_md,
            FID    = Distance,       # m
            people = Num_people,
            flock  = Num_birds)
# Flights with a known species         -> models B and C, species tests
fid_sp <- filter(fid_all, species != "Flock")

cat("\nFlights used: ", nrow(fid_all), " (", nrow(fid_sp), " with known species + ",
    nrow(fid_all) - nrow(fid_sp), " mixed-flock flights); species groups: ",
    n_distinct(fid_sp$species), "\n", sep = "")


# 2. Descriptive statistics ---------------------------------------------------
desc_overall <- fid_all |>
  summarise(n = n(), median_FID = median(FID), q25 = quantile(FID, .25),
            q75 = quantile(FID, .75), mean_FID = mean(FID), sd_FID = sd(FID),
            min_FID = min(FID), max_FID = max(FID),
            median_people = median(people), max_people = max(people),
            median_flock = median(flock), max_flock = max(flock))
# n_flock_sizes shows how much flock size varies WITHIN each species: only
# species with several different flock sizes can inform the within-species test.
desc_species <- fid_all |>
  group_by(label) |>
  summarise(n = n(), median_FID = median(FID), min_FID = min(FID),
            max_FID = max(FID), median_flock = median(flock),
            min_flock = min(flock), max_flock = max(flock),
            n_flock_sizes = n_distinct(flock),
            median_people = median(people), .groups = "drop") |>
  arrange(desc(n))

cat("\n--- FID (m), all flights ---\n");  print(desc_overall, width = Inf)
cat("\n--- Per species ---\n");           print(desc_species, n = Inf, width = Inf)
cat("\nSpearman correlation people ~ flock size: ",
    round(cor(fid_all$people, fid_all$flock, method = "spearman"), 2), "\n", sep = "")
write_csv(desc_species, file.path(out_dir, "table_descriptives_species.csv"))


# 3. Figures (descriptive) ----------------------------------------------------

## Fig 1 - FID per species (mixed flocks shown as their own group)
n_lab <- fid_all |> count(label, label_md) |>
  mutate(label_n = paste0(label_md, " (n=", n, ")"))
fid_p <- fid_all |> left_join(n_lab, by = c("label", "label_md"))
p1 <- fid_p |>
  mutate(label_n = fct_reorder(label_n, FID, median)) |>
  ggplot(aes(label_n, FID)) +
  geom_jitter(width = .15, height = 0, alpha = .6, size = 2) +
  stat_summary(fun = median, geom = "crossbar", width = .6,
               colour = "#d95f02", linewidth = .4) +
  scale_y_log10() + coord_flip() +
  labs(x = NULL, y = "Flight initiation distance (m, log scale)",
       title = "Flight initiation distance per species",
       subtitle = "Orange bar = median") +
  theme(axis.text.y = element_markdown())
save_fig(p1, "fig1_fid_per_species.png", 8, 5.5)

## Fig 2 - pooled picture: FID vs flock size (a) and vs people (b), all flights.
##         This ignores species -> shows the OVERALL association only.
p2a <- ggplot(fid_all, aes(flock, FID)) +
  geom_point(alpha = .6) +
  geom_smooth(method = "lm", colour = "#1b9e77", fill = "#1b9e77", alpha = .2) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Flock size (log scale)", y = "FID (m, log scale)", title = "a")
p2b <- ggplot(fid_all, aes(people, FID)) +
  geom_point(alpha = .6) +
  geom_smooth(method = "lm", colour = "#d95f02", fill = "#d95f02", alpha = .2) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Number of people (log scale)", y = NULL, title = "b")
save_fig(p2a | p2b, "fig2_fid_vs_flock_and_people_pooled.png", 10, 4.5)

## Fig 3 - flock size per species (shows the confounding with species)
p3 <- fid_p |>
  mutate(label_n = fct_reorder(label_n, flock, median)) |>
  ggplot(aes(label_n, flock)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90") +
  geom_jitter(width = .15, height = 0, alpha = .6) +
  scale_y_log10() + coord_flip() +
  labs(x = NULL, y = "Flock size (log scale)",
       title = "Flock size differs between species") +
  theme(axis.text.y = element_markdown())
save_fig(p3, "fig3_flock_size_per_species.png", 8, 5.5)

## Fig 4 - FID vs flock size, species by species. Coloured lines = slope WITHIN
##         each species (species with >= min_n_species flights and >= 3 flock
##         sizes); dashed black line = pooled slope across all species.
slope_sp <- fid_sp |> group_by(species) |>
  filter(n() >= min_n_species, n_distinct(flock) >= 3) |>
  pull(species) |> unique()
fid_col <- fid_sp |>
  mutate(group = if_else(species %in% slope_sp, label_md, "Other species"))
p4 <- ggplot(fid_col, aes(flock, FID, colour = group)) +
  geom_point(alpha = .7, size = 2) +
  geom_smooth(data = filter(fid_col, group != "Other species"),
              method = "lm", se = FALSE, linewidth = 1) +
  geom_smooth(data = fid_sp, aes(flock, FID), inherit.aes = FALSE,
              method = "lm", se = FALSE, colour = "black", linetype = 2) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = c(setNames(RColorBrewer::brewer.pal(max(3, length(slope_sp)), "Dark2")[seq_along(slope_sp)],
                                          sort(unique(filter(fid_col, group != "Other species")$group))),
                                 "Other species" = "grey60")) +
  labs(x = "Flock size (log scale)", y = "FID (m, log scale)", colour = NULL,
       title = "Does flock size matter within species?",
       subtitle = "Coloured lines: within-species slopes. Dashed black: pooled slope across all species") +
  theme(legend.text = element_markdown())
save_fig(p4, "fig4_flock_effect_within_species.png", 9, 5.5)

## Fig 5 - within-species relationships: every value is expressed as the
##         deviation from its species' own mean (so species differences are
##         removed). Species with a single flight are dropped here.
fid_c <- fid_sp |>
  group_by(species) |> filter(n() >= 2) |>
  mutate(fid_c    = log(FID)    - mean(log(FID)),
         flock_c  = log(flock)  - mean(log(flock)),
         people_c = log(people) - mean(log(people))) |>
  ungroup()
p5a <- ggplot(fid_c, aes(flock_c, fid_c)) +
  geom_hline(yintercept = 0, colour = "grey70") + geom_vline(xintercept = 0, colour = "grey70") +
  geom_point(alpha = .6) +
  geom_smooth(method = "lm", colour = "#1b9e77", fill = "#1b9e77", alpha = .2) +
  labs(x = "log flock size, minus species mean", y = "log FID, minus species mean", title = "a  Flock size")
p5b <- ggplot(fid_c, aes(people_c, fid_c)) +
  geom_hline(yintercept = 0, colour = "grey70") + geom_vline(xintercept = 0, colour = "grey70") +
  geom_point(alpha = .6) +
  geom_smooth(method = "lm", colour = "#d95f02", fill = "#d95f02", alpha = .2) +
  labs(x = "log people, minus species mean", y = NULL, title = "b  Number of people")
save_fig(p5a | p5b, "fig5_within_species_effects.png", 10, 4.5)


# 4. MODEL A: pooled, all birds, no species term ------------------------------
# Linear model on log(FID). FID is positive and right-skewed; predictors are
# log-transformed too, so a coefficient b means: doubling the predictor
# multiplies FID by 2^b. All 72 flights, mixed flocks included.
# CAUTION: this mixes within- and between-species differences (see header).
m_all <- lm(log(FID) ~ log(people) + log(flock), data = fid_all)

tidy_effects <- function(m) {
  ci <- confint(m)
  as_tibble(summary(m)$coefficients, rownames = "term") |>
    filter(term %in% c("log(people)", "log(flock)")) |>
    mutate(conf.low = ci[term, 1], conf.high = ci[term, 2],
           fold_per_doubling = 2^Estimate,
           fold_low = 2^conf.low, fold_high = 2^conf.high)
}
cat("\n=== MODEL A: pooled, all birds (n = ", nrow(fid_all), ") ===\n", sep = "")
print(summary(m_all))
print(tidy_effects(m_all), width = Inf)
write_csv(tidy_effects(m_all), file.path(out_dir, "table_model_A_pooled.csv"))


# 5. MODEL B: species as fixed effect (MAIN TEST) -----------------------------
# Every species gets its own intercept, so the flock-size and people effects are
# estimated only from differences WITHIN species. Species with a single flight
# or a single flock size simply contribute no information to the slopes.
m_B <- lm(log(FID) ~ species + log(people) + log(flock), data = fid_sp)

cat("\n=== MODEL B: species fixed (n = ", nrow(fid_sp), ") ===\n", sep = "")
effects_B <- tidy_effects(m_B)
print(effects_B, width = Inf)
cat("\nResidual df:", df.residual(m_B), "\n")
write_csv(effects_B, file.path(out_dir, "table_model_B_species_fixed.csv"))
cat("\n--- Collinearity ---\n"); print(check_collinearity(m_B))

cat("\n--- Does adding species improve on the pooled model (F-test)? ---\n")
print(anova(update(m_B, . ~ . - species), m_B))

cat("\n--- people x flock interaction (exploratory) ---\n")
print(summary(update(m_B, . ~ . + log(people):log(flock)))$coefficients[
  c("log(people)", "log(flock)", "log(people):log(flock)"), ])

cd <- cooks.distance(m_B)
cat("\n--- Largest Cook's distances in model B (< 0.5 = no severe outlier) ---\n")
print(fid_sp |> mutate(cooks_d = cd) |> arrange(desc(cooks_d)) |> slice_head(n = 3) |>
        select(label, FID, people, flock, cooks_d), width = Inf)


# 6. MODEL C: within- vs between-species effect of flock size -----------------
# log(flock) is split into the species mean ("between") and the deviation from
# that mean ("within"). If species were confounding the flock effect, the two
# slopes would differ. A random intercept for species is included.
fid_wb <- fid_sp |>
  group_by(species) |>
  mutate(lf = log(flock), lf_between = mean(lf), lf_within = lf - lf_between) |>
  ungroup()
m_C <- lmer(log(FID) ~ log(people) + lf_within + lf_between + (1 | species), data = fid_wb)
cat("\n=== MODEL C: within vs between species flock effect ===\n")
print(summary(m_C)$coefficients)
cat("Species variance singular (TRUE = ~0)? ", isSingular(m_C), "\n")

# Direct test of the difference: the 'lf_mean' coefficient below is the
# difference between the between- and within-species slopes.
m_C2 <- lmer(log(FID) ~ log(people) + lf + lf_mean + (1 | species),
             data = fid_wb |> group_by(species) |> mutate(lf_mean = mean(lf)) |> ungroup())
cat("\n--- Difference between-species minus within-species slope ('lf_mean') ---\n")
print(summary(m_C2)$coefficients["lf_mean", , drop = FALSE])
write_csv(as_tibble(summary(m_C)$coefficients, rownames = "term"),
          file.path(out_dir, "table_model_C_within_between.csv"))


# 7. Do species differ, and is the flock effect the same in every species? ----
fid5 <- fid_sp |> group_by(species) |> filter(n() >= min_n_species) |> ungroup()
cat("\nSpecies with >=", min_n_species, "flights:",
    paste(unique(fid5$label), collapse = ", "), "(n =", nrow(fid5), "flights)\n")

cat("\n--- (a) Raw: Kruskal-Wallis FID ~ species ---\n")
print(kruskal.test(FID ~ species, data = fid5))
cat("\n--- (b) Species effect after flock size and people ---\n")
print(anova(lm(log(FID) ~ log(people) + log(flock), data = fid5),
            lm(log(FID) ~ species + log(people) + log(flock), data = fid5)))
cat("\n--- (c) Does flock size differ between species? Kruskal-Wallis ---\n")
print(kruskal.test(flock ~ species, data = fid5))

cat("\n--- (d) Flock-size slope per species (>= ", min_n_species,
    " flights, >= 3 flock sizes) ---\n", sep = "")
slopes <- fid_sp |>
  group_by(label) |>
  filter(n() >= min_n_species, n_distinct(flock) >= 3) |>
  group_modify(~ tidy(lm(log(FID) ~ log(flock), data = .x), conf.int = TRUE) |>
                 filter(term == "log(flock)") |> mutate(n = nrow(.x))) |>
  ungroup() |> select(label, n, estimate, conf.low, conf.high, p.value)
print(slopes, width = Inf)
write_csv(slopes, file.path(out_dir, "table_flock_slopes_per_species.csv"))

cat("\n--- (e) Does the flock slope differ among species? (F-test, n>=5 species) ---\n")
print(anova(lm(log(FID) ~ species + log(people) + log(flock), data = fid5),
            lm(log(FID) ~ species * log(flock) + log(people), data = fid5)))


# 8. Robustness of the within-species flock effect (model B) ------------------
# (a) without the most extreme FID (250 m, "running people")
m_noout <- update(m_B, data = filter(fid_sp, FID < 200))
# (b) one record per approach event (keep the row with the largest flock)
fid_thin <- fid_sp |> group_by(event) |> slice_max(flock, n = 1, with_ties = FALSE) |> ungroup()
m_thin   <- update(m_B, data = fid_thin)
# (c) Gamma GLM (log link) instead of a log-normal model
m_gamma  <- glm(FID ~ species + log(people) + log(flock), family = Gamma(link = "log"),
                data = fid_sp, control = glm.control(maxit = 100))
# (d) mixed flocks added as their own group (n = 72)
m_B72    <- update(m_B, data = fid_all)

rob <- imap_dfr(list("Model B (main)"                   = m_B,
                     "Without FID >= 200 m"             = m_noout,
                     "One record per event"             = m_thin,
                     "Gamma GLM"                        = m_gamma,
                     "Mixed flocks as own group (n=72)" = m_B72),
                ~ tidy(.x, conf.int = TRUE) |>
                  filter(term %in% c("log(people)", "log(flock)")) |>
                  mutate(model = .y, n = nobs(.x)))

# (e) leave one species out (each species with >= min_n_species flights)
loo <- map_dfr(unique(fid5$species), function(sp) {
  m <- update(m_B, data = filter(fid_sp, species != sp))
  tidy(m, conf.int = TRUE) |> filter(term %in% c("log(people)", "log(flock)")) |>
    mutate(model = paste("Without", fid_sp$label[match(sp, fid_sp$species)]), n = nobs(m))
})
rob <- bind_rows(rob, loo)

cat("\n=== Robustness of the within-species effects (log scale) ===\n")
print(select(rob, model, n, term, estimate, conf.low, conf.high, p.value),
      n = Inf, width = Inf)
write_csv(rob, file.path(out_dir, "table_robustness_model_B.csv"))
if (!m_gamma$converged) cat("NOTE: Gamma GLM did not converge - treat that row with care.\n")


# 9. Diagnostics and heteroscedasticity (model B) -----------------------------
# Species with a single flight have a residual of exactly 0 by construction
# (their own intercept fits them perfectly). They are left out of the residual
# plots and the robust standard errors below; they carry no slope information.
diag_all <- fid_sp |>
  mutate(fitted = fitted(m_B), resid = resid(m_B), hat = hatvalues(m_B))
diag_df <- diag_all |>
  filter(hat < 0.999) |>
  mutate(std_resid = resid / (sigma(m_B) * sqrt(1 - hat)))
cat("\nFlights of single-flight species left out of residual plots: ",
    nrow(diag_all) - nrow(diag_df), "\n", sep = "")

cat("\n=== Heteroscedasticity tests ===\n")
cat("Breusch-Pagan, model A (pooled):\n");       print(bptest(m_all))
cat("Breusch-Pagan, model B (species fixed):\n"); print(bptest(m_B))
cat("Spearman |std. residual| ~ fitted: rho = ",
    round(cor(abs(diag_df$std_resid), diag_df$fitted, method = "spearman"), 2), ", p = ",
    signif(cor.test(abs(diag_df$std_resid), diag_df$fitted, method = "spearman",
                    exact = FALSE)$p.value, 2), "\n", sep = "")
cat("\nResidual SD by thirds of the fitted values (model B):\n")
print(diag_df |>
        mutate(fitted_third = cut(fitted, quantile(fitted, c(0, 1/3, 2/3, 1)),
                                  include.lowest = TRUE)) |>
        group_by(fitted_third) |>
        summarise(n = n(), residual_sd = sd(resid)))

# Heteroscedasticity-robust inference for model B. Even if the tests find
# nothing, these show whether the conclusions depend on equal variances.
fid_sp2  <- fid_sp |> group_by(species) |> filter(n() >= 2) |> ungroup()
m_B2     <- update(m_B, data = fid_sp2)          # same slopes as m_B
terms2   <- c("log(people)", "log(flock)")
hc3      <- coefci(m_B2, parm = terms2, vcov. = vcovHC(m_B2, type = "HC3"))
hc3_p    <- coeftest(m_B2, vcov. = vcovHC(m_B2, type = "HC3"))[terms2, 4]

set.seed(1)                                      # case bootstrap, 4000 resamples
boot <- replicate(4000, {
  b <- fid_sp[sample(nrow(fid_sp), replace = TRUE), ]
  coef(lm(log(FID) ~ species + log(people) + log(flock), data = b))[terms2]
})
boot_ci <- apply(boot, 1, quantile, c(.025, .975), na.rm = TRUE)
cl_ci   <- confint(m_B)[terms2, ]

robust_tbl <- bind_rows(
  tibble(method = "Classical (model B)",             term = terms2,
         estimate = coef(m_B)[terms2], conf.low = cl_ci[, 1], conf.high = cl_ci[, 2],
         p.value = summary(m_B)$coefficients[terms2, 4]),
  tibble(method = "HC3 robust SE",                   term = terms2,
         estimate = coef(m_B2)[terms2], conf.low = hc3[, 1], conf.high = hc3[, 2],
         p.value = hc3_p),
  tibble(method = "Case bootstrap (percentile CI)",  term = terms2,
         estimate = coef(m_B)[terms2], conf.low = boot_ci[1, ], conf.high = boot_ci[2, ],
         p.value = NA_real_))
cat("\n=== Model B: classical vs heteroscedasticity-robust inference ===\n")
print(robust_tbl, width = Inf)
write_csv(robust_tbl, file.path(out_dir, "table_robust_inference_model_B.csv"))

## Fig 6 - residual diagnostics (a: residuals vs fitted, b: scale-location, c: Q-Q)
d1 <- ggplot(diag_df, aes(fitted, resid)) +
  geom_hline(yintercept = 0, linetype = 2) + geom_point(alpha = .6) +
  geom_smooth(se = FALSE, colour = "#d95f02", method = "loess", formula = y ~ x) +
  labs(title = "a  Residuals vs fitted", x = "Fitted", y = "Residual")
d2 <- ggplot(diag_df, aes(fitted, sqrt(abs(std_resid)))) +
  geom_point(alpha = .6) +
  geom_smooth(se = FALSE, colour = "#d95f02", method = "lm", formula = y ~ x) +
  labs(title = "b  Scale-location (equal variance = flat line)",
       x = "Fitted", y = expression(sqrt("|standardised residual|")))
d3 <- ggplot(diag_df, aes(sample = resid)) + stat_qq() + stat_qq_line() +
  labs(title = "c  Normal Q-Q plot", x = "Theoretical quantiles", y = "Sample quantiles")
save_fig(d1 | d2 | d3, "fig6_diagnostics_model_B.png", 14, 4.5)

cat("\nDone. Output in: ", normalizePath(out_dir), "\n", sep = "")
cat("\n", R.version.string, "\n", sep = "")
for (p in c("lme4", "lmerTest", "ggplot2"))
  cat(p, as.character(packageVersion(p)), "\n")