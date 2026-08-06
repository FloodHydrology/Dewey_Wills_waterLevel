#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: Dewey Exploratory Analysis
#Date: 8/6/2026
#Coder: Nate Jones (cnjones7@ua.edu)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Three inundation metrics (hydroperiod, water-table recession rate, and
# inundation event count) computed per well over a common 74-day window,
# then compared across soil type, all-species mortality, oak mortality, and
# plot elevation. Figures assembled with patchwork: one figure per grouping
# variable, stacking the three metrics in a column. Export block at the end
# sizes figures for PowerPoint (centered, with margins -- not slide-filling).
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setup workspace ---------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Clear memory
remove(list = ls())

#Load packages
library(tidyverse)
library(patchwork)   #combine ggplots

#Load data
df    <- read_csv("data/dewey_wills_wL_long.csv")
sites <- read_csv("data/site_info.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Global settings ---------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Single threshold used by every metric (cm, relative to ground surface).
#-30 cm ~ effective rooting depth for the oak survivorship question.
threshold_cm <- -30

#Common window: every site has data 2025-09-23 to 2025-12-06 (74 days).
window_start <- as.POSIXct("2025-09-23 07:00:00", tz = "UTC")
window_end   <- as.POSIXct("2025-12-06 07:00:00", tz = "UTC")

#Sites to drop everywhere (broken wL math -> implausible constant saturation).
bad_sites <- c("1026", "1031", "1038")

#Mortality factor order.
mort_levels <- c("None", "Low", "High")

#Shared theme -- sized up slightly for projection legibility on slides.
theme_dewey <- theme_bw(base_size = 14) +
  theme(
    axis.title   = element_text(size = 16),
    axis.text    = element_text(size = 12),
    legend.title = element_text(size = 13),
    plot.title   = element_text(size = 15),
    plot.margin  = margin(8, 8, 8, 8)
  )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Prep data ---------------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Convert water level inches -> cm (keep raw wL), truncate to common window.
df <- df %>%
  mutate(wL_cm = wL * 2.54) %>%
  filter(datetime >= window_start, datetime <= window_end)

#Sanity check: confirm the window filter caught the intended span.
df %>% summarise(min_dt = min(datetime), max_dt = max(datetime)) %>% print()

#Daily mean water level per site (used by recession and event metrics).
daily_wL <- df %>%
  mutate(date = as.Date(datetime)) %>%
  group_by(site, date) %>%
  summarise(wL_cm_daily = mean(wL_cm, na.rm = TRUE), .groups = "drop")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Metric 1: hydroperiod (days above threshold) ---------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
site_hydro <- df %>%
  group_by(site) %>%
  summarise(
    n_obs       = sum(!is.na(wL_cm)),
    hydroperiod = sum(wL_cm > threshold_cm, na.rm = TRUE) / 24,   #hours -> days
    .groups = "drop"
  )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Metric 2: water-table recession rate (cm/day) --------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Day-to-day decline off the daily-mean series (daily averaging removes the
#hourly sensor jitter). Falling days only; median per site. Values negative.
site_rec <- daily_wL %>%
  arrange(site, date) %>%
  group_by(site) %>%
  mutate(
    dwL_cm = wL_cm_daily - lag(wL_cm_daily),
    dt_day = as.numeric(date - lag(date), units = "days"),
    rate   = dwL_cm / dt_day
  ) %>%
  filter(rate < 0) %>%
  summarise(
    n_fall         = n(),
    recession_rate = median(rate, na.rm = TRUE),
    .groups = "drop"
  )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Metric 3: inundation event count ---------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#An event = daily mean rises above threshold after being at/below it. Counted
#on the daily series to avoid noise chatter. days_wet separates flashy (many
#events) from one long flood.
site_evt <- daily_wL %>%
  arrange(site, date) %>%
  group_by(site) %>%
  mutate(
    inundated = wL_cm_daily > threshold_cm,
    onset     = inundated & !lag(inundated, default = FALSE)
  ) %>%
  summarise(
    n_events = sum(onset, na.rm = TRUE),
    days_wet = sum(inundated, na.rm = TRUE),
    .groups = "drop"
  )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Combine metrics + metadata into one site-level table -------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#One filtered site set for ALL metrics.
site_metrics <- site_hydro %>%
  left_join(site_rec, by = "site") %>%
  left_join(site_evt, by = "site") %>%
  filter(!site %in% bad_sites) %>%
  mutate(point_id = parse_number(site)) %>%
  left_join(sites, by = "point_id") %>%
  mutate(
    Mortality_all = factor(Mortality_all, levels = mort_levels),
    Mortality_oak = factor(Mortality_oak, levels = mort_levels)
  )

#Quick look
site_metrics %>%
  select(site, point_id, soil_type, Mortality_all, Mortality_oak, m,
         hydroperiod, recession_rate, n_events, days_wet) %>%
  print(n = Inf)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Helper functions --------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Boxplot for a categorical grouping.
box_panel <- function(data, metric, group, xlab, ylab) {
  data %>%
    filter(!is.na(.data[[group]])) %>%
    ggplot(aes(x = .data[[group]], y = .data[[metric]])) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, colour = "steelblue", size = 2, alpha = 0.7) +
    labs(x = xlab, y = ylab) +
    theme_dewey
}

#Scatter vs elevation with lm fit + Spearman rho in the subtitle.
elev_panel <- function(data, metric, ylab) {
  ct  <- cor.test(data$m, data[[metric]], method = "spearman")
  sub <- sprintf("Spearman rho = %.2f, p = %.3f", ct$estimate, ct$p.value)
  
  data %>%
    ggplot(aes(x = m, y = .data[[metric]])) +
    geom_smooth(method = "lm", se = TRUE, colour = "steelblue",
                fill = "grey80", lwd = 0.75) +
    geom_point(colour = "steelblue", size = 2, alpha = 0.7) +
    labs(x = "Plot Elevation (m)", y = ylab, subtitle = sub) +
    theme_dewey
}

#Axis labels.
lab_hydro <- "Hydroperiod\n[days]"
lab_rec   <- "Recession rate\n[cm/day]"
lab_evt   <- "Inundation events\n[count]"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Figure 1: all three metrics vs SOIL TYPE -------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fig_soil <-
  box_panel(site_metrics, "hydroperiod",    "soil_type", NULL, lab_hydro) /
  box_panel(site_metrics, "recession_rate", "soil_type", NULL, lab_rec)   /
  box_panel(site_metrics, "n_events",       "soil_type", "Soil type", lab_evt) +
  plot_annotation(title = "Inundation metrics by soil type")

fig_soil

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Figure 2: all three metrics vs ALL-SPECIES MORTALITY -------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fig_mort <-
  box_panel(site_metrics, "hydroperiod",    "Mortality_all", NULL, lab_hydro) /
  box_panel(site_metrics, "recession_rate", "Mortality_all", NULL, lab_rec)   /
  box_panel(site_metrics, "n_events",       "Mortality_all", "All-species mortality", lab_evt) +
  plot_annotation(title = "Inundation metrics by all-species mortality")

fig_mort

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Figure 3: all three metrics vs OAK MORTALITY --------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fig_oak <-
  box_panel(site_metrics, "hydroperiod",    "Mortality_oak", NULL, lab_hydro) /
  box_panel(site_metrics, "recession_rate", "Mortality_oak", NULL, lab_rec)   /
  box_panel(site_metrics, "n_events",       "Mortality_oak", "Oak mortality", lab_evt) +
  plot_annotation(title = "Inundation metrics by oak mortality")

fig_oak

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Figure 4: all three metrics vs PLOT ELEVATION -------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fig_elev <-
  elev_panel(site_metrics, "hydroperiod",    lab_hydro) /
  elev_panel(site_metrics, "recession_rate", lab_rec)   /
  elev_panel(site_metrics, "n_events",       lab_evt) +
  plot_annotation(title = "Inundation metrics vs plot elevation")

fig_elev

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Figure 5: inundation regime scatter -----------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#days_wet vs n_events, coloured by mortality. Bottom-left = dry;
#bottom-right = one sustained flood; top = flashy.
fig_regime <- site_metrics %>%
  filter(!is.na(Mortality_all)) %>%
  ggplot(aes(x = days_wet, y = n_events, colour = Mortality_all)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_colour_manual(values = c("None" = "grey60",
                                 "Low"  = "steelblue",
                                 "High" = "firebrick")) +
  labs(x = "Days wet [days]", y = "Inundation events\n[count]",
       colour = "Mortality",
       title  = "Inundation regime: flashy vs sustained vs dry") +
  theme_dewey

fig_regime

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Extra scatters ---------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Days wet vs recession rate
fig_wet_rec <- site_metrics %>%
  filter(!is.na(Mortality_all)) %>%
  ggplot(aes(x = days_wet, y = recession_rate, colour = Mortality_all)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_colour_manual(values = c("None" = "grey60",
                                 "Low"  = "steelblue",
                                 "High" = "firebrick")) +
  labs(x = "Days wet [days]", y = "Recession rate\n[cm/day]",
       colour = "Mortality") +
  theme_dewey

fig_wet_rec

#Inundation events vs recession rate
fig_evt_rec <- site_metrics %>%
  filter(!is.na(Mortality_all)) %>%
  ggplot(aes(x = n_events, y = recession_rate, colour = Mortality_all)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_colour_manual(values = c("None" = "grey60",
                                 "Low"  = "steelblue",
                                 "High" = "firebrick")) +
  labs(x = "Inundation events [count]", y = "Recession rate\n[cm/day]",
       colour = "Mortality") +
  theme_dewey

fig_evt_rec

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Export figures for PowerPoint -------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Slides are 13.33 x 7.5 in (16:9). Figures are sized to sit centered with
#margins -- NOT to fill the slide. 300 dpi keeps them crisp when scaled.
#bg = "white" prevents transparent PNGs that look broken on a slide master.

#Stacked columns (3 metrics tall) -> tall/narrow
ggsave("output/fig_soil.png", fig_soil, width = 5, height = 7, dpi = 300, bg = "white")
ggsave("output/fig_mort.png", fig_mort, width = 5, height = 7, dpi = 300, bg = "white")
ggsave("output/fig_oak.png",  fig_oak,  width = 5, height = 7, dpi = 300, bg = "white")
ggsave("output/fig_elev.png", fig_elev, width = 5, height = 7, dpi = 300, bg = "white")

#Single scatters -> landscape, leaves slide margins
ggsave("output/fig_regime.png",  fig_regime,  width = 7, height = 5, dpi = 300, bg = "white")
ggsave("output/fig_wet_rec.png", fig_wet_rec, width = 7, height = 5, dpi = 300, bg = "white")
ggsave("output/fig_evt_rec.png", fig_evt_rec, width = 7, height = 5, dpi = 300, bg = "white")
