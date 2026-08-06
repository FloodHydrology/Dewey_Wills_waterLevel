#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: Dewey Exploratory Analysis
#Date: 8/6/2026
#Coder: Nate Jones (cnjones7@ua.edu)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setup workspace ---------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Clear memory
remove(list=ls())

#load packages
library(tidyverse)
library(rstatix)     #pairwise_wilcox_test
library(rcompanion)  #cldList -> compact letter display

#load data 
df <- read_csv("data/dewey_wills_wL_long.csv")
sites <- read_csv("data/site_info.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Prep data ---------------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Define threshold (cm) -- surface inundation
threshold_cm <- -10

#Common window: all sites have data 2025-09-23 to 2025-12-06 (74 days)
window_start <- as.POSIXct("2025-09-23 07:00:00", tz = "UTC")
window_end   <- as.POSIXct("2025-12-06 07:00:00", tz = "UTC")

#Convert water level from inches to cm (keep raw wL for traceability)
#and truncate to the common window
df <- df %>%
  mutate(wL_cm = wL * 2.54) %>%
  filter(datetime >= window_start, datetime <= window_end)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Site-level hydroperiod above threshold ----------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Hourly readings -> days above threshold (hours above / 24)
site_duration <- df %>%
  group_by(site) %>%
  summarise(
    n_obs        = sum(!is.na(wL_cm)),
    hydroperiod  = sum(wL_cm > threshold_cm, na.rm = TRUE) / 24
  ) |> 
  filter(hydroperiod<60)

#Take a look
site_duration %>% arrange(desc(hydroperiod)) %>% print(n = Inf)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Join site metadata ------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Drop known-bad sites (broken wL math)
bad_sites <- c("1026", "1031", "1038")

site_duration <- site_duration %>%
  filter(!site %in% bad_sites) %>%
  mutate(point_id = parse_number(site)) %>%
  left_join(sites, by = "point_id")

#Check the join
site_duration %>% 
  select(site, point_id, soil_type, Mortality_all, m, hydroperiod) %>% 
  print(n = Inf)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Explore: hydroperiod by soil --------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
site_duration %>%
  ggplot(aes(x = soil_type, y = hydroperiod)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, colour = "steelblue", size = 2, alpha = 0.7) +
  labs(x = "Soil type", y = "Hydroperiod [days]") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 10)
  )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Pairwise comparison: hydroperiod by mortality --------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mort_dat <- site_duration %>%
  filter(!is.na(Mortality_all)) %>%
  mutate(Mortality_all = factor(Mortality_all, levels = c("None", "Low", "High")))

mort_dat %>% count(Mortality_all)

pw_mort <- mort_dat %>%
  wilcox_test(hydroperiod ~ Mortality_all, p.adjust.method = "BH")
pw_mort

pw_mort_letters <- pw_mort %>%
  mutate(comparison = paste(group1, group2, sep = " - ")) %>%
  select(comparison, p.adj)
cld_mort <- cldList(
  p.adj ~ comparison,
  data      = pw_mort_letters,
  threshold = 0.05
)
cld_mort

mort_dat %>%
  ggplot(aes(x = Mortality_all, y = hydroperiod)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, colour = "steelblue", size = 2, alpha = 0.7) +
  labs(x = "All-species mortality", y = "Hydroperiod [days]") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 10)
  )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Pairwise comparison: hydroperiod by oak mortality ----------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mort_oak_dat <- site_duration %>%
  filter(!is.na(Mortality_oak)) %>%
  mutate(Mortality_oak = factor(Mortality_oak, levels = c("None", "Low", "High")))

mort_oak_dat %>% count(Mortality_oak)

pw_mort_oak <- mort_oak_dat %>%
  wilcox_test(hydroperiod ~ Mortality_oak, p.adjust.method = "BH")
pw_mort_oak

pw_mort_oak_letters <- pw_mort_oak %>%
  mutate(comparison = paste(group1, group2, sep = " - ")) %>%
  select(comparison, p.adj)
cld_mort_oak <- cldList(
  p.adj ~ comparison,
  data      = pw_mort_oak_letters,
  threshold = 0.05
)
cld_mort_oak

mort_oak_dat %>%
  ggplot(aes(x = Mortality_oak, y = hydroperiod)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, colour = "steelblue", size = 2, alpha = 0.7) +
  labs(x = "Oak mortality", y = "Hydroperiod [days]") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 10)
  )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Relationship: hydroperiod vs elevation ---------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cor_elev <- site_duration %>%
  summarise(
    rho     = cor(m, hydroperiod, method = "spearman", use = "complete.obs"),
    p_value = cor.test(m, hydroperiod, method = "spearman")$p.value
  )
cor_elev

site_duration %>%
  ggplot(aes(x = m, y = hydroperiod)) +
  geom_smooth(method = "lm", se = TRUE, colour = "steelblue", 
              fill = "grey80", lwd = 0.75) +
  geom_point(colour = "steelblue", size = 2, alpha = 0.7) +
  labs(x = "Plot Elevation (m)", y = "Hydroperiod [days]") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 10)
  )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Site-level water-table recession rate ----------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Collapse to daily mean water level, then day-to-day decline (cm/day).
#Daily averaging removes hourly sensor jitter that swamped the raw-diff version.

#Step 1: daily mean water level per site
daily_wL <- df %>%
  mutate(date = as.Date(datetime)) %>%
  group_by(site, date) %>%
  summarise(wL_cm_daily = mean(wL_cm, na.rm = TRUE), .groups = "drop")

#Step 2: day-to-day change, keep falling days, one median rate per site
site_recession <- daily_wL %>%
  arrange(site, date) %>%
  group_by(site) %>%
  mutate(
    dwL_cm = wL_cm_daily - lag(wL_cm_daily),                 #change since previous day (cm)
    dt_day = as.numeric(date - lag(date), units = "days"),   #gap in days (usually 1)
    rate   = dwL_cm / dt_day                                 #cm/day
  ) %>%
  filter(rate < 0) %>%                                       #falling days only
  summarise(
    n_fall         = n(),
    recession_rate = median(rate, na.rm = TRUE)
  )

#Take a look
site_recession %>% arrange(recession_rate) %>% print(n = Inf)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Join metadata to recession rate ----------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
site_recession <- site_recession %>%
  filter(!site %in% bad_sites) %>%
  mutate(point_id = parse_number(site)) %>%
  left_join(sites, by = "point_id")

#Check the join
site_recession %>% 
  select(site, point_id, soil_type, Mortality_all, m, recession_rate) %>% 
  print(n = Inf)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Recession by soil ------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
site_recession %>%
  ggplot(aes(x = soil_type, y = recession_rate)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, colour = "steelblue", size = 2, alpha = 0.7) +
  labs(x = "Soil type", y = "Recession rate [cm/day]") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 10)
  )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Pairwise comparison: recession by mortality ----------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mort_rec <- site_recession %>%
  filter(!is.na(Mortality_all)) %>%
  mutate(Mortality_all = factor(Mortality_all, levels = c("None", "Low", "High")))

mort_rec %>% count(Mortality_all)

pw_mort_rec <- mort_rec %>%
  wilcox_test(recession_rate ~ Mortality_all, p.adjust.method = "BH")
pw_mort_rec

pw_mort_rec_letters <- pw_mort_rec %>%
  mutate(comparison = paste(group1, group2, sep = " - ")) %>%
  select(comparison, p.adj)
cld_mort_rec <- cldList(
  p.adj ~ comparison,
  data      = pw_mort_rec_letters,
  threshold = 0.05
)
cld_mort_rec

mort_rec %>%
  ggplot(aes(x = Mortality_all, y = recession_rate)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, colour = "steelblue", size = 2, alpha = 0.7) +
  labs(x = "All-species mortality", y = "Recession rate [cm/day]") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 10)
  )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Pairwise comparison: recession by oak mortality ------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mort_oak_rec <- site_recession %>%
  filter(!is.na(Mortality_oak)) %>%
  mutate(Mortality_oak = factor(Mortality_oak, levels = c("None", "Low", "High")))

mort_oak_rec %>% count(Mortality_oak)

pw_mort_oak_rec <- mort_oak_rec %>%
  wilcox_test(recession_rate ~ Mortality_oak, p.adjust.method = "BH")
pw_mort_oak_rec

pw_mort_oak_rec_letters <- pw_mort_oak_rec %>%
  mutate(comparison = paste(group1, group2, sep = " - ")) %>%
  select(comparison, p.adj)
cld_mort_oak_rec <- cldList(
  p.adj ~ comparison,
  data      = pw_mort_oak_rec_letters,
  threshold = 0.05
)
cld_mort_oak_rec

mort_oak_rec %>%
  ggplot(aes(x = Mortality_oak, y = recession_rate)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, colour = "steelblue", size = 2, alpha = 0.7) +
  labs(x = "Oak mortality", y = "Recession rate [cm/day]") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 10)
  )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Relationship: recession vs elevation -----------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cor_elev_rec <- site_recession %>%
  summarise(
    rho     = cor(m, recession_rate, method = "spearman", use = "complete.obs"),
    p_value = cor.test(m, recession_rate, method = "spearman")$p.value
  )
cor_elev_rec

site_recession %>%
  ggplot(aes(x = m, y = recession_rate)) +
  geom_smooth(method = "lm", se = TRUE, colour = "steelblue", 
              fill = "grey80", lwd = 0.75) +
  geom_point(colour = "steelblue", size = 2, alpha = 0.7) +
  labs(x = "Plot Elevation (m)", y = "Recession rate [cm/day]") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 10)
  )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Site-level inundation event count --------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#An event = water rises above surface (wL_cm > 0) after being at/below it.
#Counted on the DAILY series to avoid noise-driven threshold chatter.

site_events <- daily_wL %>%
  arrange(site, date) %>%
  group_by(site) %>%
  mutate(
    inundated = wL_cm_daily > -30,                          #wet day?
    onset     = inundated & !lag(inundated, default = FALSE)  #first wet day of a run
  ) %>%
  summarise(
    n_events    = sum(onset, na.rm = TRUE),               #number of inundation events
    days_wet    = sum(inundated, na.rm = TRUE),           #total wet days (context)
    .groups = "drop"
  )

#Take a look
site_events %>% arrange(desc(n_events)) %>% print(n = Inf)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Join metadata to event count -------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
site_events <- site_events %>%
  filter(!site %in% bad_sites) %>%
  mutate(point_id = parse_number(site)) %>%
  left_join(sites, by = "point_id")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Regime scatter: flashy vs sustained vs dry -----------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#n_events (how often it floods) vs days_wet (how long wet), coloured by mortality.
#Bottom-left = dry; bottom-right = one sustained flood; top = flashy.
site_events %>%
  filter(!is.na(Mortality_all)) %>%
  mutate(Mortality_all = factor(Mortality_all, levels = c("None", "Low", "High"))) %>%
  ggplot(aes(x = days_wet, y = n_events, colour = Mortality_all)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_colour_manual(values = c("None" = "grey60", "Low" = "steelblue", "High" = "firebrick")) +
  labs(x = "Days wet [days]", y = "Inundation events [count]", colour = "Mortality") +
  theme_bw() +
  theme(
    axis.title   = element_text(size = 14),
    axis.text    = element_text(size = 10),
    legend.title = element_text(size = 12)
  )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Events by soil ---------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
site_events %>%
  ggplot(aes(x = soil_type, y = n_events)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, colour = "steelblue", size = 2, alpha = 0.7) +
  labs(x = "Soil type", y = "Inundation events [count]") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 10)
  )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Pairwise comparison: events by mortality -------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mort_evt <- site_events %>%
  filter(!is.na(Mortality_all)) %>%
  mutate(Mortality_all = factor(Mortality_all, levels = c("None", "Low", "High")))

mort_evt %>% count(Mortality_all)

pw_mort_evt <- mort_evt %>%
  wilcox_test(n_events ~ Mortality_all, p.adjust.method = "BH")
pw_mort_evt

pw_mort_evt_letters <- pw_mort_evt %>%
  mutate(comparison = paste(group1, group2, sep = " - ")) %>%
  select(comparison, p.adj)
cld_mort_evt <- cldList(
  p.adj ~ comparison,
  data      = pw_mort_evt_letters,
  threshold = 0.05
)
cld_mort_evt

mort_evt %>%
  ggplot(aes(x = Mortality_all, y = n_events)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, colour = "steelblue", size = 2, alpha = 0.7) +
  labs(x = "All-species mortality", y = "Inundation events [count]") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 10)
  )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Pairwise comparison: events by oak mortality ---------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mort_oak_evt <- site_events %>%
  filter(!is.na(Mortality_oak)) %>%
  mutate(Mortality_oak = factor(Mortality_oak, levels = c("None", "Low", "High")))

mort_oak_evt %>% count(Mortality_oak)

pw_mort_oak_evt <- mort_oak_evt %>%
  wilcox_test(n_events ~ Mortality_oak, p.adjust.method = "BH")
pw_mort_oak_evt

pw_mort_oak_evt_letters <- pw_mort_oak_evt %>%
  mutate(comparison = paste(group1,