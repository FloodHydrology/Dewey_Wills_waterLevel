#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: Dewey Wet vs Dry Hydrographs
#Date: 8/6/2026
#Coder: Nate Jones (cnjones7@ua.edu)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Two example wetland hydrographs -- one wet site, one dry site -- drawn from
# the real (sub-daily) water-level series over each site's full record. Shaded
# band marks the rooting zone (-30 to 0 cm). Style matches the exploratory script.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setup workspace ---------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Clear memory
remove(list = ls())

#Load packages
library(tidyverse)
library(patchwork)

#Load data
df <- read_csv("data/dewey_wills_wL_long.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Global settings ---------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Rooting-depth threshold (cm).
threshold_cm <- -30

#Example sites: one wet, one dry (picked by hydroperiod from the full analysis).
#57 is genuinely wet with real dynamic range; 45 stays below rooting depth.
wet_site <- "57"
dry_site <- "45"

#Shared theme (house style, sized for slides).
theme_dewey <- theme_bw(base_size = 14) +
  theme(
    axis.title   = element_text(size = 16, lineheight = 0.9),
    axis.text    = element_text(size = 12),
    plot.title   = element_text(size = 15),
    plot.margin  = margin(8, 8, 8, 8)
  )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Prep data ---------------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Convert inches -> cm and keep the two example sites (full raw series).
wL <- df %>%
  mutate(wL_cm = wL * 2.54) %>%
  filter(site %in% c(wet_site, dry_site))

#Shared y-axis limits so the two panels are directly comparable.
y_lims <- range(wL$wL_cm, na.rm = TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Hydrograph panel function -----------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#One site's raw trace, steelblue line, shaded rooting-zone band (-30 to 0).
#show_y toggles the y-axis TITLE only (numbers stay, so panel widths match).
hydro_panel <- function(data, this_site, title, show_y = TRUE) {
  data %>%
    filter(site == this_site) %>%
    ggplot(aes(x = datetime, y = wL_cm)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = threshold_cm, ymax = 0,
             fill = "grey70", alpha = 0.25) +
    geom_line(colour = "steelblue", lwd = 0.75) +
    coord_cartesian(ylim = y_lims) +
    labs(x = NULL,
         y = if (show_y) "Water level [cm]" else NULL,
         title = title) +
    theme_dewey
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Build + combine (side by side) ------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p_wet <- hydro_panel(wL, wet_site, paste0("Wet site (", wet_site, ")"), show_y = TRUE)
p_dry <- hydro_panel(wL, dry_site, paste0("Dry site (", dry_site, ")"), show_y = FALSE)

fig_hydro <- p_wet + p_dry
fig_hydro

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Export for PowerPoint ---------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Landscape, sized to sit centered on a 16:9 slide with margins.
ggsave("output/fig_wet_dry_hydrographs.png", fig_hydro,
       width = 9, height = 4, dpi = 300, bg = "white")
