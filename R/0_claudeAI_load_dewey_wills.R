# =============================================================
# Dewey Wills PT data -> tidy long format (site, datetime, wL)
# =============================================================
# Two options below.
#  OPTION A: just read the cleaned CSV Claude already produced.
#  OPTION B: re-parse the raw .xlsx yourself (use when Campbell
#            adds sheets / you get a fresh export).
# -------------------------------------------------------------

library(tidyverse)
library(readxl)

# ---- OPTION A: read the cleaned CSV ------------------------
dw <- read_csv(
  "dewey_wills_wL_long.csv",
  col_types = cols(
    site     = col_character(),
    datetime = col_datetime(format = ""),
    wL       = col_double()
  )
)

# ---- OPTION B: re-parse the raw workbook -------------------
# Anchors wL to the column immediately right of "H gauge".
# Sheets with no "H gauge" header fall back to column P (16).

parse_dewey_wills <- function(path) {

  drop_sheets <- c(
    "Master abc", "High Mortality", "Low Mortality", "No Mortality",
    "Baro 1039", "43 Baro", "Fieldhouse Baro", "Fieldhouse Baro.1",
    "nuh uh"
  )

  well_sheets <- setdiff(excel_sheets(path), drop_sheets)

  read_one <- function(sheet) {
    # read header row to locate "H gauge"
    hdr <- read_excel(path, sheet = sheet, n_max = 1, col_names = FALSE,
                      .name_repair = "minimal")
    hdr <- as.character(unlist(hdr[1, ], use.names = FALSE))

    hg  <- which(str_to_lower(str_trim(hdr)) == "h gauge")
    wl_col <- if (length(hg) == 1) hg + 1L else 16L   # fallback: column P

    raw <- read_excel(path, sheet = sheet, col_names = FALSE, skip = 1,
                      .name_repair = "minimal")

    tibble(
      site     = sheet,
      datetime = as.POSIXct(raw[[1]]),
      wL       = suppressWarnings(as.numeric(raw[[wl_col]]))
    ) %>%
      filter(!is.na(datetime))
  }

  map_dfr(well_sheets, read_one)
}

# dw <- parse_dewey_wills("Dewey_Wills_PT_Data.xlsx")

# ---- Quick look --------------------------------------------
dw %>% count(site) %>% print(n = Inf)

# ---- Flags (see extraction_log.txt) ------------------------
# 1038 : wL positive & erratic — suspect.
# 1026, 1031 : column-P fallback, values positive — may not equal wL.
# Inspect these against raw pressure before publishing.

# ---- Example plot (your house style) -----------------------
# One site; facet or filter as needed.
dw %>%
  filter(site == "39") %>%
  ggplot(aes(datetime, wL)) +
  geom_line(colour = "steelblue", lwd = 0.75) +
  labs(x = NULL, y = "Water level") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 10)
  )
