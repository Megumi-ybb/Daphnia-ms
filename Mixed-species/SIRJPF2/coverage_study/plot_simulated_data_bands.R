library(ggplot2)
library(dplyr)
library(readxl)
library(tidyr)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
coverage_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
} else if (dir.exists("simulated_data")) {
  getwd()
} else {
  file.path("Mixed-species", "SIRJPF2", "coverage_study")
}

data_dir <- file.path(coverage_dir, "simulated_data")

real_data_files <- c(
  file.path(getwd(), "Mesocosmdata.xls"),
  file.path(getwd(), "Mesocosmdata.xlsx"),
  file.path(coverage_dir, "Mesocosmdata.xls"),
  file.path(coverage_dir, "Mesocosmdata.xlsx"),
  file.path(coverage_dir, "..", "..", "..", "Mesocosmdata.xls"),
  file.path(coverage_dir, "..", "..", "..", "Mesocosmdata.xlsx")
)
real_data_file <- real_data_files[file.exists(real_data_files)][1]

if (is.na(real_data_file)) {
  stop("Could not find Mesocosmdata.xls or Mesocosmdata.xlsx. Put it in the repo root or coverage_study folder.")
}

species_labels <- c(
  dentadult = "D. dentifera adult",
  dentinf = "D. dentifera infected",
  lumadult = "D. lumholtzi adult",
  luminf = "D. lumholtzi infected"
)

sim_files <- list.files(
  data_dir,
  pattern = "^sim_data_[0-9]+\\.rds$",
  full.names = TRUE
)

sim_data <- lapply(sim_files, function(file) {
  dataset_id <- sub("^sim_data_0*([0-9]+)\\.rds$", "\\1", basename(file))
  dataset_id <- as.integer(dataset_id)

  unit_data <- readRDS(file)

  bind_rows(lapply(names(unit_data), function(unit_name) {
    unit_data[[unit_name]] |>
      mutate(dataset = dataset_id, unit = unit_name, .before = 1)
  }))
}) |>
  bind_rows() |>
  pivot_longer(
    cols = c(dentadult, dentinf, lumadult, luminf),
    names_to = "species",
    values_to = "density"
  ) |>
  mutate(
    species = recode(species, !!!species_labels)
  )

real_data <- read_excel(real_data_file, sheet = 3) |>
  slice(91:170) |>
  select(
    unit = rep,
    day,
    dentadult = dent.adult,
    dentinf = dent.inf,
    lumadult = lum.adult,
    luminf = lum.adult.inf
  ) |>
  slice(n():1) |>
  mutate(day = (day - 1) * 5 + 7) |>
  pivot_longer(
    cols = c(dentadult, dentinf, lumadult, luminf),
    names_to = "species",
    values_to = "density"
  ) |>
  mutate(
    species = recode(species, !!!species_labels)
  )

band_data <- sim_data |>
  group_by(day, species) |>
  summarize(
    q025 = quantile(density, 0.001, na.rm = TRUE),
    q25 = quantile(density, 0.25, na.rm = TRUE),
    median = median(density, na.rm = TRUE),
    q75 = quantile(density, 0.75, na.rm = TRUE),
    q975 = quantile(density, 0.999, na.rm = TRUE),
    .groups = "drop"
  )

p <- ggplot(band_data, aes(x = day)) +
  geom_ribbon(aes(ymin = q025, ymax = q975), fill = "#7aa6c2", alpha = 0.25) +
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = "#2f78a0", alpha = 0.35) +
  geom_line(aes(y = median), color = "#123c54", linewidth = 0.8) +
  geom_line(
    data = real_data,
    aes(y = density, group = unit),
    color = "red",
    linewidth = 0.6,
    alpha = 0.75
  ) +
  facet_wrap(~ species, scales = "free_y", ncol = 2) +
  labs(
    x = "Time",
    y = "Density",
    title = "Simulated Data Bands",
    subtitle = "Dark band: 25-75%; light band: 0.1-99.9%; blue line: simulated median; red lines: real data"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey90", color = NA),
    plot.title = element_text(face = "bold")
  )

print(p)
