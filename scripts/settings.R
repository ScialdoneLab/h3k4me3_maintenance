# axis
axis_titles <- list("peak.width" = "Peak breadth", "promoter.mean" = "Peak mean", "X9hpf" = "RNA expression at 9hpf [log]", "expression.6hpf" = "RNA expression at 6hpf [log]", "expression.7hpf" = "RNA expression at 7hpf [log]", "expression.8hpf" = "RNA expression at 8hpf [log]", "expression.9hpf" = "RNA expression at 9hpf [log]", "cg.density" = "CG density", "promoter.dna.methyl" = "DNA methylation", "promoter.accessibility" = "Accessibility", "expression" = "RNA expression")

# colors
cols_chip <- setNames(c("#DDAA33", "#228833", "#BB5566", "#004488", "#DDDDDD"), c("Kept", "Recovered", "Lost", "Gained", "Absent"))
cols_chip_detail <- setNames(c("#DDCC77", "#009988", "#44BB99", "#228833", "#CC6677", "#882255", "#AA4499", "#332288", "#88CCEE", "#DDDDDD"), c("Kept", "RecoveredLate", "RecoveredEarly", "RecoveredLong", "Lost@ZGA", "Lost@Pre-ZGA", "Lost@Sperm", "GainedEarly", "GainedLate", "Absent"))
cols_chip_complete <- setNames(c("#ece0d3", "#a6afb7", "#937c65", "#625c5a", "#be8162", "#713725", "#b73a2f", "#e17e39", "#f2c237", "#6ca239", "#36753d", "#3b4ca8", "#5f91c1", "#6bc7ab", "#d0589f", "#8d55ae"), c("Y_Y_Y_Y", "Y_Y_Y_N", "Y_Y_N_Y", "Y_Y_N_N", "Y_N_Y_Y", "Y_N_Y_N", "Y_N_N_Y", "Y_N_N_N", "N_Y_Y_Y", "N_Y_Y_N", "N_Y_N_Y", "N_Y_N_N", "N_N_Y_Y", "N_N_Y_N", "N_N_N_Y", "N_N_N_N"))

cols_rna <- c(GS = "#009988", GZ = "#EE7733", ZS = "#6699CC", ND = "#666666")
cols_rna_timing <- c("Expressed@6hpf" = "#9636d1", "Expressed@7hpf" = "#f63189", "Expressed@8hpf" = "#e84c51", "Expressed@9hpf" = "#f96544", "NonZygotic" = "#b7eb34")
cols_rna_complete <- c(Monly = "#1B9E77", Ponly = "#D95F02", Zonly = "#7570B3", MPonly = "#E7298A", MZonly = "#66A61E", PZonly = "#E6AB02", MPZ = "#A6761D", ND = "#666666")

cols_dname <- cols_chip

# sperm data used
sperm_data <- "Oikawa" # Oikawa or Teperek
