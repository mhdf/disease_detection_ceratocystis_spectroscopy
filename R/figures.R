# Figures of the article, with a simpler layout.

library(patchwork)

theme_article <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(panel.grid = element_blank(), strip.background = element_rect(fill = "white"),
          plot.title = element_text(size = base_size), legend.position = "bottom")
}

wavelength_axis <- function(limits = c(400, 2450)) {
  scale_x_continuous(limits = limits, breaks = seq(400, 2400, 200), expand = expansion(mult = 0.01))
}

figure_1 <- function(meta, folds) {
  meta <- meta %>% mutate(sample = row_number())
  calibration_rows <- which(meta$dataset == "calibration")
  folds_long <- bind_rows(lapply(seq_along(folds), function(i) {
    data.frame(sample = calibration_rows[folds[[i]]], repetition = ceiling(i / N_FOLDS), fold = i %% N_FOLDS + 1)
  }))
  validation_long <- expand.grid(sample = which(meta$dataset == "validation"), repetition = seq_len(N_REPETITIONS))

  sampling <- ggplot(folds_long, aes(sample, factor(repetition))) +
    geom_tile(aes(fill = factor(fold)), colour = "grey70", width = 0.7, height = 0.9) +
    geom_tile(data = validation_long, fill = "white", colour = "black", width = 0.7, height = 0.9) +
    scale_fill_manual(values = rev(c("#3B0A45", "#7B1E5E", "#C13F60", "#E36C3F", "#F7C948")), name = "Outer fold (within repetition)") +
    scale_x_continuous(breaks = seq(1, 77, 4), limits = c(0.5, 77.5), expand = c(0, 0)) +
    labs(title = "Sample distribution: repeated k-fold cross-validation (white: validation samples)", x = NULL, y = "Repetition") +
    theme_article()

  tile_panel <- function(variable, labels, title, name) {
    colours <- if (length(labels) == 2) c("#4BC4A6", "#0A0B2E") else c("#4BC4A6", "#00909E", "#402370", "#0A0B2E")
    ggplot(meta, aes(sample, 1, fill = factor(.data[[variable]], levels = names(labels), labels = labels))) +
      geom_tile(colour = "grey70", width = 0.6, height = 0.9) +
      scale_fill_manual(values = colours, name = name) +
      scale_x_continuous(breaks = seq(1, 77, 4), limits = c(0.5, 77.5), expand = c(0, 0)) +
      labs(title = title, x = NULL, y = NULL) +
      theme_article() +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
  sampling /
    tile_panel("time", c(t1 = "12h", t2 = "24h", t3 = "48h", t4 = "60h"), "Assessment time after inoculation", "Time") /
    tile_panel("clone", c("223" = "Susceptible", "236" = "Resistant"), "Clone per sample", "Clone") /
    tile_panel("inoculation", c(no = "Non-inoculated", yes = "Inoculated"), "Inoculation per sample", NULL) +
    plot_layout(heights = c(3, 1, 1, 1))
}

figure_3 <- function(reflectance, continuum, final_sets, shap) {
  mean_spectrum <- function(table) {
    values <- table %>% filter(leaf == "apical") %>% select(all_of(BANDS)) %>% colMeans()
    data.frame(wavelength = WAVELENGTHS, value = unname(values))
  }
  reflect_mean <- mean_spectrum(reflectance)
  continuum_mean <- mean_spectrum(continuum)
  range_r <- range(reflect_mean$value)
  range_c <- range(continuum_mean$value)
  rescale_c <- function(x) (x - range_c[1]) / diff(range_c) * diff(range_r) + range_r[1]
  spectra <- ggplot() +
    geom_line(data = reflect_mean, aes(wavelength, value), linewidth = 0.8) +
    geom_line(data = continuum_mean, aes(wavelength, rescale_c(value)), linewidth = 0.8, colour = "green3") +
    scale_y_continuous(name = "Reflectance [-]", sec.axis = sec_axis(
      ~ (. - range_r[1]) / diff(range_r) * diff(range_c) + range_c[1], name = "Continuum removal [-]"
    )) +
    wavelength_axis() +
    labs(x = "Wavelength [nm]") +
    theme_article()

  plsda_bands <- shap %>%
    filter(learner == "plsda") %>%
    group_by(pretreatment) %>%
    slice_max(importance, n = ARTICLE_SIZE, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(panel = "PLS-DA") %>%
    select(panel, pretreatment, wavelength)
  selected_bands <- final_sets %>%
    filter(leaf == "apical", target_size == ARTICLE_SIZE) %>%
    mutate(panel = METHOD_LABELS[method]) %>%
    separate_rows(bands, sep = ";") %>%
    mutate(wavelength = as.numeric(sub("C", "", bands))) %>%
    select(panel, pretreatment, wavelength)
  bands <- bind_rows(plsda_bands, selected_bands) %>%
    mutate(panel = factor(panel, levels = c("PLS-DA", "RFE", "GA")),
           pretreatment = factor(pretreatment, levels = PRETREATMENTS, labels = PRETREATMENT_LABELS))
  selected <- ggplot(bands) +
    geom_vline(aes(xintercept = wavelength, colour = panel), linewidth = 0.4) +
    scale_colour_manual(values = c("PLS-DA" = "tomato", RFE = "blue", GA = "orange"), guide = "none") +
    wavelength_axis() +
    facet_grid(panel + pretreatment ~ ., switch = "y") +
    labs(x = "Wavelength [nm]") +
    theme_article() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), strip.text.y.left = element_text(angle = 0),
          panel.spacing = unit(0.15, "lines"), strip.placement = "outside")
  spectra / selected + plot_layout(heights = c(1, 1.7))
}

figure_4 <- function(shap) {
  values <- shap %>%
    mutate(panel = factor(ifelse(learner == "plsda", "PLS-DA", METHOD_LABELS[method]), levels = c("PLS-DA", "RFE", "GA")),
           pretreatment = factor(pretreatment, levels = PRETREATMENTS, labels = PRETREATMENT_LABELS))
  ggplot(values, aes(wavelength, importance, fill = panel)) +
    geom_col(width = 5) +
    scale_fill_manual(values = c("PLS-DA" = "tomato", RFE = "blue", GA = "orange"), guide = "none") +
    scale_x_reverse(breaks = seq(400, 2400, 200)) +
    coord_flip() +
    facet_wrap(~ panel + pretreatment, nrow = 1, scales = "free_x", labeller = labeller(.multi_line = FALSE)) +
    labs(x = "Wavelength [nm]", y = "Feature importance via SHAP values [% from total]") +
    theme_article(9)
}

figure_5 <- function(values) {
  plotted <- values %>%
    mutate(clone = factor(clone, levels = names(CLONE_LABELS), labels = CLONE_LABELS),
           time = factor(time, levels = names(TIME_LABELS), labels = TIME_LABELS),
           model = factor(model, levels = FIGURE_5_MODELS$label),
           upper = pmin(bacc + bacc_sd, 1))
  ggplot(plotted, aes(clone, bacc)) +
    geom_col(aes(fill = model, colour = model), position = position_dodge(0.85), width = 0.7, alpha = 0.4) +
    geom_errorbar(aes(ymin = bacc - bacc_sd, ymax = upper, colour = model), position = position_dodge(0.85), width = 0.15) +
    scale_fill_manual(values = c("PLS-DA" = "tomato", "RF+RFE" = "blue", "SVM+GA" = "orange"), name = NULL) +
    scale_colour_manual(values = c("PLS-DA" = "tomato", "RF+RFE" = "blue", "SVM+GA" = "orange"), name = NULL) +
    scale_y_continuous(limits = c(-0.05, 1.05), breaks = seq(0, 1, 0.2), labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
    facet_wrap(~ time, nrow = 1) +
    labs(x = NULL, y = "Balanced accuracy [%]") +
    theme_article()
}

figure_a1 <- function(severity) {
  panel <- function(table, title, inoculation_label, classes) {
    plotted <- table %>% mutate(class = classes$class, clone = factor(clone, levels = names(CLONE_LABELS), labels = CLONE_LABELS))
    selected <- plotted %>% filter(inoculation == inoculation_label)
    counts <- selected %>% count(clone)
    colours <- rev(dichromat::colorschemes$BrowntoBlue.10)[seq(2, 10, 2)]
    axis_breaks <- seq(min(classes$limits), max(classes$limits), length.out = 6)
    ggplot(selected, aes(severity, clone)) +
      geom_vline(xintercept = classes$limits, linetype = "longdash", linewidth = 0.3) +
      geom_jitter(aes(fill = factor(class, levels = 1:5)), shape = 21, size = 1.8, height = 0.15, width = 0) +
      geom_text(data = counts, aes(x = max(classes$limits), y = clone, label = paste0("n = ", n)), hjust = 1, vjust = -2.2, size = 3) +
      scale_fill_manual(values = colours, labels = classes$labels, name = NULL, drop = FALSE) +
      scale_x_continuous(breaks = axis_breaks, labels = sprintf("%.2f", axis_breaks)) +
      labs(title = title, subtitle = if (inoculation_label == "no") "[Non-inoculated]" else "[Inoculated]", x = "Severity [%]", y = NULL) +
      theme_article(9)
  }
  build <- function(times, title) {
    table <- severity %>% filter(time %in% times)
    classes <- severity_classes(table)
    (panel(table, title, "no", classes) / panel(table, title, "yes", classes)) + plot_layout(guides = "collect")
  }
  (build(PRE_SYMPTOM_TIMES, "Pre-visual symptoms: 12 to 60 hours after inoculation") |
    build(POST_SYMPTOM_TIMES, "Post-visual symptoms: 30 to 40 days after inoculation")) &
    theme(legend.position = "bottom")
}
