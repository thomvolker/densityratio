
dr.histogram <- function(object, sample = "both", logscale = FALSE, binwidth = NULL, bins = NULL,...) {

  # Checks
  check.object.type(object)
  check.overriden.names(vars)

  # Create data object and estimate density ratio
  data <- rbind(object$df_numerator, object$df_denominator)
  data$dr <- predict(object, newdata = data, ...)

  if (logscale) {

    if(any(data$dr <= 0)){
      # Convert negative predicted density ratios to 10e-3, so log can be computed
      count <- length(data$dr[data$dr <= 0])
      data$dr[data$dr <= 0] <- 10e-3
      warning(
        paste("Negative estimated density ratios for", count, "observation(s) converted to 10e-3 before applying logarithmic transformation"),
        call. = FALSE)
    }

    data$dr <- log(data$dr)
    x_lab <- "Log (Density Ratio)"
  } else {
  x_lab <- "Density Ratio"
  }

  # Create a sample index variable (denominator or numerator)
  obsclass <- rep(c("numerator", "denominator"),
                  c(nrow(object$df_numerator), nrow(object$df_denominator)))
  data$sample <- obsclass

  # Create a object selection variable (both, numerator, denominator)
  obsselect <- match.arg(sample, c("both", "numerator", "denominator"))

  # If not both, subset data (only num or only den)
  if (obsselect != "both") {
    data <- filter(data, obsclass == obsselect)
  }

  # Plot
  plot <-
    ggplot(data, aes(x = dr)) +
    geom_histogram(aes(fill = sample),
                   alpha = .75,
                   color = "black",
                   binwidth = if (!is.null(binwidth)) binwidth else NULL,
                   bins = if(!is.null(bins)) bins else NULL,
                   position = position_dodge2(preserve = "single",
                                              padding = 0.2,
                                              reverse = TRUE)) +
    scale_fill_manual(values = c("firebrick", "steelblue"),
                      breaks = c("numerator", "denominator"),
                      labels = c("Numerator", "Denominator")) +
    theme_bw()  +
    labs(
      x = x_lab,
      y = "Count",
      title = "Distribution of density ratio estimates",
      fill = "Sample")

  return(plot)
}

#' Title
#'
#' @param object
#' @param sample
#' @param logscale
#'
#' @return
#' @export
#'
#' @examples
plot.ulsif <- function(object, sample = "both", logscale = FALSE, binwidth = NULL,
                       bins = NULL) {
  dr.histogram(object, sample = sample, logscale = logscale, binwidth = binwidth,
               bins = bins)
}

#' Title
#'
#' @param object
#' @param sample
#' @param logscale
#'
#' @return
#' @export
#'
#' @examples
plot.kliep <- function(object, sample = "both", logscale = FALSE, binwidth = NULL,
                       bins = NULL) {
  dr.histogram(object, sample = sample, logscale = logscale, binwidth = binwidth,
               bins = bins)
}


# Write the function to make one plot, for one variable
individual_uni_plot <- function(data, var, y_lab, facet = TRUE){

  y_max <- max(1,data$dr)
  y_min <- min(-1, data$dr)
  plot <-
    ggplot(data, aes(x = .data[[var]], y = dr)) +
    geom_point(aes(col = sample),
               alpha = .6) +
    theme_bw() +
    labs(title = "Scatter plot of individual values and density ratio",
         color = "Sample",
         y = y_lab) +
    geom_hline(yintercept = 0, linetype = "dashed")+
    scale_colour_manual(values = c("firebrick", "steelblue"),
                      breaks = c("numerator", "denominator"),
                      labels = c("Numerator", "Denominator")) +
    scale_y_continuous(limits = c(y_min, y_max))

  if(facet){
    plot <- plot + facet_wrap(~sample)
  }

  return(plot)
}

#' Title
#'
#' @param object
#' @param vars
#' @param sample
#'
#' @return
#' @export
#'
#' @examples
plot_univariate <- function(object, vars, samples = "both", logscale = TRUE,
                            output = "individual", facet = FALSE,
                            nrow = NULL) {

  # Check object type
  check.object.type(object)

  # Create data object
  data <- rbind(object$df_numerator, object$df_denominator)

  # Check names in data and variable names
  check.overriden.names(data)
  check.var.names(vars, data)

  # Estimate density ratio
  data$dr <- predict(object, newdata = data)

  # Creta sample identifier
  data$sample <- rep(c("numerator", "denominator"),
                     c(nrow(object$df_numerator), nrow(object$df_denominator)))

  # Create a object selection variable (both, numerator, denominator)
  obsselect <- match.arg(samples, c("both", "numerator", "denominator"))

  if (obsselect != "both") {
    data <- filter(data, sample == obsselect)
  }


  if (logscale) {

    if(any(data$dr <= 0)){
      # Convert negative predicted density ratios to 10e-3, so log can be computed
      count <- length(data$dr[data$dr <= 0])
      data$dr[data$dr <= 0] <- 10e-3
      warning(
        paste("Negative estimated density ratios for", count, "observation(s) converted to 10e-3 before applying logarithmic transformation"),
              call. = FALSE)
      }

    data$dr <- log(data$dr)

    # Assign correct y and legend labels
    y_lab <- "Log(Density Ratio)"

  } else {
    y_lab <- "Density Ratio"
  }

  if(output == "individual"){
  # Create list storage for plots object (for iteration)
  plots <- list()
  for(var in vars){
    plots[[var]] <- individual_uni_plot(data, var, y_lab, facet)
  }
  return(plots)

  }
  if (output == "assembled"){
    data <- data %>%
      pivot_longer(cols = vars,
                   names_to = "variable",
                   values_to = "value")

    y_max <- max(1,data$dr)
    y_min <- min(-1, data$dr)


    plot <- ggplot(data) +
      geom_point(aes(x = value, y = dr, col = sample),
                 alpha = .6) +
      theme_bw() +
      labs(title = "Scatter plot of individual values and density ratio",
           color = "Sample",
           y = "Density ratio") +
      scale_color_manual(values = c("firebrick", "steelblue"),
                         breaks = c("numerator", "denominator"),
                         labels = c("Numerator", "Denominator")) +
      scale_y_continuous(limits = c(y_min, y_max))

    if(facet){
      plot <- plot +
        facet_grid(cols = vars(sample),
                   rows = vars(variable),
                   scales = "free_x")

      } else {
        plot <- plot +
          facet_wrap(~variable, scales = "free_x", nrow = nrow)
    }

  }
  return(plot)

}


# Define function to make bivariate plot
plot_twovariables <- function(data, vars, logscale, show.sample){

  dr_max <- ifelse(logscale, max(2, data$dr), max(exp(2), data$dr))
  dr_min <- ifelse(logscale, min(-2, data$dr), min(exp(-2), data$dr))

  plot <-
    ggplot(data, mapping = aes(x = .data[[vars[1]]], y = .data[[vars[2]]])) +
    geom_point(aes(colour = dr, shape = if(show.sample) sample else NULL)) +
    scale_colour_gradient2(low = "firebrick",
                           high = "steelblue",
                           mid = "lightyellow",
                           limits = c(dr_min, dr_max)) +
    theme_bw() +
    labs(title = "Scatter plot, with density ratio mapped to colour",
         colour = "Log (Density ratio)") +
    scale_shape_manual(values = c(21, 24))

  return(plot)
}

#' Title
#'
#' @param object
#' @param vars
#' @param sample
#' @param show.samples
#'
#' @return
#' @export
#'
#' @examples
plot_bivariate <- function(object, vars, samples = "both",
                           output = "assembled", logscale = TRUE, show.sample = FALSE) {

  # Check object type
  check.object.type(object)

  # Create data object and check variable names
  data <- rbind(object$df_numerator, object$df_denominator)
  check.overriden.names(data)

  # Check variable names
  check.var.names(vars, data)

  # Estimate density ratio
  data$dr <- predict(object, newdata = data)

  # Determine if DR is shown in logscale (default) or not
  if (logscale) {

    if(any(data$dr <= 0)){
      # Convert negative predicted density ratios to 10e-3, so log can be computed
      count <- length(data$dr[data$dr <= 0])
      data$dr[data$dr <= 0] <- 10e-3
      warning(
        paste("Negative estimated density ratios for", count, "observations converted to 10e-3 before applying logarithmic transformation"),
        call. = FALSE)
    }

    data$dr <- log(data$dr)

    # Assign correct y and legend labels
    colour_name <- "Log (Density ratio)"

  } else {
    colour_name <- "Density ratio"
  }


  # Create a sample index variable (denominator or numerator)
  data$sample <- rep(c("numerator", "denominator"),
                  c(nrow(object$df_numerator), nrow(object$df_denominator)))

  # Create a object selection variable (both, numerator, denominator)
  obsselect <- match.arg(samples, c("both", "numerator", "denominator"))

  # Filter data based on object selection
  if (obsselect != "both") {
    data <- filter(data, sample == obsselect)
  }

  # Create a grid of variable combinations
  var_combinations <- expand.grid(vars, vars)

  ## Remove duplicate combinations
  ## Start by sorting elements within each row
  ## This makes duplicate rows with different order of variables identical
  var_combinations <- t(apply(var_combinations, 1, sort))
  var_combinations <- unique(var_combinations) # retain unique rows only

  # Remove rows where both variables are the same
  var_combinations <- as.data.frame(apply(var_combinations, 2,  as.character))
  names(var_combinations) <- c("Var1", "Var2")

  if(output == "individual"){

  var_combinations <- var_combinations %>% filter(Var1 != Var2)
  var_combinations <- as.matrix(var_combinations)

  plots <- list()
  for(i in 1:nrow(var_combinations)){
    plots[[i]] <- plot_twovariables(data, vars = var_combinations[i,], logscale, show.sample)
  }
  return(plots)
  }

  if (output == "assembled") {

  # plot_data <- data %>%
  #   pivot_longer(cols = c("V1", "V2"), names_to = "variable", values_to = "value")
  # browser()
  combinations <- paste0(var_combinations[,1], "-", var_combinations[,2])

  plot_data <-
      inner_join(data, data, by = c("dr", "sample")) %>% # Possible error in case of duplicate DR?
      pivot_longer(cols = ends_with(".x"), names_to = "name.x", values_to = "value.x") %>%
      pivot_longer(cols = ends_with(".y"), names_to = "name.y", values_to = "value.y") %>%
      mutate(name.x = stringr::str_remove(name.x, ".x"),
             name.y = stringr::str_remove(name.y, ".y"),
             combination = paste0(name.x, "-", name.y)) %>%
      filter(combination %in% combinations)

  dr_max <- max(1, plot_data$dr)
  dr_min <- min(-1, plot_data$dr)

  plot <-
      ggplot(plot_data, mapping = aes(x = value.x, y = value.y,
                                      shape = if(show.sample) sample else NULL)) +
      geom_point(aes(colour = dr)) +
      facet_grid(rows = vars(name.y), cols = vars(name.x), scales = "free",
                 switch = "both") +
      scale_colour_gradient2(low = "firebrick",
                           high = "steelblue",
                           mid = "#ffffbf",
                           limits = c(dr_min, dr_max),
                            ) +
      scale_y_continuous(position = "left") +
      scale_x_continuous(position = "bottom") +
      theme_bw() +
      theme(strip.placement = "outside") +
      labs(title = "Scatter plots, with density ratio mapped to colour",
          x = NULL,
          y = NULL,
          colour = colour_name,
          shape = if(show.sample) "Sample" else NULL) +
      scale_shape_manual(values = c(21, 24))

  # Erase upper diagonal
  ## Create plot into a grob
  grob <- ggplotGrob(plot)
  ## Create name of empty panels in the upper diagonal
  empty_panels <- expand.grid(seq(1:length(vars)), seq(1:length(vars))) %>%
    filter(Var2 > Var1) %>%
    mutate(panel = paste0("panel-", Var1, "-", Var2)) %>%
    pull(panel)
  # Delete panels in upper diagonal, based in their index
  idx <- which(grob$layout$name %in% empty_panels)
  for (i in idx) grob$grobs[[i]] <- grid::nullGrob()

  out <- grob
  class(out) <- c("bivariateplot", class(grob))

  return(out)
  # grid::grid.newpage()
  # grid::grid.draw(grob)

  }
}
