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
      y = "Count"
    ) +
    labs(title = "Distribution of density ratio estimates",
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
individual_uni_plot <- function(data, var, y_lab){

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
                            output = "individual",
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
    plots[[var]] <- individual_uni_plot(data, var, y_lab)
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
      geom_hline(yintercept = 0, linetype = "dashed")+
      scale_color_manual(values = c("firebrick", "steelblue"),
                         breaks = c("numerator", "denominator"),
                         labels = c("Numerator", "Denominator")) +
      facet_wrap(~variable, scales = "free_x", nrow = nrow) +
      scale_y_continuous(limits = c(y_min, y_max))

  }
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
plot_bivariate <- function(object, vars, sample = "both", show.samples = TRUE,
                           output = "assembled", logscale = TRUE) {

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
  obsclass <- rep(c("numerator", "denominator"),
                  c(nrow(object$df_numerator), nrow(object$df_denominator)))

  # Create a object selection variable (both, numerator, denominator)
  obsselect <- match.arg(sample, c("both", "numerator", "denominator"))

  # Filter data based on object selection
  if (obsselect != "both") {
    data <- filter(data, obsclass == obsselect)
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
  var_combinations <- var_combinations %>% filter(Var1 != Var2)
  var_combinations <- as.matrix(var_combinations)

  if(output == "individual"){
  # Define function to make bivariate plot
  plot_twovariables <- function(data, vars, showsamples = show.samples){
    plot <-
    ggplot(data, mapping = aes(x = .data[[vars[1]]], y = .data[[vars[2]]])) +
    geom_point(aes(colour = dr, fill = dr, shape = if (showsamples) sample else NULL),
               alpha = 0.5) +
    scale_fill_viridis_c(option = "B", name = colour_name) +
    scale_colour_viridis_c(option = "B", name = colour_name) +
    theme_bw() +
    labs(title = "Density ratio estimates for combinations of values",
         shape = "Sample") +
    scale_shape_manual(values = c(21, 24))

  return(plot)
    }

  # Iterate over grid of variable combinations
  # Create a list to store plots
  plots <- list()
  for(i in 1:nrow(var_combinations)){
    plots[[i]] <- plot_twovariables(data = data, vars = var_combinations[i,])
  }
  return(plots)
  }

  if (output == "assembled") {

  # plot_data <- data %>%
  #   pivot_longer(cols = c("V1", "V2"), names_to = "variable", values_to = "value")
  # browser()
  combinations <- paste0(var_combinations[,1], "-", var_combinations[,2])

  plot_data <-
    inner_join(data, data, by = "dr") %>%
      pivot_longer(cols = ends_with(".x"), names_to = "name.x", values_to = "value.x") %>%
      pivot_longer(cols = ends_with(".y"), names_to = "name.y", values_to = "value.y") %>%
      mutate(name.x = stringr::str_remove(name.x, ".x"),
             name.y = stringr::str_remove(name.y, ".y"),
             combination = paste0(name.x, "-", name.y)) %>%
      filter(combination %in% combinations)


  plot <-
      ggplot(plot_data, mapping = aes(x = value.x, y = value.y)) +
      geom_point(aes(colour = dr, fill = dr),
                 alpha = 0.5) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      facet_grid(rows = vars(name.y), cols = vars(name.x), scales = "free",
                 switch = "both") +
      scale_fill_viridis_c(option = "B", name = colour_name) +
      scale_colour_viridis_c(option = "B", name = colour_name) +
      scale_y_continuous(position = "right") +
      scale_x_continuous(position = "top") +
      theme_bw() +
      labs(title = "Density ratio estimates for combinations of values",
          x = NULL,
          y = NULL) +
      scale_shape_manual(values = c(21, 24))

  return(suppressWarnings(print(plot)))

  }
}

plot_bivariate_heatmap <- function(object, var.x, var.y, sample = "both", show.samples = TRUE,
                                   output = "assembled", log.scale = TRUE) {

  # Check object type
  check.object.type(object)

  # Create data object and check variable names
  data <- rbind(object$df_numerator, object$df_denominator)
  check.overriden.names(data)

  # Create variables vector and check variable names
  vars <- c(var.x, var.y)
  check.var.names(vars, data)

  # Create a sample index variable (denominator or numerator)
  obsclass <- rep(c("numerator", "denominator"),
                  c(nrow(object$df_numerator), nrow(object$df_denominator)))

  # Create a object selection variable (both, numerator, denominator)
  obsselect <- match.arg(sample, c("both", "numerator", "denominator"))

  # Filter data based on object selection
  if (obsselect != "both") {
    data <- filter(data, obsclass == obsselect)
  }


  # Create a grid of variable combinations
  var_combinations <- expand.grid(var.x, var.y)

  ## Remove duplicate combinations
  ## Start by sorting elements within each row
  ## This makes duplicate rows with different order of variables identical
  var_combinations <- t(apply(var_combinations, 1, sort))
  var_combinations <- unique(var_combinations) # retain unique rows only

  # Remove rows where both variables are the same
  var_combinations <- as.data.frame(apply(var_combinations, 2,  as.character))
  names(var_combinations) <- c("Var1", "Var2")
  var_combinations <- var_combinations %>% filter(Var1 != Var2)

  var_combinations <- as.matrix(var_combinations)

  object2 <- object
  # Define function to make bivariate plot
  plot_twovariables_heatmap <- function(object = object2, data, vars, showsamples = show.samples, logscale = log.scale){

    # Create a 100x100 grid of values for the two variables, in the range of the data
    seq_var1 <- seq(min(data[[vars[1]]]), max(data[[vars[1]]]), length.out = 100)
    seq_var2 <- seq(min(data[[vars[2]]]), max(data[[vars[2]]]), length.out = 100)
    grid_data <- expand.grid(seq_var1, seq_var2)
    colnames(grid_data) <- c(vars[[1]], vars[[2]])

    # Add the rest of data variables to the grid, inputting its mean value
    # First, identify variables in data that are not in vars
    other_vars <- setdiff(names(data), vars)
    for (var in other_vars) {
      grid_data[[var]] <- mean(data[[var]])
    }

    # Predict density ratio for each combination of values
    # Assign to dataframe and reorder column in the same order as data, so that predict works
    grid_data <- as.data.frame(grid_data)
    grid_data <- grid_data[, names(data)]
    grid_data$dr <- predict(object, newdata = grid_data)

    # Estimate density ration in the original datapoints (for superposing)
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

      if(any(grid_data$dr <= 0)){
        # Convert negative predicted HEATMAP density ratios to 10e-3, so log can be computed
        grid_data$dr[grid_data$dr <= 0] <- 10e-3
        warning("Negative estimated density ratios in the heatmap converted to 10e-3 before applying logarithmic transformation",
                call. = FALSE)
      }

      grid_data$dr <- log(grid_data$dr)
      data$dr <- log(data$dr)

      # Assign correct y and legend labels
      colour_name <- "Log (Density ratio)"

    } else {
      colour_name <- "Density ratio"
    }

    # Plot
    plot <-
      ggplot(grid_data, mapping = aes(x = .data[[vars[1]]], y = .data[[vars[2]]])) +
      geom_point(data = data,
                 colour = "grey60",
                 alpha = 1,
                 shape = 1) +
      geom_raster(aes(colour = dr, fill = dr, shape = if (showsamples) sample else NULL),
                  alpha = 0.6) +
      scale_fill_viridis_c(option = "B", name = colour_name) +
      theme_bw() +
      labs(title = "Density ratio estimates for combinations of values",
           shape = "Sample") +
      scale_shape_manual(values = c(21, 24))
    return(plot)
  }

  # Iterate over grid of variable combinations
  # Create a list to store plots
  plots <- list()
  for(i in 1:nrow(var_combinations)){
    plots[[i]] <- plot_twovariables_heatmap(data = data, vars = var_combinations[i,])
  }

  # Assemble plots
  if(output == "assembled"){

    plots_assembly <- patchwork::wrap_plots(plots,
                                            guides = "collect", byrow = TRUE,
    axes = "collect",
    ncol = length(var.x), nrow = length(var.y)) & labs(title = NULL)
    return(plots_assembly)
  } else {
    return(plots)
  }
}
