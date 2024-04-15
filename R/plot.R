
#' A histogram of density ratio estimates
#'
#' Creates a histogram of the density ratio estimates. Useful to understand the
#' distribution of estimated density ratios in each sample, or compare it among
#' samples. It is the default plotting method for density ratio objects.
#'
#' @param object Density ratio object created with e.g., [kliep()], [ulsif()],
#'   or [naive()]
#' @param sample Character string indicating whether to plot the 'numerator',
#'   'denominator', or 'both' samples. Default is 'both'.
#' @param logscale Logical indicating whether to plot the density ratio
#'   estimates on a log scale. Default is TRUE.
#' @param binwidth Numeric indicating the width of the bins, passed on to
#'   `ggplot2`.
#' @param bins Numeric indicating the number of bins. Overriden by binwidth, and
#'   passed on to `ggplot2`.
#' @param ... Additional arguments passed on to `predict()`.
#'
#' @return A histogram of density ratio estimates.
#'
#'
#' @examples
dr.histogram <- function(object,
                         samples = "both",
                         logscale = TRUE,
                         binwidth = NULL,
                         bins = NULL,
                         ...) {

  # Check object type
  check.object.type(object)


  # Create data object and estimate density ratio
  data <- rbind(object$df_numerator, object$df_denominator)
  ext <- data.frame(dr = predict(object, data = data, ...),
                    sample = c(rep("numerator", nrow(object$df_numerator)),
                               rep("denominator", nrow(object$df_denominator))))

  # Check if logscale is TRUE
  if (logscale) {

    # Convert negative predicted density ratios to 10e-3, so log can be computed
    if(any(ext$dr <= 0)){
      ext$dr[ext$dr <= 0] <- 10e-3
      # Throw warning with number of converted values
      count <- length(ext$dr[ext$dr <= 0])
      warning(
        paste("Negative estimated density ratios for", count, "observation(s) converted to 10e-3 before applying logarithmic transformation"),
        call. = FALSE)
    }

    # Apply log transformation
    ext$dr <- log(ext$dr)
    x_lab <- "Log (Density Ratio)"
  } else {
    x_lab <- "Density Ratio"
  }

  # Select data
  samples <- match.arg(samples, c("both", "numerator", "denominator"))
  if(samples != "both"){
  data <- subset(data, ext$sample == samples)
  ext  <- subset(ext, ext$sample == samples)
  }

  # Plot
  plot <-
    ggplot(ext, aes(x = dr)) +
    geom_histogram(aes(fill = sample),
                   alpha = .75,
                   color = "black",
                   binwidth = binwidth,
                   bins = bins,
                   position = position_dodge2(preserve = "single",
                                              padding = 0.2,
                                              reverse = TRUE)) +

    scale_fill_viridis_d(option = "cividis",
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



#' @inheritParams dr.histogram
#' @rdname dr.histogram
#'
#' @return
#' @export
#'
#' @examples
plot.ulsif <- function(object, samples = "both", logscale = TRUE, binwidth = NULL,
                       bins = NULL) {
  dr.histogram(object, samples = samples, logscale = logscale, binwidth = binwidth,
               bins = bins)
}



#' @inheritParams dr.histogram
#' @rdname dr.histogram
#' @param logscale
#'
#' @return
#' @export
#'
#' @examples
plot.kliep <- function(object, samples = "both", logscale = TRUE, binwidth = NULL,
                       bins = NULL) {
  dr.histogram(object, samples = samples, logscale = logscale, binwidth = binwidth,
               bins = bins)
}



#' Indivual univariate plot
#'
#' Scatterplot of individual values and density ratio estimates. Used internally in [plot.univariate()]
#' @param data Data frame with the individual values and density ratio estimates
#' @param var Name of the variable to be plotted
#' @param y_lab Name of the y-axis label ("Density Ratio" or "Log Density Ratio")
#' @param sample.facet Logical indicating whether to facet the plot by sample. Default is TRUE.
#'
#' @return
#'
create_univariate_plot <- function(data, ext, var, y_lab, sample.facet = TRUE){

  y_max <- max(2, ext$dr)
  y_min <- min(-2, ext$dr)
  plot <-
    ggplot(data, aes(x = .data[[var]], y = ext$dr)) +
    geom_point(aes(col = ext$sample),
               alpha = .6) +
    theme_bw() +
    labs(title = "Scatter plot of individual values and density ratio",
         color = "Sample",
         y = y_lab) +
    geom_hline(yintercept = 0, linetype = "dashed")+
    scale_color_viridis_d(option = "cividis",
                         breaks = c("numerator", "denominator"),
                         labels = c("Numerator", "Denominator")) +
    scale_y_continuous(limits = c(y_min, y_max))

  if(sample.facet){
    plot <- plot + facet_wrap(~vars(ext$sample))
  }

  return(plot)
}

#' Scatter plot of density ratios and individual variables
#'
#' A scatter plot showing the relationship between estimated density ratios and individual variables.
#'
#' @inheritParams dr.histogram
#' @param vars Character vector of variable names to be plotted.
#' @param grid Logical indicating whether output should be a list of individual plots ("individual"), or one facetted plot with all variables ("assembled"). Defaults to "individual".
#' @param sample.facet Logical indicating whether to facet the plot by sample, i.e, showing plots separate for each sample, and side to side. Defaults to FALSE.
#' @param nrow Integer indicating the number of rows in the assembled plot. If NULL, the number of rows is automatically calculated.
#' @param ... Additional arguments passed to the predict() function.
#'
#' @return Scatter plot of density ratios and individual variables.
#' @export
#'
#' @examples
plot_univariate <- function(object, vars, samples = "both", logscale = TRUE,
                            grid = FALSE, sample.facet = FALSE,
                            nrow = NULL, ...) {
  # Check object type
  check.object.type(object)

  # Create data object and estimate density ratio
  data <- rbind(object$df_numerator, object$df_denominator)
  ext <- data.frame(dr = predict(object, data = data, ...),
                    sample = c(rep("numerator", nrow(object$df_numerator)),
                               rep("denominator", nrow(object$df_denominator))))
  # Check variable names
  check.var.names(vars, data)

  # Estimate density ratio
  ext$dr <- predict(object, newdata = data)

  # Select data
  samples <- match.arg(samples, c("both", "numerator", "denominator"))
  if(samples != "both"){
    data <- subset(data, ext$sample == samples)
    ext  <- subset(ext, ext$sample == samples)
  }


  if (logscale) {

    if(any(ext$dr <= 0)){
      # Convert negative predicted density ratios to 10e-3, so log can be computed
      count <- length(ext$dr[ext$dr <= 0])
      ext$dr[ext$dr <= 0] <- 10e-3
      warning(
        paste("Negative estimated density ratios for", count, "observation(s) converted to 10e-3 before applying logarithmic transformation"),
              call. = FALSE)
      }

    ext$dr <- log(ext$dr)

    # Assign correct y and legend labels
    y_lab <- "Log(Density Ratio)"

  } else {
    y_lab <- "Density Ratio"
  }


  if(!grid){

  plot <- lapply(vars, function(var) create_univariate_plot(data, ext, var, y_lab, sample.facet))

  } else {
    values <- data[, vars] |> unlist(use.names = FALSE)
    variable <- rep(vars, each = length(values)/length(vars))
    dr <- rep(ext$dr, length(vars))
    sample <- rep(ext$sample, length(vars))
    data <- data.frame(values = values, variable = variable)
    ext <- data.frame(dr = dr, sample = sample)

    # Maximum scale for y
    y_max <- max(1,ext$dr)
    y_min <- min(-1, ext$dr)

    plot <- ggplot(data) +
      geom_point(aes(x = values, y = ext$dr, col = ext$sample),
                 alpha = .6) +
      theme_bw() +
      labs(title = "Scatter plot of individual values and density ratio",
           color = "Sample",
           y = "Density ratio") +
      scale_color_viridis_d(option = "cividis",
                            breaks = c("numerator", "denominator"),
                            labels = c("Numerator", "Denominator")) +
      scale_y_continuous(limits = c(y_min, y_max))

    if(sample.facet){

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


#' Bivariate plot
#'
#' @inheritParams individual_uni_plot
#' @param show.sample Logical indicating whether to give different shapes to observations, depending on the sample they come from (numerator or denominator). Defaults to FALSE.
#'
#' @return Bivariate plot
#'
#' @examples
create_bivariate_plot  <- function(data, ext, vars, logscale, show.sample){


  dr_max <- ifelse(logscale, max(2, ext$dr), max(exp(2), ext$dr))
  dr_min <- ifelse(logscale, min(-2, ext$dr), min(exp(-2), ext$dr))

  plot <-
    ggplot(data, mapping = aes(x = .data[[vars[1]]], y = .data[[vars[2]]])) +
    geom_point(aes(colour = ext$dr, shape = show.sample),
               alpha = 1,
               size = 2.0) +
    # scale_colour_gradient2(low = "#00204DFF",
    #                        high = "#7D0000",
    #                        mid = "lightyellow",
    #                        midpoint = 0,
    #                        limits = c(dr_min, dr_max)) +
    scale_colour_gradient2(low = "#00204DFF",
                           high = "#FFEA46FF",
                           mid = "#7C7B78FF",
                           midpoint = 0,
                           limits = c(dr_min, dr_max)) +
    theme_bw() +
    labs(title = "Scatter plot, with density ratio mapped to colour",
         colour = "Log (Density ratio)") +
    scale_shape_manual(values = c(21, 24))

  return(plot)
}

#' Densityratio in bidimensional plot
#'
#' Plots a scatterplot of two variables, with densityratio mapped to the colour scale.
#'
#' @inheritParams plot.univariate
#' @param logscale Logical indicating whether to plot the density ratio
#'   estimates on a log scale. Default is TRUE.
#'
#'
#' @return Scatter plot of two variables, with density ratio mapped to the colour scale.
#' @export
#'
#' @examples
plot_bivariate <- function(object, vars, samples = "both",
                           grid = FALSE, logscale = TRUE, show.sample = NULL,
                           ...) {

  # Check object type
  check.object.type(object)

  # Create data object and estimate density ratio
  data <- rbind(object$df_numerator, object$df_denominator)
  ext <- data.frame(dr = predict(object, data = data, ...),
                    sample = c(rep("numerator", nrow(object$df_numerator)),
                               rep("denominator", nrow(object$df_denominator))))
  # Check variable names
  check.var.names(vars, data)

  # Determine if DR is shown in logscale (default) or not
  if (logscale) {

    if(any(ext$dr <= 0)){
      # Convert negative predicted density ratios to 10e-3, so log can be computed
      ext$dr[ext$dr <= 0] <- 10e-3
      count <- length(ext$dr[ext$dr <= 0])
      warning(
        paste("Negative estimated density ratios for", count, "observations converted to 10e-3 before applying logarithmic transformation"),
        call. = FALSE)
    }

    ext$dr <- log(ext$dr)

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

  var_combinations <- as.data.frame(apply(var_combinations, 2,  as.character))
  names(var_combinations) <- c("Var1", "Var2")

  # Remove rows where both variables are the same
  var_combinations <- var_combinations |> subset(Var1 != Var2)
  var_combinations <- as.list(as.data.frame(t(var_combinations)))

  if(!grid){


  plot <- lapply(var_combinations, function(vars) create_bivariate_plot(data, ext, vars, logscale, show.sample))

  return(plot)
  }

  if (grid) {
  browser()
  # Give variable combinations in a format we can use later

    combinations <- sapply(var_combinations, \(vars) paste0(vars, collapse = "-"))

    ext <- data.frame(data, ext$dr, ext$sample)
    datlist <- lapply(var_combinations, \(x) data.frame(values.x = ext[,x[1]],
                                              values.y = ext[,x[2]],
                                              xvar = rep(x[1], nrow(ext)),
                                              yvar = rep(x[2], nrow(ext)),
                                              sample = ext[, "ext.sample"],
                                              dr = ext[, "ext.dr"]))


    plot_data2 <- do.call(datlist, what = rbind)

  plot_data <-
      inner_join(ext, ext, by = c("dr", "sample")) %>% # Possible error in case of duplicate DR?
      pivot_longer(cols = ends_with(".x"), names_to = "name.x", values_to = "value.x") %>%
      pivot_longer(cols = ends_with(".y"), names_to = "name.y", values_to = "value.y") %>%
      mutate(name.x = substr(name.x, 1 ,nchar(name.x)-2),
             name.y = substr(name.y, 1 ,nchar(name.y)-2),
             combination = paste0(name.x, "-", name.y)) %>%
      filter(combination %in% combinations)

  dr_max <- max(1, plot_ext$dr)
  dr_min <- min(-1, plot_ext$dr)

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
 # class(out) <- c("bivariateplot", class(grob))

  grid::grid.newpage()
  grid::grid.draw(out)


  }
}
