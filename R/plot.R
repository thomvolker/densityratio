
#' A histogram of density ratio estimates
#'
#' Creates a histogram of the density ratio estimates. Useful to understand the
#' distribution of estimated density ratios in each sample, or compare it among
#' samples. It is the default plotting method for density ratio objects.
#'
#' @param object Density ratio object created with e.g., [kliep()], [ulsif()],
#'   or [naive()]
#' @param samples Character string indicating whether to plot the 'numerator',
#'   'denominator', or 'both' samples. Default is 'both'.
#' @param tol Numeric indicating the tolerance: values below this value will be set to the tolerance value, for legibility of the plots
#' @param logscale Logical indicating whether to plot the density ratio
#'   estimates on a log scale. Default is TRUE.
#' @param binwidth Numeric indicating the width of the bins, passed on to
#'   `ggplot2`.
#' @param bins Numeric indicating the number of bins. Overriden by binwidth, and
#'   passed on to `ggplot2`.
#' @param ... Additional arguments passed on to `predict()`.
#'
#' @return A histogram of density ratio estimates.
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 position_dodge2
#' @importFrom ggplot2 geom_histogram
#' @importFrom ggplot2 scale_fill_viridis_d
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 labs
#'
#'
dr.histogram <- function(object,
                         samples = "both",
                         logscale = TRUE,
                         binwidth = NULL,
                         bins = NULL,
                         tol = 10e-3,
                         ...) {

  # Check object type
  check.object.type(object)

  # Create data object and estimate density ratio
  data <- rbind(object$df_numerator, object$df_denominator)
  ext <- data.frame(dr = predict(object, newdata = data, ...),
                    sample = c(rep("numerator", nrow(object$df_numerator)),
                               rep("denominator", nrow(object$df_denominator))))

  # If logscale = TRUE, transform density ratio estimates
  ext <- check.logscale(ext, logscale, tol)

  # Assign x-axis label
  if (logscale) {
    x_lab <- "Log (Density Ratio)"
  } else {
    x_lab <- "Density Ratio"
  }

  # Filter correct subset
  samples <- match.arg(samples, c("both", "numerator", "denominator"))
  if(samples != "both"){
  data <- subset(data, ext$sample == samples)
  ext  <- subset(ext, ext$sample == samples)
  }

  # Plot
  plot <-
    ggplot2::ggplot(data, ggplot2::aes(x = ext$dr)) +
    ggplot2::geom_histogram(ggplot2::aes(fill = ext$sample),
                   alpha = .75,
                   color = "black",
                   binwidth = binwidth,
                   bins = bins,
                   position = ggplot2::position_dodge2(preserve = "single",
                                              padding = 0.2,
                                              reverse = TRUE)) +

    ggplot2::scale_fill_viridis_d(option = "cividis",
                        breaks = c("numerator", "denominator"),
                        labels = c("Numerator", "Denominator")) +
    ggplot2::theme_bw()  +
    ggplot2::labs(
      x = x_lab,
      y = "Count",
      title = "Distribution of density ratio estimates",
      fill = "Sample")

  return(plot)
}



#' @inheritParams dr.histogram
#' @rdname dr.histogram
#'
#' @return A histogram of density ratio estimates.
#' @export
#'
plot.ulsif <- function(object, samples = "both", logscale = TRUE, binwidth = NULL,
                       bins = NULL, tol = 10e-3) {
  dr.histogram(object, samples = samples, logscale = logscale, binwidth = binwidth,
               bins = bins, tol = tol)
}



#' @inheritParams dr.histogram
#' @rdname dr.histogram
#' @param logscale
#'
#' @return A histogram of density ratio estimates.
#' @export
#'

plot.kliep <- function(object, samples = "both", logscale = TRUE, binwidth = NULL,
                       bins = NULL, tol = 10e-3) {
  dr.histogram(object, samples = samples, logscale = logscale, binwidth = binwidth,
               bins = bins, tol = tol)
}

#' @inheritParams dr.histogram
#' @rdname dr.histogram
#'
#' @return A histogram of density ratio estimates.
#' @export
#'

plot.spectral <- function(object, samples = "both", logscale = TRUE, binwidth = NULL,
                       bins = NULL, tol = 10e-3) {
  dr.histogram(object, samples = samples, logscale = logscale, binwidth = binwidth,
               bins = bins, tol = tol)
}

#' @inheritParams dr.histogram
#' @rdname dr.histogram
#'
#' @return A histogram of density ratio estimates.
#' @export
#'
plot.lhss <- function(object, samples = "both", logscale = TRUE, binwidth = NULL,
                       bins = NULL, tol = 10e-3) {
  dr.histogram(object, samples = samples, logscale = logscale, binwidth = binwidth,
               bins = bins, tol = tol)
}

#' Indivual univariate plot
#'
#' Scatterplot of individual values and density ratio estimates. Used internally in [create_univariate_plot()]
#' @param data Data frame with the individual values and density ratio estimates
#' @param var Name of the variable to be plotted on the x-axis
#' @param y_lab Name of the y-axis label, typically ("Density Ratio" or "Log Density Ratio")
#' @param ext Data frame with the density ratio estimates and sample indicator
#' @param sample.facet Logical indicating whether to facet the plot by sample. Default is TRUE.
#'
#' @return A scatterplot of variable values and density ratio estimates.
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_color_viridis_d
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 geom_hline
#'
#'
create_univariate_plot <- function(data, ext, var, y_lab, sample.facet = TRUE){

  y_max <- max(2, ext$dr)
  y_min <- min(-2, ext$dr)

  plot <-
    ggplot2::ggplot(data, ggplot2::aes(x = .data[[var]], y = ext$dr)) +
    ggplot2::geom_point(ggplot2::aes(col = ext$sample),
               alpha = .6) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Scatter plot of individual values and density ratio",
         color = "Sample",
         y = y_lab) +
    ggplot2::scale_color_viridis_d(option = "cividis",
                         breaks = c("numerator", "denominator"),
                         labels = c("Numerator", "Denominator")) +
    ggplot2::scale_y_continuous(limits = c(y_min, y_max))

  if(sample.facet){
    plot <- plot + ggplot2::facet_wrap(~ext$sample)
  }

  plot <- plot +
    ggplot2::geom_hline(yintercept = ext$yintercept, linetype = "dashed")

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
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_color_viridis_d
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 facet_grid
#'

plot_univariate <- function(object, vars = NULL, samples = "both", logscale = TRUE,
                            grid = FALSE, sample.facet = FALSE,
                            nrow = NULL, tol = 10e-3, ...) {

  if(is.null(vars)){
    vars <- names(object$df_numerator)
  }
  # Check object type
  check.object.type(object)

  # Create data object, and external object with density ratio and sample indicators
  data <- rbind(object$df_numerator, object$df_denominator)
  ext <- data.frame(dr = predict(object, newdata = data, ...),
                    sample = c(rep("numerator", nrow(object$df_numerator)),
                               rep("denominator", nrow(object$df_denominator))))
  # Check variable names
  check.var.names(vars, data)

  # Filter appropriate subset
  samples <- match.arg(samples, c("both", "numerator", "denominator"))
  if(samples != "both"){
    data <- subset(data, ext$sample == samples)
    ext  <- subset(ext, ext$sample == samples)
  }


  # If logscale is TRUE, transform density ratio estimates
  ext <- check.logscale(ext, logscale, tol)

  # Set y axis label
  if(logscale) {
    y_lab <- "Log (Density Ratio)"
  }  else {
    y_lab <- "Density Ratio"
  }

  # Plot either individual plots in a list, or a grid of individual plots
  if(!grid){

  plot <- lapply(vars, function(var) create_univariate_plot(data, ext, var, y_lab, sample.facet))

  } else {

    values <- data[, vars] |> unlist(use.names = FALSE)
    variable <- rep(vars, each = length(values)/length(vars))
    dr <- rep(ext$dr, length(vars))
    sample <- rep(ext$sample, length(vars))
    ext <- data.frame(values = values, variable = variable,
                      dr = dr, sample = sample, yintercept = rep(ext$yintercept, length(vars)))

    # Maximum scale for y
    y_max <- max(1,ext$dr)
    y_min <- min(-1, ext$dr)

    plot <- ggplot2::ggplot(ext) +
      ggplot2::geom_point(ggplot2::aes(x = values, y = dr, col = sample),
                 alpha = .6) +
      ggplot2::theme_bw() +
      ggplot2::labs(title = "Scatter plot of individual values and density ratio",
           color = "Sample",
           y = "Density ratio") +
      ggplot2::scale_color_viridis_d(option = "cividis",
                            breaks = c("numerator", "denominator"),
                            labels = c("Numerator", "Denominator")) +
      ggplot2::scale_y_continuous(limits = c(y_min, y_max))

    if(sample.facet){

      plot <- plot +
        ggplot2::facet_grid(rows = ggplot2::vars(sample),
                            cols = ggplot2::vars(variable),
                   scales = "free_x")

      } else {
        plot <- plot +
          ggplot2::facet_wrap(ggplot2::vars(variable), scales = "free_x", nrow = nrow)
      }

  plot <- plot +
    ggplot2::geom_hline(yintercept = ext$yintercept, linetype = "dashed")
  }

  return(plot)
}


#' Bivariate plot
#'
#' @inheritParams create_univariate_plot
#' @param show.sample Logical indicating whether to give different shapes to observations, depending on the sample they come from (numerator or denominator). Defaults to FALSE.
#' @param vars Character vector of variable names to be plotted.
#' @param logscale Logical indicating whether the density ratio should be plotted in log scale. Defaults to TRUE.
#'
#' @return Bivariate plot
#'
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_colour_gradient2
#' @importFrom ggplot2 scale_shape_manual
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 labs
#'
create_bivariate_plot  <- function(data, ext, vars, logscale, show.sample){


  dr_max <- ifelse(logscale, max(2, ext$dr), max(exp(2), ext$dr))
  dr_min <- ifelse(logscale, min(-2, ext$dr), min(exp(-2), ext$dr))

  plot <-
    ggplot2::ggplot(data, mapping = ggplot2::aes(x = .data[[vars[1]]], y = .data[[vars[2]]])) +
    ggplot2::geom_point(ggplot2::aes(colour = ext$dr, shape = show.sample),
               alpha = 1,
               size = 2.0) +
    ggplot2::scale_colour_gradient2(low = "#00204DFF",
                           high = "#7D0000",
                           mid = "lightyellow",
                           midpoint = 0,
                           limits = c(dr_min, dr_max)) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Scatter plot, with density ratio mapped to colour",
         colour = "Log (Density ratio)") +
    ggplot2::scale_shape_manual(values = c(21, 24))

  return(plot)
}

#' Densityratio in bidimensional plot
#'
#' Plots a scatterplot of two variables, with densityratio mapped to the colour scale.
#'
#' @inheritParams plot_univariate
#' @param logscale Logical indicating whether to plot the density ratio
#'   estimates on a log scale. Default is TRUE.
#' @param show.sample Logical indicating whether to give different shapes to observations, depending on the sample they come from (numerator or denominator). Defaults to FALSE.
#' @param vars1 Character vector of variable names to be plotted on X axis
#' @param vars2 Character vector of variable names to be plotted on Y axis
#'
#' @return Scatter plot of two variables, with density ratio mapped to the colour scale.
#' @export
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 scale_colour_gradient2
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_shape_manual
#' @importFrom ggplot2 ggplotGrob
#' @importFrom grid grid.draw
#' @importFrom grid grid.newpage
#' @importFrom grid nullGrob
#'
plot_bivariate <- function(object, vars1, vars2 = NULL, samples = "both",
                           grid = FALSE, logscale = TRUE, show.sample = NULL,
                           tol = 10e-3, ...) {

  if(length(vars1) == 1 & grid == TRUE){
    grid <- FALSE
    warning("A grid cannot be displayed for one plot. Defaulting to grid = FALSE")
  }

  if(length(vars1) == 1 & is.null(vars2)){
    stop("For a bivariate plot, two variables are required as input. Please specify more variables in var1, or any variable in var2")
  }
    # Check object type
  check.object.type(object)

  # Create data object and estimate density ratio
  data <- rbind(object$df_numerator, object$df_denominator)
  ext <- data.frame(dr = predict(object, newdata = data, ...),
                    sample = c(rep("numerator", nrow(object$df_numerator)),
                               rep("denominator", nrow(object$df_denominator))))

  # Check if logscale is TRUE, then change ext
  ext <- check.logscale(ext, logscale, tol)

  # Assign correct labels depending on logscale
  if(logscale){
    colour_label <- "Log (Density Ratio)"
  } else {
    colour_label <- "Density Ratio"
  }


  # Select data
  samples <- match.arg(samples, c("both", "numerator", "denominator"))
  if(samples != "both"){
    data <- subset(data, ext$sample == samples)
    ext  <- subset(ext, ext$sample == samples)
  }

  # Create variable combinations
  var_combinations <- check.var.combinations(data, vars1, vars2)

  if(!grid){

  plot <- lapply(var_combinations, function(vars) create_bivariate_plot(data, ext, vars, logscale, show.sample))
  return(plot)

  } else {

    ext <- data.frame(data, dr = ext$dr, sample = ext$sample)
    datlist <- lapply(var_combinations, \(x) data.frame(values.x = ext[,x[1]],
                                              values.y = ext[,x[2]],
                                              xvar = rep(x[1], nrow(ext)),
                                              yvar = rep(x[2], nrow(ext)),
                                              sample = ext[, "sample"],
                                              dr = ext[, "dr"]))

    plot_data <- do.call(datlist, what = rbind)

    dr_max <- max(1, ext$dr)
    dr_min <- min(-1, ext$dr)

    plot <-
      ggplot2::ggplot(plot_data, mapping = ggplot2::aes(x = values.x, y = values.y,
                                      shape = show.sample)) +
      ggplot2::geom_point(ggplot2::aes(colour = dr),
                 size = 2.0) +
      ggplot2::facet_grid(rows = ggplot2::vars(yvar), cols = ggplot2::vars(xvar), scales = "free",
                 switch = "both") +
      ggplot2::scale_colour_gradient2(low = "#00204DFF",
                             high = "#7D0000",
                             mid = "lightyellow",
                             midpoint = 0,
                             limits = c(dr_min, dr_max)) +
      ggplot2::scale_y_continuous(position = "left") +
      ggplot2::scale_x_continuous(position = "bottom") +
      ggplot2::theme_bw() +
      ggplot2::theme(strip.placement = "outside") +
      ggplot2::labs(title = "Scatter plots, with density ratio mapped to colour",
          x = NULL,
          y = NULL,
          colour = colour_label,
          shape = show.sample) +
      ggplot2::scale_shape_manual(values = c(21, 24))

  # Erase upper diagonal
  ## Create plot into a grob
  grob <- ggplot2::ggplotGrob(plot)

  ## Create name of empty panels in the upper diagonal
  empty_panels <- expand.grid(seq(1:length(vars1)), seq(1:length(vars1))) |>
    subset(Var2 > Var1)
  empty_panels$panel <- paste0("panel-", empty_panels$Var1, "-", empty_panels$Var2)

  # Delete panels in upper diagonal, based in their index
  idx <- which(grob$layout$name %in% empty_panels$panel)
  for (i in idx) grob$grobs[[i]] <- grid::nullGrob()

  out <- grob

  grid::grid.newpage()
  grid::grid.draw(out)

  }
}
