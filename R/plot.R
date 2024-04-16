
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
                         tol = 10e-3,
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
    # Converte negative values to tol
    negdr <- ext$dr < tol
    ext$dr[negdr] <- tol

    if(any(negdr)){
    warning(
      paste("Negative estimated density ratios for", sum(negdr), "observation(s) converted to", tol, "before applying logarithmic transformation"),
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
#' @return
#' @export
#'
#' @examples
plot.ulsif <- function(object, samples = "both", logscale = TRUE, binwidth = NULL,
                       bins = NULL, tol = 10e-3) {
  dr.histogram(object, samples = samples, logscale = logscale, binwidth = binwidth,
               bins = bins, tol = tol)
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
                       bins = NULL, tol = 10e-3) {
  dr.histogram(object, samples = samples, logscale = logscale, binwidth = binwidth,
               bins = bins, tol = tol)
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
    plot <- plot + ggplot2::facet_wrap(~ext$sample) +
      ggplot2::geom_hline(yintercept = ext$yintercept, linetype = "dashed")
  } else {
    plot <- plot + ggplot2::geom_hline(yintercept = ext$yintercept, linetype = "dashed")
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
plot_univariate <- function(object, vars = NULL, samples = "both", logscale = TRUE,
                            grid = FALSE, sample.facet = FALSE,
                            nrow = NULL, tol = 10e-3, ...) {

  if(is.null(vars)){
    vars <- names(object$df_numerator)
  }
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


  # Check if logscale is TRUE
  if (logscale) {
    negdr <- ext$dr < tol
    ext$dr[negdr] <- tol
    if(any(negdr)){
      warning(
        paste("Negative estimated density ratios for", sum(negdr), "observation(s) converted to", tol, "before applying logarithmic transformation"),
        call. = FALSE)
    }
    # For the plots
    ext$yintercept <- 0
    # Apply log transformation
    ext$dr <- log(ext$dr)
    y_lab <- "Log (Density Ratio)"
  }
     else {
    y_lab <- "Density Ratio"
    ext$yintercept <- 1
  }


  if(!grid){

  plot <- lapply(vars, function(var) create_univariate_plot(data, ext, var, y_lab, sample.facet))

  } else {

    values <- data[, vars] |> unlist(use.names = FALSE)
    variable <- rep(vars, each = length(values)/length(vars))
    dr <- rep(ext$dr, length(vars))
    sample <- rep(ext$sample, length(vars))
    data <- data.frame(values = values, variable = variable)
    ext <- data.frame(dr = dr, sample = sample, yintercept = rep(ext$yintercept, length(vars)))

    # Maximum scale for y
    y_max <- max(1,ext$dr)
    y_min <- min(-1, ext$dr)

    plot <- ggplot2::ggplot(data) +
      ggplot2::geom_point(ggplot2::aes(x = values, y = ext$dr, col = ext$sample),
                 alpha = .6) +
      ggplot2::theme_bw() +
      ggplot2::geom_hline(ggplot2::aes(yintercept = ext$yintercept), linetype = "dashed") +
      ggplot2::labs(title = "Scatter plot of individual values and density ratio",
           color = "Sample",
           y = "Density ratio") +
      ggplot2::scale_color_viridis_d(option = "cividis",
                            breaks = c("numerator", "denominator"),
                            labels = c("Numerator", "Denominator")) +
      ggplot2::scale_y_continuous(limits = c(y_min, y_max))

    if(sample.facet){

      plot <- plot +
        ggplot2::facet_grid(cols = vars(sample),
                   rows = vars(variable),
                   scales = "free_x")

      } else {
        plot <- plot +
          ggplot2::facet_wrap(~variable, scales = "free_x", nrow = nrow)
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
#' @inheritParams plot.univariate
#' @param logscale Logical indicating whether to plot the density ratio
#'   estimates on a log scale. Default is TRUE.
#'
#'
#' @return Scatter plot of two variables, with density ratio mapped to the colour scale.
#' @export
#'
#' @examples
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
  ext <- data.frame(dr = predict(object, data = data, ...),
                    sample = c(rep("numerator", nrow(object$df_numerator)),
                               rep("denominator", nrow(object$df_denominator))))

  # Determine if DR is shown in logscale (default) or not
  # Check if logscale is TRUE
  if (logscale) {
    negdr <- ext$dr < tol
    ext$dr[negdr] <- tol
    if(any(negdr)){
      warning(
        paste("Negative estimated density ratios for", sum(negdr), "observation(s) converted to", tol, "before applying logarithmic transformation"),
        call. = FALSE)
    }

    # Apply log transformation
    ext$dr <- log(ext$dr)
    colour_name <- "Log (Density Ratio)"
  } else {
    colour_name <- "Density ratio"
  }


  # Select data
  samples <- match.arg(samples, c("both", "numerator", "denominator"))
  if(samples != "both"){
    data <- subset(data, ext$sample == samples)
    ext  <- subset(ext, ext$sample == samples)
  }


  if (is.null(vars2)) {
    # Check variable names
    check.var.names(vars1, data)
    var_combinations <- expand.grid(vars1, vars1)
  } else {
    if(length(vars1) != length(vars2)){
      stop("The number of variables in vars1 and vars2 must be the same")
    }
    # Check variable names
    check.var.names(vars1, data)
    check.var.names(vars2, data)

    var_combinations <- expand.grid(vars1, vars2)
  }


  ## Remove duplicate combinations
  ## Start by sorting elements within each row
  ## This makes duplicate rows with different order of variables identical
  var_combinations <- t(apply(var_combinations, 1, sort))
  var_combinations <- unique(var_combinations) # retain unique rows only

  if(length(vars1) != 1){
  var_combinations <- as.data.frame(apply(var_combinations, 2,  as.character))
  } else {
  var_combinations <- as.data.frame(var_combinations)
  }

  names(var_combinations) <- c("Var1", "Var2")

  # Remove rows where both variables are the same
  var_combinations <- var_combinations |> subset(Var1 != Var2)
  var_combinations <- as.list(as.data.frame(t(var_combinations)))

  if(!grid){


  plot <- lapply(var_combinations, function(vars) create_bivariate_plot(data, ext, vars, logscale, show.sample))

  return(plot)
  }

  if (grid) {

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
      ggplot2::geom_point(aes(colour = dr),
                 size = 2.0) +
      ggplot2::facet_grid(rows = vars(yvar), cols = vars(xvar), scales = "free",
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
          colour = colour_name,
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
