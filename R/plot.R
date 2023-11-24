plot.histogram <- function(object, sample = "both", logscale = FALSE, binwidth = NULL) {

  check.object.type(object)

  data <- rbind(object$df_numerator, object$df_denominator)

  check.overriden.names(vars)

  data$dr <- predict(object, newdata = data)

  if(logscale){
  data$dr[data$dr < 0] <- 10e-8
  data$dr <- log(data$dr)
  warning("Negative estimated density ratios converted to 10e-0.8 before applying logarithmic transformation")
  }

  obsclass <- rep(c("numerator", "denominator"),
                  c(nrow(object$df_numerator), nrow(object$df_denominator)))

  data$sample <- obsclass

  obsselect <- match.arg(sample, c("both", "numerator", "denominator"))

  if (obsselect != "both") {
    data <- filter(data, obsclass == obsselect)
  }

  plot <-
    ggplot(data, aes(x = dr)) +
    geom_histogram(aes(fill = sample),
                   alpha = .75,
                   color = "black",
                   binwidth = if (!is.null(binwidth)) binwidth else NULL,
                   position = position_dodge2(preserve = "single",
                                              padding = 0.2)) +
    scale_fill_manual(values = c("firebrick", "steelblue")) +
    theme_bw()  +
    labs(
      x = "Density ratio",
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
plot.ulsif <- function(object, sample = "both", logscale = FALSE, binwidth = NULL) {
  plot.histogram(object, sample = sample, logscale = logscale, binwidth = binwidth)
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
plot.kliep <- function(object, sample = "both", logscale = FALSE, binwidth = NULL) {
  plot.histogram(object, sample = sample, logscale = logscale, binwidth = binwidth)
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
plot_univariate <- function(object, vars, sample = "both", logscale = TRUE) {

  check.object.type(object)

  # Data handling
  data <- rbind(object$df_numerator, object$df_denominator)

  check.overriden.names(data)
  check.var.names(vars, data)

  data$dr <- predict(object, newdata = data)

  obsclass <- rep(c("numerator", "denominator"),
                  c(nrow(object$df_numerator), nrow(object$df_denominator)))

  obsselect <- match.arg(sample, c("both", "numerator", "denominator"))

  if (obsselect != "both") {
    data <- filter(data, obsclass == obsselect)
  }

  if (logscale) {
    if(any(data$dr < 0)){
      warning("Negative estimated density ratios converted to 10e-8 before applying logarithmic transformation")
      data$dr[data$dr < 0] <- 10e-8
    }

    data$dr <- log(data$dr)
    y_lab <- "Log(Density Ratio)"

  } else {
    y_lab <- "Density Ratio"
  }

  plots <- list()

  plot_onevariable <- function(var, shape = "sample"){
    plot <-
      ggplot(data, aes(x = .data[[var]], y = dr )) +
      geom_point(aes(col = dr, shape = sample)) +
      theme_bw() +
      labs(title = "Scatter plot of individual values and density ratio",
           shape = "Sample",
           y = y_lab) +
      geom_hline(yintercept = 0, linetype = "dashed")+
      scale_colour_viridis_c(option = "B", name ="Density ratio")  +
      scale_shape_manual(values = c(16, 3))  +
      scale_y_continuous(breaks = c(-15,-10,-5,0,1,2,3,4, 8))


    return(plot)
  }


  for(var in vars){
    plots[[var]] <- plot_onevariable(var)
  }
  return(plots)
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
plot_bivariate <- function(object, var.x,var.y, sample = "both", show.samples = TRUE, output = "assembled") {

  check.object.type(object)

  data <- rbind(object$df_numerator, object$df_denominator)

  check.overriden.names(data)

  vars <- c(var.x, var.y)

  check.var.names(vars, data)
  data$dr <- predict(object, newdata = data)

  obsclass <- rep(c("numerator", "denominator"),
                  c(nrow(object$df_numerator), nrow(object$df_denominator)))

  data$sample <- obsclass

  obsselect <- match.arg(sample, c("both", "numerator", "denominator"))

  if (obsselect != "both") {
    data <- filter(data, obsclass == obsselect)
  }

  plots <- list()



  var_combinations <- expand.grid(var.x, var.y)
  var_combinations <- as.data.frame(apply(var_combinations, 2,  as.character))
  var_combinations <- var_combinations %>% filter(Var1 != Var2)
  var_combinations <- as.matrix(var_combinations)

  # Write function to remove duplicate rows
  remove_duplicate_rows <- function(data_matrix) {
    # Convert matrix to data frame
    data_frame <- as.data.frame(data_matrix)

    # Sort elements within each row
    sorted_data_frame <- t(apply(data_frame, 1, function(row) sort(row)))

    # Convert back to matrix
    sorted_matrix <- as.matrix(sorted_data_frame)

    # Remove duplicate rows
    unique_matrix <- unique(sorted_matrix)

    return(unique_matrix)
  }

  var_combinations <- remove_duplicate_rows(var_combinations)

  # helper function to plot two variables
  plot_twovariables <- function(data, vars, showsamples = show.samples){
    plot <-
    ggplot(data, mapping = aes(x = .data[[vars[1]]], y = .data[[vars[2]]])) +
    geom_point(aes(colour = dr, fill = dr, shape = if (showsamples) sample else NULL),
               alpha = 0.5) +
    scale_fill_viridis_c(option = "B", name ="Density ratio") +
    scale_colour_viridis_c(option = "B", name ="Density ratio") +
    theme_bw() +
    labs(title = "Density ratio estimates for combinations of values",
         shape = "Sample") +
    scale_shape_manual(values = c(21, 24))

  return(plot)
  }

  for(i in 1:nrow(var_combinations)){
    plots[[i]] <- plot_twovariables(data = data, vars = var_combinations[i,])
  }

  if(output == "assembled"){
  plots_assembly <- patchwork::wrap_plots(plots, guides = "collect", byrow = TRUE,ncol = length(var.x), nrow = length(var.y)) & labs(title = NULL)
  plots_assembly <- plots_assembly + plot_annotation(title = "Density ratio estimates for combinations of values")
  return(plots_assembly)
  } else {
  return(plots)
  }
}

plot_bivariate_heatmap <- function(object, var.x,var.y, sample = "both", show.samples = TRUE, output = "assembled") {

  check.object.type(object)

  data <- rbind(object$df_numerator, object$df_denominator)

  check.overriden.names(data)

  vars <- c(var.x, var.y)

  check.var.names(vars, data)
  data$dr <- predict(object, newdata = data)

  obsclass <- rep(c("numerator", "denominator"),
                  c(nrow(object$df_numerator), nrow(object$df_denominator)))

  data$sample <- obsclass

  obsselect <- match.arg(sample, c("both", "numerator", "denominator"))

  if (obsselect != "both") {
    data <- filter(data, obsclass == obsselect)
  }

  plots <- list()

  var_combinations <- expand.grid(var.x, var.y)
  var_combinations <- as.data.frame(apply(var_combinations, 2,  as.character))
  var_combinations <- var_combinations %>% filter(Var1 != Var2)
  var_combinations <- as.matrix(var_combinations)


  # helper function to plot two variables
  plot_twovariables <- function(object, data, vars, showsamples = show.samples){

    grid_var1 <- seq(min(data[[vars[1]]]), max(data[[vars[1]]]))
    grid_var2 <- seq(min(data[[vars[2]]]), max(data[[vars[2]]]))

    grid_data <- expand.grid(grid_var1, grid_var2)

    colnames(grid_data) <- c(vars[[1]], vars[[2]])
    #browser()
    # Identify variables in data that are not in vars
    other_vars <- setdiff(names(data), vars)

    # Assign the mean value of each other variable to the corresponding variable in grid_data
    for (var in other_vars) {
      grid_data[[var]] <- mean(data[[var]])
    }

    grid_data <- as.data.frame(grid_data)

    grid_data <- grid_data[, names(data)] # reorder names


    grid_data$dr <- predict(object, newdata = grid_data)

    plot <-
      ggplot(grid_data, mapping = aes(x = .data[[vars[1]]], y = .data[[vars[2]]])) +
      geom_raster(aes(colour = dr, fill = dr, shape = if (showsamples) sample else NULL),
                  alpha = 0.5) +
      scale_fill_viridis_c(option = "B", name ="Density ratio") +
      scale_colour_viridis_c(option = "B", name ="Density ratio") +
      theme_bw() +
      labs(title = "Density ratio estimates for combinations of values",
           shape = "Sample") +
      scale_shape_manual(values = c(21, 24))
    return(plot)
  }
  for (i in 1:nrow(var_combinations)) {
    plots[[i]] <- plot_twovariables(object = object, data = data, vars = as.character(var_combinations[i,]))
  }

  if(output == "assembled"){
    plots_assembly <- patchwork::wrap_plots(plots, guides = "collect", byrow = TRUE,ncol = length(var.x), nrow = length(var.y))  & labs(title = NULL)
    plots_assembly <- plots_assembly + plot_annotation(title = "Density ratio estimates for combinations of values")
    return(plots_assembly)
  } else {
    return(plots)
  }
}
