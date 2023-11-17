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
plot_univariate <- function(object, vars, sample = "both") {

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


  plots <- list()

  plot_onevariable <- function(var, shape = "sample"){
    plot <-
      ggplot(data, aes(x = .data[[var]], y = dr )) +
      geom_point(aes(col = dr, shape = sample)) +
      theme_bw() +
      labs(y = "Density ratio") +
      labs(title = "Scatter plot of individual values and density ratio",
           shape = "Sample") +
      scale_colour_viridis_c(option = "B", name ="Density ratio")  +
      scale_shape_manual(values = c(16, 3))


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
plot_bivariate <- function(object, var.x,var.y, sample = "both", show.samples = TRUE) {

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
  return(plots)
}
