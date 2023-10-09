
plot.ulsif <- function(object, sample = "both", logscale = FALSE) {

  data <- rbind(object$df_numerator, object$df_denominator)

  if("dr" %in% names(data) | "sample" %in% names(data)) {
    stop("Variables in your dataframe cannot have name 'dr' or 'sample'. Please rename your variable(s)")
  }

  data$dr <- predict(object, newdata = data)

  if(logscale){
  data$dr <- log(data$dr)
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
    geom_histogram(aes(fill = sample), alpha = .75,
                   color = "black",
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


plot_univariate <- function(object, vars, sample = "both") {

  # Data handling
  data <- rbind(object$df_numerator, object$df_denominator)

  if("dr" %in% names(data)) {
    stop("Variables in your dataframe cannot have name 'dr'. Please rename your variable")
  }

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

plot_bivariate <- function(object, vars, sample = "both", show.samples = FALSE) {

  data <- rbind(object$df_numerator, object$df_denominator)

  if("dr" %in% names(data) | "sample" %in% names(data)) {
    stop("Variables in your dataframe cannot have name 'dr' or 'sample'. Please rename your variable(s)")
  }

  data$dr <- predict(object, newdata = data)

  obsclass <- rep(c("numerator", "denominator"),
                  c(nrow(object$df_numerator), nrow(object$df_denominator)))

  data$sample <- obsclass

  obsselect <- match.arg(sample, c("both", "numerator", "denominator"))

  if (obsselect != "both") {
    data <- filter(data, obsclass == obsselect)
  }

  plot <-
    ggplot(data, mapping = aes(x = .data[[vars[1]]], y = .data[[vars[2]]])) +
    geom_point(aes(col = dr, shape = if (show.samples) sample else NULL)) +
    scale_colour_viridis_c(option = "B", name ="Density ratio") +
    theme_bw() +
    labs(title = "Density ratio estimates for combinations of values",
         shape = "Sample") +
    scale_shape_manual(values = c(16, 3))

  return(plot)
}
