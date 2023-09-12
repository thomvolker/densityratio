
plot.ulsif <- function(object, type = "univariate", var1 = NULL, var2, ...) {

  # Data handling
  data <- rbind(object$df_numerator, object$df_denominator)
  data$dr <- predict(object, newdata = data)

  # Univariate plot

  if(type == "univariate") {
    if(is.null(var1)){
      cor <- cor(data)[1:ncol(data)-1,ncol(data)]
      var <- names(cor[which.max(cor)])
    }
    else {
      var <- var1
    }

    plot <-
      ggplot(data, aes(x = .data[[var]], y = .data[["dr"]])) +
      geom_point(aes(col = dr)) +
      theme_bw() +
      labs(y = "Density ratio") +
      labs(title = "Scatter plot of individual values and density ratio") +
      scale_colour_viridis_c(option = "B", name ="Density ratio")

  }

  # Histogram
  if(type == "histogram") {
    plot <-
      ggplot(data, aes(x = dr)) +
      geom_histogram(fill = "grey25", color = "black",
                     bins = max(data$dr)) +
      scale_x_continuous(breaks = seq(round(min(data$dr), 0),round(max(data$dr), 0))) +
      theme_bw()  +
      labs(
        x = "Density ratio",
        y = "Count"
      ) +
      labs(title = "Distribution of density ratio estimates")
  }

  # Bivariate
  if(type == "bivariate") {
    plot <-
      ggplot(data, mapping = aes(x = .data[[var1]], y = .data[[var2]])) +
      geom_point(aes(col = dr)) +
      scale_colour_viridis_c(option = "B", name ="Density ratio") +
      theme_bw() +
      labs(x = var1,
           y = var2) +
      labs(title = "Density ratio estimates for combinations of values")

  }

  return(plot)
}
