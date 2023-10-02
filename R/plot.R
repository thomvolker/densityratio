
plot.ulsif <- function(object, sample = "both") {

  data <- rbind(object$df_numerator, object$df_denominator, make.row.names = TRUE)
  data$dr <- predict(object, newdata = data)

  data$sample[1:nrow(object$df_numerator)] <- "Numerator"
  data$sample[(nrow(object$df_numerator) + 1):nrow(data)] <- "Denominator"

  if(sample == "numerator"){
    data <- data %>% filter(sample == "Numerator")
  }
  if(sample == "denominator"){
    data <- data %>% filter(sample == "Denominator")
  }
  if(sample == "both"){
    data <- data
  }

  plot <-
    ggplot(data, aes(x = dr)) +
    geom_histogram(aes(fill = sample), alpha = .75,
                   color = "black",
                   position = "identity",
                   bins = max(data$dr)) +
    scale_x_continuous(breaks = seq(round(min(data$dr), 0),round(max(data$dr), 0))) +
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
  data$dr <- predict(object, newdata = data)

  data$sample[1:nrow(object$df_numerator)] <- "Numerator"
  data$sample[(nrow(object$df_numerator) + 1):nrow(data)] <- "Denominator"

  if(sample == "numerator"){
    data <- data %>% filter(sample == "Numerator")
  }
  if(sample == "denominator"){
    data <- data %>% filter(sample == "Denominator")
  }
  if(sample == "both"){
    data <- data
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

plot_bivariate <- function(object, vars, sample = "both") {

  data <- rbind(object$df_numerator, object$df_denominator)
  data$dr <- predict(object, newdata = data)

  data$sample[1:nrow(object$df_numerator)] <- "Numerator"
  data$sample[(nrow(object$df_numerator) + 1):nrow(data)] <- "Denominator"

  if(sample == "numerator"){
    data <- data %>% filter(sample == "Numerator")
  }
  if(sample == "denominator"){
    data <- data %>% filter(sample == "Denominator")
  }
  if(sample == "both"){
    data <- data
  }


  plot <-
    ggplot(data, mapping = aes(x = .data[[vars[1]]], y = .data[[vars[2]]])) +
    geom_point(aes(col = dr, shape = sample)) +
    scale_colour_viridis_c(option = "B", name ="Density ratio") +
    theme_bw() +
    labs(title = "Density ratio estimates for combinations of values",
         shape = "Sample") +
    scale_shape_manual(values = c(16, 3))

  return(plot)
}
