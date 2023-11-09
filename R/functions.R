#' descriptive_stats
#'
#' @param data
#'
#' @return A data.frame/tibble.

descriptive_stats <- function(data) {
  data %>%
    dplyr::group_by(metabolite) %>%
    dplyr::summarise(dplyr::across(value, list(mean = mean, sd = sd))) %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ round(.x, digits = 1)))
}

#' plot_distributions
#'
#' @param data
#'
#' @return A plot object.

plot_distributions <- function(data) {
  metabolite_distribution_plot <- ggplot2::ggplot(data, aes(x = value)) +
    geom_histogram() +
    facet_wrap(vars(metabolite), scales = "free")
  metabolite_distribution_plot
}

#' column_to_snakecase
#'
#' @param data Data with the string columns
#' @param cols The column you wish to convert to snakecase
#'
#' @return A data frame
column_values_to_snakecase <- function(data, cols) {
  data %>%
    dplyr::mutate(dplyr::across({{ cols }}, snakecase::to_snake_case))
}

#' metabolites_to_wider
#'
#' @param data
#'
#' @return A dataframe in long format
metabolites_to_wider <- function(data) {
  data %>%
    tidyr::pivot_wider(
      names_from = metabolite,
      values_from = value,
      values_fn = mean,
      names_prefix = "metabolite_"
    )
}

#' A transformation recipe to pre-process the data
#'
#' @param data The lipidomics dataset
#' @param metabolite_variable
#'
#' @return
create_recipe_spec <- function(data, metabolite_variable) {
  recipes::recipe(data) %>%
    recipes::update_role({{ metabolite_variable }}, age, gender,
      new_role = "predictor"
    ) %>%
    recipes::update_role(class, new_role = "outcome") %>%
    recipes::step_normalize(tidyselect::starts_with("metabolite_"))
}
