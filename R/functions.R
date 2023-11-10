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

#' Create a workflow object of the model and the transformations
#'
#' @param model_specs the model specifications
#' @param recipe_specs the recipe specifications
#'
#' @return a workflow object
create_model_workflow <- function(model_specs, recipe_specs) {
  workflows::workflow() %>%
    workflows::add_model(model_specs) %>%
    workflows::add_recipe(recipe_specs)
}

#' Create a tidy output of the model results
#'
#' @param workflow_fitted_model The model workflow opbject that has been fitted
#'
#' @return A dataframe
tidy_model_output <- function(workflow_fitted_model) {
  workflow_fitted_model %>%
    workflows::extract_fit_parsnip() %>%
    broom::tidy(exponentiate = TRUE)
}

#' Convert the long form data set into a list of wide form dataframes
#'
#' @param data Lipidomics data set
#'
#' @return A list of dataframes
split_by_metabolite <- function(data) {
  data %>%
    column_values_to_snakecase(metabolite) %>%
    dplyr::group_split(metabolite) %>%
    purrr::map(metabolites_to_wider)
}

#' Generate the results of the model
#'
#' @param data lipidomics data set
#'
#' @return A dataframe
generate_model_results <- function(data) {
    create_model_workflow(
        parsnip::logistic_reg() %>%
            parsnip::set_engine("glm"),
        data %>%
            create_recipe_spec(tidyselect::starts_with("metabolite_"))
    ) %>%
        parsnip::fit(data) %>%
        tidy_model_output()
}

#' The data frame with the model results
#'
#' @param model_results the data frame with the model results
#' @param data the original unprocessed lipidomics data set
#'
#' @return A data frame
add_original_metabolite_names <- function(model_results, data) {
    data %>%
        dplyr::select(metabolite) %>%
        dplyr::mutate(term = metabolite) %>%
        column_values_to_snakecase(term) %>%
        dplyr::mutate(term = stringr::str_c("metabolite_", term)) %>%
        dplyr::distinct(term, metabolite) %>%
        dplyr::right_join(model_results, by = "term")
}

#' Calculate the estimates for the model for each metabolite.
#'
#' @param data The lipidomics dataset.
#'
#' @return A data frame.
calculate_estimates <- function(data) {
    data %>%
        column_values_to_snakecase(metabolite) %>%
        dplyr::group_split(metabolite) %>%
        purrr::map(metabolites_to_wider) %>%
        purrr::map(generate_model_results) %>%
        purrr::list_rbind() %>%
        dplyr::filter(stringr::str_detect(term, "metabolite_")) %>%
        add_original_metabolite_names(data)
}

