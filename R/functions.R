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
