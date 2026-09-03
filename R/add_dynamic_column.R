#' @import dplyr
#' @import tidyr
#' @export
add_dynamic_column <- function(df, param, col_name) {
  df %>% mutate(!!sym(col_name) := param[[col_name]])
}
