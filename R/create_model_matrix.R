#' Helper Function for Training Data Matrix Creation
#' 
#' Creates a model matrix from training data with proper factor handling
#' and stores factor information for consistent prediction.
#' 
#' @param data Training data frame
#' @param features Character vector of feature names
#' 
#' @return A list containing:
#'   \item{matrix}{Numeric matrix for XGBoost training}
#'   \item{factor_info}{List of factor information (xlevels, contrasts, assign)}
#'   \item{formula}{The formula used for model matrix creation}
#'   \item{feature_names}{Column names of the matrix}
#'   
#' @keywords internal
#' @noRd
create_model_matrix <- function(data, features) {
  # Validate inputs
  if (!all(features %in% colnames(data))) {
    missing <- setdiff(features, colnames(data))
    stop(paste("Missing features in data:", paste(missing, collapse = ", ")))
  }
  
  # Create formula
  formula <- as.formula(paste("~", paste(features, collapse = " + ")))
  
  # Create model matrix
  mm <- model.matrix(formula, data = data)
  
  # Store factor information before modifying
  factor_info <- list(
    xlevels = .getXlevels(terms(formula), data),
    contrasts = attr(mm, "contrasts"),
    assign = attr(mm, "assign")
  )
  
  # Remove intercept
  feature_matrix <- mm[, -1, drop = FALSE]
  feature_names <- colnames(feature_matrix)
  
  # CRAN-safe conversion to properly aligned matrix
  aligned_matrix <- matrix(
    as.numeric(feature_matrix),
    nrow = nrow(feature_matrix),
    ncol = ncol(feature_matrix)
  )
  colnames(aligned_matrix) <- feature_names
  storage.mode(aligned_matrix) <- "double"
  
  return(list(
    matrix = aligned_matrix,
    factor_info = factor_info,
    formula = formula,
    feature_names = feature_names
  ))
}
