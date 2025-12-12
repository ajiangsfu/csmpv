#' Create Consistent Model Matrix for XGBoost Predictions
#' 
#' This internal function creates a model matrix from new data that matches the 
#' encoding used during model training. It ensures that factor levels, contrasts, 
#' and column ordering are consistent between training and prediction.
#' 
#' @param newdata A data frame containing the new data for prediction.
#' @param stored_formula The formula object stored from the training phase.
#' @param factor_info A list containing factor information from training, 
#'   including xlevels and contrasts.
#' @param stored_feature_names Character vector of feature names in the order 
#'   used during training.
#'   
#' @return A numeric matrix with proper memory alignment for XGBoost, with 
#'   columns in the same order as the training data.
#'   
#' @details This function is critical for ensuring that predictions are made 
#'   using the same feature encoding as training. It handles:
#'   \itemize{
#'     \item Factor level consistency
#'     \item Missing factor levels in new data
#'     \item Column ordering
#'     \item Memory alignment for XGBoost compatibility
#'   }
#'   
#' @keywords internal
#' @noRd
#' 
#' @importFrom stats terms
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats .getXlevels
create_consistent_matrix <- function(newdata, stored_formula, factor_info, stored_feature_names) {
  
  # Validate inputs
  if (is.null(newdata) || !is.data.frame(newdata) || nrow(newdata) == 0) {
    stop("newdata must be a non-empty data frame")
  }
  
  if (is.null(stored_formula)) {
    stop("stored_formula cannot be NULL")
  }
  
  if (is.null(stored_feature_names) || length(stored_feature_names) == 0) {
    stop("stored_feature_names cannot be NULL or empty")
  }
  
  # Check if all original features exist in newdata
  original_features <- all.vars(stored_formula)
  
  # Remove response variable if present (formula with response has length 2)
  if (length(stored_formula) == 3) {
    original_features <- original_features[-1]  # Remove response variable
  }
  
  missing_features <- setdiff(original_features, colnames(newdata))
  if (length(missing_features) > 0) {
    stop(paste("Missing features in newdata:", paste(missing_features, collapse = ", ")))
  }
  
  # Create model matrix with same encoding as training
  if (!is.null(factor_info) && !is.null(factor_info$xlevels)) {
    # Use stored factor levels and contrasts
    tt <- terms(stored_formula, data = newdata, xlev = factor_info$xlevels)
    mf <- model.frame(tt, data = newdata, xlev = factor_info$xlevels)
    
    # Use stored contrasts if available
    if (!is.null(factor_info$contrasts)) {
      new_matrix <- model.matrix(tt, mf, contrasts.arg = factor_info$contrasts)
    } else {
      new_matrix <- model.matrix(tt, mf)
    }
  } else {
    # Fallback: create model matrix without stored factor info
    new_matrix <- model.matrix(stored_formula, data = newdata)
  }
  
  # Remove intercept column if present
  intercept_col <- which(colnames(new_matrix) == "(Intercept)")
  if (length(intercept_col) > 0) {
    new_matrix <- new_matrix[, -intercept_col, drop = FALSE]
  }
  
  # Ensure column order matches training
  if (!all(stored_feature_names %in% colnames(new_matrix))) {
    # Some columns might be missing (e.g., factor levels not in new data)
    missing_cols <- setdiff(stored_feature_names, colnames(new_matrix))
    
    if (length(missing_cols) > 0) {
      # Create a matrix with all expected columns, initialized to 0
      aligned_matrix <- matrix(0, nrow = nrow(newdata), ncol = length(stored_feature_names))
      colnames(aligned_matrix) <- stored_feature_names
      
      # Fill in available columns
      common_cols <- intersect(stored_feature_names, colnames(new_matrix))
      if (length(common_cols) > 0) {
        aligned_matrix[, common_cols] <- new_matrix[, common_cols]
      }
      
      new_matrix <- aligned_matrix
    }
  } else {
    # Reorder to match training
    new_matrix <- new_matrix[, stored_feature_names, drop = FALSE]
  }
  
  # CRAN-safe conversion: create new matrix with explicit memory allocation
  # This avoids pointer misalignment issues in XGBoost 2.x
  aligned_matrix <- matrix(
    as.numeric(new_matrix),
    nrow = nrow(new_matrix),
    ncol = ncol(new_matrix)
  )
  
  # Ensure numeric type and proper storage mode
  storage.mode(aligned_matrix) <- "double"
  
  # Preserve column names for debugging
  colnames(aligned_matrix) <- colnames(new_matrix)
  
  return(aligned_matrix)
}


