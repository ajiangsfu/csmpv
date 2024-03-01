#' A function rather aimed at developers
#' @import ggpubr
#' @import ggplot2
#' @importFrom stats lm
#' @importFrom stats coef
#' @noRd

validate_continuous = function(yhat, yobs){
  compdat = data.frame(cbind(yhat, yobs))
  #rmse = sqrt(mean((yhat - yobs)**2, na.rm = T))
  # Calculate Mean Squared Error (MSE)
  mse = mean((yhat - yobs)^2, na.rm = T)
  f1 = lm(yobs ~ yhat, data = compdat)
  p1 <- ggplot(data = compdat, aes(x = yobs, y = yhat)) +
    geom_point(color = "grey", size = 0.5) +
    geom_smooth(method = "lm", se = TRUE,  fill = "lightgray") +
    ylab("Actual values") +
    xlab("Predicted values") +
    ggtitle(paste0("Mean Squared Error (MSE): ", sprintf("%.4f", mse))) +
    #geom_abline(intercept = 0, slope = 1) +  ## removed the idential line
    geom_text(aes(x = min(yobs, na.rm = T), y = max(yhat, na.rm = T),
                  label = paste0("y = ", sprintf("%.2f", coef(f1)[1]), "+", sprintf("%.2f", coef(f1)[2]), "x")),
              hjust = 0, vjust = 1) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 6))+
    theme(plot.title = element_text(size = 11)) 
  
  residuals <- compdat$yobs - compdat$yhat
  
  p2 <- ggplot(data = compdat, aes(x = yhat, y = residuals)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    labs(x = "Predicted Values", y = "Residuals", title = "Residuals vs. Predicted Values") +
    theme(plot.title = element_text(size = 11)) 
  
  
  p3 <- ggplot(compdat, aes(sample = residuals)) +
    geom_qq() +
    geom_qq_line() +
    labs(title = "Q-Q Plot of Residuals",
         x = "Theoretical Quantiles",
         y = "Sample Quantiles") +
    theme(plot.title = element_text(size = 11)) 


  p4 =  ggplot(compdat, aes(x = residuals)) +
    geom_histogram(binwidth = 0.5, fill = "lightblue", color = "black") +
    labs(title = "Histogram of Residuals",
         x = "Residuals",
         y = "Frequency") +
    theme(plot.title = element_text(size = 11)) 
  
  # Use ggarrange to arrange the plots
  pall <- ggarrange(p1, p2, p3,p4, ncol = 2, nrow = 2)
  return(pall)
}







