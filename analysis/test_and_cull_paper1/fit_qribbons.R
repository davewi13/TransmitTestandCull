library(quantreg)

fit_qribbons <- function(df, xvar = "R0", yvar = "Time", lambda = 1) {
  
  taus <- c(0.025, 0.5, 0.975)
  
  fits <- lapply(taus, function(tau) {
    rqss(
      formula = reformulate(
        paste0("qss(", xvar, ", lambda=", lambda, ")"),
        response = yvar
      ),
      tau = tau,
      data = df
    )
  })
  
  xseq <- seq(
    min(df[[xvar]], na.rm = TRUE),
    max(df[[xvar]], na.rm = TRUE),
    length.out = 1000
  )
  
  newdata <- setNames(data.frame(xseq), xvar)
  
  preds <- lapply(
    fits,
    function(fit) predict(fit, newdata = newdata)
  )
  
  tibble(
    x = xseq,
    lwr = preds[[1]],
    fit = preds[[2]],
    upr = preds[[3]]
  )
}
