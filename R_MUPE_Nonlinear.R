# General nonlinear instantiation of the Minimum Unbiased Percent Error technique (MUPE)
# for multiplicative error models, which utilizes Iteratively Re-weighted Least Squares
# (IRLS) with weights equal to the squared inverse predictions from the prior iteration.
# Utilizes the MINPACK library implementation of the Levenberg-Marquardt algorithm. To
# install this, run 'install.packages("minpack.lm")' from the R console.
# Usage Example:
#   mlist = mupe_nonlinear(formula_str = "y ~ a * x1^b * c^x2", data = df, 
#                          start = c(a=10, b=1, c=1))
#     - 'formula_str' must be a character string that resembles an R 'nls' formula object
#     - 'data' must be a dataframe containing the variables listed in formula_str
#     - 'start' is a named numeric vector of initial guesses (parameter names must 
#       match those in formula_str). Whenever possible, provide values of the correct
#       sign and order of magnitude. For log-linear model forms, use the LOLS or PING
#       solution as the initial guess.
# Returns a list containing a standard R 'nls' object and accompanying details.
# 
mupe_nonlinear = function(formula_str, data, start) {
  library(minpack.lm)                 # load Levenberg-Marquardt algorithm
  f = as.formula(formula_str)         # convert string to R formula object
  wt = rep(1, nrow(data))             # initialize weights
  pbeta = start;  conv = 1.0;  i = 0  # initialize other variables
  while (conv > 1e-5) {
    model = suppressWarnings(nlsLM(f, data, pbeta, weights=wt, control=list(maxiter=10)))
    wt = 1 / model$m$pred()^2           # calculate weights
    beta = model$m$getAllPars()         # solution of current iteration
    conv = max(abs((beta-pbeta)/beta))  # maximum fractional change in any parameter
    pbeta = beta                        # reset prior beta
    i = i + 1;  if (i == 200) break     # force stop, if necessary
  }
  return(list(model=model, start=start, mupe_iters=i))
}