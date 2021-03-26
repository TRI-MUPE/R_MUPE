# Linear instantiation of the Minimum Unbiased Percent Error technique (MUPE) for
# multiplicative error models, which utilizes Iteratively Re-weighted Least Squares (IRLS)
# with weights equal to the squared inverse predictions from the prior iteration.
# Usage Example:
#   mupe = mupe_linear(formula_str = "y ~ 0 + x", data = df)
#     - 'formula_str' must be a character string that resembles an R 'lm' formula object
#       (default is no intercept, i.e. a simple factor model)
#     - 'data' must be a dataframe containing the variables listed in formula_str
# Returns a list containing a standard R 'lm' object and the number of iterations.
# 
mupe_linear = function(formula_str = "y ~ 0 + x", data) {
  f = as.formula(formula_str)           # convert string to R formula object
  model = lm(f, data)                   # 1st iteration (Ordinary Least Squares)
  conv = 1.0;  i = 1                    # initialize convergence and counter
  while (conv > 1e-5) {
    wt = 1 / model$fitted^2             # calculate weights
    pbeta = model$coef                  # solution of prior iteration
    model = lm(f, data, weights=wt)     # weighted least squares
    beta = model$coef                   # solution of current iteration
    conv = max(abs((beta-pbeta)/beta))  # maximum fractional change in any parameter
    i = i + 1;  if (i == 200) break     # force stop, if necessary
  }
  return(list(model=model, mupe_iters=i))
}