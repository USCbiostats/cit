// Test 4 resampling engine for the causal inference test.
//
// For each column b of `perm`, build G* = Ghat + Gres[perm[, b], ] (fitted values
// from G ~ L + CG plus row-permuted residuals) and return the statistic for L in
//     T ~ A + G* + L   versus   T ~ A + G*
// where A holds the intercept and outcome covariates. The statistic is the
// F statistic (linear outcome) or the likelihood-ratio chi-square (binary outcome).
//
// All randomness lives in R: `perm` is generated with R's RNG, so set.seed()
// controls the results and this code is deterministic. Missing data are removed
// in R before this is called. A failed fit returns NA.

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// residual sum of squares of y on X; false if X is rank deficient
static bool rss_ls(const mat& X, const vec& y, double& rss) {
  vec b;
  bool ok = solve(b, X, y, solve_opts::no_approx);
  if (!ok || !b.is_finite()) return false;
  vec r = y - X * b;
  rss = dot(r, r);
  return true;
}

// binomial deviance of a logistic regression fitted by IRLS, following glm.fit:
// start at mu = (y + 0.5) / 2, stop when |dev - dev_old| / (|dev| + 0.1) < tol
static bool dev_logit(const mat& X, const vec& y, int maxiter, double tol, double& dev) {
  const double eps = 2.220446e-16;
  vec mu = (y + 0.5) / 2.0;
  vec eta = log(mu / (1.0 - mu));
  double dev_old = -2.0 * accu(y % log(mu) + (1.0 - y) % log(1.0 - mu));
  for (int it = 0; it < maxiter; ++it) {
    vec w = mu % (1.0 - mu);
    vec z = eta + (y - mu) / w;
    vec sw = sqrt(w);
    mat Xw = X.each_col() % sw;
    vec b;
    if (!solve(b, Xw, z % sw, solve_opts::no_approx) || !b.is_finite()) return false;
    eta = X * b;
    mu = 1.0 / (1.0 + exp(-eta));
    mu.clamp(eps, 1.0 - eps);
    dev = -2.0 * accu(y % log(mu) + (1.0 - y) % log(1.0 - mu));
    if (std::fabs(dev - dev_old) / (std::fabs(dev) + 0.1) < tol) return true;
    dev_old = dev;
  }
  return false;   // not converged
}

// [[Rcpp::export(.boot_test4)]]
Rcpp::NumericVector boot_test4(const arma::vec& y, const arma::mat& A, const arma::mat& L,
                               const arma::mat& Ghat, const arma::mat& Gres,
                               const arma::imat& perm, int family,
                               int maxiter = 25, double tol = 1e-8) {
  const uword n = y.n_elem, a = A.n_cols, m = Ghat.n_cols, k = L.n_cols;
  const uword B = perm.n_cols;
  const double df2 = (double) n - a - m - k;
  Rcpp::NumericVector out(B, NA_REAL);
  mat X1(n, a + m + k);
  X1.cols(0, a - 1) = A;
  X1.cols(a + m, a + m + k - 1) = L;

  for (uword b = 0; b < B; ++b) {
    if (b % 256 == 0) Rcpp::checkUserInterrupt();
    for (uword i = 0; i < n; ++i) {
      uword src = (uword) perm(i, b);          // 0-based row index
      for (uword j = 0; j < m; ++j) X1(i, a + j) = Ghat(i, j) + Gres(src, j);
    }
    mat X0 = X1.cols(0, a + m - 1);
    if (family == 0) {
      double rss0, rss1;
      if (rss_ls(X0, y, rss0) && rss_ls(X1, y, rss1) && rss1 > 0) {
        double ss = rss0 - rss1;
        if (ss < 0) ss = 0;                    // rounding
        out[b] = (ss / k) / (rss1 / df2);
      }
    } else {
      double dev0, dev1;
      if (dev_logit(X0, y, maxiter, tol, dev0) && dev_logit(X1, y, maxiter, tol, dev1)) {
        double chi = dev0 - dev1;
        out[b] = (chi < 0) ? 0.0 : chi;
      }
    }
  }
  return out;
}
