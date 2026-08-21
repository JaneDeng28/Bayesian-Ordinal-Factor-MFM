// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// =============================================================================
// C++ backend for the Bayesian ordinal factor MFM model
// =============================================================================
// This file contains the computationally intensive components called from
// ordinal_mfm.R. It is compiled at runtime with Rcpp::sourceCpp().
//
// Data convention:
//   Y: p x n integer matrix (p features, n cells), categories 1,...,5
//   B: p x q factor-loading matrix
// eta: q x n cell-specific latent factor scores
//
// Main R-exported functions:
//   ordinal_factor_mfm_mcmc(): one MCMC chain
//   dahl_partition_cpp():      posterior similarity and Dahl partition
//   to_ordinal():              continuous-to-ordinal utility
// =============================================================================
#include <RcppArmadillo.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>
#include <limits>

using namespace Rcpp;
using namespace arma;

// -----------------------------------------------------------------------------
// Truncated-normal random-number generators
// -----------------------------------------------------------------------------
// Separate algorithms are used for central, one-sided-tail, and finite
// intervals to avoid numerical problems and inefficient rejection in the tails.
double rtruncnorm_std_lower(double a) {
  if (!std::isfinite(a)) stop("Invalid lower truncation point.");

  if (a <= 0.0) {
    double z;
    do { z = R::rnorm(0.0, 1.0); } while (z <= a);
    return z;
  }

  const double alpha = 0.5 * (a + std::sqrt(a * a + 4.0));
  while (true) {
    const double z = a + R::rexp(1.0 / alpha);
    const double log_accept = -0.5 * (z - alpha) * (z - alpha);
    if (std::log(R::runif(0.0, 1.0)) <= log_accept) return z;
  }
}

double rtruncnorm_std(double a, double b) {
  if (std::isnan(a) || std::isnan(b) || !(a < b))
    stop("Invalid truncation interval.");

  const bool a_inf = !std::isfinite(a);
  const bool b_inf = !std::isfinite(b);

  if (a_inf && b_inf) return R::rnorm(0.0, 1.0);

  if (b_inf) return rtruncnorm_std_lower(a);

  if (a_inf) return -rtruncnorm_std_lower(-b);

  if (b < 0.0) return -rtruncnorm_std(-b, -a);

  const double width = b - a;

  if (width <= 2.0) {
    double mode = 0.0;
    if (a > 0.0) mode = a;
    const double mode_sq = mode * mode;

    while (true) {
      const double z = R::runif(a, b);
      const double log_accept = -0.5 * (z * z - mode_sq);
      if (std::log(R::runif(0.0, 1.0)) <= log_accept) return z;
    }
  }

  if (a > 0.0) {
    double z;
    do { z = rtruncnorm_std_lower(a); } while (z >= b);
    return z;
  }

  double z;
  do { z = R::rnorm(0.0, 1.0); } while (z <= a || z >= b);
  return z;
}

double rtruncnorm(double mu, double sd, double lo, double hi) {
  if (!std::isfinite(mu) || !std::isfinite(sd) || sd <= 0.0)
    stop("Invalid normal mean or standard deviation.");
  if (std::isnan(lo) || std::isnan(hi) || !(lo < hi))
    stop("Invalid truncation bounds.");

  const double a = (lo - mu) / sd;
  const double b = (hi - mu) / sd;
  const double z = rtruncnorm_std(a, b);
  double x = mu + sd * z;

  if (std::isfinite(lo) && x <= lo) x = std::nextafter(lo, hi);
  if (std::isfinite(hi) && x >= hi) x = std::nextafter(hi, lo);
  return x;
}

double rtruncnorm(double mu, double lo, double hi) {
  return rtruncnorm(mu, 1.0, lo, hi);
}

mat ensure_spd(const mat& A, const std::string& context,
               double rel_floor = 1e-10, double abs_floor = 1e-12) {
  if (A.n_rows == 0 || A.n_rows != A.n_cols)
    stop(context + ": matrix must be nonempty and square.");
  if (!A.is_finite())
    stop(context + ": matrix contains NA, NaN, or Inf.");

  const int d = A.n_rows;
  mat S = 0.5 * (A + A.t());
  mat L;
  if (chol(L, S, "lower")) return S;

  double scale = norm(S, "fro") / std::max(1.0, (double)d);
  scale = std::max(scale, mean(abs(S.diag())));
  scale = std::max(scale, 1e-12);

  double jitter = std::max(abs_floor, rel_floor * scale);
  for (int attempt = 0; attempt < 8; attempt++) {
    mat S_try = S + jitter * eye(d, d);
    if (chol(L, S_try, "lower")) return S_try;
    jitter *= 10.0;
  }

  vec eigval;
  mat eigvec;
  if (!eig_sym(eigval, eigvec, S))
    stop(context + ": eigendecomposition failed while repairing covariance matrix.");

  const double eig_floor = std::max(abs_floor, rel_floor * scale);
  for (uword i = 0; i < eigval.n_elem; i++)
    if (!std::isfinite(eigval(i)) || eigval(i) < eig_floor) eigval(i) = eig_floor;

  mat S_fixed = eigvec * diagmat(eigval) * eigvec.t();
  S_fixed = 0.5 * (S_fixed + S_fixed.t());
  if (!chol(L, S_fixed, "lower"))
    stop(context + ": covariance matrix could not be made positive definite.");
  return S_fixed;
}

// -----------------------------------------------------------------------------
// Stable positive-definite matrix and multivariate-distribution utilities
// -----------------------------------------------------------------------------
// Covariance matrices are symmetrized and, only when required, repaired with
// scale-aware diagonal jitter/eigenvalue flooring before Cholesky operations.
mat chol_spd_safe(const mat& A, const std::string& context) {
  mat S = ensure_spd(A, context);
  mat L;
  if (!chol(L, S, "lower"))
    stop(context + ": Cholesky decomposition failed after SPD repair.");
  return L;
}

// use Cholesky instead of inv_sympd
mat inv_spd_safe(const mat& A, const std::string& context) {
  mat L = chol_spd_safe(A, context);
  mat I = eye(L.n_rows, L.n_cols);
  mat L_inv = solve(trimatl(L), I);
  mat A_inv = L_inv.t() * L_inv;
  return 0.5 * (A_inv + A_inv.t());
}

double dmvnorm_log(const vec& x, const vec& mu, const mat& Sigma) {
  int K = x.n_elem;
  mat L = chol_spd_safe(Sigma, "dmvnorm_log Sigma");
  vec z = solve(trimatl(L), x - mu);
  double logdet = 2.0 * accu(log(L.diag()));
  return -0.5 * (K * log(2.0 * M_PI) + logdet + dot(z, z));
}

double dmvt_log(const vec& x, const vec& mu, const mat& Sigma, double nu) {
  int K = x.n_elem;
  if (!std::isfinite(nu) || nu <= 0.0) stop("dmvt_log: degrees of freedom must be positive.");
  mat L = chol_spd_safe(Sigma, "dmvt_log Sigma");
  vec z = solve(trimatl(L), x - mu);
  double logdet = 2.0 * accu(log(L.diag()));
  double quad = dot(z, z);
  double logc = lgamma((nu + K) / 2.0) - lgamma(nu / 2.0)
                - 0.5 * K * log(nu * M_PI) - 0.5 * logdet;
  return logc - 0.5 * (nu + K) * log(1.0 + quad / nu);
}

mat rwish(double nu, const mat& V) {
  int K = V.n_rows;
  if (!std::isfinite(nu) || nu <= K - 1.0)
    stop("rwish: degrees of freedom must exceed dimension minus one.");
  mat L = chol_spd_safe(V, "rwish scale V");
  mat A(K, K, fill::zeros);
  for (int i = 0; i < K; i++) {
    double df_i = nu - i;
    if (df_i <= 0.0) stop("rwish: invalid Bartlett chi-square degrees of freedom.");
    A(i, i) = sqrt(R::rchisq(df_i));
    for (int j = 0; j < i; j++) A(i, j) = R::rnorm(0.0, 1.0);
  }
  mat LA = L * A;
  mat W = LA * LA.t();
  return ensure_spd(W, "rwish sampled precision");
}

mat riwish(double nu, const mat& S) {
  int K = S.n_rows;
  if (!std::isfinite(nu) || nu <= K - 1.0)
    stop("riwish: degrees of freedom must exceed dimension minus one.");
  mat S_spd = ensure_spd(S, "riwish scale S");
  mat S_inv = inv_spd_safe(S_spd, "riwish inverse scale");
  mat W = rwish(nu, S_inv);
  return ensure_spd(inv_spd_safe(W, "riwish sampled Wishart precision"),
                    "riwish covariance draw");
}

vec rmvnorm(const vec& mu, const mat& Sigma) {
  int K = mu.n_elem;
  mat L = chol_spd_safe(Sigma, "rmvnorm covariance");
  vec z(K);
  for (int k = 0; k < K; k++) z(k) = R::rnorm(0.0, 1.0);
  return mu + L * z;
}

// Convert proportions to the five right-closed intervals used in the paper:
// [0,.2], (.2,.4], (.4,.6], (.6,.8], and (.8,1].
// [[Rcpp::export]]
arma::imat to_ordinal(const arma::mat& W) {
  int p = W.n_rows, n = W.n_cols;
  imat Y(p, n);
  for (int j = 0; j < n; j++)
    for (int i = 0; i < p; i++) {
      double v = W(i, j);
      if (v <= 0.2) Y(i, j) = 1;
      else if (v <= 0.4) Y(i, j) = 2;
      else if (v <= 0.6) Y(i, j) = 3;
      else if (v <= 0.8) Y(i, j) = 4;
      else Y(i, j) = 5;
    }
  return Y;
}

// -----------------------------------------------------------------------------
// MFM normalizing terms
// -----------------------------------------------------------------------------
// Compute log V_n(t) under a Poisson(lambda) prior on K-1. Calculations use
// log-sum-exp stabilization because the factorial/rising-factorial terms can be
// extremely small for realistic sample sizes.
vec compute_log_Vn_poisson(int T_max, int n, double alpha_mfm,
                           double lambda_pois, int K_sum_max) {
  vec log_Vn(T_max + 1);
  for (int t = 0; t <= T_max; t++) {
    std::vector<double> log_terms;
    int K_start = std::max(t, 1);
    for (int Kv = K_start; Kv <= K_sum_max; Kv++) {
      double log_fall = 0.0;
      for (int s = 0; s < t; s++) log_fall += log((double)(Kv - s));
      double Ka = alpha_mfm * (double)Kv;
      double log_rise = lgamma(Ka + (double)n) - lgamma(Ka);
      double log_pK = -lambda_pois + (double)(Kv - 1) * log(lambda_pois)
                      - lgamma((double)Kv);
      log_terms.push_back(log_fall - log_rise + log_pK);
    }
    double max_lt = *std::max_element(log_terms.begin(), log_terms.end());
    double sum_exp = 0.0;
    for (size_t i = 0; i < log_terms.size(); i++)
      sum_exp += exp(log_terms[i] - max_lt);
    log_Vn(t) = max_lt + log(sum_exp);
  }
  return log_Vn;
}

// -----------------------------------------------------------------------------
// One Gibbs-sampling chain for the ordinal factor MFM model
// -----------------------------------------------------------------------------
// Arguments use K for the latent factor dimension (called q in the R interface
// and manuscript). K_init initializes the partition but does not fix the final
// number of occupied clusters. Cluster parameters follow a conjugate
// Normal-Inverse-Wishart prior. The function returns post-burn-in allocation,
// cluster-count, likelihood, mean, and variance traces plus final latent state.
// [[Rcpp::export]]
List ordinal_factor_mfm_mcmc(
    const arma::imat& Y,
    int K,
    int n_iter,
    int n_burn,
    int K_init = 15,             
    double sigma_B = 5.0,
    double alpha_mfm = 1.0,
    int K_max_cluster = 40,
    double kappa0 = 0.1,
    double nu0_add = 2.0,
    double lambda_pois = 10.0,   
    double S0_scale_mult = 1.0,  
    int K_sum_max = 200,
    int loglik_every = 10,      
    bool verbose = true,
    int print_every = 100
) {
  // Dimensions and input checks are repeated here so direct C++ calls fail
  // safely even when the user bypasses the R wrapper.
  int p = Y.n_rows, n = Y.n_cols, C = 5;
  if (p < 1 || n < 2) stop("Y must have at least one row and two columns.");
  if (K < 1 || K > p) stop("K must satisfy 1 <= K <= nrow(Y).");
  if (n_iter <= 0 || n_burn < 0 || n_burn >= n_iter)
    stop("Require n_iter > 0 and 0 <= n_burn < n_iter.");
  if (K_init < 1 || K_max_cluster < 1)
    stop("K_init and K_max_cluster must be positive.");
  if (K_sum_max < K_max_cluster + 1)
    stop("K_sum_max must be at least K_max_cluster + 1.");
  if (sigma_B <= 0.0 || alpha_mfm <= 0.0 || kappa0 <= 0.0 ||
      lambda_pois <= 0.0 || S0_scale_mult <= 0.0)
    stop("All scale and concentration parameters must be positive.");
  if (nu0_add <= 1.0)
    stop("nu0_add must exceed 1 so the inverse-Wishart prior has a finite mean.");
  if (loglik_every < 1 || print_every < 1)
    stop("loglik_every and print_every must be at least 1.");
  for (int j = 0; j < n; j++)
    for (int i = 0; i < p; i++)
      if (Y(i, j) < 1 || Y(i, j) > C)
        stop("Y contains a category outside 1,...,5.");

  if (verbose) Rcout << "p=" << p << ", n=" << n << ", K=" << K
                     << ", K_init=" << K_init
                     << ", lambda=" << lambda_pois
                     << ", alpha=" << alpha_mfm
                     << ", loglik_every=" << loglik_every << std::endl;

  // Initial latent-Gaussian cutpoints. Interior cutpoints are updated below.
  vec gamma = {-datum::inf, 0.0, 0.5, 1.0, 1.5, datum::inf};

  // Initialize B under the positive lower-triangular (PLT) identification
  // constraint on its leading K x K block.
  mat B(p, K, fill::zeros);
  for (int i = 0; i < p; i++)
    for (int k = 0; k < K; k++) {
      if (i < K) {
        if (k < i) B(i, k) = R::rnorm(0.0, 0.5);
        else if (k == i) B(i, k) = std::abs(R::rnorm(1.0, 0.3));
      } else B(i, k) = R::rnorm(0.0, 0.5);
    }

  // Initialize latent Gaussian observations within category-specific intervals.
  mat Z(p, n, fill::zeros);
  for (int j = 0; j < n; j++)
    for (int i = 0; i < p; i++) {
      int c = Y(i, j);
      Z(i, j) = rtruncnorm(0.0, gamma(c - 1), gamma(c));
    }

  // Initialize cell factor scores with a regularized least-squares projection.
  mat eta = solve(B.t() * B + 0.01 * eye(K, K), B.t() * Z,
                  solve_opts::likely_sympd);
  if (!eta.is_finite()) stop("Initial eta contains non-finite values.");

  double nu0 = K + nu0_add;
  vec m0(K, fill::zeros);

  // Initialize the MFM partition and cluster parameters with K-means in factor
  // space. This affects the starting state only, not the inferred cluster count.
  ivec z_clust(n, fill::zeros);
  int n_clust = 1;
  std::vector<vec> mu_clust(K_max_cluster, zeros(K));
  std::vector<mat> Sig_clust(K_max_cluster, eye(K, K));

  mat eta_cov_total; 
  {
    vec em = mean(eta, 1);
    mat ec = eta.each_col() - em;
    eta_cov_total = ec * ec.t() / double(n - 1);
  }

  double s0_scale = std::max(nu0 - K - 1.0, 1.0);
  double eta_scale = std::max(std::abs(trace(eta_cov_total)) / K, 1e-8);
  double eta_jitter = std::max(1e-8, 1e-4 * eta_scale);
  mat S0 = ensure_spd(eta_cov_total * s0_scale + eta_jitter * eye(K, K),
                      "initial S0 from total eta covariance");

  {
    mat centers;
    int K_use = std::max(1, std::min(std::min(K_init, K_max_cluster), n));
    bool ok = false;
    if (K_use > 1) ok = kmeans(centers, eta, K_use, random_subset, 30, false);
    if (ok) {
      n_clust = K_use;
      for (int j = 0; j < n; j++) {
        double best_d = datum::inf; int best_k = 0;
        for (int k = 0; k < K_use; k++) {
          double d = accu(square(eta.col(j) - centers.col(k)));
          if (d < best_d) { best_d = d; best_k = k; }
        }
        z_clust(j) = best_k;
      }

      // Estimate initial component means/covariances and pooled within-cluster
      // covariance used to calibrate the NIW scale matrix S0.
      mat S_pool(K, K, fill::zeros);
      int pool_df = 0;
      for (int k = 0; k < n_clust; k++) {
        std::vector<int> mem;
        for (int j = 0; j < n; j++) if (z_clust(j) == k) mem.push_back(j);
        int nk = mem.size();
        if (nk >= 2) {
          vec ebar(K, fill::zeros);
          for (size_t ii = 0; ii < mem.size(); ii++) ebar += eta.col(mem[ii]);
          ebar /= nk;
          mat Sk(K, K, fill::zeros);
          for (size_t ii = 0; ii < mem.size(); ii++) {
            vec d = eta.col(mem[ii]) - ebar;
            Sk += d * d.t();
          }
          mu_clust[k]  = ebar;
          Sig_clust[k] = ensure_spd(
            Sk / std::max(nk - 1, 1) + 1e-4 * eye(K, K),
            "initial cluster covariance from K-means");

          S_pool += Sk;
          pool_df += nk - 1;
        } else if (nk == 1) {
          mu_clust[k]  = eta.col(mem[0]);
          Sig_clust[k] = ensure_spd(
            eta_cov_total * 0.1 + eta_jitter * eye(K, K),
            "initial singleton cluster covariance");
        } else {
          mu_clust[k]  = centers.col(k);
          Sig_clust[k] = ensure_spd(
            eta_cov_total * 0.1 + eta_jitter * eye(K, K),
            "initial singleton cluster covariance");
        }
      }

      if (pool_df > 0) {
        mat within_cov = S_pool / (double)pool_df;
        double within_scale = std::max(std::abs(trace(within_cov)) / K, 1e-8);
        double within_jitter = std::max(1e-8, 1e-4 * within_scale);
        S0 = ensure_spd(
          within_cov * s0_scale * S0_scale_mult + within_jitter * eye(K, K),
          "S0 from pooled within-cluster covariance");
        if (verbose) {
          Rcout << "K-means init: " << n_clust << " clusters, sizes=";
          ivec cs(n_clust, fill::zeros);
          for (int j = 0; j < n; j++) cs(z_clust(j))++;
          for (int k = 0; k < n_clust; k++) Rcout << cs(k) << " ";
          Rcout << std::endl;
          Rcout << "S0 diag (from within-cluster cov): "
                << S0.diag().t();
          Rcout << "eta_total_cov diag: " << eta_cov_total.diag().t();
          Rcout << "Ratio (within/total): "
                << (within_cov.diag() / eta_cov_total.diag()).t();
        }
      }
    }
  }

  // Precompute MFM V_n ratios once; allocation updates reuse them each iteration.
  int T_max = K_max_cluster + 1;
  vec log_Vn = compute_log_Vn_poisson(T_max, n, alpha_mfm, lambda_pois, K_sum_max);
  vec log_Vn_ratio(T_max + 2, fill::zeros);
  for (int t = 1; t <= T_max; t++)
    log_Vn_ratio(t) = log_Vn(t) - log_Vn(t - 1);

  if (verbose) {
    Rcout << "log V_n(t)/V_n(t-1) (t=2..8): ";
    for (int t = 2; t <= std::min(8, T_max); t++) Rcout << log_Vn_ratio(t) << " ";
    Rcout << std::endl;
  }

  // Allocate post-burn-in storage. z_trace columns are complete cell partitions.
  int n_save = n_iter - n_burn;
  ivec nclust_trace(n_save, fill::zeros);
  vec loglik_trace(n_save, fill::zeros);
  imat z_trace(n, n_save, fill::zeros);

  mat mu_trace(K_max_cluster * K, n_save);
  mat variance_trace(K_max_cluster * K, n_save);
  mu_trace.fill(NA_REAL);
  variance_trace.fill(NA_REAL);

  mat BtB(K, K);
  double sig2B_inv = 1.0 / (sigma_B * sigma_B);
  int save_idx = 0;

  for (int iter = 0; iter < n_iter; iter++) {

    // Step A: sample latent Gaussian Z subject to the observed ordinal category.
    mat M = B * eta;
    for (int j = 0; j < n; j++)
      for (int i = 0; i < p; i++) {
        int c = Y(i, j);
        Z(i, j) = rtruncnorm(M(i, j), gamma(c - 1), gamma(c));
      }

    // Step B: sample cell factor scores conditional on their current component.
    BtB = B.t() * B;
    std::vector<mat> Sig_inv(n_clust), V_eta(n_clust);
    for (int k = 0; k < n_clust; k++) {
      Sig_inv[k] = inv_spd_safe(
        Sig_clust[k], "eta update: inverse of Sig_clust, cluster " + std::to_string(k));
      V_eta[k] = inv_spd_safe(
        BtB + Sig_inv[k], "eta update: posterior covariance, cluster " + std::to_string(k));
    }
    mat BtZ = B.t() * Z;
    for (int j = 0; j < n; j++) {
      int k = z_clust(j);
      vec m_star = V_eta[k] * (BtZ.col(j) + Sig_inv[k] * mu_clust[k]);
      eta.col(j) = rmvnorm(m_star, V_eta[k]);
    }

    // Step C: sample factor loadings B while enforcing PLT identification.
    mat eta_etaT = eta * eta.t();
    mat EZt = eta * Z.t();                                      // K×p
    mat V_full = inv_spd_safe(
      eta_etaT + sig2B_inv * eye(K, K), "B update: V_full precision");
    mat L_full = chol_spd_safe(V_full, "B update: V_full covariance");
    for (int i = 0; i < p; i++) {
      if (i < K) {
        int dim_i = i + 1;
        mat V_i = inv_spd_safe(
          eta_etaT.submat(0, 0, i, i) + sig2B_inv * eye(dim_i, dim_i),
          "B update: constrained row covariance " + std::to_string(i));
        vec rhs_i = EZt.col(i).subvec(0, i);
        vec m_i = V_i * rhs_i;

        double m_ii = m_i(i), v_ii = V_i(i, i);
        if (!std::isfinite(v_ii) || v_ii <= 0.0)
          stop("Non-positive or non-finite PLT diagonal variance.");
        double sd_ii = std::sqrt(v_ii);
        double b_ii = rtruncnorm(m_ii, sd_ii, 0.0, datum::inf);
        if (!std::isfinite(b_ii) || b_ii <= 0.0)
          stop("Invalid positive PLT diagonal draw.");
        vec b_new(dim_i, fill::zeros);
        b_new(i) = b_ii;
        if (i > 0) {
          vec V_ri   = V_i.submat(0, i, i - 1, i);              // Cov(b_rest, b_i)
          vec m_rest = m_i.subvec(0, i - 1);
          vec cond_mean = m_rest + V_ri * ((b_ii - m_ii) / v_ii);
          mat cond_cov = V_i.submat(0, 0, i - 1, i - 1)
                         - (V_ri * V_ri.t()) / v_ii;
          cond_cov = ensure_spd(cond_cov,
            "PLT conditional covariance, row " + std::to_string(i));
          b_new.subvec(0, i - 1) = rmvnorm(cond_mean, cond_cov);
        }
        B.row(i).zeros();
        for (int k = 0; k <= i; k++) B(i, k) = b_new(k);
      } else {
        vec rhs_i = EZt.col(i);
        vec zr(K); for (int k = 0; k < K; k++) zr(k) = R::rnorm(0.0, 1.0);
        B.row(i) = (V_full * rhs_i + L_full * zr).t();        
      }
    }

    // Step D: sample free ordered cutpoints from their feasible intervals.
    for (int c = 2; c <= C - 1; c++) {
      double lo_data = -datum::inf, hi_data = datum::inf;
      for (int j = 0; j < n; j++)
        for (int i = 0; i < p; i++) {
          if (Y(i, j) == c     && Z(i, j) > lo_data) lo_data = Z(i, j);
          if (Y(i, j) == c + 1 && Z(i, j) < hi_data) hi_data = Z(i, j);
        }
      double lo = std::max(gamma(c - 1), lo_data);
      double hi = std::min(gamma(c + 1), hi_data);
      if (lo < hi) gamma(c) = R::runif(lo, hi);
    }

    // Step E: update cell allocations under the MFM predictive probabilities.
    // Existing clusters use Gaussian component densities; a new cluster uses
    // the NIW-integrated multivariate-t prior predictive density.
    double t_scale = (kappa0 + 1.0) / (kappa0 * std::max(nu0 - K + 1.0, 1.0));
    mat t_Sigma = ensure_spd(t_scale * S0, "MFM prior predictive covariance");
    double t_df = std::max(nu0 - K + 1.0, 1.0);
    double log2pi_e = log(2.0 * M_PI);

    ivec n_k_global(K_max_cluster, fill::zeros);
    for (int jj = 0; jj < n; jj++) n_k_global(z_clust(jj))++;
    int t_all = 0;
    for (int k = 0; k < n_clust; k++) if (n_k_global(k) > 0) t_all++;

    mat logdens(n, K_max_cluster); logdens.fill(-1e300);
    for (int k = 0; k < n_clust; k++) {
      mat Lk = chol_spd_safe(
        Sig_clust[k], "allocation density covariance, cluster " + std::to_string(k));
      double logdet = 2.0 * accu(log(Lk.diag()));
      mat D  = eta.each_col() - mu_clust[k];      // K×n
      mat Zs = solve(trimatl(Lk), D);           
      logdens.col(k) = -0.5 * (K * log2pi_e + logdet) - 0.5 * (sum(Zs % Zs, 0).t());
    }

    for (int j = 0; j < n; j++) {
      vec eta_j = eta.col(j);
      int old_k = z_clust(j);

      n_k_global(old_k)--;
      if (n_k_global(old_k) == 0) t_all--;
      int t_minus_j = t_all;

      vec log_prob(n_clust + 1);
      for (int k = 0; k < n_clust; k++) {
        if (n_k_global(k) > 0)
          log_prob(k) = log((double)n_k_global(k) + alpha_mfm) + logdens(j, k);
        else
          log_prob(k) = -1e300;
      }

      int new_idx = n_clust;
      if (n_clust < K_max_cluster) {
        int Vn_idx = std::min(t_minus_j + 1, T_max);
        log_prob(new_idx) = log_Vn_ratio(Vn_idx)
                          + log(alpha_mfm)
                          + dmvt_log(eta_j, m0, t_Sigma, t_df);
      } else {
        log_prob(new_idx) = -1e300;
      }

      double max_lp = log_prob.subvec(0, new_idx).max();
      vec prob(new_idx + 1);
      for (int k = 0; k <= new_idx; k++) prob(k) = exp(log_prob(k) - max_lp);
      prob /= accu(prob);

      double u = R::runif(0.0, 1.0), cum = 0.0;
      int chosen = new_idx;
      for (int k = 0; k <= new_idx; k++) {
        cum += prob(k);
        if (u <= cum) { chosen = k; break; }
      }

      if (chosen == new_idx && n_clust < K_max_cluster) {
        double kappa_n = kappa0 + 1.0, nu_n = nu0 + 1.0;
        vec m_n = (kappa0 * m0 + eta_j) / kappa_n;
        mat S_n = ensure_spd(
          S0 + (kappa0 / kappa_n) * (eta_j - m0) * (eta_j - m0).t(),
          "new cluster inverse-Wishart posterior scale");
        Sig_clust[n_clust] = riwish(nu_n, S_n);
        mu_clust[n_clust]  = rmvnorm(m_n, Sig_clust[n_clust] / kappa_n);

        {
          mat Lk = chol_spd_safe(
            Sig_clust[n_clust], "new cluster density covariance");
          double logdet = 2.0 * accu(log(Lk.diag()));
          mat D  = eta.each_col() - mu_clust[n_clust];
          mat Zs = solve(trimatl(Lk), D);
          logdens.col(n_clust) = -0.5 * (K * log2pi_e + logdet) - 0.5 * (sum(Zs % Zs, 0).t());
        }

        z_clust(j) = n_clust;
        n_k_global(n_clust) = 1;
        t_all++;
        n_clust++;
      } else {
        z_clust(j) = chosen;
        n_k_global(chosen)++;
        if (n_k_global(chosen) == 1) t_all++;
      }
    }

    // Remove empty components and relabel occupied clusters consecutively.
    // Cluster numbers themselves have no scientific meaning.
    ivec cluster_count(n_clust, fill::zeros);
    for (int j = 0; j < n; j++) cluster_count(z_clust(j))++;
    ivec label_map(n_clust, fill::ones); label_map *= -1;
    int new_label = 0;
    std::vector<vec> mu_new; std::vector<mat> Sig_new;
    for (int k = 0; k < n_clust; k++) {
      if (cluster_count(k) > 0) {
        label_map(k) = new_label++;
        mu_new.push_back(mu_clust[k]);
        Sig_new.push_back(Sig_clust[k]);
      }
    }
    for (int j = 0; j < n; j++) z_clust(j) = label_map(z_clust(j));
    n_clust = new_label;
    for (int k = 0; k < n_clust; k++) {
      mu_clust[k]  = mu_new[k];
      Sig_clust[k] = Sig_new[k];
    }

    // Step F: update component means/covariances from the NIW posterior.
    double max_diag_this_iter = 0.0;
    for (int k = 0; k < n_clust; k++) {
      std::vector<int> members;
      for (int j = 0; j < n; j++) if (z_clust(j) == k) members.push_back(j);
      int nk = members.size();
      if (nk == 0) continue;

      vec eta_bar(K, fill::zeros);
      for (size_t idx = 0; idx < members.size(); idx++) eta_bar += eta.col(members[idx]);
      eta_bar /= nk;

      mat S_k(K, K, fill::zeros);
      for (size_t idx = 0; idx < members.size(); idx++) {
        vec diff = eta.col(members[idx]) - eta_bar;
        S_k += diff * diff.t();
      }
      mat R_k = (kappa0 * nk / (kappa0 + nk)) * (eta_bar - m0) * (eta_bar - m0).t();
      double kappa_n = kappa0 + nk;
      double nu_n    = nu0 + nk;
      vec m_n  = (kappa0 * m0 + nk * eta_bar) / kappa_n;
      mat S_n = ensure_spd(
        S0 + S_k + R_k,
        "cluster inverse-Wishart posterior scale, cluster " + std::to_string(k));

      Sig_clust[k] = riwish(nu_n, S_n);
      mu_clust[k]  = rmvnorm(m_n, Sig_clust[k] / kappa_n);

      double md = Sig_clust[k].diag().max();
      if (md > max_diag_this_iter) max_diag_this_iter = md;
    }
    // Save retained draws after burn-in. The observed-data log likelihood is
    // evaluated only at the requested interval because it requires a p x n loop.
    if (iter >= n_burn && save_idx < n_save) {
      nclust_trace(save_idx) = n_clust;
      z_trace.col(save_idx) = z_clust + 1;

      for (int k = 0; k < n_clust; k++) {
        for (int d = 0; d < K; d++) {
          mu_trace(k * K + d, save_idx) = mu_clust[k](d);
          variance_trace(k * K + d, save_idx) = Sig_clust[k](d, d);
        }
      }

      if (save_idx % loglik_every == 0) {
        double ll = 0.0;
        mat Mu = B * eta;
        for (int j = 0; j < n; j++)
          for (int i = 0; i < p; i++) {
            int c = Y(i, j);
            double p_lo = R::pnorm(gamma(c - 1) - Mu(i, j), 0, 1, 1, 0);
            double p_hi = R::pnorm(gamma(c)     - Mu(i, j), 0, 1, 1, 0);
            ll += log(std::max(p_hi - p_lo, 1e-300));
          }
        loglik_trace(save_idx) = ll;
      } else {
        loglik_trace(save_idx) = NA_REAL;
      }
      save_idx++;
    }

    if (verbose && (iter + 1) % print_every == 0) {
      ivec cs(n_clust, fill::zeros);
      for (int j = 0; j < n; j++) cs(z_clust(j))++;
      Rcout << "Iter " << iter + 1 << "/" << n_iter
            << " | K=" << n_clust << " | max Sig_k diag="
            << max_diag_this_iter << " | sizes: ";
      for (int k = 0; k < std::min(n_clust, 12); k++) Rcout << cs(k) << " ";
      Rcout << std::endl;
    }
    if ((iter + 1) % 50 == 0) Rcpp::checkUserInterrupt();
  }

  // Return both posterior traces and the final latent state for diagnostics.
  return List::create(
    Named("B")            = B,
    Named("eta")          = eta,
    Named("gamma")        = gamma,
    Named("z_clust")      = z_clust + 1,
    Named("n_clust")      = n_clust,
    Named("mu_trace")       = mu_trace,
    Named("variance_trace") = variance_trace,
    Named("nclust_trace")   = nclust_trace,
    Named("loglik_trace")   = loglik_trace,
    Named("z_trace")        = z_trace,
    Named("log_Vn")       = log_Vn,
    Named("log_Vn_ratio") = log_Vn_ratio,
    Named("S0")           = S0,
    Named("n_save")       = save_idx
  );
}

// -----------------------------------------------------------------------------
// Dahl least-squares representative partition
// -----------------------------------------------------------------------------
// Z has cells in rows and retained allocation draws in columns. Numeric labels
// may switch across iterations/chains; pairwise co-clustering indicators are
// label-invariant. The selected partition is the sampled draw closest (squared
// Frobenius loss) to the posterior similarity matrix.
// [[Rcpp::export]]
Rcpp::List dahl_partition_cpp(const arma::imat& Z) {
  const int n = Z.n_rows;
  const int M = Z.n_cols;

  if (n <= 0 || M <= 0)
    Rcpp::stop("dahl_partition_cpp: Z must have at least one row and one column.");

  // Posterior similarity matrix:
  // psm(i,j) = posterior proportion of saved draws in which i and j
  // belong to the same cluster.
  arma::mat psm(n, n, arma::fill::zeros);
  psm.diag().ones();

  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      int same_count = 0;
      for (int m = 0; m < M; ++m)
        same_count += (Z(i, m) == Z(j, m));

      const double value = static_cast<double>(same_count) /
                           static_cast<double>(M);
      psm(i, j) = value;
      psm(j, i) = value;
    }
  }

  // Dahl least-squares loss.
  arma::vec loss(M, arma::fill::zeros);
  int best_index = 0;
  double best_loss = std::numeric_limits<double>::infinity();

  for (int m = 0; m < M; ++m) {
    double current_loss = 0.0;

    for (int i = 0; i < n; ++i) {
      for (int j = i + 1; j < n; ++j) {
        const double a_ij = (Z(i, m) == Z(j, m)) ? 1.0 : 0.0;
        const double diff = a_ij - psm(i, j);
        current_loss += 2.0 * diff * diff;
      }
    }

    loss(m) = current_loss;
    if (current_loss < best_loss) {
      best_loss = current_loss;
      best_index = m;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("z_dahl") = Z.col(best_index),
    Rcpp::Named("best_m") = best_index + 1,
    Rcpp::Named("psm") = psm,
    Rcpp::Named("loss") = loss
  );
}

