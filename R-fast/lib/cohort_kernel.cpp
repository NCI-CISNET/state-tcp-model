// Compiled kernel for cohort_prevs() and the SAD sum in mort_outcomes() (tcp_fast.R).
// Each value is computed with the same IEEE double operations, in the same order, as the R
// expressions in R/02_model_functions.R, so results are bit-identical. Fused multiply-add must stay disabled for that to hold (load_kernel() also
// compiles with -ffp-contract=off), and the per-cohort former-smoker sum uses long double like
// R's rowSums().
//
// Former smokers are held one cohort at a time in an age x years-since-quit buffer (years since
// quit fastest, 32 KB, stays in L1 cache), and sY / dY arrive with years since quit as their
// first dimension, so the inner loops read contiguous memory.
#ifdef __clang__
#pragma clang fp contract(off)
#endif
#include <Rcpp.h>
#include <algorithm>
#include <vector>
using namespace Rcpp;

// ip, cp: initiation / cessation probabilities, age x cohort (100 x nc)
// sNS, sCS: 1 - death probability for never / current smokers, age x cohort
// sYt: 1 - death probability for former smokers, years since quit x age x cohort (40 x 100 x nc)
// popdist, dCS: age x calendar year (100 x nc); dYt: years since quit x age x calendar year
// [[Rcpp::export]]
List cohort_kernel(NumericMatrix ip, NumericMatrix cp, NumericMatrix sNS, NumericMatrix sCS,
                   NumericVector sYt, NumericVector births, NumericMatrix popdist,
                   NumericMatrix dCS, NumericVector dYt, bool keep_fs) {
  const int na = 100, nc = ip.ncol(), nj = 40;
  const R_xlen_t plane = (R_xlen_t)na * nc;
  NumericMatrix pop(na, nc), NSprev(na, nc), CSprev(na, nc), SAD(na, nc);
  NumericVector FSprev(keep_fs ? plane * nj : 0);
  std::vector<double> F((size_t)na * nj);
  double NS[100], CS[100];

  // SAD cells before any cohort's birth (calendar year t < age r): prevalence is 0.
  for (int t = 0; t < nc; t++) for (int r = t + 1; r < na; r++) {
    const R_xlen_t ap = (R_xlen_t)t * na + r;
    const double pd = popdist[ap];
    const double *dy = &dYt[ap * nj];
    double sad = pd * (0.0 * dCS[ap]);
    for (int j = 0; j < nj; j++) sad = sad + pd * 0.0 * dy[j];
    SAD[ap] = sad;
  }

  for (int b = 0; b < nc; b++) {
    const R_xlen_t c0 = (R_xlen_t)b * na;
    std::fill(F.begin(), F.end(), 0.0);
    NS[0] = births[b]; CS[0] = 0.0;
    for (int i = 1; i < na; i++) {
      const R_xlen_t k1 = c0 + i - 1;
      const double ns = NS[i - 1], cs = CS[i - 1], p_i = ip[k1], p_c = cp[k1];
      NS[i] = ns * (1 - p_i) * sNS[k1];
      CS[i] = ns * p_i * sNS[k1] + cs * (1 - p_c) * sCS[k1];
      double *f = &F[(size_t)i * nj];
      const double *f1 = &F[(size_t)(i - 1) * nj];
      const double *sy = &sYt[k1 * nj];
      f[0] = cs * (p_c) * sCS[k1];
      for (int j = 1; j < nj - 1; j++) f[j] = f1[j - 1] * sy[j - 1];
      f[nj - 1] = f1[nj - 2] * sy[nj - 2] + f1[nj - 1] * sy[nj - 1];
    }
    for (int r = 0; r < na; r++) {
      const R_xlen_t k = c0 + r;
      const double *f = &F[(size_t)r * nj];
      long double fs_sum = 0.0;
      for (int j = 0; j < nj; j++) fs_sum += f[j];
      const double p = NS[r] + CS[r] + (double)fs_sum;
      pop[k] = p;
      NSprev[k] = NS[r] / p;
      CSprev[k] = CS[r] / p;
      if (keep_fs) for (int j = 0; j < nj; j++) FSprev[k + j * plane] = f[j] / p;
      // SAD in age x calendar-year layout: this cell is calendar year t = b + r.
      const int t = b + r;
      if (t < nc) {
        const R_xlen_t ap = (R_xlen_t)t * na + r;
        const double pd = popdist[ap];
        const double *dy = &dYt[ap * nj];
        double sad = pd * (CSprev[k] * dCS[ap]);
        for (int j = 0; j < nj; j++) sad = sad + pd * (f[j] / p) * dy[j];
        SAD[ap] = sad;
      }
    }
  }
  if (keep_fs) FSprev.attr("dim") = IntegerVector::create(na, nc, nj);

  return List::create(_["NSprev"] = NSprev, _["CSprev"] = CSprev, _["pop"] = pop,
                      _["FSprev"] = keep_fs ? (SEXP)FSprev : R_NilValue, _["SAD"] = SAD);
}
