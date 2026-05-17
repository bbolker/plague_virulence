// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
using namespace Rcpp;

// smooth-minimum: approaches min(a,b) as both grow, with a == 0 guard
inline double phi_minfun(double a, double b) {
    if (a <= 0.0) return 0.0;
    return a * (1.0 - std::exp(-b / a));
}

inline double plogis_cpp(double x) {
    return 1.0 / (1.0 + std::exp(-x));
}

// Evaluate raw degree-2 polynomial in 4 variables.
// Variables: x1 = R01, x2 = R02, x3 = log(I10), x4 = log(I20).
// Coefficient order matches R poly(..., degree = 2, raw = TRUE):
//   [0]  intercept
//   [1]  1.0.0.0  x1
//   [2]  2.0.0.0  x1^2
//   [3]  0.1.0.0  x2
//   [4]  1.1.0.0  x1*x2
//   [5]  0.2.0.0  x2^2
//   [6]  0.0.1.0  x3
//   [7]  1.0.1.0  x1*x3
//   [8]  0.1.1.0  x2*x3
//   [9]  0.0.2.0  x3^2
//   [10] 0.0.0.1  x4
//   [11] 1.0.0.1  x1*x4
//   [12] 0.1.0.1  x2*x4
//   [13] 0.0.1.1  x3*x4
//   [14] 0.0.0.2  x4^2
inline double eval_poly(double x1, double x2, double x3, double x4,
                        const NumericVector& b) {
    return b[0]
        + b[1]*x1      + b[2]*x1*x1
        + b[3]*x2      + b[4]*x1*x2   + b[5]*x2*x2
        + b[6]*x3      + b[7]*x1*x3   + b[8]*x2*x3  + b[9]*x3*x3
        + b[10]*x4     + b[11]*x1*x4  + b[12]*x2*x4 + b[13]*x3*x4 + b[14]*x4*x4;
}

//' Two-strain metapopulation simulation (C++ core)
//'
//' Called by simulate_metapopulation_2strain_cpp(); not intended for direct use.
//'
//' @param n_patches Number of patches
//' @param n_years Number of years
//' @param K Carrying capacity per patch
//' @param r Host growth rate
//' @param c0 Migration proportion
//' @param nu Dead-host-to-flea scaling constant
//' @param rho Flea carrying capacity per migrating host
//' @param alpha Transmission efficiency per flea
//' @param D Disease-induced mortality factor
//' @param R01 Basic reproduction number, strain 1
//' @param R02 Basic reproduction number, strain 2
//' @param invade_year Year of strain 2 introduction (1-indexed, matches R convention)
//' @param initial_inf_ratio_1 Initial fraction of patches infected with strain 1
//' @param initial_inf_ratio_2 Fraction of patches seeded with strain 2 at invasion
//' @param initial_pop_ratio Initial population as fraction of K
//' @param early_stop Stop early when both strains go extinct after invasion
//' @param z_vec Pre-computed final sizes: c(z1, z2, z_mean)
//' @param coef_finalsize 15 raw-polynomial coefficients for logit(final size)
//' @param coef_ratio     15 raw-polynomial coefficients for logit(I1 fraction)
//' @param coinf_approx Coinfection approximation: "polyfit" or "yy"
//' @return List with N (3-D array), I1, I2, S1, S2, total_inf1, total_inf2
// [[Rcpp::export]]
List sim_metapop_2strain_cpp(
    int    n_patches,
    int    n_years,
    double K,
    double r,
    double c0,
    double nu,
    double rho,
    double alpha,
    double D,
    double R01,
    double R02,
    int    invade_year,
    double initial_inf_ratio_1,
    double initial_inf_ratio_2,
    double initial_pop_ratio,
    bool   early_stop,
    NumericVector z_vec,
    NumericVector coef_finalsize,
    NumericVector coef_ratio,
    std::string   coinf_approx
) {
    const int np = n_patches, ny = n_years;
    const bool use_polyfit = (coinf_approx == "polyfit");
    const double z1     = z_vec[0];
    const double z2     = z_vec[1];
    const double z_mean = z_vec[2];
    const double fac1   = std::max(0.0, 1.0 - 1.0 / R01);
    const double fac2   = std::max(0.0, 1.0 - 1.0 / R02);
    // regression data have R01 >= R02; swap inputs when R01 < R02
    const bool swap_strains = (R01 < R02);

    // N: (n_patches x n_years x 4) column-major, stages 0..3 =
    //   begin / after_growth / after_colonization / end
    NumericVector N(np * ny * 4, 0.0);

    // patch x year matrices
    NumericMatrix I1(np, ny), I2(np, ny), S1(np, ny), S2(np, ny);
    NumericVector total_inf1(ny, 0.0), total_inf2(ny, 0.0);

    // N[patch, year, stage]  (column-major)
    auto Nref = [&](int p, int y, int s) -> double& {
        return N[p + y * np + s * np * ny];
    };

    // ---- Year 0 initialisation ----
    for (int p = 0; p < np; p++)
        Nref(p, 0, 0) = initial_pop_ratio * K;

    int n_init1 = std::max(1, (int)(initial_inf_ratio_1 * np));
    IntegerVector init1_idx = Rcpp::sample(np, n_init1, false);   // 1-indexed
    for (int i = 0; i < n_init1; i++) {
        int p = init1_idx[i] - 1;
        S1(p, 0)     = z1 * D * Nref(p, 0, 0);
        I1(p, 0)     = 1.0;
        total_inf1[0] += S1(p, 0);
    }

    NumericVector b_vec(np, 0.0);   // strain-2 seeds; held at zero before invasion

    // ---- Main loop: C++ year index k = 1 .. ny-1
    //      corresponds to R loop k = 2 .. n_years                         ----
    for (int k = 1; k < ny; k++) {

        // Stage 1 -> 2: logistic growth from previous year-end population.
        // For k == 1 (R k == 2), use N_begin[0] - S1[0] because N_end[0]
        // was never stored (matches the same special-case in the R code).
        NumericVector N_ag(np);
        for (int p = 0; p < np; p++) {
            double n_prev = (k == 1)
                ? Nref(p, 0, 0) - S1(p, 0)
                : Nref(p, k - 1, 3);
            N_ag[p]       = n_prev + r * n_prev * (1.0 - n_prev / K);
            Nref(p, k, 1) = N_ag[p];
        }

        // Stage 2 -> 3: mean-field migration
        double mean_pop = 0.0;
        for (int p = 0; p < np; p++) mean_pop += N_ag[p];
        mean_pop /= np;

        NumericVector N_ac(np);
        for (int p = 0; p < np; p++) {
            N_ac[p]       = (1.0 - c0) * N_ag[p] + c0 * mean_pop;
            Nref(p, k, 2) = N_ac[p];
        }

        // Force of infection: accumulate dispersed flea pools
        NumericVector q1(np, 0.0), q2(np, 0.0);
        NumericVector E_tot(np, 0.0), R_tot(np, 0.0);
        double sum_q1E = 0.0, sum_q2E = 0.0;

        for (int p = 0; p < np; p++) {
            double s1 = S1(p, k - 1), s2 = S2(p, k - 1);
            double s_tot = s1 + s2;
            if (s_tot > 0.0) {
                q1[p] = s1 / s_tot;
                q2[p] = s2 / s_tot;
            }
            double nu_s  = nu * s_tot;
            double cap   = rho * c0 * N_ag[p];
            E_tot[p]     = phi_minfun(nu_s, cap);
            R_tot[p]     = nu_s - E_tot[p];
            sum_q1E     += q1[p] * E_tot[p];
            sum_q2E     += q2[p] * E_tot[p];
        }

        // Poisson seeds for strain 1 (every year from k >= 1)
        NumericVector a_vec(np);
        for (int p = 0; p < np; p++) {
            double lam1 = alpha * (q1[p] * R_tot[p] + sum_q1E / np) * fac1;
            a_vec[p] = (double) R::rpois(lam1);
        }

        // Poisson seeds for strain 2 (only from invade_year onward).
        // R uses 1-indexed k; invade_year in R corresponds to k = invade_year - 1 here.
        if (k >= invade_year - 1) {
            for (int p = 0; p < np; p++) {
                double lam2 = alpha * (q2[p] * R_tot[p] + sum_q2E / np) * fac2;
                b_vec[p] = (double) R::rpois(lam2);
            }
        }
        // Force at least one seed in a random subset of patches at the invasion year
        if (k == invade_year - 1) {
            int n_inv = std::max(1, (int)(initial_inf_ratio_2 * np));
            IntegerVector inv_idx = Rcpp::sample(np, n_inv, false);
            for (int i = 0; i < n_inv; i++)
                b_vec[inv_idx[i] - 1] += 1.0;
        }

        // Per-patch epidemic outcomes
        for (int p = 0; p < np; p++) {
            bool ind1 = (a_vec[p] > 0.0);
            bool ind2 = (b_vec[p] > 0.0);
            I1(p, k)  = ind1 ? 1.0 : 0.0;
            I2(p, k)  = ind2 ? 1.0 : 0.0;

            double fs = 0.0, i1f = 0.0;   // final size fraction, I1 fraction

            if (ind1 && !ind2) {
                fs  = z1;
                i1f = 1.0;
            } else if (ind2 && !ind1) {
                fs  = z2;
                i1f = 0.0;
            } else if (ind1 && ind2) {
                if (use_polyfit && N_ac[p] > 0.0) {
                    // Polynomial approximation of co-infection outcome.
                    // Regression data have R01 >= R02; swap when reversed.
                    double I10 = a_vec[p] / N_ac[p];
                    double I20 = b_vec[p] / N_ac[p];
                    double x1_, x2_, x3_, x4_;
                    if (swap_strains) {
                        x1_ = R02; x2_ = R01;
                        x3_ = std::log(I20); x4_ = std::log(I10);
                    } else {
                        x1_ = R01; x2_ = R02;
                        x3_ = std::log(I10); x4_ = std::log(I20);
                    }
                    fs  = plogis_cpp(eval_poly(x1_, x2_, x3_, x4_, coef_finalsize));
                    i1f = plogis_cpp(eval_poly(x1_, x2_, x3_, x4_, coef_ratio));
                    if (swap_strains) i1f = 1.0 - i1f;
                } else {
                    // YY approximation: average R0 for final size, seed-ratio partition
                    fs  = z_mean;
                    i1f = a_vec[p] / (a_vec[p] + b_vec[p]);
                }
            }

            double Z      = fs * D * N_ac[p];
            S1(p, k)      = (ind1 || ind2) ? Z * i1f         : 0.0;
            S2(p, k)      = (ind1 || ind2) ? Z * (1.0 - i1f) : 0.0;
            Nref(p, k, 3) = N_ac[p] - S1(p, k) - S2(p, k);
            total_inf1[k] += S1(p, k);
            total_inf2[k] += S2(p, k);
        }

        // Carry end-of-year population to begin of next year
        if (k < ny - 1) {
            for (int p = 0; p < np; p++)
                Nref(p, k + 1, 0) = Nref(p, k, 3);
        }

        // Early stopping: both strains absent after invasion
        if (early_stop && k >= invade_year) {
            bool any_inf = false;
            for (int p = 0; p < np; p++) {
                if (I1(p, k) > 0.0 || I2(p, k) > 0.0) { any_inf = true; break; }
            }
            if (!any_inf) break;
        }
    }

    N.attr("dim") = IntegerVector::create(np, ny, 4);

    return List::create(
        Named("N")          = N,
        Named("I1")         = I1,
        Named("I2")         = I2,
        Named("S1")         = S1,
        Named("S2")         = S2,
        Named("total_inf1") = total_inf1,
        Named("total_inf2") = total_inf2
    );
}
