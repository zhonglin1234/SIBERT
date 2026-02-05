#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP _TRRident_calculate_sum_cpp(SEXP, SEXP);
extern SEXP _TRRident_compute_prob_cpp(SEXP, SEXP);
extern SEXP _TRRident_cpp_update_beta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TRRident_cpp_update_beta_mat(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TRRident_cpp_update_phi_MCMC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TRRident_cpp_update_phi_MCMC_vec(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TRRident_MC_sampling_x_cpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TRRident_rpg_R(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_TRRident_calculate_sum_cpp", (DL_FUNC) &_TRRident_calculate_sum_cpp, 2},
    {"_TRRident_compute_prob_cpp", (DL_FUNC) &_TRRident_compute_prob_cpp, 2},
    {"_TRRident_cpp_update_beta", (DL_FUNC) &_TRRident_cpp_update_beta, 7},
    {"_TRRident_cpp_update_beta_mat", (DL_FUNC) &_TRRident_cpp_update_beta_mat, 7},
    {"_TRRident_cpp_update_phi_MCMC", (DL_FUNC) &_TRRident_cpp_update_phi_MCMC, 8},
    {"_TRRident_cpp_update_phi_MCMC_vec", (DL_FUNC) &_TRRident_cpp_update_phi_MCMC_vec, 8},
    {"_TRRident_MC_sampling_x_cpp", (DL_FUNC) &_TRRident_MC_sampling_x_cpp, 5},
    {"_TRRident_rpg_R", (DL_FUNC) &_TRRident_rpg_R, 3},
    {NULL, NULL, 0}
};

void R_init_TRRident(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
