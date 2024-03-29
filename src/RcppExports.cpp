// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// arma_onestage
arma::mat arma_onestage(const arma::mat& Y, const arma::colvec& samp_unit_ids, const arma::colvec& strata_ids, const arma::colvec& strata_samp_sizes, const arma::colvec& strata_pop_sizes, const Rcpp::CharacterVector& singleton_method, const Rcpp::LogicalVector& use_singleton_method_for_domains, const int& stage);
RcppExport SEXP _survey_arma_onestage(SEXP YSEXP, SEXP samp_unit_idsSEXP, SEXP strata_idsSEXP, SEXP strata_samp_sizesSEXP, SEXP strata_pop_sizesSEXP, SEXP singleton_methodSEXP, SEXP use_singleton_method_for_domainsSEXP, SEXP stageSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type samp_unit_ids(samp_unit_idsSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type strata_ids(strata_idsSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type strata_samp_sizes(strata_samp_sizesSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type strata_pop_sizes(strata_pop_sizesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type singleton_method(singleton_methodSEXP);
    Rcpp::traits::input_parameter< const Rcpp::LogicalVector& >::type use_singleton_method_for_domains(use_singleton_method_for_domainsSEXP);
    Rcpp::traits::input_parameter< const int& >::type stage(stageSEXP);
    rcpp_result_gen = Rcpp::wrap(arma_onestage(Y, samp_unit_ids, strata_ids, strata_samp_sizes, strata_pop_sizes, singleton_method, use_singleton_method_for_domains, stage));
    return rcpp_result_gen;
END_RCPP
}
// arma_multistage
arma::mat arma_multistage(arma::mat Y, arma::mat samp_unit_ids, arma::mat strata_ids, arma::mat strata_samp_sizes, arma::mat strata_pop_sizes, Rcpp::CharacterVector singleton_method, Rcpp::LogicalVector use_singleton_method_for_domains, Rcpp::LogicalVector use_only_first_stage, int stage);
RcppExport SEXP _survey_arma_multistage(SEXP YSEXP, SEXP samp_unit_idsSEXP, SEXP strata_idsSEXP, SEXP strata_samp_sizesSEXP, SEXP strata_pop_sizesSEXP, SEXP singleton_methodSEXP, SEXP use_singleton_method_for_domainsSEXP, SEXP use_only_first_stageSEXP, SEXP stageSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type samp_unit_ids(samp_unit_idsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type strata_ids(strata_idsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type strata_samp_sizes(strata_samp_sizesSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type strata_pop_sizes(strata_pop_sizesSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type singleton_method(singleton_methodSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type use_singleton_method_for_domains(use_singleton_method_for_domainsSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type use_only_first_stage(use_only_first_stageSEXP);
    Rcpp::traits::input_parameter< int >::type stage(stageSEXP);
    rcpp_result_gen = Rcpp::wrap(arma_multistage(Y, samp_unit_ids, strata_ids, strata_samp_sizes, strata_pop_sizes, singleton_method, use_singleton_method_for_domains, use_only_first_stage, stage));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_survey_arma_onestage", (DL_FUNC) &_survey_arma_onestage, 8},
    {"_survey_arma_multistage", (DL_FUNC) &_survey_arma_multistage, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_survey(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
