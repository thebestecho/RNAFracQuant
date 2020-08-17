getmixingratiosTSP <- function(ctdata,iter=1000,
                                control=list(adapt_delta=0.85),...) {
    # head(ctdata) %>% print()
    stansummary <- fit_mixingTSP(ctdata=ctdata,iter=iter,control=control,...) %>%
        summary()
    # return medians
    data.frame(
        mixing.Sup=stansummary$summary["mixing_sup","50%"],
        mixing.P100=stansummary$summary["mixing_p100","50%"],
        lp.n_eff  =stansummary$summary["lp__","n_eff"],
        lp.Rhat   =stansummary$summary["lp__","Rhat"])
}
