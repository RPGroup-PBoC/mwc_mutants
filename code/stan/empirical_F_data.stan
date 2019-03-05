data {
    int<lower=1> N_mut;
    int<lower=1> N_wt;
    real ref_bohr;
    real fc_mut[N_mut];
    real fc_wt[N_wt];
}
transformed data {
    f
}

parameters {
    real<lower=0, upper=1> fc_mu_wt;
    real<lower=0> fc_sigma_wt;
    real<lower=0, upper=1> fc_mu_mut;
    real<lower=0> fc_sigma_mut;

}

model {
fc_mu_wt ~ uniform(0, 1);
fc_mu_mut ~ uniform(0, 1);
fc_sigma_mut ~ normal(0, 0.01);
fc_sigma_wt ~ normal(0, 0.01);
fc_mut ~ normal(fc_mu_mut, fc_sigma_mut);
fc_wt ~ normal(fc_mu_wt, fc_sigma_wt);
}
generated quantities {
    real bohr_wt = log((1/fc_mu_wt) - 1);
    real bohr_mut = log((1/fc_mu_mut) - 1);
    real delta_bohr = bohr_wt - bohr_mut;
}

