data {
    int<lower=1> N_mut;
    int<lower=1> N_wt;
    real ref_bohr;
    real fc_mut[N_mut];
    real fc_wt[N_wt];
}

parameters {
    real<lower=0, upper=1> fc_mu_wt;
    real log_fc_sigma_wt;
    real<lower=0, upper=1> fc_mu_mut;
    real log_fc_sigma_mut;
}

transformed parameters {
    real fc_sigma_wt = exp(log_fc_sigma_wt);
    real fc_sigma_mut = exp(log_fc_sigma_mut);
}

model {
fc_mu_mut ~ uniform(0, 1);
log_fc_sigma_mut ~ normal(0, 3);
log_fc_sigma_wt ~ normal(0, 3);
fc_mut ~ normal(fc_mu_mut, fc_sigma_mut);
fc_wt ~ normal(fc_mu_wt, fc_sigma_wt);
}

