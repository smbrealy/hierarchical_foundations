data { 
    int<lower=0> N; // number of data
    vector[N] wn; // variates
    array[N] int I; // index
    int<lower=0> K; // number of domains
    vector[6] c; // coefficients

    // hyperpriors
     
    vector[2] Emu0; 
    vector[2] Vmu0; 
    //vector<lower=0>[2] Esig0; // not to be learned...
    vector[2] Esig0;  
    vector[2] Vsig0; 
    //vector<lower=0>[2] Vsig0; 
    vector[2] gamma0;
} 

parameters {
    vector[N] snorm; // vector of standard normally distributed values
    
    // used to draw local means
    real<lower=0> Emu;
    real<lower=0> Vmu;
    
    // used to draw local stds
    real<lower=0> Esig;
    real<lower=0> Vsig;
    
    // local vals
    //vector[K] Es;
    vector<lower=0>[K] Es;
    vector<lower=0>[K] Vs;

    real<lower=0> gamma; // noise
} 

transformed parameters {
    vector[N] s; // stiffness // putting <lower=0> here results in samples below 0 errors
    for (n in 1:N) { // for all data
        s[n] = Es[I[n]] + Vs[I[n]]*snorm[n]; // draw local stiffness using standard normal
    }
}
model {
    vector[N] Ewn; // expected natural frequency (without noise)

    Emu ~ normal(Emu0[1],Emu0[2]);
    Vmu ~ cauchy(Vmu0[1],Vmu0[2]);
    Esig ~ gamma(Esig0[1],Esig0[2]); // this param is causing some errors e.g. 0 or inf but must be + finite
    Vsig ~ gamma(Vsig0[1],Vsig0[2]); // one or two similar errors with this param too.
    //Esig ~ cauchy(Esig0[1],Esig0[2]); // this param is causing some errors e.g. 0 or inf but must be + finite
    //Vsig ~ cauchy(Vsig0[1],Vsig0[2]);
    
    gamma ~ cauchy(gamma0[1],gamma0[2]);

    for (k in 1:K) {
        Es[k] ~ normal(Emu,Vmu);
        Vs[k] ~ normal(Esig,Vsig);
    }

    for (n in 1:N) {
        snorm[n] ~ std_normal();
        Ewn[n] = c[1]*pow(s[n],5) + c[2]*pow(s[n],4) + c[3]*pow(s[n],3) + c[4]*pow(s[n],2) + c[5]*pow(s[n],1) + c[6];
    }
        
    wn ~ normal(Ewn, gamma); // likelihood  
}



