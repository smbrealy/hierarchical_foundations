data { 
    int<lower=0> N; // number of data points
    vector[N] wn; // natural frequency observations
    array[N] int I; // structure number for each data point
    int<lower=0> K; // number of structures
    vector[B] C; // coefficients

    // hyperpriors
    vector[2] Emu0; 
    vector[2] Vmu0; 
    vector[2] Esig0;  
    vector[2] Vsig0; 
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
    
    // global vals
    real<lower=0> gamma; // meaurement noise

    // local vals
    vector<lower=0>[K] Es; // local mean
    vector<lower=0>[K] Vs; // local std
} 

transformed parameters {
    vector[N] s; // stiffness values
    for (n in 1:N) { // for all data points
        s[n] = Es[I[n]] + Vs[I[n]]*snorm[n]; // draw local stiffness using standard normal
    }
}

model {
    vector[N] Ewn; // expected natural frequency (without noise)

    Emu ~ gamma(Emu0[1],Emu0[2]);
    Vmu ~ gamma(Vmu0[1],Vmu0[2]);
    Esig ~ gamma(Esig0[1],Esig0[2]); 
    Vsig ~ gamma(Vsig0[1],Vsig0[2]); 
    gamma ~ gamma(gamma0[1],gamma0[2]);

    for (k in 1:K) {
        Es[k] ~ normal(Emu,Vmu);
        Vs[k] ~ normal(Esig,Vsig);
    }

    for (n in 1:N) {
        snorm[n] ~ std_normal();
        Ewn[n] = C[1]*pow(s[n],5) + C[2]*pow(s[n],4) + C[3]*pow(s[n],3) + C[4]*pow(s[n],2) + C[5]*pow(s[n],1) + C[6]; // surrogate fe model
    }
        
    wn ~ normal(Ewn, gamma); // likelihood  
}



