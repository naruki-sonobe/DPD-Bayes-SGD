// funcitons
functions {
  real DPD_normal_lpdf(real x, real mu, real sigma, real beta, real w){
    real density;
    density = (1.0/beta)*exp(normal_lpdf(x| mu,sigma))^(beta)-(1.0/((2.0*pi())^(beta/2.0)*(1.0 + beta)^1.5*sigma^beta));
    return w*density; 
  }
}


// The input data is a vector 'x' of length 'N'.
data {
  int<lower=0> N;
  vector[N] x;
  real beta;
  vector[2] theta0;
  real<lower=0> w;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real mu;
  real<lower=0> sigma;
}

// The model to be estimated. We model the output
// 'x' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  
  for(n in 1:N) {
    target += DPD_normal_lpdf(x[n]|mu, sigma, beta, w);
  }
}
