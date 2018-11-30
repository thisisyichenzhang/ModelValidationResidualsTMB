// Estimate and validate a Ricker model based on data simulated from the logistic map
//
// Compare Thygesen et al (submitted, 2016): Validation of state space models
// fitted as mixed effects models
//
// This file implements the "Theta logistic population model" from
// Pedersen et al 2012, Ecol. Modelling. With theta=1, this is the Ricker
// model.
//
// Uffe Høgsbro Thygesen and Kasper Kristensen, 2016

#include <TMB.hpp>

template <class Type>
Type objective_function<Type>::operator()()
{
  /* Data section */
  DATA_VECTOR(Y);                  // Counted abundance
  DATA_VECTOR_INDICATOR(keep, Y);  // For one-step predictions

  /* Parameter section */
  PARAMETER_VECTOR(X);  // Latent states. As last as long as Y;
                        // extra elements are not used
  PARAMETER(logr);      // Growth rate
  PARAMETER(logK);      // Carrying capacity
  PARAMETER(logQ);      // Process noise
  PARAMETER(logS);		// Sample size controlling measurement noise
  PARAMETER(logitp);    // Probability of having extra zeros (zero-inflated Poisson)
  
  /* Procedure section */

  Type r = exp(logr);
  Type K = exp(logK);
  Type Q = exp(logQ);
  Type S = exp(logS);
  Type p = 1.0 /(1.0 + exp(-logitp));


  int timeSteps = Y.size();
  Type nll = 0;

  // Contributions from state transitions
  for (int i = 1; i < timeSteps; i++) {
    Type m = X[i - 1] + r * (1.0 - exp(X[i - 1]) / K);
    nll -= dnorm(X[i], m, sqrt(Q), true);
  }

  // Contributions from observations
  for (int i = 0; i < timeSteps; i++) {
    nll -= keep(i) * dzipois(Y[i], S * exp(X[i]), p ,true);  // keep(i) for one-step predictions
  }

  // Make parameters available to R (for report) 
  ADREPORT(r);
  ADREPORT(K);
  ADREPORT(Q);
  ADREPORT(p);
  return nll;
}