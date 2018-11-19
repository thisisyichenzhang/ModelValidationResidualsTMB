// This file serves as the C++ code for implementation of the beaver example in TMB 
// Model and Data: Reynolds, P.S. (1994) Time-series analyses of beaver body temperatures.

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator()()
{

// Data Declarations
  DATA_VECTOR(Y); // Temperatue Observations
  DATA_VECTOR_INDICATOR(keep, Y);  // For one-step predictions
  DATA_FACTOR(act);  // Indicator for activities, 
                    // 1 for beaver outside of the retreat (high-intensity activity)
                   // 0 otherwise 

// Parameter Declarations
  PARAMETER(logitPhi);
  PARAMETER(logSdState);
  PARAMETER(logSdObs);
  PARAMETER_VECTOR(mu);
  PARAMETER_VECTOR(X);

 // Parameter Transformations
  Type phi = Type(2.0) / (Type(1.0) + exp(-logitPhi)) - 1.0;
  Type sdState = exp(logSdState);
  Type sdObs = exp(logSdObs);

 // Initialization of Negative Log Likelihood (nll)
  Type nll = 0.0;

// Update nll 
  nll -= dnorm(X(0), mu(act(0)), sdState/sqrt(1.0-pow(phi,2.0)), true);
  // Update nll looping over  states X 
  for(int i = 1; i < X.size(); ++i)
    nll -= dnorm(X(i),mu(act(i)) + phi*(X(i-1) - mu(act(i-1))), sdState, true);
  // Update nll looping over observations Y 
  for(int i = 0; i < Y.size(); ++i)
    nll -= keep(i) * dnorm(Y(i), X(i), sdObs, true); // For one-step prediction

// Make parameters available to R (for report) 
  ADREPORT(phi);
  ADREPORT(sdState);
  ADREPORT(sdObs);

  return nll;
}