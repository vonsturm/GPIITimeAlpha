// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "TimeAlpha.h"

#include <BAT/BCMath.h>

// ---------------------------------------------------------
TimeAlpha::TimeAlpha(const char * name) : BCModel(name) {
	// constructor
	// define parameters here. For example:
	// AddParameter("mu",-2,1,"#mu");
	// and set priors, if using built-in priors. For example:
	// SetPriorGauss("mu",-1,0.25);
}

// ---------------------------------------------------------
TimeAlpha::~TimeAlpha() {
	// destructor
}

// ---------------------------------------------------------
double TimeAlpha::LogLikelihood(const std::vector<double> & parameters) {
	// This methods returns the logarithm of the conditional probability
	// p(data|parameters). This is where you have to define your model.

	// access parameters from vector by remembering their positions, e.g.
	// double mu = parameters[0];
	// or by looking up their indicies, e.g.
	// double mu = parameters[fParameters.Index("mu")];

	// Calculate your likelihood according to your model. You may find
	// the built in functions such as BCMath::LogPoisson helpful.
	// Return the logarithm of this likelood

	return -1;
}

// ---------------------------------------------------------
// double TimeAlpha::LogAPrioriProbability(const std::vector<double> & parameters) {
// 	// This method returns the logarithm of the prior probability for the
// 	// parameters p(parameters).

// 	// You need not overload this function, if you are using built-in
// 	// priors through the function SetPriorGauss, SetPriorConstant, etc.
// }
