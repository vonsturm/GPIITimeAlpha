// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__TIMEALPHA__H
#define __BAT__TIMEALPHA__H

#include <BAT/BCModel.h>

// This is a TimeAlpha header file.
// Model source code is located in file TimeAlpha/TimeAlpha.cxx

// ---------------------------------------------------------
class TimeAlpha : public BCModel {
 public:

	// Constructor and destructor
	TimeAlpha(const char * name = "TimeAlpha");
	~TimeAlpha();

	// Methods to overload, see file TimeAlpha.cxx
	double LogLikelihood(const std::vector<double> & parameters);
	// double LogAPrioriProbability(const std::vector<double> & parameters);

};
// ---------------------------------------------------------

#endif
