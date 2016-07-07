// ***************************************************************
// This file was created using the bat-project script
// for project TimeAlpha.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCSummaryTool.h>

#include "TimeAlpha.h"

int main()
{
	// set nicer style for drawing than the ROOT default
	BCAux::SetStyle();

	// open log file
	BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

	// create new TimeAlpha object
	TimeAlpha* m = new TimeAlpha("TimeAlpha");

	// set precision
	m -> MCMCSetPrecision(BCEngineMCMC::kMedium);

	BCLog::OutSummary("Test model created");

	//////////////////////////////
	// perform your analysis here

	// normalize the posterior, i.e. integrate posterior over the full
	// parameter space
	// m -> SetIntegrationMethod(BCIntegrate::kIntDefault);
	// m -> Normalize();

	// run MCMC and marginalize posterior w/r/t all parameters and all
	// combinations of two parameters
	// m -> MarginalizeAll(BCIntegrate::kMargMetropolis);

	// run mode finding; by default using Minuit
	// m -> FindMode( m->GetBestFitParameters() );

	// draw all marginalized distributions into a PDF file
	// m -> PrintAllMarginalized("TimeAlpha_plots.pdf");

	// print all summary plots
	// m -> PrintParameterPlot("TimeAlpha_parameters.pdf");
	// m -> PrintCorrelationPlot("TimeAlpha_correlation.pdf");
	// m -> PrintCorrelationMaxtrix("TimeAlpha_correlationMatrix.pdf");

	// create a new summary tool object, to print change from prior -> posterior
	// BCSummaryTool * summary = new BCSummaryTool(m);
	// summary -> PrintKnowledgeUpdatePlots("TimeAlpha_update.pdf");

	// calculate p-value
	// m -> CalculatePValue( m->GetBestFitParameters() );

	// print results of the analysis into a text file
	//  m -> PrintResults("TimeAlpha_results.txt");

	delete m;
	// delete summary;

	BCLog::OutSummary("Exiting");

	// close log file
	BCLog::CloseLog();

	return 0;
}
