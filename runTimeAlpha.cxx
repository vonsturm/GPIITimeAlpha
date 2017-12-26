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

using namespace std;

void Usage();

int main( int argc, char* argv[]  )
{
	string keylist;
	string precision = "kMedium";

	if( argc == 2 ) keylist = argv[1];
	else if( argc == 3 )
	{
		keylist = argv[1];
		precision = argv[2];
	}
	else if( argc > 3 )
	{
		cout << "Too many parameters" << endl;
		Usage();
	}
	else
	{
		cout << "Not enough parameters" << endl;
		Usage();
		return -1;
	}

	// set nicer style for drawing than the ROOT default
	BCAux::SetStyle();

	// open log file
	BCLog::OpenLog( "log.txt", BCLog::detail, BCLog::detail );

	// create new TimeAlpha object
	TimeAlpha * m = new TimeAlpha( "TimeAlpha" );
	m -> SetFittingDataSet( "enrBEGe" ); // enrBEGe, enrCoax, natCoax, all, single detector
	m -> SetBinningInDays( 10 );

	// NBinsHistograms( 180, 0., 180. );
	m -> ReadDataPhaseII( keylist );

	// set precision
	if( precision == "kLow" )			m -> MCMCSetPrecision( BCEngineMCMC::kLow );
	else if( precision == "kMedium" ) 	m -> MCMCSetPrecision( BCEngineMCMC::kMedium );
	else if( precision == "kHigh" ) 	m -> MCMCSetPrecision( BCEngineMCMC::kHigh );
	else if( precision == "kVeryHigh" ) m -> MCMCSetPrecision( BCEngineMCMC::kVeryHigh );

	BCLog::OutSummary("Test model created");

	//////////////////////////////
	// perform your analysis here

	// normalize the posterior, i.e. integrate posterior over the full
	// parameter space
	m -> SetIntegrationMethod(BCIntegrate::kIntDefault);
	m -> Normalize();

	// run MCMC and marginalize posterior w/r/t all parameters and all
	// combinations of two parameters
	m -> MarginalizeAll(BCIntegrate::kMargMetropolis);

	// run mode finding; by default using Minuit
	m -> FindMode( m->GetBestFitParameters() );

	// draw all marginalized distributions into a PDF file
//	m -> PrintAllMarginalized("TimeAlpha_plots.pdf");

	m -> WriteOutput( "./out/TimeAlpha_model.root" );

	double pValue = m -> EstimatePValue();
	cout << "P-Value estimated with MCMC algorithm: " << pValue << endl;

//	m -> WriteDistributions();

	// create a new summary tool object, to print change from prior -> posterior
	BCSummaryTool * summary = new BCSummaryTool(m);
	summary -> PrintKnowledgeUpdatePlots("./out/TimeAlpha_update.pdf");
	summary -> PrintParameterPlot("./out/TimeAlpha_parameters.pdf");
	summary -> PrintCorrelationPlot("./out/TimeAlpha_correlation.pdf");
	summary -> PrintCorrelationMatrix("./out/TimeAlpha_correlationMatrix.pdf");

	// calculate p-value
//	m -> CalculatePValue( m->GetBestFitParameters() );

	// print results of the analysis into a text file
	m -> PrintResults("./out/TimeAlpha_results.txt");

	delete m;

	// delete summary;

	BCLog::OutSummary("Exiting");

	// close log file
	BCLog::CloseLog();

	return 0;
}


void Usage()
{



	return;
}
