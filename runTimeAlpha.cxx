// ***************************************************************
// This file was created using the bat-project script
// for project TimeAlpha.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "BAT/BCLog.h"
#include "BAT/BCAux.h"
#include "BAT/BCSummaryTool.h"
#include "BAT/BCH2D.h"

#include "TimeAlpha.h"
#include "TH2D.h"

using namespace std;

void Usage();

int main( int argc, char* argv[]  )
{
	string keylist;
	string data_set = "enrBEGe"; // enrBEGe, enrCoax, natCoax, all, single detector
	string precision = "kMedium";

	if( argc == 2 )	keylist = argv[1];
	else if( argc == 3 )
	{
		keylist = argv[1];
		data_set = argv[2];
	}
	else if( argc == 4 )
	{
		keylist = argv[1];
		data_set = argv[2];
		precision = argv[3];
	}
	else if( argc > 4 )
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
	m -> SetFittingDataSet( data_set );
	m -> SetBinningInDays( 20 );

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

	double pValue = m -> EstimatePValue();
	cout << "P-Value estimated with MCMC algorithm: " << pValue << endl;

//	m -> WriteDistributions();

	// create a new summary tool object, to print change from prior -> posterior
	BCSummaryTool * summary = new BCSummaryTool(m);
	summary -> PrintKnowledgeUpdatePlots( Form( "./out/%s_TimeAlpha_update.pdf", m->GetDataSet().c_str() ) );
	summary -> PrintParameterPlot( Form( "./out/%s_TimeAlpha_parameters.pdf", m->GetDataSet().c_str() ) );
	summary -> PrintCorrelationPlot( Form( "./out/%s_TimeAlpha_correlation.pdf", m->GetDataSet().c_str() ) );
	summary -> PrintCorrelationMatrix( Form( "./out/%s_TimeAlpha_correlationMatrix.pdf", m->GetDataSet().c_str() ) );

	// int index = 0 * npar + 1;
    // double corr = summary -> fCorrCoeff.at(index);
	double corr = m -> GetMarginalized(0, 1) -> GetHistogram() -> GetCorrelationFactor();

	cout << "Correlation factor between constant and amplitude parameters" << endl;

	m -> WriteOutput( Form("./out/%s_TimeAlpha_model.root", m->GetDataSet().c_str() ), corr );

	// calculate p-value
//	m -> CalculatePValue( m->GetBestFitParameters() );

	// print results of the analysis into a text file
	m -> PrintResults( Form( "./out/%s_TimeAlpha_results.txt", m->GetDataSet().c_str() ) );

	delete m;

	// delete summary;

	BCLog::OutSummary("Exiting");

	// close log file
	BCLog::CloseLog();

	return 0;
}


void Usage()
{
	cout << "Use for Time Alpha fitting " << endl;
	cout << "\t runTimeAlpha keylist [data-set|enrBEGe] [precision|kMedium]" << endl;
	cout << "\t data-set: enrBEGe, enrCoax, natCoax, all, single-detector-name " << endl;
	cout << "\t precision: kLow, kMedium, kHigh, kVeryHigh" << endl;

	return;
}
