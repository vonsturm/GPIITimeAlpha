// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include <string>
#include <climits>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <cmath>
#include <algorithm>

#include "TChain.h"
#include "TH1D.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "BAT/BCMath.h"
#include "BAT/BCH1D.h"

// gerda-ada includes
#include "FileMap.h"
#include "DataLoader.h"

#include "GETRunConfiguration.hh"
#include "GERunConfigurationManager.hh"

#include "ProgressBar.h"

#include "TimeAlpha.h"

using namespace std;
using namespace gada;

// ---------------------------------------------------------
TimeAlpha::TimeAlpha( const char * name ) : BCModel(name) {

	fDetectorType = kUNKNOWN; // all detectors
	fDetectorEnriched = kUNKNOWN; // all detectors
	fModel = LinAndExp; // OnlyLin, OnlyExp

	fNBins = 100;
	fHMaximum = 0.;
	fHMinimum = 100.;

	fHTimeAlpha = NULL;
	fHLiveTimeFraction = NULL;

	fVTimeAlpha = vector<double>(0);
	fVLiveTimeFraction = vector<double>(0);

	DefineParametersAndPriors();

	cout << "Created model" << endl;
}

int TimeAlpha::DefineParametersAndPriors()
{
	// define parameters
	AddParameter( "constant", 0., 2. );
	AddParameter( "amplitude", 0., 5. );
	AddParameter( "halflife", 137., 140. );

	// set priors
	switch ( fModel )
	{
		case LinAndExp:
			fModelName = "linexp";
			SetPriorConstant( 0 );
			SetPriorConstant( 1 );
//			SetPriorConstant( 2 );
			SetPriorGauss( 2, 138.4, 0.2 );
//			SetPriorDelta( 2, 138.3763 );
			break;

		case OnlyExp:
			fModelName = "exp";
			SetPriorDelta( 0, 0. );
			SetPriorConstant( 1 );
			SetPriorConstant( 2 );
//			SetPriorGauss( 2, 138.4, 0.2 );
			break;

		case OnlyLin:
			fModelName = "lin";
			SetPriorConstant( 0 );
			SetPriorDelta( 1, 0. );
			SetPriorDelta( 2, 0. );
			break;
		default:
			cout << "Model not set!" << endl;
			return -1;
	}

	return 0;
}


// ---------------------------------------------------------
TimeAlpha::~TimeAlpha() {
	// destructor
	if(fHTimeAlpha) fHTimeAlpha->Delete();
	if(fHLiveTimeFraction) fHLiveTimeFraction->Delete();
}

// ---------------------------------------------------------
int TimeAlpha::InitializeHistograms()
{
	fHTimeAlpha = new TH1D( "HTimeAlpha", "HTimeAlpha", fNBins, fHMinimum, fHMaximum );
	fHLiveTimeFraction = new TH1D( "HLiveTimeFraction", "HLiveTimeFraction", fNBins, fHMinimum, fHMaximum );

	fHTimeAlpha_fine = new TH1D( "HTimeAlpha_fine", "HTimeAlpha_fine", (int)fHMaximum, fHMinimum, fHMaximum );
	fHLiveTimeFraction_fine = new TH1D( "HLiveTimeFraction_fine", "HLiveTimeFraction_fine", (int)fHMaximum, fHMinimum, fHMaximum );

	if( fHTimeAlpha && fHLiveTimeFraction && fHTimeAlpha_fine && fHLiveTimeFraction_fine )
		return 0;
	else
		return -1;
}

// ---------------------------------------------------------
// Read Data from Runs key list tier3
// ---------------------------------------------------------
int TimeAlpha::ReadDataPhaseII( string keylist )
{
	string MU_CAL = getenv("MU_CAL");

	// check if $MU_CAL is set
	if( MU_CAL == "" )
	{
		cout << "Environment variable MU_CAL not set. Have it point to the directory where runconfiguration.db can be found.\n"
			 << "On LNGS that would be e.g. /nfs/gerda5/gerda-data/blind/active/meta/config/_aux/geruncfg " << endl;
		return -1;
	}

	cout << "MU_CAL = " << MU_CAL << endl;

	// initialize run configuration manager
	GERunConfigurationManager * RunConfManager = new GERunConfigurationManager();
	GETRunConfiguration * RunConf = 0;
	RunConfManager -> AllowRunConfigurationSwitch(true);
	RunConfManager -> SetVerbosity(1);

	// tell the data map where to find everything
	string GERDA_PHASEII_DATA = getenv("GERDA_PHASEII_DATA");
	GERDA_PHASEII_DATA += "/gen";
	string GERDA_DATA_SETS = getenv("GERDA_DATA_SETS"); GERDA_DATA_SETS += "/";
	string data_set = GERDA_DATA_SETS; data_set += keylist;

	cout << "Reading from data dir: " << GERDA_PHASEII_DATA << endl;
	cout << "\t data set dir: " << GERDA_DATA_SETS << endl;
	cout << "\t from keylist: " << keylist << endl;
	cout << "This may take a while..." << endl;

	// setting FileMap and DataLoader
	gada::FileMap myMap;
	myMap.SetRootDir( GERDA_PHASEII_DATA );
	myMap.BuildFromListOfKeys( data_set );

	gada::DataLoader l;
	l.AddFileMap(&myMap);
	l.BuildTier3();

	TChain * chain = l.GetSharedMasterChain();
	int nentries = chain->GetEntries();

	cout << "There are " << nentries << " events in the chain!" <<endl;

	// fill the data in histograms
	int eventChannelNumber;
	unsigned long long timestamp;
	//unsigned int decimalTimestamp;
	int multiplicity;
	int isTP;
	int isVetoedInTime;
	vector<int> * firedFlag = new vector<int>(fNDetectors);
	vector<double> * energy = new vector<double>(fNDetectors);
	vector<double> * risetime = new vector<double>(fNDetectors);
	vector<int> * failedFlag = new vector<int>(fNDetectors);
	vector<int> * failedFlag_isPhysical = new vector<int>(fNDetectors);
	vector<int> * failedFlag_isSaturated = new vector<int>(fNDetectors);

	chain -> SetBranchAddress("eventChannelNumber", &eventChannelNumber);
	chain -> SetBranchAddress("timestamp",&timestamp);
	//chain -> SetBranchAddress("decimalTimestamp",&decimalTimestamp);
	chain -> SetBranchAddress("multiplicity",&multiplicity);
	chain -> SetBranchAddress("isTP",&isTP);
	chain -> SetBranchAddress("isVetoedInTime", &isVetoedInTime);
	chain -> SetBranchAddress("firedFlag", &firedFlag);
	chain -> SetBranchAddress("rawEnergyGauss",&energy);
	chain -> SetBranchAddress("risetime",&risetime);
	chain -> SetBranchAddress("failedFlag_isPhysical",&failedFlag_isPhysical);
	chain -> SetBranchAddress("failedFlag_isSaturated",&failedFlag_isSaturated);

	cout << "Branches are set" << endl;

	chain -> GetEntry( 0 ); // Get timestamp of first event in chain
	unsigned long long time0 = timestamp; // time in hours

	chain -> GetEntry( nentries - 1 ); // Get timestamp of last event in chain
	unsigned long long timeN = timestamp; // time in hours

	// calculate fit interval in days
	int timeInDays = (int)( ( timeN - time0 ) / ( 60. * 60. * 24. ) );

	int bins = timeInDays/fBinning + 1;
	double min = 0., max = ( timeInDays - timeInDays%fBinning + fBinning ) * ( 60. * 60. * 24. );
	SetNBinsHistograms( bins, min, max );
	InitializeHistograms();

	cout << "Fitting " << timeInDays << " days of data." << endl;
	cout << "Starting loop over event chain..." << endl;

	// ProgressBar
	ProgressBar bar( nentries, '#', false );

	for( int e = 0; e < nentries; e++ )
	{
		bar.Update();

		chain->GetEntry(e);

		unsigned long long time = timestamp - time0;

		if( isTP )
		{
			fHLiveTimeFraction -> Fill( time );
			fHLiveTimeFraction_fine -> Fill( time );
			continue;
		}

		// Apply cuts
		if( isVetoedInTime ) continue;
		if( multiplicity > 1 ) continue;

		for( int d = 0; d < fNDetectors; d++ )
		{
			if( !firedFlag -> at( d ) ) continue;
			if( failedFlag -> at( d ) ) continue;

			// use kUNKNOWN to select all detectors
			RunConf = RunConfManager -> GetRunConfiguration( timestamp );
			const GETDetector * det = RunConf -> GetDetector( d );
			DetectorType_t dType = det -> GetDetectorType(); // kIsBEGe=1, kIsCOAX=2, kUNKNOWN=3
			Bool_t dIsEnriched = det -> IsEnriched();
			TString dDetName = det -> GetDetectorTable() -> at( det->GetID() ) -> GetName();

			if( !RunConf -> IsOn( d ) ) continue;

			// fit of single detector: accept only this detector namespace
			// fit of enrBEGe detectors: match Type and Enriched and exclude GD91C
			// fit of enrCoax detectors: -"- and exclude ANG4
			// fit of natCoax det: -"- exlude GTF45
			if( fSingleDetectorFit ){ if( dDetName != fDataSet ) continue; }
			else if( fDetectorType != kUNKNOWN && fDetectorEnriched != kUNKNOWN )
			{
				if( dType != fDetectorType ) continue;
				if( dIsEnriched != fDetectorEnriched ) continue;
				if( fDataSet == "enrBEGe" && dDetName == "GD91C" ) continue;
				if( fDataSet == "enrCoax" && dDetName == "ANG4" )  continue;
				if( fDataSet == "natCoax" && dDetName == "GTF45_2" ) continue;
			}

			if( ( multiplicity == 1 && !failedFlag_isPhysical -> at( d ) && energy -> at( d ) > 3500. )
				|| !failedFlag_isSaturated -> at( d ) )
			{
				fHTimeAlpha -> Fill( time );
				fHTimeAlpha_fine -> Fill( time );

				if(  failedFlag_isSaturated -> at( d ) )
					cout << " Saturated Event: " << e << "(" << d << ")" << endl;
			}
		}
	}

	cout << "...DONE" << endl;

	// Test Pulser rate is 0.05Hz = 1/20s
	double TPFrequency = 1. / 20.;

	// Scale live time fraction with TP rate
	fHLiveTimeFraction -> Scale( 1./( TPFrequency * fHLiveTimeFraction -> GetBinWidth( 1 ) ) );
	fHLiveTimeFraction_fine -> Scale( 1./( TPFrequency * fHLiveTimeFraction_fine -> GetBinWidth( 1 ) ) );

	const string outfilename = Form( "./out/%s_TimeAlpha_Data.root", fDataSet.c_str() );
	TFile * out = new TFile( outfilename.c_str(), "RECREATE" );
	fHTimeAlpha -> Write();
	fHTimeAlpha_fine -> Write();
	fHLiveTimeFraction -> Write();
	fHLiveTimeFraction_fine -> Write();
	out -> Close();

	FillDataArrays();

	return 0;
}


// ---------------------------------------------------------
int TimeAlpha::FillDataArrays()
{
	cout << "Loading data in vectors for faster access." << endl;

	for( int b = 1; b <= fNBins; b++ )
	{
		fVTimeAlpha.push_back( fHTimeAlpha->GetBinContent(b) );
		fVLiveTimeFraction.push_back( fHLiveTimeFraction->GetBinContent(b) );
	}

	return 0;
}


// ---------------------------------------------------------
double TimeAlpha::LogLikelihood(const std::vector<double> & parameters)
{
	// This methods returns the logarithm of the conditional probability
	// p(data|parameters). This is where you have to define your model.

	// access parameters from vector by remembering their positions, e.g.
	// double mu = parameters[0];
	// or by looking up their indicies, e.g.
	// double mu = parameters[fParameters.Index("mu")];

	// Calculate your likelihood according to your model. You may find
	// the built in functions such as BCMath::LogPoisson helpful.
	// Return the logarithm of this likelood

	// translate parameters given in days in seconds
	double secondsInOneDay = 24.*60.*60.;
	double constant = parameters[0] / secondsInOneDay;
	double amplitude = parameters[1] / secondsInOneDay;
	double halflife = parameters[2] * secondsInOneDay;

	double BinWidth = (fHMaximum - fHMinimum) / fNBins;

	double logprob = 0.;

	for( int b = 0; b < fNBins; b++)
	{
		if( fVLiveTimeFraction.at(b) <= fLTLimit ) continue;

		double t1 = fHMinimum + b * BinWidth;
		double t2 = t1 + BinWidth;

		double tau = halflife / log(2.);

		double lambda = constant * BinWidth - amplitude * tau * ( exp( -t2/tau ) - exp( -t1/tau ) );
		lambda *= fVLiveTimeFraction.at(b);

		if( lambda <= 0. ) continue;

	    double data = fVTimeAlpha.at(b);

	    double sum = data * log(lambda) - lambda - BCMath::LogFact( (int)data );

		logprob += sum;
	  }

	  return logprob;
}


// ---------------------------------------------------------
double TimeAlpha::EstimatePValue()
{
	//Allen's routine for the evaluation of p-value
  	//This is derived from PRD 83 (2011) 012004, appendix
  	//taken from Luciano

	vector<double> parameters = GetBestFitParameters();
	double logp0 = LogLikelihood( parameters ); //likelihood at the mode

  	/*
    	Now we initialize the Markov Chain by setting the number of entries in a bin to
    	the integer
    	part of the mean in that bin.  This will give the maximum possible likelihood.
  	*/

	double sumlog = 0;

	/* mean is the array where you have the expected mean in each bin (calculated
     	from the parameter values at the mode. Nom is the nominal value of entries in
     	each bin (integer part)
  	*/

  	vector<double> mean( fNBins, 0 );
	vector<int> nom( fNBins, 0 );

	double secondsInOneDay = 24.*60.*60.;

	double constant = parameters[0] / secondsInOneDay;
	double amplitude = parameters[1] / secondsInOneDay;
	double halflife = parameters[2] * secondsInOneDay;

	double BinWidth = (fHMaximum - fHMinimum) / fNBins;

  	for( int ibin = 0; ibin < fNBins; ibin++ )
	{
		//Require a live time fraction of at least 10% otherwise skip data point
		if( fVLiveTimeFraction.at(ibin) <= fLTLimit ) continue;

                double t1 = fHMinimum + ibin * BinWidth;
                double t2 = t1 + BinWidth;

                double tau = halflife / log(2.);

                double lambda = constant * BinWidth - amplitude * tau * ( exp( -t2/tau ) - exp( -t1/tau ) );
                lambda *= fVLiveTimeFraction.at(ibin);

		mean.at( ibin ) = std::max( lambda, 1e-8 );
      		nom.at( ibin ) = int( mean.at( ibin ) );
      		sumlog += BCMath::LogPoisson( nom[ ibin ], mean[ ibin ] );
  	}

  	cout << "Logprob for best: " << sumlog << endl;

  	/*
  	Now we run the Markov chain to generate new configurations.  One iteration means
  	a loop over all the bins, with an attempt to vary each bin up or down by one unit.  We
  	accept/reject at each step  and compare the data logprob to the simulated at the end of each iteration.
  	*/

  	const int nloops = 100000;
  	int Pgood = 0;

  	for( int iloop = 0; iloop < nloops; iloop++ )
  	{
      		for( int ibin = 0; ibin < fNBins; ibin++)
      		{
			if( fVLiveTimeFraction.at(ibin) <= fLTLimit ) continue;

    	  		if ( rand() > RAND_MAX/2 ) // Try to increase the bin content by 1
    	  		{
    		  		double r = mean[ ibin ]/( nom[ ibin ] + 1 );
    		  		double rtest = double( rand() )/RAND_MAX;
    		  		if( rtest < r ) //Accept
    		  		{
    			  		nom[ ibin ] = nom[ ibin ] + 1;
    			  		sumlog += log(r);
    		  		}
    	  		}
    	  		else // Try to decrease the bin content by 1
    	  		{
    		  		double r = nom[ ibin ]/mean[ ibin ];
    		  		double rtest = double( rand() )/RAND_MAX;
    		  		if ( rtest < r ) //Accept
    		  		{
    			  		nom[ ibin ] = nom[ ibin ] - 1;
    			  		sumlog += log(r);
    		  		}
    	  		}
      		}

		if ( sumlog < logp0 ) Pgood++;
  	}

  	double pvalue = double(Pgood)/double(nloops);

  	cout << "p-value is " << pvalue << endl;

  	return pvalue;
}



// ---------------------------------------------------------
void TimeAlpha::WriteOutput( string outputfilename, double corr, string timeformat )
{
	double secondsInOneDay = 24.*60.*60.;

	const vector<double> BestFitMarginalized = GetBestFitParametersMarginalized();

	const std::vector<double> BestFitGlobal = GetBestFitParameters();
	const std::vector<double> BestFitGlobalErrors = GetBestFitParameterErrors();

	// only uncertainties on contant and amplitude parameters are taken into account
	TF1 * MF = new TF1( "MF", "[0]+[1]*exp(-x*log(2.)/[2])", fHMinimum, fHMaximum );
	MF -> SetParameters( BestFitGlobal.at(0), BestFitGlobal.at(1), BestFitGlobal.at(2) * secondsInOneDay);
	MF -> SetLineColor(kAzure-1); MF -> SetLineWidth(1);

	TF1 * MF_up = new TF1( "MF_up", "[0]+[1]*exp(-x*log(2.)/[2]) + sqrt( [3]*[3] + exp(-x*log(2.)/[2])*exp(-x*log(2.)/[2])*[4]*[4] + 2*exp(-x*log(2.)/[2])*[3]*[4]*[5] )", fHMinimum, fHMaximum );
	MF_up -> SetParameters( BestFitGlobal.at(0), BestFitGlobal.at(1), BestFitGlobal.at(2) * secondsInOneDay,
		BestFitGlobalErrors.at(0), BestFitGlobalErrors.at(1), corr );
	MF_up -> SetLineColor(kAzure-1); MF_up -> SetLineStyle(2); MF_up -> SetLineWidth(1);

	TF1 * MF_low = new TF1( "MF_low", "[0]+[1]*exp(-x*log(2.)/[2]) - sqrt( [3]*[3] + exp(-x*log(2.)/[2])*exp(-x*log(2.)/[2])*[4]*[4] + 2*exp(-x*log(2.)/[2])*[3]*[4]*[5] )", fHMinimum, fHMaximum );
	MF_low -> SetParameters( BestFitGlobal.at(0), BestFitGlobal.at(1), BestFitGlobal.at(2) * secondsInOneDay,
		BestFitGlobalErrors.at(0), BestFitGlobalErrors.at(1), corr );
	MF_low -> SetLineColor(kAzure-1); MF_low -> SetLineStyle(2); MF_low -> SetLineWidth(1);

	TH1D * copyHTimeAlpha = (TH1D*)fHTimeAlpha -> Clone( "copyHTimeAlpha" );
	copyHTimeAlpha -> Divide( fHLiveTimeFraction );
	copyHTimeAlpha -> Sumw2();
	copyHTimeAlpha -> Scale( 1./fBinning );
	copyHTimeAlpha -> GetXaxis() -> SetTimeDisplay(1);
	copyHTimeAlpha -> GetXaxis() -> SetTimeFormat( timeformat.c_str() );
	copyHTimeAlpha -> GetXaxis() -> SetTitle( "date" );
	copyHTimeAlpha -> GetYaxis() -> SetTitle( "alphas / live-day" );
	copyHTimeAlpha -> GetYaxis() -> SetNdivisions( 7 + 100*5 + 10000*0 );
	copyHTimeAlpha -> GetYaxis() -> SetTitleOffset( 0.7 );

	TLegend * l = new TLegend( 0.65, 0.7, 0.85, 0.9 );
	l -> AddEntry( copyHTimeAlpha, "data", "pl" );
	l -> AddEntry( MF, "model", "l" );
	l -> AddEntry( MF, "f(t) = C + N exp#left(#frac{-ln(2) t}{T_{1/2}}#right)", "" );
	l -> AddEntry( MF, Form("C = (%.1f +- %.1f) cts/live-day", BestFitGlobal.at(0), BestFitGlobalErrors.at(0)), "" );
	l -> AddEntry( MF, Form("A = (%.1f +- %.1f) cts/live-day", BestFitGlobal.at(1), BestFitGlobalErrors.at(1)), "" );
	l -> SetLineColor( kWhite );

	const string title = Form( "Time dependence of alpha events: %s", fDataSet.c_str() );
	TCanvas * c = new TCanvas( "TimeAlpha", title.c_str(), 1000, 500 );
	copyHTimeAlpha -> Draw();
	MF -> Draw("same");
	MF_up -> Draw("same");
	MF_low -> Draw("same");
	l -> Draw();

	TFile * out = new TFile( outputfilename.c_str(), "RECREATE" );
	c -> Write();
	MF -> Write();
	MF_up -> Write();
	MF_low -> Write();
	out -> Close();

	return;
}


// ---------------------------------------------------------
void TimeAlpha::WriteDistributions()
{
	const string filename = Form("./out/%s_TimeAlphas_posteriors.root", fDataSet.c_str() );
	TFile * test = new TFile( filename.c_str(), "RECREATE" );

	for( uint i = 0; i < GetNParameters(); i++ )
	{
		BCH1D * parHisto = MCMCGetH1Marginalized( i );

		TH1D * histo = parHisto -> GetHistogram();

		histo -> Write();
	}

	test -> Close();

	return;
}


// ---------------------------------------------------------
// double TimeAlpha::LogAPrioriProbability(const std::vector<double> & parameters) {
// 	// This method returns the logarithm of the prior probability for the
// 	// parameters p(parameters).

// 	// You need not overload this function, if you are using built-in
// 	// priors through the function SetPriorGauss, SetPriorConstant, etc.
// }
