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

#include "TimeAlpha.h"

#include <BAT/BCMath.h>
#include <BAT/BCH1D.h>

#include "TChain.h"
#include "TH1D.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"

// gerda-ada includes
#include "FileMap.h"
#include "DataLoader.h"

using namespace std;
using namespace gada;

// ---------------------------------------------------------
TimeAlpha::TimeAlpha( const char * name ) : BCModel(name) {

	fDetectorType = all;
//	fModel = OnlyLin;
//	fModel = OnlyExp;
	fModel = LinAndExp;

	fNBins = 100;
	fHMaximum = 0.;
	fHMinimum = 100.;

	fHTimeAlpha = NULL;
	fHLiveTimeFraction = NULL;

	fVTimeAlpha = vector<double>(0);
	fVLiveTimeFraction = vector<double>(0);

	// constructor
	// define parameters here. For example:
	AddParameter( "constant", 0., 15. );
	AddParameter( "amplitude", 0., 25. );
	AddParameter( "halflife", 137., 140. );

	// and set priors, if using built-in priors. For example:
	switch ( fModel )
	{
		case LinAndExp:
			fModelName = "linexp";
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
	}

	cout << "Created model" << endl;
}

// ---------------------------------------------------------
TimeAlpha::~TimeAlpha() {
	// destructor
	if(fHTimeAlpha) fHTimeAlpha->Delete();
	if(fHLiveTimeFraction) fHLiveTimeFraction->Delete();
}

// ---------------------------------------------------------
// Read Data from Runs key list tier3
// ---------------------------------------------------------
int TimeAlpha::ReadData( string keylist )
{
	string GERDA_PHASEII_DATA = getenv("GERDA_PHASEII_DATA");
	string GERDA_DATA_SETS = getenv("GERDA_DATA_SETS"); GERDA_DATA_SETS += "/";
	string data_set = GERDA_DATA_SETS; data_set += keylist;

	cout << "Reading data from keylist: " << data_set << endl;
	cout << "This may take a while..." << endl;

	//
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
	unsigned int decimalTimestamp;
	vector<int> * firedFlag = new vector<int>(fNDetectors);
	int multiplicity;
	vector<double> * energy = new vector<double>(fNDetectors);
	int isTP;
	int isVetoed;
	int isVetoedInTime;
	vector<int> * failedFlag = new vector<int>(fNDetectors);

	chain -> SetBranchAddress("eventChannelNumber", &eventChannelNumber);
	chain -> SetBranchAddress("timestamp",&timestamp);
	chain -> SetBranchAddress("decimalTimestamp",&decimalTimestamp);
	chain -> SetBranchAddress("firedFlag", &firedFlag);
	chain -> SetBranchAddress("multiplicity",&multiplicity);
	chain -> SetBranchAddress("rawEnergyGauss",&energy);
	chain -> SetBranchAddress("isVetoed", &isVetoed);
	chain -> SetBranchAddress("isTP",&isTP);
	chain -> SetBranchAddress("isVetoedInTime", &isVetoedInTime);
	chain -> SetBranchAddress("failedFlag",&failedFlag);

	// Prepare histograms cut the last piece of data if it doen't fit in the 20 days scheme
	unsigned long long secondsInAnHour= 60 * 60;

	chain->GetEntry(0);

	double time0 = (double)timestamp / (double)secondsInAnHour; // time in hours

	fHTimeAlpha = new TH1D( "HTimeAlpha", "HTimeAlpha", fNBins, fHMinimum, fHMaximum );
	fHLiveTimeFraction = new TH1D( "HLiveTimeFraction", "HLiveTimeFraction", fNBins, fHMinimum, fHMaximum );

	TH1D * HTimeAlpha_fine = new TH1D( "HTimeAlpha_fine", "HTimeAlpha_fine", (int)fHMaximum, fHMinimum, fHMaximum );
	TH1D * HLiveTimeFraction_fine = new TH1D( "HLiveTimeFraction_fine", "HLiveTimeFraction_fine", (int)fHMaximum, fHMinimum, fHMaximum );

	cout << "Loop over event chain." << endl;

	vector<int> enrBEGeChannels = {	0,1,2,3,4,5,6,7,
									11,12,13,14,15,16,17,18,
									19,20,21,22,23,24,25,26,
									30,31,32,33,34,35 };
	vector<int> enrCoaxChannels = {8,9,10,27,28,29,36};
	vector<int> natCoaxChannels = {37,38,39};

	for( int e = 0; e < nentries; e++ )
	{
		chain->GetEntry(e);

		double time = (double)timestamp / (double)secondsInAnHour;
		time -= time0;

		if( e%10000 == 0 ) cout << "h: " << time << " d: " << time/24. << endl;

		if( isTP )
		{
			fHLiveTimeFraction->Fill( time / 24. );
			HLiveTimeFraction_fine->Fill( time / 24. );
			continue;
		}

		// Apply cuts
		if( multiplicity != 1 ) continue;
		if( isVetoed ) 	continue;
		if( isVetoedInTime ) 	continue;

		for( int i = 0; i < fNDetectors; i++ )
		{
			if( failedFlag->at(i) ) continue;

			if( fDetectorType == enrBEGe )
			{
				if( find( enrBEGeChannels.begin(), enrBEGeChannels.end(), i) == enrBEGeChannels.end() )
					continue;
			}
			if( fDetectorType == enrCoax )
			{
				if( find( enrCoaxChannels.begin(), enrCoaxChannels.end(), i) == enrCoaxChannels.end() )
					continue;
			}
			if( fDetectorType == natCoax )
			{
				if( find( natCoaxChannels.begin(), natCoaxChannels.end(), i) == natCoaxChannels.end() )
					continue;
			}

			if( energy->at(i) > 3500. && energy->at(i) < 5300. )
			{
				fHTimeAlpha->Fill( time / 24. );
				HTimeAlpha_fine->Fill( time / 24. );
			}
		}
	}

	// Expected number of TP in 1 day: Test Pulser rate is 0.05Hz = 1/20s
	double TPExpected = 24. * 60. * 60. / 20.;

	fHLiveTimeFraction->Scale( 1./( TPExpected * fHLiveTimeFraction -> GetBinWidth( 1 ) ) );
	HLiveTimeFraction_fine->Scale( 1./( TPExpected * HLiveTimeFraction_fine -> GetBinWidth( 1 ) ) );

	TFile * out = new TFile( "TimeAlpha_Data.root", "RECREATE" );
	HTimeAlpha_fine -> Write();
	fHTimeAlpha -> Write();
	HLiveTimeFraction_fine -> Write();
	fHLiveTimeFraction -> Write();
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

	double constant = parameters[0];
	double amplitude = parameters[1];
	double halflife = parameters[2];

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

	    	double data = fVTimeAlpha.at(b);

	    	double sum = data * log(lambda) - lambda - BCMath::LogFact( (int)data );

		// require at least 10% livetime fraction otherwise ignore data point
		if( lambda > 0 ) logprob += sum;
	  }

	  return logprob;
}


// ---------------------------------------------------------
double TimeAlpha::EstimatePValue()
{
	//Allen's routine for the evaluation of p-value
  	//This is derived from PRD 83 (2011) 012004, appendix
  	//taken from Luciano
  	
	double logp0 = LogLikelihood( GetBestFitParameters() ); //likelihood at the mode

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

  	vector<double> parameters = GetBestFitParameters();
	double constant = parameters[0];
	double amplitude = parameters[1];
	double halflife = parameters[2];

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
void TimeAlpha::WriteOutput( string outputfilename )
{
	// Write Livetime Fraction

	// Write Data

	// Write Model
	const vector<double> BestFit = GetBestFitParametersMarginalized();

	vector<double> * BestFitParameters = new vector<double>;
	vector<double> * BestFitParameterErrors1s = new vector<double>;
//	vector<double> * BestFitParameterErrors2s = new vector<double>;
//	vector<double> * BestFitParameterErrors3s = new vector<double>;

	for( uint i = 0; i < GetNParameters(); i++ )
	{
		BCH1D * parHisto = MCMCGetH1Marginalized( i );

		double value =  parHisto->GetMode();

		BestFitParameters -> push_back( value );

		double xmin, xmax;
	    parHisto->GetSmallestInterval( xmin, xmax, 0.683 );
	    BestFitParameterErrors1s -> push_back( (xmax - xmin) / 2. );

	    cout << "Par0 " << value << "(" << xmin << "," << xmax << ")" << endl;

	    /*
	    parHisto->GetSmallestInterval( xmin, xmax, 0.954 );
	    BestFitParameterErrors2s -> push_back( (xmax - xmin) / 2. );
	    parHisto->GetSmallestInterval( xmin, xmax, 0.997 );
	    BestFitParameterErrors3s -> push_back( (xmax - xmin) / 2. );
	    */
	}

	TF1 * modelfunction = new TF1( "modelfunction", "[0]+[1]*exp(-x*log(2.)/[2])", fHMinimum, fHMaximum );
	modelfunction -> SetParameters( BestFitParameters->at(0), BestFitParameters->at(1), BestFitParameters->at(2) );
	modelfunction -> SetParError( 0, BestFitParameterErrors1s->at(0) );
	modelfunction -> SetParError( 1, BestFitParameterErrors1s->at(1) );
	modelfunction -> SetParError( 2, BestFitParameterErrors1s->at(2) );

	cout << fHMinimum << " " << fHMaximum << endl;
	cout << BestFitParameters->at(0) << " " <<	BestFitParameters->at(1) << " "<<	BestFitParameters->at(2) << " " << endl;
	cout << BestFitParameterErrors1s->at(0) << " " <<	BestFitParameterErrors1s->at(1) << " "<<	BestFitParameterErrors1s->at(2) << " " << endl;

	TCanvas * c = new TCanvas("c1");
	modelfunction->Draw();

	TH1D * expected = new TH1D( "histo_model", "histo_model", fNBins,  fHMinimum, fHMaximum );

	double BinWidth = (fHMaximum - fHMinimum) / fNBins;

	for( int b = 0; b < fNBins; b++)
	{
		double t1 = fHMinimum + b * BinWidth;
		double t2 = t1 + BinWidth;

		double tau = BestFitParameters->at(2) / log(2.);

		double lambda = BestFitParameters->at(0) * BinWidth
				- BestFitParameters->at(1) * tau
				* ( exp( -t2/tau ) - exp( -t1/tau ) );
		lambda *= fVLiveTimeFraction.at(b);

		expected -> SetBinContent( b+1, lambda );
	}

	TFile * out = new TFile( outputfilename.c_str(), "RECREATE" );
	c->Write();
	modelfunction -> Write();
	fHTimeAlpha -> Write();
	fHLiveTimeFraction -> Write();
	expected -> Write();
	out -> Close();

	return;
}


// ---------------------------------------------------------
void TimeAlpha::WriteDistributions()
{	
	TFile * test = new TFile( "posteriors.root", "RECREATE" );

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
