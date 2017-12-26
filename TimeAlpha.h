// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__TIMEALPHA__H
#define __BAT__TIMEALPHA__H

#include <cstdlib>
#include <iostream>
#include <vector>

#include "TH1D.h"

#include "GETDetector.hh"

#include <BAT/BCModel.h>

// Enumerate variables for detector and model selection
enum tModel { LinAndExp, OnlyExp, OnlyLin };

// This is a TimeAlpha header file.
// Model source code is located in file TimeAlpha/TimeAlpha.cxx

// ---------------------------------------------------------
class TimeAlpha : public BCModel {
 public:

	// Constructor and destructor
	TimeAlpha( const char * name = "TimeAlpha" );
	~TimeAlpha();

    int DefineParametersAndPriors();

    // FIXME read a list of keylists instead
	int ReadDataPhaseII( std::string keylist );
    int ReadDataPhaseI();

    void SetBinningInDays( int b ){ fBinning = b; };
    int GetBinningInDays(){ return fBinning; };

    // FIXME make me private
    void SetNBinsHistograms( int n, double min, double max )
    {
	fNBins = n;
	fHMinimum = min;
	fHMaximum = max;
    };

    int InitializeHistograms();

    // FIXME implement single detector reading
    void SetFittingDataSet( std::string set )
    {
        fDataSet = set;

        if( set == "enrBEGe" || set.find("GD") == 0 )
            { fDetectorType = kIsBEGe; fDetectorEnriched = true; }
        else if( set == "enrCoax" || set.find("ANG") == 0 || set.find("RG") == 0 )
            { fDetectorType = kIsCOAX; fDetectorEnriched = true; }
        else if( set == "natCoax" || set.find("GTF") == 0 )
            { fDetectorType = kIsCOAX; fDetectorEnriched = false; }
        else if( set == "all"){ fDetectorType = kUNKNOWN; fDetectorType = kUNKNOWN; }
        else if( set.find("GD") == 0 || set.find("ANG") == 0 ||
                 set.find("RG") == 0 || set.find("GTF") == 0 )
        {
            fSingleDetectorFit = true;
        }
        else
            std::cout << "Data set unknown " << set << std::endl;
    }

	void SetDetectorType( DetectorType_t det ){ fDetectorType = det; };
    DetectorType_t GetDetectorType(){ return fDetectorType; };

    void SetDetectorEnriched( bool b ){ fDetectorEnriched = b; };
    bool GetDetectorEnriched(){ return fDetectorEnriched; };

	void SetModel( tModel model ){ fModel = model; };

	void SetModelName( std::string name ){ fModelName = name; };
	std::string GetModelName( tModel model ){ return fModelName; };

	void SetNDetectors( int n ){ fNDetectors = n; };
	int GetNDetectors(){ return fNDetectors; };

	// Methods to overload, see file TimeAlpha.cxx
	double LogLikelihood(const std::vector<double> & parameters);
	double EstimatePValue();

	// double LogAPrioriProbability(const std::vector<double> & parameters);

	void WriteOutput( std::string outputfilename );
	void WriteDistributions();

 private:

	DetectorType_t fDetectorType;
    bool fDetectorEnriched;
    bool fSingleDetectorFit;
    std::string fDataSet;
	tModel fModel;

	std::string fModelName;

    // FIXME read nDet from RunConfiguration
	int fNDetectors = 40;

    int fBinning;
	int fNBins;
	double fHMaximum, fHMinimum;

	// Data arrays with time distribution of events and livetime fraction in time (counting pulser events)
	TH1D * fHTimeAlpha;
    TH1D * fHTimeAlpha_fine;
	TH1D * fHLiveTimeFraction;
    TH1D * fHLiveTimeFraction_fine;

	double fLTLimit = 0.;

	std::vector<double> fVTimeAlpha;
	std::vector<double> fVLiveTimeFraction;

	int FillDataArrays();

};
// ---------------------------------------------------------

#endif
