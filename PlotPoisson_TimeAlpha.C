////////////////////////////////////
/////  N.Becerici-Schmidt   ////////
//////// 10 October 2012 ///////////
////////////////////////////////////

#include <TROOT.h>
#include <TUnixSystem.h>

#include <iostream>
#include <iomanip>

#include <TFile.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>

#include <RA_PoiStat.h>
#include <fstream>

using namespace std;

void PlotPoisson_TimeAlpha( string filename, string outfilename = "default.root", double Xlow = 0., double Xup = 180. )
{

// style settings
//.............................................................................//
  gROOT->SetStyle("Plain");
  gStyle->SetLabelFont(42,"xy");
  gStyle->SetTitleFont(42,"xy");
  gStyle->SetOptFit(111111);
  gStyle->SetOptStat(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLabelSize(0.07,"xy");
  gStyle->SetTitleSize(0.07,"xy");
  gStyle->SetLabelOffset(0.02,"x");
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(1);
  gStyle->SetTitleOffset(1.1,"x");
  //   gStyle->SetTitleOffset(1.0,"x");
  gStyle->SetTitleOffset(0.5,"y");

  gStyle->SetPadTopMargin(0.03);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.11);

  gROOT->ForceStyle();
//.............................................................................//

  TFile* fitfile = new TFile( filename.c_str() );

  // read in the histograms

  TH1D* hdata=(TH1D*)fitfile->Get("HTimeAlpha");
  hdata->SetFillColor(kGray+2);
  hdata->SetLineColor(kGray+2);
  hdata->GetYaxis()->SetTitle("counts/day");
  hdata->GetXaxis()->SetTitle("time [days]");

  TF1 * fModelFunction = (TF1*) fitfile->Get("modelfunction");
  TH1D * hLiveTimeFraction = (TH1D*) fitfile->Get("HLiveTimeFraction");

  vector<TGraph*> NotConsidered;
  for( int i = 1; i <= hLiveTimeFraction->GetNbinsX(); i++ )
  {
	if( hLiveTimeFraction->GetBinContent(i) <= 0.15 )
	{
		double up = 100., down = 0.;
		double left = hLiveTimeFraction->GetBinLowEdge(i);
		double right = left + hLiveTimeFraction->GetBinWidth(1);

		TGraph * g = new TGraph( 4 );
		g -> SetPoint(0, left, down);
		g -> SetPoint(1, right, down);
		g -> SetPoint(2, right, up);
		g -> SetPoint(3, left, up);
		g -> SetFillStyle( 3002 );
		g -> SetFillColor( kBlack );

		NotConsidered.push_back(g);
	}
  }

  TH1D * hExpectedEvents = (TH1D*) fitfile->Get("histo_model");
  hExpectedEvents->SetLineColor(1);
  hExpectedEvents->SetLineWidth(2);
  hExpectedEvents->SetLineStyle(1);

  TCanvas * c1 = new TCanvas("c1", "c1", 5, 5, 600, 500);

  TPad * pad1 = new TPad( "pad1" , "pad1" , 0.02, 0.40, 0.98, 0.98);
  TPad * pad2 = new TPad("pad2", "The pad with the histogram",0.02,0.02,0.98,0.40);
  pad1 -> SetGrid();
  pad1 -> Draw();
  pad2 -> Draw();

  pad1->cd();

  fModelFunction -> GetYaxis() -> SetRangeUser( 0, 20 );
  fModelFunction -> GetYaxis() -> SetTitle("cts / day");
  fModelFunction -> GetYaxis() -> SetTitleOffset(0.5);
  fModelFunction -> Draw();
  c1->Update();

  cout << "Uymax: " << pad1 -> GetUymax() << endl;

  double rightmax = 1.1 * hLiveTimeFraction -> GetBinContent( hLiveTimeFraction -> GetMaximumBin() );
  double scale = pad1 -> GetUymax() / rightmax;
  hLiveTimeFraction -> SetLineColor( kAzure );
  hLiveTimeFraction -> SetFillColor( kAzure );
  hLiveTimeFraction -> SetFillStyle( 3005 );
  hLiveTimeFraction -> Scale( scale );
  hLiveTimeFraction -> Draw("same");

  TGaxis *axis = new TGaxis( pad1 -> GetUxmax(), pad1 -> GetUymin(), pad1 -> GetUxmax(), pad1 -> GetUymax(), 0, rightmax, 510, "+L" );
  axis->SetLineColor( kAzure );
  axis->SetTextColor( kAzure );
  axis->SetLabelColor( kAzure );
  axis->SetTitle("f_{LT}");
  axis->SetLabelSize( fModelFunction -> GetYaxis() -> GetLabelSize() );
  axis->SetTitleSize( fModelFunction -> GetYaxis() -> GetTitleSize() );
  axis->SetLabelFont( fModelFunction -> GetYaxis() -> GetLabelFont() );
  axis->SetTitleFont( fModelFunction -> GetYaxis() -> GetTitleFont() );
  axis->Draw();
/*
  for( int n = 0; n < NotConsidered.size(); n++ )
  {
	NotConsidered.at(n)->Draw( "fsame" );
  }
*/
//  TLine * fLTThreshold = new TLine( 0, 0.15 * scale, 180., 0.15 * scale );
//  fLTThreshold -> Draw( "same" );

  TLegend * myLeg1 = new TLegend( 0.31,0.75,0.42,0.91,"" );
  myLeg1 -> SetBorderSize(0);
  myLeg1 -> SetFillColor(kWhite);
  myLeg1 -> SetTextSize(0.07);
  myLeg1 -> SetTextFont(42);
  myLeg1 -> AddEntry( fModelFunction,"model","L");
  myLeg1 -> AddEntry( hLiveTimeFraction,"livetime fraction","L");
//  myLeg1 -> AddEntry( NotConsidered.at(0),"Skipped Bins","F");
//  myLeg1 -> AddEntry( fLTThreshold,"f_{LT} > 15% threshold","F");
  myLeg1 -> Draw();


// draw color bands
  pad2 -> cd();
  pad2 -> SetTickx(1);
  pad2 -> SetBottomMargin(0.25);
  pad2 -> SetTopMargin(0);

  TH1D * hmodel1 = new TH1D();
  hExpectedEvents->Copy(*hmodel1);
  hmodel1->GetXaxis()->SetRangeUser(Xlow, Xup);

  TH1D * hdata1 = new TH1D();
  hdata->Copy(*hdata1);
  hdata1->GetXaxis()->SetRangeUser(Xlow, Xup);

  hdata1->GetYaxis()->SetTitle("data/model ratio ");
  hdata1->GetXaxis()->SetTitle("time [days]");
  hdata1->SetTitleOffset(0.30,"y");
  hdata1->SetLabelSize(0.11,"xy");
  hdata1->SetTitleSize(0.11,"xy");
  hdata1->GetYaxis()->SetDecimals(1);
  hdata1->GetYaxis()->SetNdivisions(505);

  hdata1->GetXaxis()->SetNdivisions(105);
  hdata1->SetMarkerStyle(24);
  hdata1->SetMarkerSize(0.7);

  TLegend * myLeg11 = new TLegend(0.15,0.83,0.25,0.96,"");
  myLeg11->SetBorderSize(0);
  myLeg11->SetFillColor(kWhite);
  myLeg11->SetTextSize(0.095);
  myLeg11->SetTextFont(42);
  myLeg11->AddEntry(hdata1,"data","P");
  myLeg11->AddEntry(hmodel1,"model","l");

  TLegend * myLeg12 = new TLegend(0.38,0.83,0.48,0.96,"");
  myLeg12->SetBorderSize(0);
  myLeg12->SetFillColor(kWhite);
  myLeg12->SetTextSize(0.095);
  myLeg12->SetTextFont(42);
  TH1D*h1=new TH1D();
  h1->SetFillColor(kGreen-4);
  myLeg12->AddEntry(h1,"68%","F");

  TLegend * myLeg13 = new TLegend(0.50,0.83,0.60,0.96,"");
  myLeg13->SetBorderSize(0);
  myLeg13->SetFillColor(kWhite);
  myLeg13->SetTextSize(0.095);
  myLeg13->SetTextFont(42);
  TH1D*h2=new TH1D();
  h2->SetFillColor(kYellow-4);
  myLeg13->AddEntry(h2,"95%","F");

  TLegend * myLeg14 = new TLegend(0.62,0.83,0.72,0.96,"");
  myLeg14->SetBorderSize(0);
  myLeg14->SetFillColor(kWhite);
  myLeg14->SetTextSize(0.095);
  myLeg14->SetTextFont(42);
  TH1D*h3=new TH1D();
  h3->SetFillColor(kRed-7);
  myLeg14->AddEntry(h3,"99.7%","F");

  RA_PoiStat * h2_=  new RA_PoiStat();
  h2_->RA_PoiStat::Plot_w3ProbLines_ratio( hmodel1, hdata1, 1., 0., "Smallest", 0.68, 0.95, 0.997, "infinite" );
//  h2_->RA_PoiStat::Plot_w3ProbLines_lin( hmodel1, hdata1, 1., 0., "Smallest", 0.68, 0.95, 0.997, "lin", "infinite" );
/*
  for( int n = 0; n < NotConsidered.size(); n++ )
  {
	NotConsidered.at(n)->Draw( "fsame" );
  }
*/
  myLeg11->Draw();
  myLeg12->Draw();
  myLeg13->Draw();
  myLeg14->Draw();

  TFile * file = new TFile( outfilename.c_str() , "RECREATE");

  c1->Write();

  file->Close();

  outfilename.replace( outfilename.end()-4, outfilename.end(), "svg" );
  c1->SaveAs( outfilename.c_str() );
}
