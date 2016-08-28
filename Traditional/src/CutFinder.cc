#include "UTIL.h"
#include "Cutter.h"
#include "CutFinder.h"

//OXSX libaries

// #include <BinnedPdf.h>        
// #include <Rand.h>
// #include <ROOTNtuple.h>        
// #include <BinnedNLLH.h>        
// #include <GridSearch.h>
// #include <MetropolisHastings.h> 
// #include <Minuit.h> 
// #include <Histogram.h>
// #include <PdfConverter.h>        
// #include <CompositePdf.h>        
// #include <Convolution.h>        
// #include <Gaussian.h>        
// #include <DataSetGenerator.h>        
// #include <OXSXDataSet.h>        
// #include <BoolCut.h>        
// #include <BoxCut.h>        

//ROOT libaries

#include <TString.h>
#include <TH1D.h>        
#include <TH2D.h>        
#include <TStyle.h> 
#include <TPaveStats.h> 
#include <TAttFill.h>
#include <TGraph2D.h>
#include <TMath.h>
#include <TMultiGraph.h>
#include <TFrame.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TPaveText.h>
#include <TCanvas.h>
#include <TSystemDirectory.h>
#include <TList.h>
#include <THStack.h>
#include <TLegend.h>
#include <TLine.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TPad.h>

//Standard libaries

#include <iomanip>
#include <sstream>
#include <vector>        
#include <string>        
#include <math.h>	
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <numeric>


void CutFinder::FindBoundary(){

				double startingCut=-500;
				double step=0.1;
				double cutValue=startingCut;
				double accept=0;
				double mistagged=0;
				double rejection=0;
				double minBin=acceptor->->GetXaxis()->GetXmin();
				double maxBin=compareBi210->GetXaxis()->GetXmax();
				double sliceWidth=compareBi210->GetXaxis()->GetBinWidth(1);
				cout<<"minBin = "<<minBin<<endl;
				cout<<"maxBin = "<<maxBin<<endl;
				cout<<"sliceWidth = "<<sliceWidth<<endl;
				double threshold=0.99;
				TAxis* bi_x=compareBi210->GetXaxis();
				TAxis* bi_y=compareBi210->GetYaxis();
				TAxis* po_x=comparePo210->GetXaxis();
				TAxis* po_y=comparePo210->GetYaxis();

				TH1D* complete_bi=compareBi210->ProjectionY("complete_bi",bi_x->FindBin(0.),bi_x->FindBin(4.));
				TH1D* complete_po=comparePo210->ProjectionY("complete_po",po_x->FindBin(0.),po_x->FindBin(4.));
				double bi_numEV_complete= compareBi210->GetEntries();
				double po_numEV_complete= comparePo210->GetEntries();
				cout<<"Number of entries in complete bi = "<< bi_numEV_complete<<endl;
				cout<<"Number of entries in complete po = "<< po_numEV_complete<<endl;
				double totalMistaging = 0;

				for(double energy =minBin; energy<maxBin-sliceWidth;energy+=sliceWidth){
			

}

void CutFinder::FindCutValue(){

			

}

//You have to think about this now. You have two cutter objects which contains stuff for alphas/betas. 
//But you know want to combine these to because they are both needed for the findCutsEnergy().
// void Cutter::findCutsEnergy(){
// 								#<{(| This function should take a 2d histrogram split it into energy 
// 								 * strips and then find a cut value which retain ~99% of the signal.
// 								 |)}>#
// 								//for(double energy =0; energy<3.5;energy+=0.1){
//
// 								double startingCut=-500;
// 								double step=0.1;
// 								double cutValue=startingCut;
// 								double accept=0;
// 								double mistagged=0;
// 								double rejection=0;
// 								double minBin=compareBi210->GetXaxis()->GetXmin();
// 								double maxBin=compareBi210->GetXaxis()->GetXmax();
// 								double sliceWidth=compareBi210->GetXaxis()->GetBinWidth(1);
// 								cout<<"minBin = "<<minBin<<endl;
// 								cout<<"maxBin = "<<maxBin<<endl;
// 								cout<<"sliceWidth = "<<sliceWidth<<endl;
// 								double threshold=0.99;
// 								TAxis* bi_x=compareBi210->GetXaxis();
// 								TAxis* bi_y=compareBi210->GetYaxis();
// 								TAxis* po_x=comparePo210->GetXaxis();
// 								TAxis* po_y=comparePo210->GetYaxis();
//
// 								TH1D* complete_bi=compareBi210->ProjectionY("complete_bi",bi_x->FindBin(0.),bi_x->FindBin(4.));
// 								TH1D* complete_po=comparePo210->ProjectionY("complete_po",po_x->FindBin(0.),po_x->FindBin(4.));
// 								double bi_numEV_complete= compareBi210->GetEntries();
// 								double po_numEV_complete= comparePo210->GetEntries();
// 								cout<<"Number of entries in complete bi = "<< bi_numEV_complete<<endl;
// 								cout<<"Number of entries in complete po = "<< po_numEV_complete<<endl;
// 								double totalMistaging = 0;
//
// 								for(double energy =minBin; energy<maxBin-sliceWidth;energy+=sliceWidth){
// 												#<{(| for(double energy =0; energy<2.5;energy+=0.1){ |)}>#
// 												cout<<"Energy = "<< energy<<endl;
//
// 												//These histograms contain the sliced projections
// 												TH1D* slice_bi=compareBi210->ProjectionY(SSTR(energy).c_str(),bi_x->FindBin(energy),bi_x->FindBin(energy+sliceWidth));
// 												TH1D* slice_po=comparePo210->ProjectionY(("Po"+SSTR(energy)).c_str(),po_x->FindBin(energy),po_x->FindBin(energy+sliceWidth));
//
// 												//These two numbers are the num of events before cuts
// 												double bi_numEV=slice_bi->GetEntries();
// 												double po_numEV=slice_po->GetEntries();
//
// 												TAxis* bi_slice_x=slice_bi->GetXaxis();
// 												TAxis* po_slice_x=slice_po->GetXaxis();
//
// 												while(accept<threshold){
//
// 																//This loop finds the cut value to a certain threshold is achevied, it does so by incrementing by 'step'
// 																double bi_remain=slice_bi->Integral(bi_slice_x->FindBin(startingCut),bi_slice_x->FindBin(cutValue+=step));
// 																//double po_remain=slice_po->Integral(po_slice_x->FindBin(startingCut),po_slice_x->FindBin(startingCut+step));
// 																// cout<<"Number of entries remaining in bi slice = "<< bi_remain<<endl;
// 																//cout<<"Number of entries remaining in po slice = "<< po_remain<<endl;
// 																accept=bi_remain/bi_numEV;	
// 																// cout<<"Accepted fraction = "<<accept<<endl;
// 												}
//
// 												mistagged=slice_po->Integral(po_slice_x->FindBin(startingCut),po_slice_x->FindBin(cutValue));
// 												totalMistaging+= mistagged;
// 												rejection= po_numEV/(mistagged+0.00001);
// 												// rejection= po_numEV/(mistagged);
// 												// cout<<mistagged<<" events were mistagged. Rejection= "<<rejection<<"."<<endl;
// 												Rejection_values.push_back(rejection);
//
// 												if(mistagged>0){
// 																Rejection_errors.push_back(rejection/sqrt(mistagged));
// 												}else{
// 																Rejection_errors.push_back(0);
// 												}
//
// 												cutValues.push_back(cutValue);
// 												energyValues.push_back(energy+ sliceWidth/2);
// 												energy_errors.push_back(sliceWidth/2);
//
// 												if(false){
//
// 																slice_bi->SetLineColorAlpha(kBlue,0.5);
// 																slice_bi->SetMarkerColorAlpha(kBlue,0.5);
// 																slice_bi->SetFillColorAlpha(kBlue,0.5);
//
// 																slice_po->SetLineColorAlpha(kRed,0.5);
// 																slice_po->SetMarkerColorAlpha(kRed,0.5);
// 																slice_po->SetFillColorAlpha(kRed,0.5);
//
// 																TCanvas* c1 = new TCanvas();
// 																c1->cd();
// 																double maximum;
// 																// if(energy== minBin) maximum = ( slice_bi->GetMaximum() > slice_po->GetMaximum() ) ? slice_bi->GetMaximum() : slice_po->GetMaximum();
// 																if(energy== minBin) maximum = 2000;
// 																std::cout << "maximum = "<< maximum << std::endl;
// 																TLegend* leg = new TLegend( 0.70, 0.70, 0.9, 0.9 );
// 																TLine* cutLine = new TLine( cutValue, slice_bi->GetMinimum(), cutValue, maximum );
// 																cutLine->SetLineWidth(2);
// 																cutLine->SetLineStyle(4);
// 																slice_bi->SetMaximum(maximum);
// 																slice_bi->SetTitle("BerkeleyAlphaBeta {0.4<E<0.5 MeV}");
// 																slice_bi->Draw();
// 																slice_po->Draw("same");
// 																cutLine->Draw("same");
// 																// c1->Print(("plots/energySlices/Energy_slice_"+SSTR(energy)+"_E_"+SSTR(energy+sliceWidth)+".png").c_str());
// 																leg->AddEntry(slice_bi,"#beta","f");
// 																leg->AddEntry(slice_po,"#alpha","f");
// 																leg->AddEntry(cutLine,"99%","l");
// 																leg->Draw();
// 																std::cout << "energy = "<<energy<<" energy+sliceWidth = "<<energy+sliceWidth << std::endl;
// 																c1->Print(Form("plots/energySlices/forTalk/Energy_slice_%.3f_E_%.3f.png",energy,energy+sliceWidth));
//
// 												}
// 												cutValue=startingCut;
// 												accept=0;
// 								}
//
// 								cout<<"Totel number of mistagged alphas= "<< totalMistaging<<endl;
// 								cout<<"Total number of po_numEV_complete = "<<po_numEV_complete<<endl;
// 								cout<<"Percentage of beta leaf after E-dep cut = "<<totalMistaging*100/po_numEV_complete<<endl;
//
// 								}
//
