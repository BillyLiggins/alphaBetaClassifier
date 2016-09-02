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

#define SSTR( x ) static_cast< std::ostringstream & >( \
								( std::ostringstream() << std::dec << x ) ).str()


void CutFinder::FindBoundary(){
				std::cout<<"Inside FindBoundary"<<std::endl;

				double startingCut=0;

				// if( acceptor->GetPID().compare("beta")==0){
				if( acceptor->GetPID()=="beta"){
						startingCut=-1000;
				}else if( rejector->GetPID()=="alpha"){

						startingCut=4000;
				}else{
							std::cout<<"Cutter needs a PID"<< std::endl; 
				}

				double step=0.1;
				double cutValue=startingCut;
				double accept=0;
				double mistagged=0;
				double rejection=0;
				double minBin=acceptor->GetHist()->GetXaxis()->GetXmin();
				std::cout<<"Stopper"<<std::endl;
				double maxBin=acceptor->GetHist()->GetXaxis()->GetXmax();
				double sliceWidth=acceptor->GetHist()->GetXaxis()->GetBinWidth(1);
				std::cout<<"minBin = "<<minBin<<std::endl;
				std::cout<<"maxBin = "<<maxBin<<std::endl;
				std::cout<<"sliceWidth = "<<sliceWidth<<std::endl;
				TAxis* accept_x=acceptor->GetHist()->GetXaxis();
				TAxis* reject_x=rejector->GetHist()->GetXaxis();

				double acceptor_numEV_complete= acceptor->GetHist()->GetEntries();
				double rejector_numEV_complete= rejector->GetHist()->GetEntries();
				std::cout<<"Number of entries in acceptor = "<< acceptor_numEV_complete<<std::endl;
				std::cout<<"Number of entries in rejector = "<< rejector_numEV_complete<<std::endl;
				double totalMistaging = 0;
				

				std::cout<<"threshold = "<<this->GetThreshold()<<std::endl;
				for(double energy =minBin; energy<maxBin-sliceWidth;energy+=sliceWidth){

								// std::cout<<"Finding cut in energy slice "<<energy<<"< E < "<<energy+sliceWidth<<std::endl;

								//These histograms contain the sliced projections
								TH1D* slice_acceptor=acceptor->GetHist()->ProjectionY(SSTR(energy).c_str(),accept_x->FindBin(energy),accept_x->FindBin(energy+sliceWidth));
								TH1D* slice_rejector=rejector->GetHist()->ProjectionY(("Po"+SSTR(energy)).c_str(),reject_x->FindBin(energy),reject_x->FindBin(energy+sliceWidth));

								//These two numbers are the num of events before cuts
								double acceptor_numEV=slice_acceptor->GetEntries();
								double rejector_numEV=slice_rejector->GetEntries();

								TAxis* acceptor_slice_x=slice_acceptor->GetXaxis();
								TAxis* rejector_slice_x=slice_rejector->GetXaxis();

								while(accept<threshold){

												//This loop finds the cut value to a certain threshold is achevied, it does so by incrementing by 'step'

												double acceptor_remain=slice_acceptor->Integral(acceptor_slice_x->FindBin(startingCut),acceptor_slice_x->FindBin(cutValue+=step));
												accept= acceptor_remain/acceptor_numEV;	
								}

								mistagged=slice_rejector->Integral(rejector_slice_x->FindBin(startingCut),rejector_slice_x->FindBin(cutValue));
								this->EnterMistaggedValue(mistagged);
								this->EnterMistaggedError(sqrt(mistagged));
								totalMistaging+= mistagged;
								rejection= rejector_numEV/(mistagged+0.00001);

								this->EnterRejectionValue(rejection);

								if(mistagged>0){
												this->EnterRejectionError(rejection/sqrt(mistagged));
								}else{
												this->EnterRejectionError(0);
								}

								this->EnterCutValues(cutValue);
								this->EnterEnergyValues(energy + sliceWidth/2);
								this->EnterEnergyError(sliceWidth/2);

								cutValue=startingCut;
								accept=0;
				}

								TGraph* cutBoundary= new TGraph(this->GetEnergyValuesVector().size(),&(this->GetEnergyValuesVector())[0],&(this->GetCutValuesVector())[0]);
								TF1 *f_E = new TF1("f_E", "[1]*x +[0]",0,2.5);
								cutBoundary->Fit(f_E);

								TFile* function = new TFile("cutBoundary.root","recreate");
								f_E->Write();
								function->Close();

								double Intercept = f_E->GetParameter(0);
								double Gradient = f_E->GetParameter(1);

								rejector->SetIntercept(Intercept);
								rejector->SetGradient(Gradient);

								acceptor->SetIntercept(Intercept);
								acceptor->SetGradient(Gradient);

								std::cout<<"Total number of mistagged alphas= "<< totalMistaging<<std::endl;
								std::cout<<"Total number of rejector_numEV_complete (alpha) = "<<rejector_numEV_complete<<std::endl;
								std::cout<<"Percentage of alpha left after E-dep cut = "<<totalMistaging*100/rejector_numEV_complete<<std::endl;

				}


void CutFinder::Holder(){

				ofstream File_bi;
				ofstream File_po;

				File_bi.open("Classifier_data_bi.dat");
				File_po.open("Classifier_data_po.dat");

				File_bi << "mcEdepQuenched,"<<"posr," << "BerekelyAlphaBeta" <<  std::endl;
				File_po << "mcEdepQuenched,"<<"posr," << "BerekelyAlphaBeta" << std::endl;

				UTIL* util = new UTIL();


				// Cutter* rejector = new Cutter("alpha");
				// Cutter* acceptor = new Cutter("beta");

				rejector->SetHistLimits(100,0,2.5,260,-160,100);
				acceptor->SetHistLimits(100,0,2.5,260,-160,100);

				acceptor->FillCutter("/data/snoplus/liggins/year1/fitting/fitting/alphaSims/output_electron/ntuple","electron",File_bi);
				acceptor->PrintHist();

				rejector->FillCutter("/data/snoplus/liggins/year1/fitting/fitting/alphaSims/output/ntuple","alpha",File_po);
				rejector->PrintHist();

				File_bi.close();
				File_po.close();


				// SetThreshold(0.9999);
				this->FindBoundary();


}

void CutFinder::Plotting(){

				std::vector<double> energyValues = this->GetEnergyValuesVector();
				std::vector<double> energyErrors = this->GetEnergyErrorVector();

				std::vector<double> rejection = this->GetRejectionValuesVector();
				std::vector<double> rejectionErrors = this->GetRejectionErrorVector();

				std::vector<double> mistagged = this->GetMistaggedValueVector();
				std::vector<double> mistaggedErrors = this->GetMistaggedErrorVector();

				std::vector<double> cutValues = this->GetCutValuesVector();


				std::cout<<"Number of acceptor entries = "<< acceptor->GetNumberOfEntries()	<<std::endl;
				std::cout<<"Number of  rejector entries = "<< rejector->GetNumberOfEntries()	<<std::endl;

				/************************************************************
				 *******************Plotting*********************************
				 */


				{
								TGraphErrors* cutBoundary = new TGraphErrors(energyValues.size(),&energyValues[0],&rejection[0],&energyErrors[0],&rejectionErrors[0]);
								TCanvas * c1 = new TCanvas();
								c1->cd();
								c1->SetLogy();
								cutBoundary->SetTitle(Form("Rejection across energy with %.0f\% #beta retention",this->GetThreshold()*100));

								cutBoundary->GetXaxis()->SetTitle("mcEdepQuenched (MeV)");
								cutBoundary->GetYaxis()->SetTitle("Rejection Power");
								cutBoundary->SetMaximum(10e5);
								cutBoundary->Draw("ap");
								c1->Print("plots/RejectionAcrossEnergy.png");
				}
				
				{
								TGraphErrors* mistagged_events= new TGraphErrors(energyValues.size(),&energyValues[0],&mistagged[0],&energyErrors[0],&mistaggedErrors[0]);
								TCanvas * c1 = new TCanvas();
								c1->cd();
								c1->SetLogy();
								mistagged_events->SetTitle(Form("Mistagged events across energy with %.0f\% #beta retention",this->GetThreshold()*100));
								mistagged_events->GetXaxis()->SetTitle("mcEdepQuenched (MeV)");
								mistagged_events->GetYaxis()->SetTitle("Mistagged events");
								mistagged_events->SetMaximum(10e5);
								mistagged_events->Draw("ap");
								c1->Print("plots/MistaggedAcrossEnergy.png");
				}

				{
								TGraph* cutBoundary= new TGraph(energyValues.size(),&energyValues[0],&cutValues[0]);
								TF1 *f_E = new TF1("f_E", "[1]*x +[0]",0,2.5);
								f_E->SetLineColor(1);
								f_E->SetLineStyle(2);
								f_E->SetLineWidth(3);
								cutBoundary->Fit(f_E);


								TCanvas * c1= new TCanvas();
								c1->cd();

								rejector->GetHist()->Draw();
								acceptor->GetHist()->Draw("same");
								// cutBoundary->Draw("same .");
								f_E->Draw("same");

								TLegend* t2 = new TLegend( 0.11, 0.11, 0.31, 0.31 );
								t2->AddEntry(acceptor->GetHist(), "#beta 's","f");
								t2->AddEntry(rejector->GetHist(), "#alpha 's","f");
								t2->AddEntry( f_E, "99%","l");
								t2->Draw();

								c1->Print("plots/cutBoundary.png");
				}
}


