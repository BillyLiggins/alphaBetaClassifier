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

				double startingCut=0;

				// if( acceptor->GetPID().compare("beta")==0){
				std::cout<<"got here"<<std::endl;
				if( acceptor->GetPID()=="beta"){
						startingCut=-1000;
				}else if( rejector->GetPID()=="alpha"){

						startingCut=4000;
				}else{
							std::cout<<"Cutter needs a PID"<< std::endl; 
				}
				std::cout<<"after here"<<std::endl;

				double step=0.1;
				double cutValue=startingCut;
				double accept=0;
				double mistagged=0;
				double rejection=0;
				double minBin=acceptor->GetHist()->GetXaxis()->GetXmin();
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

								// std::cout<<"energy= "<<energy<<std::endl;
								// std::cout<<"threshold = "<<threshold<<std::endl;
								
								while(accept<threshold){

												//This loop finds the cut value to a certain threshold is achevied, it does so by incrementing by 'step'

												double acceptor_remain=slice_acceptor->Integral(acceptor_slice_x->FindBin(startingCut),acceptor_slice_x->FindBin(cutValue+=step));
												accept= acceptor_remain/acceptor_numEV;	
								}

								mistagged=slice_rejector->Integral(rejector_slice_x->FindBin(startingCut),rejector_slice_x->FindBin(cutValue));
								rejector->EnterMistaggedValue(mistagged);
								rejector->EnterMistaggedError(sqrt(mistagged));
								totalMistaging+= mistagged;
								rejection= rejector_numEV/(mistagged+0.00001);

								rejector->EnterRejectionValue(rejection);

								if(mistagged>0){
												rejector->EnterRejectionError(rejection/sqrt(mistagged));
								}else{
												rejector->EnterRejectionError(0);
								}

								rejector->EnterCutValues(cutValue);
								rejector->EnterEnergyValues(energy + sliceWidth/2);
								rejector->EnterEnergyError(sliceWidth/2);

								cutValue=startingCut;
								accept=0;
				}
								TGraph* cutBoundary= new TGraph(rejector->GetEnergyValuesVector().size(),&(rejector->GetEnergyValuesVector())[0],&(rejector->GetCutValuesVector())[0]);
								TF1 *f_E = new TF1("f_E", "[1]*x +[0]",0,2.5);
								cutBoundary->Fit(f_E);

								{
												TCanvas* c1 = new TCanvas();

												f_E->SetLineColor(kBlack);
												acceptor->GetHist()->Draw();
												rejector->GetHist()->Draw("same");
												f_E->Draw("same");

												c1->Print("testLine.png");
								}



								//This is for passing the gradient and intercept to the cutters.
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


