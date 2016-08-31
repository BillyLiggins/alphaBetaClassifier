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


// void CutFinder::FindCutValue(){
// }

void CutFinder::FindBoundary(){

				double startingCut=0;
				std::cout<<"You got in."<<std::endl;

				std::cout<<"acceptor->GetPID() = "<<acceptor->GetPID()<<std::endl;
				// if( acceptor->GetPID().compare("beta")==0){
				if( acceptor->GetPID()=="beta"){
						startingCut=-500;
				}else if( rejector->GetPID()=="alpha"){

						startingCut=1000;
				}else{
							std::cout<<"Cutter needs a PID"<< std::endl; 
				}

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

								std::cout<<"Finding cut in energy slice "<<energy<<"< E < "<<energy+sliceWidth<<std::endl;

								//These histograms contain the sliced projections
								TH1D* slice_acceptor=acceptor->GetHist()->ProjectionY(SSTR(energy).c_str(),accept_x->FindBin(energy),accept_x->FindBin(energy+sliceWidth));
								TH1D* slice_rejector=rejector->GetHist()->ProjectionY(("Po"+SSTR(energy)).c_str(),reject_x->FindBin(energy),reject_x->FindBin(energy+sliceWidth));

								//These two numbers are the num of events before cuts
								double acceptor_numEV=slice_acceptor->GetEntries();
								double rejector_numEV=slice_rejector->GetEntries();

								TAxis* acceptor_slice_x=slice_acceptor->GetXaxis();
								TAxis* rejector_slice_x=slice_rejector->GetXaxis();

								while(accept<threshold){

				// std::cout<<"You got in."<<std::endl;
												//This loop finds the cut value to a certain threshold is achevied, it does so by incrementing by 'step'

												double acceptor_remain=slice_acceptor->Integral(acceptor_slice_x->FindBin(startingCut),acceptor_slice_x->FindBin(cutValue+=step));

												// std::cout<<"Number of entries remaining in bi slice = "<< bi_remain<<std::endl;
												//std::cout<<"Number of entries remaining in po slice = "<< po_remain<<std::endl;
												accept= acceptor_remain/acceptor_numEV;	
												// std::cout<<"Accepted fraction = "<<accept<<std::endl;
								}

								mistagged=slice_rejector->Integral(rejector_slice_x->FindBin(startingCut),rejector_slice_x->FindBin(cutValue));
								rejector->EnterMistaggedValue(mistagged);
								rejector->EnterMistaggedError(sqrt(mistagged));
								totalMistaging+= mistagged;
								rejection= rejector_numEV/(mistagged+0.00001);

								// std::cout<<mistagged<<" events were mistagged. Rejection= "<<rejection<<"."<<std::endl;
								rejector->EnterRejectionValue(rejection);

								if(mistagged>0){
												rejector->EnterRejectionError(rejection/sqrt(mistagged));
								}else{
												rejector->EnterRejectionError(0);
								}

								rejector->EnterCutValues(cutValue);
								rejector->EnterEnergyValues(energy + sliceWidth/2);
								rejector->EnterEnergyError(sliceWidth/2);

								if(false){

												slice_acceptor->SetLineColorAlpha(kBlue,0.5);
												slice_acceptor->SetMarkerColorAlpha(kBlue,0.5);
												slice_acceptor->SetFillColorAlpha(kBlue,0.5);

												slice_rejector->SetLineColorAlpha(kRed,0.5);
												slice_rejector->SetMarkerColorAlpha(kRed,0.5);
												slice_rejector->SetFillColorAlpha(kRed,0.5);

												TCanvas* c1 = new TCanvas();
												c1->cd();
												double maximum;
												// if(energy== minBin) maximum = ( slice_acceptor->GetMaximum() > slice_rejector->GetMaximum() ) ? slice_acceptor->GetMaximum() : slice_rejector->GetMaximum();
												if(energy== minBin) maximum = 2000;
												std::cout << "maximum = "<< maximum << std::endl;
												TLegend* leg = new TLegend( 0.70, 0.70, 0.9, 0.9 );
												TLine* cutLine = new TLine( cutValue, slice_acceptor->GetMinimum(), cutValue, maximum );
												cutLine->SetLineWidth(2);
												cutLine->SetLineStyle(4);
												slice_acceptor->SetMaximum(maximum);
												slice_acceptor->SetTitle("BerkeleyAlphaBeta {0.4<E<0.5 MeV}");
												slice_acceptor->Draw();
												slice_rejector->Draw("same");
												cutLine->Draw("same");
												// c1->Print(("plots/energySlices/Energy_slice_"+SSTR(energy)+"_E_"+SSTR(energy+sliceWidth)+".png").c_str());
												leg->AddEntry(slice_acceptor,"#beta","f");
												leg->AddEntry(slice_rejector,"#alpha","f");
												leg->AddEntry(cutLine,"99%","l");
												leg->Draw();
												std::cout << "energy = "<<energy<<" energy+sliceWidth = "<<energy+sliceWidth << std::endl;
												c1->Print(Form("plots/energySlices/forTalk/Energy_slice_%.3f_E_%.3f.png",energy,energy+sliceWidth));

								}
								cutValue=startingCut;
								accept=0;
				}

				std::cout<<"Totel number of mistagged alphas= "<< totalMistaging<<std::endl;
				std::cout<<"Total number of rejector_numEV_complete = "<<rejector_numEV_complete<<std::endl;
				std::cout<<"Percentage of beta leaf after E-dep cut = "<<totalMistaging*100/rejector_numEV_complete<<std::endl;


				}




//You have to think about this now. You have two cutter objects which contains stuff for alphas/betas. 
//But you know want to combine these to because they are both needed for the findCutsEnergy().
