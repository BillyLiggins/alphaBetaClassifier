#include "UTIL.h"
#include "Cutter.h"

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



void Cutter::FillHist(TFile* file,ofstream& outputfile)
{

				TTree* Tree = (TTree*) file->Get("output");
				Double_t para, mcEdepQuenched,posr,mcPosr;
				Bool_t  Qfit;
				Int_t pdg1, pdg2, evIndex;
				Int_t parentpdg1,parentpdg2;

				Tree->SetBranchAddress("pdg1",&pdg1);
				Tree->SetBranchAddress("pdg2",&pdg2);
				Tree->SetBranchAddress("parentpdg1",&parentpdg1);
				Tree->SetBranchAddress("parentpdg2",&parentpdg2);

				Tree->SetBranchAddress("berkeleyAlphaBeta",&para);
				Tree->SetBranchAddress("mcEdepQuenched",&mcEdepQuenched);
				Tree->SetBranchAddress("posr",&posr);
				Tree->SetBranchAddress("mcPosr",&mcPosr);
				Tree->SetBranchAddress("fitValid",&Qfit);
				Tree->SetBranchAddress("evIndex",&evIndex);
				Int_t n = (Int_t)Tree->GetEntries();

				for( Int_t i =0;i<n;i++){
								Tree->GetEntry(i);

								if( Qfit && evIndex==0 && mcPosr<4000  ){

												BabVsEnergy->Fill(mcEdepQuenched,para);

												outputfile<< mcEdepQuenched<<","<<posr<<","<< para << std::endl;
								}
				}

				file->Close();
}

void Cutter::SetHist(){

				BabVsEnergy = new TH2D("BabVsEnergy","berkeleyAlphaBeta",25,0,2.5,260,-60,40);
				if(PID=="alpha"){
								BabVsEnergy->SetLineColorAlpha(kRed,0.2);
								BabVsEnergy->SetLineWidth(3);
								BabVsEnergy->SetMarkerColorAlpha(kRed,0.2);
								BabVsEnergy->SetFillColorAlpha(kRed,0.2);
				}else if (PID=="beta"){
								BabVsEnergy->SetLineColorAlpha(kBlue,0.2);
								BabVsEnergy->SetLineWidth(3);
								BabVsEnergy->SetMarkerColorAlpha(kBlue,0.2);
								BabVsEnergy->SetFillColorAlpha(kBlue,0.2);
				}else{
								std::cout<<"You need to supply a PID"<<std::endl;
				}
				BabVsEnergy->SetTitle("BerekelyAlphaBeta across energy");
				BabVsEnergy->GetXaxis()->SetTitle("mcEdepQuenched (Mev)");
}

void Cutter::PrintHist(){

				TCanvas * c1= new TCanvas();
				c1->cd();
				BabVsEnergy->Draw();
				c1->Print(Form("plots/BabVsMCEdepQuenched_%s.png",PID.c_str()));

}

