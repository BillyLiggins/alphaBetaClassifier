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




void Cutter::FindNEntries(std::string folder, std::string fileStart){

				UTIL* util = new UTIL();

				std::vector<std::string> FileList= util->glob(folder ,fileStart);


				for( int i=0; i<FileList.size(); i++ ){
								TFile * file= TFile::Open(FileList[i].c_str());	
								// std::cout<<"Loaded file "<<FileList[i]<<std::endl;
								TTree* Tree = (TTree*) file->Get("output");
								Double_t berkeleyAlphaBeta, mcEdepQuenched,posr,mcPosr;
								Bool_t  Qfit;
								Int_t evIndex;
								Int_t parentpdg1,parentpdg2;

								Tree->SetBranchAddress("berkeleyAlphaBeta",&berkeleyAlphaBeta);
								Tree->SetBranchAddress("mcEdepQuenched",&mcEdepQuenched);
								Tree->SetBranchAddress("posr",&posr);
								Tree->SetBranchAddress("mcPosr",&mcPosr);
								Tree->SetBranchAddress("fitValid",&Qfit);
								Tree->SetBranchAddress("evIndex",&evIndex);
								Int_t n = (Int_t)Tree->GetEntries();

								for( Int_t i =0;i<n;i++){
												Tree->GetEntry(i);

												if( Qfit && evIndex==0 && mcPosr<radialCut  ){
																numberOfEntries++;
													}
												}
								file->Close();

				}
}

void Cutter::FindTotalNEntries(std::string folder, std::string fileStart){

				UTIL* util = new UTIL();

				std::vector<std::string> FileList= util->glob(folder ,fileStart);

				for( int i=0; i<FileList.size(); i++ ){
								TFile * file= TFile::Open(FileList[i].c_str());	
								TTree* Tree = (TTree*) file->Get("output");
								Int_t n = (Int_t)Tree->GetEntries();

								TotalNumberOfEntries+=n;
								file->Close();

				}
}

void Cutter::FillHist(TFile* file,ofstream& outputfile)
{

				TTree* Tree = (TTree*) file->Get("output");
				Double_t berkeleyAlphaBeta, mcEdepQuenched,posr,mcPosr;
				Bool_t  Qfit;
				Int_t evIndex;
				Int_t parentpdg1,parentpdg2;


				Tree->SetBranchAddress("berkeleyAlphaBeta",&berkeleyAlphaBeta);
				Tree->SetBranchAddress("mcEdepQuenched",&mcEdepQuenched);
				Tree->SetBranchAddress("posr",&posr);
				Tree->SetBranchAddress("mcPosr",&mcPosr);
				Tree->SetBranchAddress("fitValid",&Qfit);
				Tree->SetBranchAddress("evIndex",&evIndex);
				Int_t n = (Int_t)Tree->GetEntries();

				for( Int_t i =0;i<n;i++){
								Tree->GetEntry(i);

								if( Qfit && evIndex==0 && mcPosr<radialCut  ){
												BabVsEnergy->Fill(mcEdepQuenched,berkeleyAlphaBeta);
												numberOfEntries++;
												outputfile<< mcEdepQuenched<<","<<posr<<","<< berkeleyAlphaBeta<<","<< Qfit <<","<< evIndex << std::endl;
								}
				}

				file->Close();
}

void Cutter::FillHist(TFile* file)
{

				TTree* Tree = (TTree*) file->Get("output");
				Double_t berkeleyAlphaBeta, mcEdepQuenched,posr,mcPosr;
				Bool_t  Qfit;
				Int_t evIndex;

				Tree->SetBranchAddress("berkeleyAlphaBeta",&berkeleyAlphaBeta);
				Tree->SetBranchAddress("mcEdepQuenched",&mcEdepQuenched);
				Tree->SetBranchAddress("posr",&posr);
				Tree->SetBranchAddress("mcPosr",&mcPosr);
				Tree->SetBranchAddress("fitValid",&Qfit);
				Tree->SetBranchAddress("evIndex",&evIndex);
				Int_t n = (Int_t)Tree->GetEntries();

				for( Int_t i =0;i<n;i++){
								Tree->GetEntry(i);

								if( Qfit && evIndex==0 && mcPosr<radialCut){

												BabVsEnergy->Fill(mcEdepQuenched,berkeleyAlphaBeta);
												numberOfEntries++;
								}
				}

				file->Close();
}




void Cutter::SetHistLimits(double Ebins,double ELow,double EHigh,double BabBins, double BabLow, double BabHigh){

				BabVsEnergy = new TH2D("BabVsEnergy","berkeleyAlphaBeta",Ebins,ELow,EHigh,BabBins,BabLow,BabHigh);
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

void Cutter::ApplyCut(TFile * file){
				TTree* Tree = (TTree*) file->Get("output");
				Double_t berkeleyAlphaBeta, mcEdepQuenched,posr,mcPosr;
				Bool_t  Qfit;
				Int_t  evIndex;

				Tree->SetBranchAddress("berkeleyAlphaBeta",&berkeleyAlphaBeta);
				Tree->SetBranchAddress("mcEdepQuenched",&mcEdepQuenched);
				Tree->SetBranchAddress("posr",&posr);
				Tree->SetBranchAddress("mcPosr",&mcPosr);
				Tree->SetBranchAddress("fitValid",&Qfit);
				Tree->SetBranchAddress("evIndex",&evIndex);
				Int_t n = (Int_t)Tree->GetEntries();

				for( Int_t i =0;i<n;i++){
								Tree->GetEntry(i);

								if( Qfit && evIndex==0 && mcPosr<radialCut){
												// NumberOfEntries is given in the FillHist method.
												if(berkeleyAlphaBeta<(gradient*mcEdepQuenched+intercept)) remainingAfterCut++;
								}
				}
				file->Close();
}

void Cutter::FillCutter(std::string folder, std::string fileStart){

				UTIL* util = new UTIL();

				// std::vector<std::string> FileList= util->glob("/data/snoplus/liggins/year1/fitting/fitting/alphaSims/output_electron/ntuple","electron");
				std::vector<std::string> FileList= util->glob(folder ,fileStart);

				for( int i=0; i<FileList.size(); i++ ){
								TFile * file= TFile::Open(FileList[i].c_str());	
								// std::cout<<"Loaded file "<<FileList[i]<<std::endl;
								this->FillHist(file);
								file->Close();
				}


}

void Cutter::FillCutter(std::string folder, std::string fileStart,ofstream&  File_bi){

				UTIL* util = new UTIL();

				std::vector<std::string> FileList= util->glob(folder ,fileStart);


				for( int i=0; i<FileList.size(); i++ ){
								TFile * file= TFile::Open(FileList[i].c_str());	
								// std::cout<<"Loaded file "<<FileList[i]<<std::endl;
								this->FillHist(file,File_bi);
								file->Close();
				}

}

void Cutter::ApplyBoundary(std::string folder,std::string fileStart){


				UTIL* util = new UTIL();

				std::vector<std::string> FileList= util->glob(folder,fileStart);

				for( int i=0; i<FileList.size(); i++ ){
								TFile * file= TFile::Open(FileList[i].c_str());	
								this->ApplyCut(file);
								file->Close();
				}

				this->SetRemainingPercentage( (this->GetRemainingAfterCut())*100/(this->GetNumberOfEntries()) );
}

void Cutter::FindRejection(std::string folder,std::string fileStart){


				UTIL* util = new UTIL();

				std::vector<std::string> FileList= util->glob(folder,fileStart);
				std::vector<double>  E, Eerror;

				double minBin=this->GetHist()->GetXaxis()->GetXmin();
				double maxBin=this->GetHist()->GetXaxis()->GetXmax();
				double sliceWidth=this->GetHist()->GetXaxis()->GetBinWidth(1);
				for (double energy = minBin; energy <maxBin; energy+=sliceWidth) {
								E.push_back(energy+sliceWidth/2);
								Eerror.push_back(sliceWidth/2);
								std::cout<<"energy = "<< energy <<std::endl;
								double N=0;
								double N_remain=0;
								for( int i=0; i<FileList.size(); i++ ){
											
							
								TFile * file= TFile::Open(FileList[i].c_str());	
								// this->ApplyCut(file);
								TTree* Tree = (TTree*) file->Get("output");
								Double_t berkeleyAlphaBeta, mcEdepQuenched,posr,mcPosr;
								Bool_t  Qfit;
								Int_t evIndex;

								Tree->SetBranchAddress("berkeleyAlphaBeta",&berkeleyAlphaBeta);
								Tree->SetBranchAddress("mcEdepQuenched",&mcEdepQuenched);
								Tree->SetBranchAddress("posr",&posr);
								Tree->SetBranchAddress("mcPosr",&mcPosr);
								Tree->SetBranchAddress("fitValid",&Qfit);
								Tree->SetBranchAddress("evIndex",&evIndex);
								Int_t n = (Int_t)Tree->GetEntries();

								for( Int_t i =0;i<n;i++){
												Tree->GetEntry(i);

												if( Qfit && evIndex==0 && mcPosr<radialCut){
																if(energy<mcEdepQuenched && mcEdepQuenched<energy+sliceWidth){
																				N++;
																				//This cut is if below the line then you pass. So beta should be high and alpha should be low.
																				if(berkeleyAlphaBeta<(gradient*mcEdepQuenched+intercept)) N_remain++;

																				}
																}
												}


								file->Close();
												}//end of file n loop.

								
								rejection.push_back(N/N_remain);
								rejection_errors.push_back(N/N_remain*sqrt((N+N_remain)/(N*N_remain)));

								}//energy of filelist loop.

				for (int i = 0; i < rejection.size(); ++i) {
								std::cout<<"rejection = "<<rejection[i]<<" +/- "<<rejection_errors[i]<<std::endl;
				}

				{
								TCanvas* c1 =new TCanvas();
								TGraphErrors* rejectionGraph = new TGraphErrors(rejection.size(),&E[0],&rejection[0],&Eerror[0],&rejection_errors[0]);
								rejectionGraph->SetTitle(Form("%s rejection across energy {mcPosr < %0.f }",this->GetPID().c_str(),this->GetRadialCut()));
								rejectionGraph->GetXaxis()->SetTitle("mcEdepQuenched (MeV) ");
								rejectionGraph->GetYaxis()->SetTitle("Rejection");
								rejectionGraph->GetYaxis()->SetTitleOffset(1.4);
								rejectionGraph->Draw("ap");
								c1->SetLogy();
								c1->SetGrid();

								c1->Print(Form("plots/%s_rejection_across_energy_%f.png",this->GetPID().c_str(),this->GetRadialCut()));
								c1->Print(Form("plots/%s_rejection_across_energy_%f.tex",this->GetPID().c_str(),this->GetRadialCut()));


								TFile outFile("plots/Rejection.root","update");
								// fractionGraph.SetName(Form("radCut=%f",radCut));
								rejectionGraph->Write();
								outFile.Close();
				}
}
