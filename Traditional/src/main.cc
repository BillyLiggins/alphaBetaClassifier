#include "UTIL.h"
#include "Cutter.h"
#include "CutFinder.h"

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
#include <TGraphErrors.h>
#include <TF1.h>

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

void findBoundary();
void scan();
void applyBoundary();

int main(){

				
				findBoundary();
				applyBoundary();
				// scan();

				return 0;
}


void findBoundary(){

				ofstream File_bi;
				ofstream File_po;

				File_bi.open("Classifier_data_bi.dat");
				File_po.open("Classifier_data_po.dat");

				File_bi << "mcEdepQuenched,"<<"posr," << "BerekelyAlphaBeta" <<  std::endl;
				File_po << "mcEdepQuenched,"<<"posr," << "BerekelyAlphaBeta" << std::endl;

				UTIL* util = new UTIL();


				Cutter* alpha = new Cutter("alpha");
				Cutter* beta = new Cutter("beta");

				alpha->SetHistLimits(100,0,2.5,260,-160,100);
				beta->SetHistLimits(100,0,2.5,260,-160,100);

				beta->FillCutter("/data/snoplus/liggins/year1/fitting/fitting/alphaSims/output_electron/ntuple","electron",File_bi);
				beta->PrintHist();

				alpha->FillCutter("/data/snoplus/liggins/year1/fitting/fitting/alphaSims/output/ntuple","alpha",File_po);
				alpha->PrintHist();

				File_bi.close();
				File_po.close();


				CutFinder* cutFinder= new CutFinder(beta,alpha);
				cutFinder->SetThreshold(0.9999);
				cutFinder->FindBoundary();


				std::vector<double> energyValues = alpha->GetEnergyValuesVector();
				std::vector<double> energyErrors = alpha->GetEnergyErrorVector();

				std::vector<double> rejection = alpha->GetRejectionValuesVector();
				std::vector<double> rejectionErrors = alpha->GetRejectionErrorVector();

				std::vector<double> mistagged = alpha->GetMistaggedValueVector();
				std::vector<double> mistaggedErrors = alpha->GetMistaggedErrorVector();

				std::vector<double> cutValues = alpha->GetCutValuesVector();


				std::cout<<"Number of beta entries = "<< beta->GetNumberOfEntries()	<<std::endl;

				/************************************************************
				 *******************Plotting*********************************
				 */


				{
								TGraphErrors* cutBoundary = new TGraphErrors(energyValues.size(),&energyValues[0],&rejection[0],&energyErrors[0],&rejectionErrors[0]);
								TCanvas * c1 = new TCanvas();
								c1->cd();
								c1->SetLogy();
								cutBoundary->SetTitle(Form("Rejection across energy with %.0f\% #beta retention",cutFinder->GetThreshold()*100));

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
								mistagged_events->SetTitle(Form("Mistagged events across energy with %.0f\% #beta retention",cutFinder->GetThreshold()*100));
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

								alpha->GetHist()->Draw();
								beta->GetHist()->Draw("same");
								// cutBoundary->Draw("same .");
								f_E->Draw("same");

								TLegend* t2 = new TLegend( 0.11, 0.11, 0.31, 0.31 );
								t2->AddEntry(beta->GetHist(), "#beta 's","f");
								t2->AddEntry(alpha->GetHist(), "#alpha 's","f");
								t2->AddEntry( f_E, "99%","l");
								t2->Draw();

								c1->Print("plots/cutBoundary.png");
				}
}


void applyBoundary(){

				UTIL* util = new UTIL();

				std::vector<std::string> betaFileList= util->glob("/data/snoplus/liggins/year1/fitting/fitting/alphaSims/output_electron/ntuple","electron");
				std::vector<std::string> alphaFileList= util->glob("/data/snoplus/liggins/year1/fitting/fitting/alphaSims/output/ntuple","alpha");

				// std::vector<std::string> betaFileList= util->glob("/home/billy/workspace/PhD/testData/beta/output/ntuple","electron");
				// std::vector<std::string> alphaFileList= util->glob("/home/billy/workspace/PhD/testData/alpha/output/ntuple","alpha");

				Cutter* alpha = new Cutter("alpha");
				Cutter* beta = new Cutter("beta");

				TFile *cutBoundary = TFile::Open("cutBoundary.root");
				TF1 * f_E =(TF1*) cutBoundary->Get("f_E");	

				double Intercept = f_E->GetParameter(0);
				double Gradient = f_E->GetParameter(1);

				alpha ->SetIntercept(Intercept);
				alpha->SetGradient(Gradient);

				beta->SetIntercept(Intercept);
				beta->SetGradient(Gradient);

				alpha->SetHistLimits(25,0,2.5,260,-160,100);
				beta->SetHistLimits(25,0,2.5,260,-160,100);

				for( int i=0; i<betaFileList.size(); i++ ){
								TFile * file= TFile::Open(betaFileList[i].c_str());	
								beta->ApplyCut(file);
				}

				for( int i=0; i<alphaFileList.size(); i++ ){
								TFile * file= TFile::Open(alphaFileList[i].c_str());	
								alpha->ApplyCut(file);
				}
				
				beta->SetRemainingPercentage((beta->GetRemainingAfterCut())*100/(beta->GetNumberOfEntries()));
				std::cout<<"======================================================="<<std::endl;
				std::cout<<"For the beta Cutter"<<std::endl;
				std::cout<<"Before and after the linear cut:"<<std::endl;
				std::cout<<"Before cut "<< beta->GetNumberOfEntries()<<" entries existed"<<std::endl;
				std::cout<<"After cut "<< beta->GetRemainingAfterCut()<<" entries existed"<<std::endl;
				std::cout<<"Remaining after the being cut away as a percentage = "<<(beta->GetRemainingAfterCut())*100/(beta->GetNumberOfEntries())<< std::endl;
				std::cout<<"======================================================="<<std::endl;

				alpha->SetRemainingPercentage((alpha->GetRemainingAfterCut())*100/(alpha->GetNumberOfEntries()));
				std::cout<<"======================================================="<<std::endl;
				std::cout<<"For the alpha Cutter"<<std::endl;
				std::cout<<"Before and after the linear cut:"<<std::endl;
				std::cout<<"Before cut "<< alpha->GetNumberOfEntries()<<" entries existed"<<std::endl;
				std::cout<<"After cut "<< alpha->GetRemainingAfterCut()<<" entries existed"<<std::endl;
				std::cout<<"Remaining after the being cut away as a percentage = "<<(alpha->GetRemainingAfterCut())*100/(alpha->GetNumberOfEntries())<< std::endl;
				std::cout<<"======================================================="<<std::endl;
}




void scan(){

				UTIL* util = new UTIL();

				Cutter* alpha = new Cutter("alpha");
				Cutter* beta = new Cutter("beta");

				alpha->SetHistLimits(100,0,2.5,260,-160,100);
				beta->SetHistLimits(100,0,2.5,260,-160,100);

				beta->FillCutter("/data/snoplus/liggins/year1/fitting/fitting/alphaSims/output_electron/ntuple","electron");

				alpha->FillCutter("/data/snoplus/liggins/year1/fitting/fitting/alphaSims/output/ntuple","alpha");

				std::vector<double>  effListError, effList;
				std::vector<double>  remainingListError, remainingList;

				for( double eff= 0.99; eff<1;eff+=0.0001){
								CutFinder* cutFinder= new CutFinder(beta,alpha);
								cutFinder->SetThreshold(eff);
								cutFinder->FindBoundary();

								beta->ApplyBoundary("/data/snoplus/liggins/year1/fitting/fitting/alphaSims/output_electron/ntuple","electron");

								alpha->ApplyBoundary("/data/snoplus/liggins/year1/fitting/fitting/alphaSims/output/ntuple","alpha");

								beta->SetRemainingPercentage( (beta->GetRemainingAfterCut())*100/(beta->GetNumberOfEntries()) );
								std::cout<<"======================================================="<<std::endl;
								std::cout<<"For the beta Cutter"<<std::endl;
								std::cout<<"Before and after the linear cut:"<<std::endl;
								std::cout<<"Before cut "<< beta->GetNumberOfEntries()<<" entries existed"<<std::endl;
								std::cout<<"After cut "<< beta->GetRemainingAfterCut()<<" entries existed"<<std::endl;
								std::cout<<"Remaining after the being cut away as a percentage = "<<(beta->GetRemainingAfterCut())*100/(beta->GetNumberOfEntries())<< std::endl;
								std::cout<<"======================================================="<<std::endl;

								alpha->SetRemainingPercentage( (alpha->GetRemainingAfterCut())*100/(alpha->GetNumberOfEntries()) );
								std::cout<<"======================================================="<<std::endl;
								std::cout<<"For the alpha Cutter"<<std::endl;
								std::cout<<"Before and after the linear cut:"<<std::endl;
								std::cout<<"Before cut "<< alpha->GetNumberOfEntries()<<" entries existed"<<std::endl;
								std::cout<<"After cut "<< alpha->GetRemainingAfterCut()<<" entries existed"<<std::endl;
								std::cout<<"Remaining after the being cut away as a percentage = "<<(alpha->GetRemainingAfterCut())*100/(alpha->GetNumberOfEntries())<< std::endl;
								std::cout<<"======================================================="<<std::endl;

								effList.push_back(beta->GetRemainingPercentage());
								effListError.push_back(beta->GetRemainingPercentageError());
								remainingList.push_back(alpha->GetRemainingPercentage());
								remainingListError.push_back(alpha->GetRemainingPercentageError());
				}

				for(int i =0 ; i< effList.size();i++){

								std::cout<<"beta eff = "<<effList[i]<<" gives "<<remainingList[i]<<" alpha particles remaining."<<std::endl;
				}


				
				{
								TCanvas* c1 =new TCanvas();
								TGraphErrors* fractionGraph = new TGraphErrors(effList.size(),&effList[0],&remainingList[0],&effListError[0],&remainingListError[0]);
								fractionGraph ->SetTitle("Surviving #alpha 's over accepted #beta 's {mcPosr < 6000 }");
								fractionGraph ->GetXaxis()->SetTitle("Remaining #beta 's % ");
								fractionGraph ->GetYaxis()->SetTitle("Remaining #alpha 's % ");
								fractionGraph ->GetYaxis()->SetTitleOffset(1.4);
								fractionGraph ->Draw("ap");
								c1->SetGrid();

								c1->Print("plots/sig_sack_alpha_highRad_6000.png");
								c1->Print("plots/sig_sack_alpha_highRad_6000.tex");
				}
}
