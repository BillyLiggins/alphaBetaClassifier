// This file has a main function that goes and find cuts for the berkeleyAlphaBeta across energy and radius (independently). 
// You were trying to get a radial fit done.
// TODOs:
//
// 	1) combine the radial and energy varibles in to one TH3D 
// 	2) Clean it 
//

#include <RAT/DU/DSReader.hh>
#include <RAT/DS/Entry.hh>

#include <TGraph2D.h>
#include <TMath.h>
#include <string>
#include <TMultiGraph.h>
#include <TH2D.h>
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
#include <TStyle.h>
#include <TPad.h>


#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <iostream>
#include <fstream>
#include <numeric>

#define SSTR( x ) static_cast< std::ostringstream & >( \
								( std::ostringstream() << std::dec << x ) ).str()

using namespace std;
std::string partflag="";

vector<std::string> glob( const string& path, const string& start )
{
				vector<string> result;
				TSystemDirectory dir(path.c_str(), path.c_str());
				TList *files = dir.GetListOfFiles();
				if (files) {
								TSystemFile *file;
								TString fname;
								TIter next(files);
								while ((file=(TSystemFile*)next())) {
												fname = file->GetName();
												if (!file->IsDirectory() && fname.EndsWith(".root") && fname.BeginsWith( start ) ) {
																stringstream fullPath; fullPath << path << "/" << fname.Data();
																result.push_back(fullPath.str());
												}
								}
				}
				return result;
}


void FillHist(TFile* file,TH1D * hist,TH2D * hist2,TH1D *radial,TH2D *compareRadial,ofstream& outputfile){

				TTree* Tree = (TTree*) file->Get("output");
				Double_t Bab, mcEdepQuenched,posr,mcPosr;
				Bool_t  Qfit;
				Int_t pdg1, pdg2, evIndex;
				Int_t parentpdg1,parentpdg2;

				Tree->SetBranchAddress("pdg1",&pdg1);
				Tree->SetBranchAddress("pdg2",&pdg2);
				Tree->SetBranchAddress("parentpdg1",&parentpdg1);
				Tree->SetBranchAddress("parentpdg2",&parentpdg2);

				Tree->SetBranchAddress("berkeleyAlphaBeta",&Bab);
				Tree->SetBranchAddress("mcEdepQuenched",&mcEdepQuenched);
				Tree->SetBranchAddress("posr",&posr);
				Tree->SetBranchAddress("mcPosr",&mcPosr);
				Tree->SetBranchAddress("fitValid",&Qfit);
				Tree->SetBranchAddress("evIndex",&evIndex);
				Int_t n = (Int_t)Tree->GetEntries();
				Int_t code;

				for( Int_t i =0;i<n;i++){
								Tree->GetEntry(i);
								if( Qfit && evIndex==0 && mcPosr<4000  ){
												// if( Qfit && evIndex==0 ){
												/* if( Qfit){ */
												hist->Fill(Bab);
												hist2->Fill(mcEdepQuenched,Bab);
												radial->Fill(posr);
												compareRadial->Fill(posr,Bab);
												outputfile<< mcEdepQuenched<<","<<posr<<","<< Bab << endl;
								}
								}

								file->Close();
								}


								// void findCutsEnergy(TH2D *compareBi210,TH2D *comparePo210,vector<double>& energyValues, vector<double>& cutValues,vector<double>& Rejection_values){

double findCutsEnergy(TH2D *compareBi210,TH2D *comparePo210,vector<double>& energyValues, vector<double>& cutValues,vector<double>& Rejection_values,vector<double>& Rejection_errors,vector<double>& energy_errors){
				/* This function should take a 2d histrogram split it into energy 
				 * strips and then find a cut value which retain ~99% of the signal.
				 */

				double startingCut=4000;
				double step=1;
				double cutValue=startingCut;
				double accept=0;
				double mistagged=0;
				double rejection=0;
				double minBin=compareBi210->GetXaxis()->GetXmin();
				double maxBin=compareBi210->GetXaxis()->GetXmax();
				double sliceWidth=compareBi210->GetXaxis()->GetBinWidth(1);
				cout<<"minBin = "<<minBin<<endl;
				cout<<"maxBin = "<<maxBin<<endl;
				cout<<"sliceWidth = "<<sliceWidth<<endl;
				double threshold=0.995;
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

				double TotalMistag=0;

				for(double energy =minBin; energy<maxBin-sliceWidth;energy+=sliceWidth){
								/* for(double energy =0; energy<2.5;energy+=0.1){ */
								cout<<"Energy = "<< energy<<endl;

								//These histograms contain the sliced projections
								TH1D* slice_bi=compareBi210->ProjectionY(SSTR(energy).c_str(),bi_x->FindBin(energy),bi_x->FindBin(energy+sliceWidth));
								TH1D* slice_po=comparePo210->ProjectionY(("Po"+SSTR(energy)).c_str(),po_x->FindBin(energy),po_x->FindBin(energy+sliceWidth));

								//These two numbers are the num of events before cuts
								double bi_numEV=slice_bi->GetEntries();
								double po_numEV=slice_po->GetEntries();

								TAxis* bi_slice_x=slice_bi->GetXaxis();
								TAxis* po_slice_x=slice_po->GetXaxis();

								while(accept<threshold){

												//This loop finds the cut value to a certain threshold is achevied, it does so by incrementing by 'step'
												double po_remain=slice_po->Integral(po_slice_x->FindBin(cutValue-=step),po_slice_x->FindBin(startingCut));
												//double po_remain=slice_po->Integral(po_slice_x->FindBin(startingCut),po_slice_x->FindBin(startingCut+step));
												// cout<<"cutvalue = "<<cutValue<<endl;
												// cout<<"Number of entries remaining in bi slice = "<< bi_remain<<endl;
												// cout<<"Number of entries remaining in po slice = "<< po_remain<<endl;
												accept=po_remain/po_numEV;	
												// cout<<"Accepted fraction = "<<accept<<endl;
								}

								mistagged=slice_bi->Integral(bi_slice_x->FindBin(cutValue),bi_slice_x->FindBin(startingCut));

								rejection= bi_numEV/(mistagged+0.00001);
								// rejection= po_numEV/(mistagged);
								// cout<<mistagged<<" events were mistagged. Rejection= "<<rejection<<"."<<endl;
								Rejection_values.push_back(rejection);

								if(mistagged>0){
												Rejection_errors.push_back(rejection/sqrt(mistagged));
								}else{
												Rejection_errors.push_back(0);
								}

								cutValues.push_back(cutValue);
								energyValues.push_back(energy+ sliceWidth/2);
								energy_errors.push_back(sliceWidth/2);

								if(false){

												TCanvas* c1 = new TCanvas();
												c1->cd();
												double maximum;
												// if(energy== minBin) maximum = ( slice_bi->GetMaximum() > slice_po->GetMaximum() ) ? slice_bi->GetMaximum() : slice_po->GetMaximum();
												if(energy== minBin) maximum = 4000;
												std::cout << "maximum = "<< maximum << std::endl;
												TLegend* leg = new TLegend( 0.61, 0.61, 0.9, 0.9 );
												TLine* cutLine = new TLine( cutValue, slice_bi->GetMinimum(), cutValue, maximum );
												cutLine->SetLineWidth(2);
												cutLine->SetLineStyle(4);
												slice_bi->SetMaximum(maximum);
												slice_bi->Draw();
												slice_po->Draw("same");
												cutLine->Draw("same");
												// c1->Print(("plots/energySlices/Energy_slice_"+SSTR(energy)+"_E_"+SSTR(energy+sliceWidth)+"_alpha.png").c_str());
												leg->AddEntry(slice_bi,"#beta","f");
												leg->AddEntry(slice_po,"#alpha","f");
												leg->AddEntry(cutLine,"99%","l");
												leg->Draw();
												std::cout << "energy = "<<energy<<" energy+sliceWidth = "<<energy+sliceWidth << std::endl;
												c1->Print(Form("plots/energySlices/Energy_slice_%.3f_E_%.3f_alpha.png",energy,energy+sliceWidth));
												c1->Print(Form("plots/energySlices/forTalk/Energy_slice_%.3f_E_%.3f_alpha.png",energy,energy+sliceWidth));

								}
								cutValue=startingCut;
								accept=0;
								TotalMistag=TotalMistag+mistagged;
								cout<<"Number of mistagged betas = "<< mistagged<<endl;
				}

				cout<<"Totel number of mistagged betas = "<<TotalMistag<<endl;
				cout<<"Total number of bi_numEV_complete = "<<bi_numEV_complete<<endl;
				cout<<"Percentage of beta leaf after E-dep cut = "<< TotalMistag*100/bi_numEV_complete<<endl;

				return TotalMistag;
				}

void forStephie( std::vector<string> biFileList, std::vector<string> poFileList, std::vector<double> energyValues, std::vector<double> bi_cuts){

				TGraph* cutGraph = new TGraph(bi_cuts.size(),&energyValues[0],&bi_cuts[0]);
				TF1 *f = new TF1("f_Rad", "[1]*x +[0]");
				cutGraph->Fit(f);

				TCanvas* c1 = new TCanvas();
				cutGraph->Draw();
				c1->Print("forStephie/cutLine.png");


				//Extracting fit Babmeters.
				double P0 = f->GetParameter(0);
				double P1 = f->GetParameter(1);
				std::cout << "Intercept P0 = "<< f->GetParameter(0) << std::endl;
				std::cout << "Gradient P1 = "<< f->GetParameter(1) << std::endl;

				//Now look at the number of betas that get through the cut.
				double TotalNumberOfEntries, TotalNumberOfMistags;

				for ( int i=0; i<biFileList.size();i++){


								TFile * file= TFile::Open(biFileList[i].c_str());	
								TTree* Tree = (TTree*) file->Get("output");
								Double_t Bab, mcEdepQuenched,posr,mcPosr;
								Bool_t  Qfit;
								Int_t evIndex;

								Tree->SetBranchAddress("berkeleyAlphaBeta",&Bab);
								Tree->SetBranchAddress("mcEdepQuenched",&mcEdepQuenched);
								Tree->SetBranchAddress("posr",&posr);
								Tree->SetBranchAddress("mcPosr",&mcPosr);
								Tree->SetBranchAddress("fitValid",&Qfit);
								Tree->SetBranchAddress("evIndex",&evIndex);
								Int_t n = (Int_t)Tree->GetEntries();

								for( Int_t i =0;i<n;i++){
												Tree->GetEntry(i);
												if( Qfit && evIndex==0 && mcPosr<4000  ){
																// if( Qfit && evIndex==0 ){
																/* if( Qfit){ */

																TotalNumberOfEntries++;			
																if(Bab> (P0+P1*mcEdepQuenched)) TotalNumberOfMistags++;


												}
												}

												file->Close();
												}

												cout<<"TotalNumberOfEntries = "<<TotalNumberOfEntries<<endl;
												cout<<"TotalNumberOfMistags= "<<TotalNumberOfMistags<<endl;
												double betaPercentRemaining=100*TotalNumberOfMistags/TotalNumberOfEntries;
												double percentError= betaPercentRemaining * sqrt(1/TotalNumberOfEntries+1/TotalNumberOfMistags);

												cout<<"Percent of beta remaining = "<<betaPercentRemaining<<"+/-"<<percentError<<endl;
												//Now look at the number of alpha that get through the cut.
												TotalNumberOfEntries=0;
												TotalNumberOfMistags=0;

												for ( int i=0; i<poFileList.size();i++){


																TFile * file= TFile::Open(poFileList[i].c_str());	
																TTree* Tree = (TTree*) file->Get("output");
																Double_t Bab, mcEdepQuenched,posr,mcPosr;
																Bool_t  Qfit;
																Int_t evIndex;

																Tree->SetBranchAddress("berkeleyAlphaBeta",&Bab);
																Tree->SetBranchAddress("mcEdepQuenched",&mcEdepQuenched);
																Tree->SetBranchAddress("posr",&posr);
																Tree->SetBranchAddress("mcPosr",&mcPosr);
																Tree->SetBranchAddress("fitValid",&Qfit);
																Tree->SetBranchAddress("evIndex",&evIndex);
																Int_t n = (Int_t)Tree->GetEntries();

																for( Int_t i =0;i<n;i++){
																				Tree->GetEntry(i);
																				if( Qfit && evIndex==0 && mcPosr<4000  ){
																								// if( Qfit && evIndex==0 ){
																								/* if( Qfit){ */

																								TotalNumberOfEntries++;			
																								if(Bab> (P0+P1*mcEdepQuenched)) TotalNumberOfMistags++;


																				}
																				}

																				file->Close();
																				}

																				cout<<"TotalNumberOfEntries = "<<TotalNumberOfEntries<<endl;
																				cout<<"TotalNumberOfMistags= "<<TotalNumberOfMistags<<endl;
																				double alphaPercentRemaining=100*TotalNumberOfMistags/TotalNumberOfEntries;
																				double percentError_alpha= alphaPercentRemaining * sqrt(1/TotalNumberOfEntries+1/TotalNumberOfMistags);

																				cout<<"Percent of alpha remaining = "<<alphaPercentRemaining<<"+/-"<<percentError_alpha<<endl;
																}


int cutter_alpha(){
				//====================================================================	
				// double binNum=8;
				// double binNum=25;
				double binNum=100;
				gStyle->SetOptStat(0);
				TH1D *hBi210   = new TH1D("hBi210","berkeleyAlphaBeta",100,-100,100);
				hBi210->SetLineColor(kBlue);hBi210->SetLineWidth(3);
				hBi210->SetFillColor(kBlue);
				TH1D *hPo210   = new TH1D("hPo210","berkeleyAlphaBeta",100,-1000,1000);
				hPo210->SetLineColor(kRed);hPo210->SetLineWidth(3);
				TH1D *hBlank   = new TH1D("hBlank","berkeleyAlphaBeta",100,-100,100);

				TH2D* compareBi210   = new TH2D("compareBi210","berkeleyAlphaBeta",binNum,0,2.5,260,-160,100);
				compareBi210->SetLineColorAlpha(kBlue,0.5);
				compareBi210->SetLineWidth(3);
				compareBi210->SetMarkerColorAlpha(kBlue,0.5);
				compareBi210->SetFillColorAlpha(kBlue,0.5);
				TH2D* comparePo210   = new TH2D("comparePo212","berkeleyAlphaBeta",binNum,0,2.5,260,-160,100);
				comparePo210->SetLineColorAlpha(kRed,0.5);
				comparePo210->SetLineWidth(3);
				comparePo210->SetMarkerColorAlpha(kRed,0.5);
				comparePo210->SetFillColorAlpha(kRed,0.5);

				//==========================Radial Histrograms ========================
				TH1D *hEle_radial   = new TH1D("hEle_radial","berkeleyAlphaBeta",6000,0,6000);
				hEle_radial->SetLineColor(kBlue);hEle_radial->SetLineWidth(3);
				hEle_radial->SetFillColor(kBlue);
				TH1D *hAlpha_radial   = new TH1D("hAlpha_radial","berkeleyAlphaBeta",6000,0,6000);
				hAlpha_radial->SetLineColor(kRed);hAlpha_radial->SetLineWidth(3);

				TH2D* compareEle_radial  = new TH2D("compareEle_radial","berkeleyAlphaBeta",6000,0,6000,460,-260,200);
				compareEle_radial->SetLineColor(kBlue);compareEle_radial->SetLineWidth(3);
				compareEle_radial->SetMarkerColor(kBlue);
				compareEle_radial->SetFillColor(kBlue);
				TH2D* compareAlpha_radial   = new TH2D("compareAlpha_radial","berkeleyAlphaBeta",6000,0,6000,460,-260,200);
				compareAlpha_radial->SetLineColor(kRed);compareAlpha_radial->SetLineWidth(3);
				compareAlpha_radial->SetMarkerColor(kRed);
				compareAlpha_radial->SetFillColor(kRed);
				//TH2D *hBlank2   = new TH2D("hBlank2","berkeleyAlphaBeta",100,-100,100);

				TLegend* t1 = new TLegend( 0.6, 0.7, 0.89, 0.88 );
				TLegend* t2 = new TLegend( 0.11, 0.11, 0.31, 0.31 );
				TLegend* t3 = new TLegend( 0.11, 0.11, 0.21, 0.21 );

				ofstream File_bi;
				ofstream File_po;

				File_bi.open("Classifier_data_bi.dat");
				File_po.open("Classifier_data_po.dat");

				File_bi << "mcEdepQuenched,"<<"posr," << "BerekelyAlphaBeta" <<  endl;
				File_po << "mcEdepQuenched,"<<"posr," << "BerekelyAlphaBeta" << endl;

				vector<string> biFileList= glob("/data/snoplus/liggins/year1/fitting/fitting/alphaSims/output_electron/ntuple","electron");
				vector<string> poFileList= glob("/data/snoplus/liggins/year1/fitting/fitting/alphaSims/output/ntuple","alpha");

				for( int i=0; i<biFileList.size(); i++ ){
								TFile * file= TFile::Open(biFileList[i].c_str());	
								partflag="Bi";
								FillHist( file, hBi210,compareBi210,hEle_radial,compareEle_radial,File_bi);
				}

				for( int i=0; i<poFileList.size(); i++ ){
								TFile * file= TFile::Open(poFileList[i].c_str());	
								std::cout<<"Loaded file "<<poFileList[i]<<std::endl;
								partflag="Po";
								FillHist( file, hPo210,comparePo210,hEle_radial,compareAlpha_radial,File_po);
				}
				File_bi.close();
				File_po.close();
				/* TCanvas * tester = new TCanvas(); */
				/* compareBi210->Draw(); */
				/* comparePo210->Scale(0.1); */
				/* comparePo210->Draw(); */
				/* hPo210->Draw(); */

				/* return 0; */

				vector<double> bi_cuts, energyValues,Rejection_energy;
				vector<double> radial_cuts, RadiusValues,Rejection_radius;
				vector<double> Rejection_errors, energy_errors;

				double TotalMistagged = findCutsEnergy(compareBi210,comparePo210,energyValues,bi_cuts,Rejection_energy,Rejection_errors,energy_errors);
				// forStephie( biFileList, poFileList,energyValues,bi_cuts);
				// return 0;
				// double TotalRej=std::accumulate(Rejection_energy.begin(), Rejection_energy.end(), 0);

				// double percentLeft= 100/TotalRej;

				// double RejError=TotalRej/sqrt(TotalMistagged);

				// double percentError= percentLeft*TotalRej/TotalMistagged;
				// percentLeft*rejection_Error[i]/rej[i];

				double thres=0.9;
				int counter=1;
				vector< vector<double> > listOfCuts;

				Bool_t QcutScan=false;
				// Bool_t QcutScan=true;
				vector<TF1*> listOfFunction;


				//------Drawing and saving---------------


				TGraph* cutGraph = new TGraph(bi_cuts.size(),&energyValues[0],&bi_cuts[0]);
				// TF1 *f = new TF1("f_Rad", "[2]*x*x+[1]*x +[0]");
				/* cutGraph->Fit(f); */

				TGraph* cutGraph_radial = new TGraph(radial_cuts.size(),&RadiusValues[0],&radial_cuts[0]);
				TF1 *f_Rad = new TF1("f_Rad", "[1]*x +[0]");
				cutGraph_radial->Fit(f_Rad);



				TCanvas * c2 = new TCanvas();
				c2->cd();
				cout<<"# of entries in hBi210 = "<< hBi210->GetEntries()<<endl;
				cout<<"# of entries in hPo210 = "<< hPo210->GetEntries()<<endl;
				TF1 *f_E = new TF1("f_E", "[1]*x +[0]",0,2.5);
				f_E->SetLineStyle(2);
				f_E->SetLineColor(kBlack);
				f_E->SetLineWidth(2);
				cutGraph->Fit(f_E);
				std::cout << "+++++++Fit+++++" << std::endl;
				std::cout << "P0 = "<< f_E->GetParameter(0) << std::endl;
				std::cout << "P1 = "<< f_E->GetParameter(1) << std::endl;

				compareBi210->SetTitle("Classifiers values across truth quenched energy");
				compareBi210->GetXaxis()->SetTitle("mcEdepQuenched (MeV)");
				compareBi210->GetYaxis()->SetTitle("BerkeleyAlphaBeta");
				compareBi210->Draw();
				comparePo210->Draw("same");
				/* compareBi210->Draw("box"); */
				/* comparePo210->Draw("same box"); */
				/* compareBi210->Draw("surf"); */
				/* comparePo210->Draw("same surf"); */
				/* cutGraph->Draw("AP*1"); */
				// cutGraph->Draw("same");
				f_E->Draw("same");

				t2->AddEntry( compareBi210, "#beta 's","f");
				t2->AddEntry( comparePo210, "#alpha 's","f");
				t2->AddEntry( f_E, "99%","l");
				t2->Draw();
				c2->Print("plots/BerkeleyAlphaBetaVsMCEnergy_withCuts_alpha.png");
				TFile stephieROOT("AlphaRejection.root","RECREATE");
				stephieROOT.cd();
				
				compareBi210->Write();
				comparePo210->Write();
				f_E->Write();
				c2->Write();
				
				stephieROOT.Close();
				// c2->Print("plots/BerkeleyAlphaBetaVsMCEnergy_withCuts_alpha.tex");

				if(QcutScan){
								TCanvas* c_test= new TCanvas();

								TLegend* scanLeg= new TLegend( 0.11, 0.11, 0.31, 0.41 );
								compareBi210->Draw();
								comparePo210->Draw("same");
								scanLeg->AddEntry( compareBi210, "#beta 's","f");
								scanLeg->AddEntry( comparePo210, "#alpha 's","f");

								double percentLeft=90;
								for(int i=0;i<listOfFunction.size();i++){
												listOfFunction[i]->Draw("same");

												// scanLeg->AddEntry(listOfFunction[i], (""+SSTR(0.9/pow(10,i))+"").c_str(),"l");
												scanLeg->AddEntry(listOfFunction[i], (""+SSTR(percentLeft)+"%").c_str(),"l");
												percentLeft=percentLeft+9/pow(10,i);
								}
								scanLeg->Draw();

								c_test->Print("plots/functions_alpha.png");

				}//end of QcutScan


				TCanvas* RejectionEnergy_can = new TCanvas();
				// TGraph* RejectionEnergyGraph = new TGraph(Rejection_energy.size(),&energyValues[0],&Rejection_energy[0]);
				TGraphErrors* RejectionEnergyGraph = new TGraphErrors(Rejection_energy.size(),&energyValues[0],&Rejection_energy[0],&energy_errors[0],&Rejection_errors[0]);
				//TF1 *f = new TF1("f", "[1]*x +[0]");
				//cutGraph->Fit(f);
				//
				// for (int i = 0; i < Rejection_energy.size(); i++) {
				// 				std::cout << "ith rejection value = "<<Rejection_energy[i]<<"+/-"<<Rejection_errors[i] << std::endl;
				// 				std::cout << "ith energy value = "<<energyValues[i]<<"+/-"<<energy_errors[i] << std::endl;
				// }

				RejectionEnergyGraph->SetTitle("Rejection across energy with 99% #beta retention");
				RejectionEnergyGraph->GetXaxis()->SetTitle("mcEdepQuenched (MeV)");
				RejectionEnergyGraph->GetYaxis()->SetTitle("Rejection Power");
				RejectionEnergyGraph->SetMaximum(10e5);
				// RejectionEnergyGraph->SetMaximum(0);
				RejectionEnergyGraph->Draw("ap");
				// RejectionEnergyGraph->Draw("a");
				RejectionEnergy_can->SetLogy();
				RejectionEnergy_can->Print("plots/RejectionAcrossEnergy_alpha.png");
				RejectionEnergy_can->Print("plots/RejectionAcrossEnergy_alpha.tex");



				TFile fileout("BerkeleyAlphaBetaTraditionalCutPlots.root","RECREATE");
				fileout.cd();
				// 	-----Energy-----
				compareBi210->Write();
				comparePo210->Write();
				RejectionEnergyGraph->Write();
				RejectionEnergy_can->Write();
				// 	-----Radial-----
				/* compareEle_radial->Write(); */
				/* compareAlpha_radial->Write(); */
				/* RejectionRadialGraph->Write(); */
				/* RejectionRadial_can->Write(); */
				fileout.Close();

				// cout<<"Total Mistagged = "<<TotalMistagged<<endl;
				// cout<<"99\% alpha retention gives  beta percent left = "<< percentLeft<< "+/-"<< percentError<<endl;	

				return 0;

}


