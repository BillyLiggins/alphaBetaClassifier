// This file will go and find a cut boundary for BerekelyAlphaBeta classifier and then computes the full rejection of that boundary. Then the plan is to find a rejection against effiency.
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
// #include "util.h"

#define SSTR( x ) static_cast< std::ostringstream & >( \
								( std::ostringstream() << std::dec << x ) ).str()


std::string partflag="";

vector<string> glob( const string& path, const string& start )
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


void FillHist(TFile* file,TH2D * hist2, const double highRad=6000){

				TTree* Tree = (TTree*) file->Get("output");
				Double_t berkeleyAlphaBeta, mcEdepQuenched,mcPosr;
				Bool_t  Qfit;
				Int_t pdg1, pdg2, evIndex;
				Int_t parentpdg1,parentpdg2;

				Tree->SetBranchAddress("pdg1",&pdg1);
				Tree->SetBranchAddress("pdg2",&pdg2);
				Tree->SetBranchAddress("parentpdg1",&parentpdg1);
				Tree->SetBranchAddress("parentpdg2",&parentpdg2);

				Tree->SetBranchAddress("berkeleyAlphaBeta",&berkeleyAlphaBeta);
				Tree->SetBranchAddress("mcEdepQuenched",&mcEdepQuenched);
				Tree->SetBranchAddress("mcPosr",&mcPosr);
				Tree->SetBranchAddress("fitValid",&Qfit);
				Tree->SetBranchAddress("evIndex",&evIndex);
				Int_t n = (Int_t)Tree->GetEntries();
				Int_t code;

				for( Int_t i =0;i<n;i++){
								Tree->GetEntry(i);
								// if( Qfit && evIndex==0 && mcPosr<4000){
								if( Qfit && evIndex==0 && mcPosr<highRad){
												hist2->Fill(mcEdepQuenched,berkeleyAlphaBeta);
								}
				}

				}

				double findCutsEnergy(TH2D *compareBi210,TH2D *comparePo210, double threshold,std::vector<double>& rej_error, TH1D* Mistaggedhist){
								/* This function returns the total rejection for a threshold.
								*/

								vector<double> Rejection_values;
								double startingCut=-1000;
								double step=0.001;
								double cutValue=startingCut;
								double accept=0;
								double mistagged=0;//number of alphas in the beta accpeted window.
								double rejection=0;
								double minBin=compareBi210->GetXaxis()->GetXmin();
								double maxBin=compareBi210->GetXaxis()->GetXmax();
								double sliceWidth=compareBi210->GetXaxis()->GetBinWidth(1);
								TAxis* bi_x=compareBi210->GetXaxis();
								TAxis* bi_y=compareBi210->GetYaxis();
								TAxis* po_x=comparePo210->GetXaxis();
								TAxis* po_y=comparePo210->GetYaxis();

								TH1D* complete_bi=compareBi210->ProjectionY("complete_bi",bi_x->FindBin(0.),bi_x->FindBin(4.));
								TH1D* complete_po=comparePo210->ProjectionY("complete_po",po_x->FindBin(0.),po_x->FindBin(4.));
								double bi_numEV_complete= compareBi210->GetEntries();
								double po_numEV_complete= comparePo210->GetEntries();
								// cout<<"Number of entries in complete bi = "<< bi_numEV_complete<<endl;
								// cout<<"Number of entries in complete po = "<< po_numEV_complete<<endl;

								double TotalMistagged=0;
								double TotalNumOfPo=0;

								for(double energy =minBin; energy<maxBin;energy+=sliceWidth){
												// cout<<"Energy = "<< energy<<endl;

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
																double bi_remain=slice_bi->Integral(bi_slice_x->FindBin(startingCut),bi_slice_x->FindBin(cutValue+=step));
																accept=bi_remain/bi_numEV;	
												}
												mistagged=slice_po->Integral(po_slice_x->FindBin(startingCut),po_slice_x->FindBin(cutValue));
												Mistaggedhist->Fill(mistagged);
												TotalMistagged+=mistagged;
												TotalNumOfPo+=po_numEV;

												rejection= po_numEV/(mistagged);
												// if(po_numEV>0 && mistagged>0){
												// 				rejection= po_numEV/(mistagged);
												// }else{
												// 				rejection= 0;
												// }

												std::cout << "acceptance = "<<threshold<< std::endl;
												std::cout <<energy<< " <E< "<<energy+sliceWidth<< std::endl;
												std::cout << "po_numEV = "<<po_numEV << std::endl;
												std::cout << "mistagged = "<<mistagged<< std::endl;

												Rejection_values.push_back(rejection);

												if(false){
																TCanvas* c1 = new TCanvas();
																c1->cd();
																//TLine* cutLine = new TLine( cutValue, slice_bi->GetYaxis()->GetXmin(), cutValue, slice_bi->GetYaxis()->GetXmax() );
																TLine* cutLine = new TLine( cutValue, slice_bi->GetMinimum(), cutValue, slice_bi->GetMaximum() );
																cutLine->SetLineWidth(2);
																cutLine->SetLineStyle(4);
																slice_bi->SetTitle(("Energy slice "+SSTR(energy)+" < E < "+SSTR(energy+sliceWidth)+" { #beta efficiency = "+SSTR(threshold)+"}").c_str());
																// slice_bi->Draw();
																// slice_po->Draw("same");
																slice_po->Draw();
																cutLine->Draw("same");
																c1->Print(("plots/energySlice/Energy_slice_betaEfficiency_"+SSTR(threshold)+"_"+SSTR(energy)+"_E_"+SSTR(energy+sliceWidth)+".png").c_str());
												}
												cutValue=startingCut;
												accept=0;
								}

								double TotalRej=std::accumulate(Rejection_values.begin(), Rejection_values.end(), 0);

								// rej_error.push_back(rejection*sqrt((TotalNumOfPo+TotalMistagged)/(TotalNumOfPo*TotalMistagged)));
								rej_error.push_back(TotalRej/sqrt(TotalMistagged));

								return TotalRej;
								Rejection_values.clear();
				}

				int scan(const double highRad=6000){
								//====================================================================	
								// gStyle->SetOptStat(0);
								double eBins=26;
								double eLow=0;
								double eHigh=2.6;
								double BabBins=260;
								double BabLow=-160;
								double BabHigh=100;
								double ebinWidth=(eHigh-eLow)/eBins;

								TH2D* compareBi210   = new TH2D("compareBi210","berkeleyAlphaBeta",eBins,eLow,eHigh,BabBins,BabLow,BabHigh);
								compareBi210->SetLineColor(kBlue);compareBi210->SetLineWidth(3);
								compareBi210->SetMarkerColor(kBlue);
								compareBi210->SetFillColor(kBlue);

								TH2D* comparePo210   = new TH2D("comparePo210","berkeleyAlphaBeta",eBins,eLow,eHigh,BabBins,BabLow,BabHigh);
								comparePo210->SetLineColor(kRed);comparePo210->SetLineWidth(3);
								comparePo210->SetMarkerColor(kRed);
								comparePo210->SetFillColor(kRed);

								TH1D* mistaggedHist= new TH1D("mistaggedHist","berkeleyAlphaBeta",1000,0,10000);
								mistaggedHist->SetLineColorAlpha(kBlue,0.2);
								mistaggedHist->SetLineWidth(3);
								mistaggedHist->SetFillColorAlpha(kBlue,0.2);
								//==========================Radial Histrograms ========================


								TLegend* t1 = new TLegend( 0.6, 0.7, 0.89, 0.88 );
								TLegend* t2 = new TLegend( 0.11, 0.11, 0.21, 0.21 );
								TLegend* t3 = new TLegend( 0.11, 0.11, 0.21, 0.21 );

								vector<string> biFileList= glob("/data/snoplus/liggins/year1/fitting/fitting/alphaSims/output_electron/ntuple","electron");
								vector<string> poFileList= glob("/data/snoplus/liggins/year1/fitting/fitting/alphaSims/output/ntuple","alpha");

								for( int i=0; i<biFileList.size(); i++ ){
												TFile * file= TFile::Open(biFileList[i].c_str());	
												partflag="Bi";
												FillHist(file,compareBi210,highRad);
												file->Close();
								}

								for( int i=0; i<poFileList.size(); i++ ){
												TFile * file= TFile::Open(poFileList[i].c_str());	
												partflag="Po";
												FillHist(file,comparePo210,highRad);
												file->Close();
								}

								vector<double> bi_cuts, energyValues,Rejection_energy;
								vector<double> radial_cuts, RadiusValues,Rejection_radius;

								vector<double> eff,rej;
								std::vector<double> rejection_Error;

								for(double i=0.80;i<1;i+=0.01){
												// cout<<"Finding Cut for eff = "<<i<<endl;
												double TotRej=findCutsEnergy(compareBi210,comparePo210,i,rejection_Error,mistaggedHist);
												eff.push_back(i*100);
												rej.push_back(TotRej);
								}

								vector<double> energyValues_rev,cuts_rev;
								vector<double> eff_rev,rej_rev;
								vector<TGraph*> graphs;


								//------Drawing and saving---------------


								std::vector<double> percent_Error,eff_Error,percent;

								for (int i = 0; i < eff.size(); i++)
								{
												percent.push_back(100/rej[i]);	
												percent_Error.push_back(percent[i]*rejection_Error[i]/rej[i]);	
												eff_Error.push_back(0);
												std::cout << "percent = "<<percent[i]<<"+/-"<<percent_Error[i] << std::endl;
												std::cout << "rej = "<<rej[i]<<"+/-"<<rejection_Error[i] << std::endl;
								}
								std::cout << "FLTMIN = "<<FLT_MIN << std::endl;

								TCanvas* c_mistagged =new TCanvas();
								c_mistagged->SetGrid();
								mistaggedHist->SetTitle("Number of mistagged events");
								mistaggedHist->GetXaxis()->SetTitle("Mistagged events");
								mistaggedHist->SetMaximum(100);
								mistaggedHist->Draw();

								c_mistagged->Print(Form("plots/mistagged_total_highRad_%.f.png",highRad));
								c_mistagged->Print(Form("plots/mistagged_total_highRad_%.f.tex",highRad));

								TCanvas* c_new =new TCanvas();
								TGraphErrors* fractionGraph = new TGraphErrors(percent.size(),&eff[0],&percent[0],&eff_Error[0],&percent_Error[0]);
								// TF1 *f_percent_alpha = new TF1("f", "[3]*x*x*x +[2]*x*x +[1]*x +[0]");
								// fractionGraph ->Fit(f_percent_alpha);
								fractionGraph ->SetTitle(Form("Surviving #alpha 's over accepted #beta 's {mcPosr < %f }",highRad));
								fractionGraph ->GetXaxis()->SetTitle("Remaining #beta 's ");
								fractionGraph ->GetYaxis()->SetTitle("Remaining #alpha 's ");
								// fractionGraph ->GetXaxis()->SetTitle("Remaining #beta 's %");
								// fractionGraph ->GetYaxis()->SetTitle("Remaining #alpha 's %");
								fractionGraph ->Draw("ap");
								c_new->SetGrid();

								c_new->Print(Form("plots/sig_sack_alpha_highRad_%.f.png",highRad));
								c_new->Print(Form("plots/sig_sack_alpha_highRad_%.f.tex",highRad));

								c_new->SetLogy();
								fractionGraph ->GetYaxis()->SetMoreLogLabels();
								fractionGraph ->Draw("ap");

								c_new->Print(Form("plots/sig_sack_alpha_log_highRad_%.f.png",highRad));
								c_new->Print(Form("plots/sig_sack_alpha_log_highRad_%.f.tex",highRad));

								TCanvas* c_complete = new TCanvas();
								fractionGraph->SetMarkerColor(highRad/1000);
								fractionGraph->SetLineColor(highRad/1000);
								fractionGraph->SetFillColor(highRad/1000);
								TLegend* leg = new TLegend(0.11,0.11,0.22,0.22);
								leg->AddEntry(fractionGraph,SSTR(highRad).c_str(),"f");
								// if(highRad==1000) fractionGraph ->Draw("ap");
								fractionGraph ->Draw("ap same");
								leg->Draw();
								c_new->SetGrid();
								c_complete->Print("plots/All_sig_sack_alpha.png");
								c_complete->Print("plots/All_sig_sack_alpha.tex");

								c_complete->SetLogy();
								fractionGraph ->Draw("ap same");
								c_complete->Print("plots/All_sig_sack_alpha_log.png");
								c_complete->Print("plots/All_sig_sack_alpha_log.tex");

								TFile fileout("RejectionsAndEfficiency.root","UPDATE");
								fileout.cd();
								// 	-----Energy-----
								fractionGraph->Write("ap");
								// RejectionGraph_rev->Write();
								// 	-----Radial-----
								/* compareEle_radial->Write(); */
								/* compareAlpha_radial->Write(); */
								/* RejectionRadialGraph->Write(); */
								/* RejectionRadial_can->Write(); */
								fileout.Close();
								std::cout << "beta number of entries = "<<compareBi210->GetEntries() << std::endl;
								std::cout << "alpha number of entries = "<<comparePo210->GetEntries() << std::endl;

								return 0;

				}

				void scanRad(){

								for (double i = 1000; i < 6001; i+=1000) {
												scan(i);	
								}

				}
