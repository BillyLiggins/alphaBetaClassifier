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


void findCuts(TH2D *compareBi210,TH2D *comparePo210,vector<double>& energyValues, vector<double>& cutValues,vector<double>& Rejection_values,double threshold);
string partflag="";

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


void FillHist(TFile* file,TH1D * hist,TH2D * hist2,TH1D *radial,TH2D *compareRadial,ofstream& outputfile){

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
				Int_t code;

				for( Int_t i =0;i<n;i++){
								Tree->GetEntry(i);
								if (partflag=="Bi"){
												code=11;
								}else if( partflag=="Po"){
												code=1000020040;
								}
								if( Qfit && evIndex==0 && mcPosr<4000  ){
												/* if( Qfit){ */
												hist->Fill(para);
												hist2->Fill(mcEdepQuenched,para);
												radial->Fill(posr);
												compareRadial->Fill(posr,para);
												outputfile<< mcEdepQuenched<<","<<posr<<","<< para << endl;
								}
								}

								file->Close();
				}


				void findCutsEnergy(TH2D *compareBi210,TH2D *comparePo210,vector<double>& energyValues, vector<double>& cutValues,vector<double>& Rejection_values){
								/* This function should take a 2d histrogram split it into energy 
								 * strips and then find a cut value which retain ~99% of the signal.
								 */
								//for(double energy =0; energy<3.5;energy+=0.1){

								double startingCut=-500;
								double step=0.1;
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

								for(double energy =minBin; energy<maxBin;energy+=sliceWidth){
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
																double bi_remain=slice_bi->Integral(bi_slice_x->FindBin(startingCut),bi_slice_x->FindBin(cutValue+=step));
																//double po_remain=slice_po->Integral(po_slice_x->FindBin(startingCut),po_slice_x->FindBin(startingCut+step));
																cout<<"Number of entries remaining in bi slice = "<< bi_remain<<endl;
																//cout<<"Number of entries remaining in po slice = "<< po_remain<<endl;
																accept=bi_remain/bi_numEV;	
																cout<<"Accepted fraction = "<<accept<<endl;
												}
												mistagged=slice_po->Integral(po_slice_x->FindBin(startingCut),po_slice_x->FindBin(cutValue));

												rejection= po_numEV/(mistagged+0.00001);
												cout<<mistagged<<" events were mistagged. Rejection= "<<rejection<<"."<<endl;
												Rejection_values.push_back(rejection);
												cutValues.push_back(cutValue);
												energyValues.push_back(energy+ sliceWidth/2);

												if(false){
																TCanvas* c1 = new TCanvas();
																c1->cd();
																//TLine* cutLine = new TLine( cutValue, slice_bi->GetYaxis()->GetXmin(), cutValue, slice_bi->GetYaxis()->GetXmax() );
																TLine* cutLine = new TLine( cutValue, slice_bi->GetMinimum(), cutValue, slice_bi->GetMaximum() );
																cutLine->SetLineWidth(2);cutLine->SetLineStyle(4);
																slice_bi->Draw();
																slice_po->Draw("same");
																cutLine->Draw("same");
												}
												cutValue=startingCut;
												accept=0;
								}

								}


								int cutter(){
												//====================================================================	
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
												compareBi210->SetLineColor(kBlue);compareBi210->SetLineWidth(3);
												compareBi210->SetMarkerColorAlpha(kBlue,0.2);
												compareBi210->SetFillColor(kBlue);
												TH2D* comparePo210   = new TH2D("comparePo210","berkeleyAlphaBeta",binNum,0,2.5,260,-160,100);
												comparePo210->SetLineColor(kRed);comparePo210->SetLineWidth(3);
												comparePo210->SetMarkerColorAlpha(kRed,0.2);
												comparePo210->SetFillColor(kRed);

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

												findCutsEnergy(compareBi210,comparePo210,energyValues,bi_cuts,Rejection_energy);
												double thres=0.9;
												int counter=1;
												vector< vector<double> > listOfCuts;

												Bool_t QcutScan=false;
																vector<TF1*> listOfFunction;
												if (QcutScan) {


																while(thres<0.999999){
																				vector<double> cuts ;

																				findCuts(compareBi210,comparePo210,energyValues,cuts,Rejection_energy,thres);
																				listOfCuts.push_back(cuts);
																				thres+=0.9/pow(10,counter);
																				counter++;
																}
																vector<TGraph*> listOfGraphs;
																for(int i=0;i<listOfCuts.size();i++){
																				vector<double> loader=listOfCuts[i];
																				listOfGraphs.push_back(	new TGraph(loader.size(),&energyValues[0],&loader[0]));
																				TF1 *f_E = new TF1("f_E", "[1]*x +[0]",0,2.5);
																				// f_E->SetLineColor(kBlack);
																				f_E->SetLineColor(i*2+40);
																				f_E->SetLineStyle(i+1);
																				f_E->SetLineWidth(3);
																				listOfGraphs[i]->Fit(f_E);
																				listOfFunction.push_back(f_E);

																}



												}//end of QcutScan


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
												c2->Print("plots/BerkeleyAlphaBetaVsMCEnergy_withCuts.png");
												c2->Print("plots/BerkeleyAlphaBetaVsMCEnergy_withCuts.tex");

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

																c_test->Print("functions.png");

												}//end of QcutScan


												TCanvas* RejectionEnergy_can = new TCanvas();
												TGraph* RejectionEnergyGraph = new TGraph(Rejection_energy.size(),&energyValues[0],&Rejection_energy[0]);
												//TF1 *f = new TF1("f", "[1]*x +[0]");
												//cutGraph->Fit(f);
												//
												for (int i = 0; i < Rejection_energy.size(); i++) {
																std::cout << "ith rejection value = "<<Rejection_energy[i] << std::endl;
												}

												RejectionEnergyGraph->SetTitle("Rejection across energy with 99% #beta retention");
												RejectionEnergyGraph->GetXaxis()->SetTitle("Energy (MeV)");
												RejectionEnergyGraph->GetYaxis()->SetTitle("Rejection");
												// RejectionEnergyGraph->SetMaximum(100000);
												RejectionEnergyGraph->Draw("apl");
												// RejectionEnergyGraph->Draw("a");
												RejectionEnergy_can->SetLogy();
												RejectionEnergy_can->Print("plots/RejectionAcrossEnergy.png");
												RejectionEnergy_can->Print("plots/RejectionAcrossEnergy.tex");

												TCanvas * c3 = new TCanvas();
												cout<<"# of entries in hBi210 = "<< hBi210->GetEntries()<<endl;
												cout<<"# of entries in hPo210 = "<< hPo210->GetEntries()<<endl;
												compareEle_radial->SetTitle("Classifier values across reconstructed radius");
												compareEle_radial->GetXaxis()->SetTitle("posr (mm)");
												compareEle_radial->GetYaxis()->SetTitle("BerkeleyAlphaBeta");
												compareEle_radial->Draw();
												compareAlpha_radial->Draw("same");
												//cutGraph->Draw("AP*1");
												cutGraph_radial->Draw("* same");

												t3->AddEntry( compareBi210, "#beta 's","f");
												t3->AddEntry( comparePo210, "#alpha 's","f");
												t3->Draw();
												c3->Print("plots/compareBerkeleyAlphaBetaVsPosr_withCuts.png");
												c3->Print("plots/compareBerkeleyAlphaBetaVsPosr_withCuts.tex");

												// TCanvas* RejectionRadial_can= new TCanvas();
												// TGraph* RejectionRadialGraph = new TGraph(radial_cuts.size(),&RadiusValues[0],&Rejection_radius[0]);
												// double TotalRejectionInPosr=std::accumulate(Rejection_radius.begin(), Rejection_radius.end(), 0);
												// cout<<"Total Rejection = "<<TotalRejectionInPosr<<endl;
												//
												// TPaveText *pt = new TPaveText(500,500,2500,900,"NB");
												// pt->AddText(("Total Rejection = "+SSTR(TotalRejectionInPosr)+".").c_str());
												// RejectionRadialGraph->SetTitle("Rejection across posr");
												// RejectionRadialGraph->GetXaxis()->SetTitle("posr (mm)");
												// RejectionRadialGraph->GetYaxis()->SetTitle("Rejection");
												// RejectionRadialGraph->Draw();
												// pt->Draw();
												// RejectionRadialGraph->SetMaximum(1000);
												// RejectionRadial_can->Print("RejectionAcrossRadius.png");


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

												return 0;

								}

								void findCuts(TH2D *compareBi210,TH2D *comparePo210,vector<double>& energyValues, vector<double>& cutValues,vector<double>& Rejection_values,double threshold){
												/* This function should take a 2d histrogram split it into energy 
												 * strips and then find a cut value which retain ~99% of the signal.
												 */
												//for(double energy =0; energy<3.5;energy+=0.1){

												double startingCut=-500;
												double step=0.1;
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
												// double threshold=0.99;
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

												for(double energy =minBin; energy<maxBin;energy+=sliceWidth){
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
																				double bi_remain=slice_bi->Integral(bi_slice_x->FindBin(startingCut),bi_slice_x->FindBin(cutValue+=step));
																				//double po_remain=slice_po->Integral(po_slice_x->FindBin(startingCut),po_slice_x->FindBin(startingCut+step));
																				cout<<"Number of entries remaining in bi slice = "<< bi_remain<<endl;
																				//cout<<"Number of entries remaining in po slice = "<< po_remain<<endl;
																				accept=bi_remain/bi_numEV;	
																				cout<<"Accepted fraction = "<<accept<<endl;
																}
																mistagged=slice_po->Integral(po_slice_x->FindBin(startingCut),po_slice_x->FindBin(cutValue));

																rejection= po_numEV/(mistagged+0.00001);
																cout<<mistagged<<" events were mistagged. Rejection= "<<rejection<<"."<<endl;
																Rejection_values.push_back(rejection);
																cutValues.push_back(cutValue);
																energyValues.push_back(energy+ sliceWidth/2);

																if(false){
																				TCanvas* c1 = new TCanvas();
																				c1->cd();
																				//TLine* cutLine = new TLine( cutValue, slice_bi->GetYaxis()->GetXmin(), cutValue, slice_bi->GetYaxis()->GetXmax() );
																				TLine* cutLine = new TLine( cutValue, slice_bi->GetMinimum(), cutValue, slice_bi->GetMaximum() );
																				cutLine->SetLineWidth(2);cutLine->SetLineStyle(4);
																				slice_bi->Draw();
																				slice_po->Draw("same");
																				cutLine->Draw("same");
																}
																cutValue=startingCut;
																accept=0;
												}

												}

