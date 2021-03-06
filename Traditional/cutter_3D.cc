#include <RAT/DU/DSReader.hh>
#include <RAT/DS/Entry.hh>

#include <TGraph2D.h>
#include <TMath.h>
#include <string>
#include <TMultiGraph.h>
#include <TH2D.h>
#include <TH3D.h>
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
#include <fstream>

#define SSTR( x ) static_cast< std::ostringstream & >( \
								( std::ostringstream() << std::dec << x ) ).str()


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


void FillHist(TFile* file,TH3D* hist){

				TTree* Tree = (TTree*) file->Get("output");
				Double_t berkeleyAlphaBeta, mcEdepQuenched,mcPosr;
				Bool_t  Qfit;
				Int_t pdg1, pdg2;
				Int_t parentpdg1,parentpdg2;

				Tree->SetBranchAddress("pdg1",&pdg1);
				Tree->SetBranchAddress("pdg2",&pdg2);
				Tree->SetBranchAddress("parentpdg1",&parentpdg1);
				Tree->SetBranchAddress("parentpdg2",&parentpdg2);

				Tree->SetBranchAddress("berkeleyAlphaBeta",&berkeleyAlphaBeta);
				Tree->SetBranchAddress("mcEdepQuenched",&mcEdepQuenched);
				Tree->SetBranchAddress("mcPosr",&mcPosr);
				// Tree->SetBranchAddress("posr",&posr);
				Tree->SetBranchAddress("fitValid",&Qfit);
				Int_t n = (Int_t)Tree->GetEntries();
				Int_t code;

				for( Int_t i =0;i<n;i++){
								Tree->GetEntry(i);
								if (partflag=="Bi"){
												code=11;
								}else if( partflag=="Po"){
												code=1000020040;
								}
								if( Qfit && pdg1==code  ){
												// x-axis- MC energy
												// y-axis- posr
												// z-axis- berkeleyAlphaBeta
												// hist->Fill(mcEdepQuenched,posr,berkeleyAlphaBeta);
												hist->Fill(mcEdepQuenched,mcPosr,berkeleyAlphaBeta);
								}
				}
				file->Close();
}


/* void findCutsEnergy(TH3D *ele_hist,TH3D *alpha_hist,vector<double>& energyValues, vector<double>& cutValues,vector<double>& Rejection_values){ */
/* 	* This function should take a 2d histrogram split it into energy */ 
/* 	 * strips and then find a cut value which retain ~99% of the signal. */
/* 	 / */
/* 	//for(double energy =0; energy<3.5;energy+=0.1){ */

/* 	double startingCut=-1000; */
/* 	double step=0.1; */
/* 	double cutValue=startingCut; */
/* 	double accept=0; */
/* 	double mistagged=0; */
/* 	double rejection=0; */
/* 	double sliceWidth=0.1; */
/* 	double threshold=0.995; */
/* 	TAxis* bi_x=compareBi210->GetXaxis(); */
/* 	TAxis* bi_y=compareBi210->GetYaxis(); */
/* 	TAxis* po_x=comparePo210->GetXaxis(); */
/* 	TAxis* po_y=comparePo210->GetYaxis(); */

/* 	TH1D* complete_bi=compareBi210->ProjectionY("complete_bi",bi_x->FindBin(0.),bi_x->FindBin(4.)); */
/* 	TH1D* complete_po=comparePo210->ProjectionY("complete_po",po_x->FindBin(0.),po_x->FindBin(4.)); */
/* 	double bi_numEV_complete= compareBi210->GetEntries(); */
/* 	double po_numEV_complete= comparePo210->GetEntries(); */
/* 	cout<<"Number of entries in complete bi = "<< bi_numEV_complete<<endl; */
/* 	cout<<"Number of entries in complete po = "<< po_numEV_complete<<endl; */

/* 	for(double energy =0; energy<2.5;energy+=0.1){ */
/* 		cout<<"Energy = "<< energy<<endl; */

/* 		//These histograms contain the sliced projections */
/* 		TH1D* slice_bi=compareBi210->ProjectionY(SSTR(energy).c_str(),bi_x->FindBin(energy),bi_x->FindBin(energy+sliceWidth)); */
/* 		TH1D* slice_po=comparePo210->ProjectionY(("Po"+SSTR(energy)).c_str(),po_x->FindBin(energy),po_x->FindBin(energy+sliceWidth)); */

/* 		//These two numbers are the num of events before cuts */
/* 		double bi_numEV=slice_bi->GetEntries(); */
/* 		double po_numEV=slice_po->GetEntries(); */

/* 		TAxis* bi_slice_x=slice_bi->GetXaxis(); */
/* 		TAxis* po_slice_x=slice_po->GetXaxis(); */

/* 		while(accept<threshold){ */

/* 			//This loop finds the cut value to a certain threshold is achevied, it does so by incrementing by 'step' */
/* 			double bi_remain=slice_bi->Integral(bi_slice_x->FindBin(startingCut),bi_slice_x->FindBin(cutValue+=step)); */
/* 			//double po_remain=slice_po->Integral(po_slice_x->FindBin(startingCut),po_slice_x->FindBin(startingCut+step)); */
/* 			cout<<"Number of entries remaining in bi slice = "<< bi_remain<<endl; */
/* 			//cout<<"Number of entries remaining in po slice = "<< po_remain<<endl; */
/* 			accept=bi_remain/bi_numEV; */	
/* 			cout<<"Accepted fraction = "<<accept<<endl; */
/* 		} */
/* 		mistagged=slice_po->Integral(po_slice_x->FindBin(startingCut),po_slice_x->FindBin(cutValue)); */

/* 		rejection= po_numEV/mistagged; */
/* 		cout<<mistagged<<" events were mistagged. Rejection= "<<rejection<<"."<<endl; */
/* 		Rejection_values.push_back(rejection); */
/* 		cutValues.push_back(cutValue); */
/* 		energyValues.push_back(energy+0.05); */

/* 		if(false){ */
/* 		TCanvas* c1 = new TCanvas(); */
/* 		c1->cd(); */
/* 		//TLine* cutLine = new TLine( cutValue, slice_bi->GetYaxis()->GetXmin(), cutValue, slice_bi->GetYaxis()->GetXmax() ); */
/* 		TLine* cutLine = new TLine( cutValue, slice_bi->GetMinimum(), cutValue, slice_bi->GetMaximum() ); */
/* 		cutLine->SetLineWidth(2);cutLine->SetLineStyle(4); */
/* 		slice_bi->Draw(); */
/* 		slice_po->Draw("same"); */
/* 		cutLine->Draw("same"); */
/* 		} */
/* 		cutValue=startingCut; */
/* 		accept=0; */
/* 	} */

/* } */

void findCuts(TH3D* ele_hist,TH3D* alpha_hist,vector<double>& energyValues,vector<double>& rValues, vector<double>& cutValues,vector<double>& Rejection_values){

				double threshold=0.99;

				double startingCut= -1000;
				double step=0.05;
				double cutValue=startingCut;

				double accept=0;
				double mistagged=0;//number of alphas in the beta accpeted window.
				double rejection=0;

				double minBin_x=ele_hist->GetXaxis()->GetXmin();
				double maxBin_x=ele_hist->GetXaxis()->GetXmax();
				double sliceWidth_x=ele_hist->GetXaxis()->GetBinWidth(1);

				double minBin_y=alpha_hist->GetYaxis()->GetXmin();
				double maxBin_y=alpha_hist->GetYaxis()->GetXmax();
				double sliceWidth_y=alpha_hist->GetYaxis()->GetBinWidth(1);

				TAxis* ele_x=ele_hist->GetXaxis();
				TAxis* ele_y=ele_hist->GetYaxis();
				TAxis* alpha_x=alpha_hist->GetXaxis();
				TAxis* alpha_y=alpha_hist->GetYaxis();

				double ele_numEV_complete= ele_hist->GetEntries();
				double alpha_numEV_complete= alpha_hist->GetEntries();

				double TotalMistagged=0;
				double TotalNumOfPo=0;

				for(double E =minBin_x; E<maxBin_x; E+=sliceWidth_x){
								for (double r = 0; r < maxBin_y; r+=sliceWidth_y) {

												//These histograms contain the sliced projections
												TH1D* slice_ele=ele_hist->ProjectionZ("",ele_x->FindBin(E),ele_x->FindBin(E+sliceWidth_x),ele_y->FindBin(r),ele_y->FindBin(r+sliceWidth_y));
												TH1D* slice_alpha=alpha_hist->ProjectionZ("holder",alpha_x->FindBin(E),alpha_x->FindBin(E+sliceWidth_x),alpha_y->FindBin(r),alpha_y->FindBin(r+sliceWidth_y));

												//These two numbers are the num of events before cuts
												double ele_numEV=slice_ele->GetEntries();
												double alpha_numEV=slice_alpha->GetEntries();

												TAxis* ele_slice_x=slice_ele->GetXaxis();
												TAxis* alpha_slice_x=slice_alpha->GetXaxis();
												if(ele_numEV>100 && alpha_numEV>100){
																while(accept<threshold){

																				//This loop finds the cut value to a certain threshold is achevied, it does so by incrementing by 'step'
																				double ele_remain=slice_ele->Integral(ele_slice_x->FindBin(startingCut),ele_slice_x->FindBin(cutValue+=step));

																				//double alpha_remain=slice_alpha->Integral(alpha_slice_x->FindBin(startingCut),alpha_slice_x->FindBin(startingCut+step));
																				cout<<"Number of entries remaining in ele slice = "<< ele_remain<<endl;
																				// cout<<"Number of entries remaining in alpha slice = "<< alpha_remain<<endl;
																				accept=ele_remain/ele_numEV; 
																				cout<<"Accepted fraction = "<<accept<<endl;
																}
																mistagged=slice_alpha->Integral(alpha_slice_x->FindBin(startingCut),alpha_slice_x->FindBin(cutValue));

																rejection= alpha_numEV/mistagged;
																// cout<<mistagged<<" events were mistagged. Rejection= "<<rejection<<"."<<endl;
																Rejection_values.push_back(rejection);
																cutValues.push_back(cutValue);
																energyValues.push_back(E+sliceWidth_x/2);
																rValues.push_back(r+sliceWidth_y/2);
																// std::vector<double> ent;
																// ent.push_back(E);
																// ent.push_back(r);
																// ent.push_back(cutValues);
																// data.push_back(ent);
												}else{
																Rejection_values.push_back(0);
																cutValues.push_back(0);
																energyValues.push_back(E+sliceWidth_x/2);
																rValues.push_back(r+sliceWidth_y/2);
																// std::vector<double> ent;
																// ent.push_back(E);
																// ent.push_back(r);
																// ent.push_back(cutValues);
																// data.push_back(ent);

												}

												cutValue=startingCut;
												accept=0;
								}//End of r loop.

				}//End of E loop.
}

int main(){
				//====================================================================	
				gStyle->SetOptStat(0);

				// TH3D* ele_hist= new TH3D("ele_hist","berkeleyAlphaBeta",25,0,2.5,60,0,6000,80,-400,400);
				TH3D* ele_hist= new TH3D("ele_hist","berkeleyAlphaBeta",5,0,2.5,6,0,6000,8,-400,400);
				ele_hist->SetLineColor(kBlue);ele_hist->SetLineWidth(3);
				ele_hist->SetMarkerColor(kBlue);
				ele_hist->SetFillColor(kBlue);
				// TH3D* alpha_hist= new TH3D("alpha_hist","berkeleyAlphaBeta",25,0,2.5,60,0,6000,40,-400,400);
				TH3D* alpha_hist= new TH3D("alpha_hist","berkeleyAlphaBeta",5,0,2.5,6,0,6000,8,-400,400);
				alpha_hist->SetLineColor(kRed);alpha_hist->SetLineWidth(3);
				alpha_hist->SetMarkerColor(kRed);
				alpha_hist->SetFillColor(kRed);


				TLegend* t1 = new TLegend( 0.11, 0.11, 0.21, 0.21 );

				vector<string>  biFileList= glob("/data/snoplus/liggins/year1/fitting/fitting/alphaSims/output_electron/ntuple","electron");
				vector<string> poFileList= glob("/data/snoplus/liggins/year1/fitting/fitting/alphaSims/output/ntuple","alpha");

				for( int i=0; i<biFileList.size(); i++ ){
								TFile * file= TFile::Open(biFileList[i].c_str());	
								partflag="Bi";
								FillHist( file, ele_hist);
				}

				for( int i=0; i<poFileList.size(); i++ ){
								TFile * file= TFile::Open(poFileList[i].c_str());	
								partflag="Po";
								FillHist( file, alpha_hist);
				}

				vector<double> cuts,  energyValues, radiusValues,Rejection_energy;
				vector<double> radial_cuts, RadiusValues,Rejection_radius;

				vector< vector<double> > data;

				TGraph2D * bound = new TGraph2D();
				findCuts(ele_hist,alpha_hist,energyValues,radiusValues,cuts,Rejection_energy);
				std::ofstream output_file("./example.txt");
				output_file<< "E,r,cut" << std::endl;					

				double minBin_x=ele_hist->GetXaxis()->GetXmin();
				double maxBin_x=ele_hist->GetXaxis()->GetXmax();
				double sliceWidth_x=ele_hist->GetXaxis()->GetBinWidth(1);

				double minBin_y=alpha_hist->GetYaxis()->GetXmin();
				double maxBin_y=alpha_hist->GetYaxis()->GetXmax();
				double sliceWidth_y=alpha_hist->GetYaxis()->GetBinWidth(1);

				double numofbin_x=(maxBin_x-minBin_x)/sliceWidth_x;
				double numofbin_y=(maxBin_y-minBin_y)/sliceWidth_y;

				double start=0;
				double point=0;
				for (int i = 0; i < numofbin_x; i++) {
								// for (int j = start*numofbin_y; j < start*numofbin_y+numofbin_y; j++) {
								for (int j = start*numofbin_y; j < start*numofbin_y+numofbin_y; j++) {
									
								output_file<<energyValues[j]<<","<< radiusValues[j] <<","<<cuts[j]<<std::endl;					
								bound->SetPoint(point,energyValues[j],radiusValues[j],cuts[j]);
								point++;
								}
								start++;
				}

				output_file.close();

				/* findCutsRadial(compareEle_radial,compareAlpha_radial,RadiusValues,radial_cuts,Rejection_radius); */
				/* // ----Finding Rejection Powers--------- */


				/* //------Drawing and saving--------------- */


				/* TGraph* cutGraph = new TGraph(bi_cuts.size(),&energyValues[0],&bi_cuts[0]); */
				/* TF1 *f = new TF1("f", "[3]*x*x*x +[2]*x*x +[1]*x +[0]"); */
				/* cutGraph->Fit(f); */

				/* TGraph* cutGraph_radial = new TGraph(radial_cuts.size(),&RadiusValues[0],&radial_cuts[0]); */
				/* TF1 *f_Rad = new TF1("f_Rad", "[1]*x +[0]"); */
				/* cutGraph_radial->Fit(f_Rad); */



				TCanvas * c1 = new TCanvas();
				bound->Draw("ap");
				
				return 0;

				TCanvas * c2 = new TCanvas();
				cout<<"# of entries in ele_hist= "<< ele_hist->GetEntries()<<endl;
				cout<<"# of entries in alpha_hist = "<< alpha_hist->GetEntries()<<endl;
				ele_hist->SetTitle("Classifiers values across truth quenched energy and reconstruted radius");
				ele_hist->GetXaxis()->SetTitle("mcEdepQuenched (MeV)");
				ele_hist->GetYaxis()->SetTitle("posr(mm)");
				ele_hist->GetZaxis()->SetTitle("BerkeleyAlphaBeta");
				ele_hist->SetFillColor(kRed);
				ele_hist->SetMarkerColor(kRed);
				ele_hist->SetLineColor(kRed);
				ele_hist->Draw();
				alpha_hist->SetFillColor(kBlue);
				alpha_hist->SetMarkerColor(kBlue);
				alpha_hist->SetLineColor(kBlue);
				alpha_hist->Draw("same");

				t1->AddEntry(ele_hist, "beta","f");
				t1->AddEntry(alpha_hist, "alpha","f");
				t1->Draw();

				c2->Print("compareBerkeleyAlphaBetaVsMCEnergy_withCuts.png");
				TFile fileout("3D_Plot.root","RECREATE");
				fileout.cd();
				alpha_hist->Write();
				ele_hist->Write();
				fileout.Close();

				return 0;

}

