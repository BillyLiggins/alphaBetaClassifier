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


void FillHist(TFile* file,TH2D * hist2){

	TTree* Tree = (TTree*) file->Get("output");
	Double_t para, mcEdepQuenched,posr;
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
		if( Qfit && evIndex==0 && posr<5900  ){
			hist2->Fill(mcEdepQuenched,para);
		}
	}

}

double findCutsEnergy(TH2D *compareBi210,TH2D *comparePo210, double threshold,std::vector<double>& rej_vect){
	/* This function returns the total rejection for a threshold.
	 */

	vector<double> Rejection_values;
	double startingCut=-1000;
	double step=0.05;
	double cutValue=startingCut;
	double accept=0;
	double mistagged=0;
	double rejection=0;
	double minBin=compareBi210->GetXaxis()->GetXmin();
	double maxBin=compareBi210->GetXaxis()->GetXmax();
	double sliceWidth=compareBi210->GetXaxis()->GetBinWidth(1);
	/* cout<<"minBin = "<<minBin<<endl; */
	/* cout<<"maxBin = "<<maxBin<<endl; */
	/* cout<<"sliceWidth = "<<sliceWidth<<endl; */
	//double threshold=0.995;
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

	double TotalMistagged=0;
	double TotalNumOfPo=0;
	for(double energy =minBin; energy<maxBin;energy+=sliceWidth){
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
			accept=bi_remain/bi_numEV;	
		}
		mistagged=slice_po->Integral(po_slice_x->FindBin(startingCut),po_slice_x->FindBin(cutValue));
		TotalMistagged=mistagged+TotalMistagged;
		TotalNumOfPo=po_numEV+TotalNumOfPo;

		if(po_numEV>1 && mistagged>1){
		rejection= po_numEV/(mistagged);
		}else{
		rejection= 0;
		}
		Rejection_values.push_back(rejection);

		if(false){
			TCanvas* c1 = new TCanvas();
			c1->cd();
			//TLine* cutLine = new TLine( cutValue, slice_bi->GetYaxis()->GetXmin(), cutValue, slice_bi->GetYaxis()->GetXmax() );
			TLine* cutLine = new TLine( cutValue, slice_bi->GetMinimum(), cutValue, slice_bi->GetMaximum() );
			cutLine->SetLineWidth(2);
			cutLine->SetLineStyle(4);
			slice_bi->SetTitle(("Energy slice "+SSTR(energy)+" < E < "+SSTR(energy+sliceWidth)+" { #beta efficiency = "+SSTR(threshold)+"}").c_str());
			slice_bi->Draw();
			slice_po->Draw("same");
			cutLine->Draw("same");
			c1->Print(("plots/Energy_slice_betaEfficiency_"+SSTR(threshold)+"_"+SSTR(energy)+"_E_"+SSTR(energy+sliceWidth)+".png").c_str());
		}
		cutValue=startingCut;
		accept=0;
	}

	rej_vect.push_back(rejection*sqrt((TotalNumOfPo+TotalMistagged)/(TotalNumOfPo*TotalMistagged)));
	return std::accumulate(Rejection_values.begin(), Rejection_values.end(), 0);
	Rejection_values.clear();
}

double findCutsEnergy_reverse(TH2D *compareBi210,TH2D *comparePo210, double threshold,vector<double> &energyValues, vector<double> &cutsValues){
	/* This function returns the total rejection for a threshold.
	 */

	vector<double> Rejection_values;
	double startingCut=1000;
	double step=0.01;
	double cutValue=startingCut;
	double accept=0;
	double mistagged=0;
	double rejection=0;
	double minBin=comparePo210->GetXaxis()->GetXmin();
	double maxBin=comparePo210->GetXaxis()->GetXmax();
	double sliceWidth=comparePo210->GetXaxis()->GetBinWidth(1);
	cout<<"minBin = "<<minBin<<endl;
	cout<<"maxBin = "<<maxBin<<endl;
	cout<<"sliceWidth = "<<sliceWidth<<endl;
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
		cout<<"Energy = "<< energy<<endl;

		energyValues.push_back(energy);
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
			accept=po_remain/po_numEV;	
		}
		cout<<"cutValue = "<<cutValue<<endl;
		cutsValues.push_back(cutValue);
		mistagged=slice_bi->Integral(bi_slice_x->FindBin(cutValue),bi_slice_x->FindBin(startingCut));

		if(bi_numEV>1 && mistagged>1){
		rejection= bi_numEV/(mistagged);
		}else{
		rejection= 0;
		}
		Rejection_values.push_back(rejection);

		if(false){
			TCanvas* c1 = new TCanvas();
			c1->cd();
			//TLine* cutLine = new TLine( cutValue, slice_bi->GetYaxis()->GetXmin(), cutValue, slice_bi->GetYaxis()->GetXmax() );
			TLine* cutLine = new TLine( cutValue, slice_po->GetMinimum(), cutValue, slice_po->GetMaximum() );
			cutLine->SetLineWidth(2);cutLine->SetLineStyle(4);
			slice_bi->SetTitle(("Rev energy "+SSTR(energy)+" ").c_str());
			slice_bi->Draw();
			slice_po->Draw("same");
			cutLine->Draw("same");
		}
		cutValue=startingCut;
		accept=0;
	}


	return std::accumulate(Rejection_values.begin(), Rejection_values.end(), 0);
}


int scan(){
//====================================================================	
	gStyle->SetOptStat(0);
	TH2D* compareBi210   = new TH2D("compareBi210","berkeleyAlphaBeta",26,0,2.6,260,-160,100);
	compareBi210->SetLineColor(kBlue);compareBi210->SetLineWidth(3);
	compareBi210->SetMarkerColor(kBlue);
	compareBi210->SetFillColor(kBlue);
	TH2D* comparePo210   = new TH2D("comparePo210","berkeleyAlphaBeta",26,0,2.6,260,-160,100);
	comparePo210->SetLineColor(kRed);comparePo210->SetLineWidth(3);
	comparePo210->SetMarkerColor(kRed);
	comparePo210->SetFillColor(kRed);

//==========================Radial Histrograms ========================

	
	TLegend* t1 = new TLegend( 0.6, 0.7, 0.89, 0.88 );
	TLegend* t2 = new TLegend( 0.11, 0.11, 0.21, 0.21 );
	TLegend* t3 = new TLegend( 0.11, 0.11, 0.21, 0.21 );
	
	vector<string> biFileList= glob("/data/snoplus/liggins/year1/fitting/fitting/alphaSims/output_electron/ntuple","electron");
	vector<string> poFileList= glob("/data/snoplus/liggins/year1/fitting/fitting/alphaSims/output/ntuple","alpha");

	for( int i=0; i<biFileList.size(); i++ ){
		TFile * file= TFile::Open(biFileList[i].c_str());	
		partflag="Bi";
		FillHist( file,compareBi210);
		file->Close();
	  }

	for( int i=0; i<poFileList.size(); i++ ){
		TFile * file= TFile::Open(poFileList[i].c_str());	
		partflag="Po";
		FillHist( file,comparePo210);
		file->Close();
	  }

	vector<double> bi_cuts, energyValues,Rejection_energy;
	vector<double> radial_cuts, RadiusValues,Rejection_radius;

	vector<double> eff,rej;
	std::vector<double> rejection_Error;

	for(double i=0.90;i<1;i+=0.01){
		cout<<"Finding Cut for eff = "<<i<<endl;
		double TotRej=findCutsEnergy(compareBi210,comparePo210,i,rejection_Error);
		eff.push_back(i);
		rej.push_back(TotRej);
		}
		
	vector<double> energyValues_rev,cuts_rev;
	vector<double> eff_rev,rej_rev;
	vector<TGraph*> graphs;

	for(double i=0.90;i<1;i+=0.01){
		cout<<"Finding Cut for eff_rev = "<<i<<endl;
		double TotRej=findCutsEnergy_reverse(compareBi210,comparePo210,i,energyValues_rev,cuts_rev);
		graphs.push_back(new TGraph(cuts_rev.size(),&energyValues_rev[0],&cuts_rev[0]));

		eff_rev.push_back(i);
		rej_rev.push_back(TotRej);
	}

//------Drawing and saving---------------

	TCanvas * c1= new TCanvas();
	TGraph* RejectionGraph = new TGraph(eff.size(),&eff[0],&rej[0]);
	TF1 *f = new TF1("f", "[3]*x*x*x +[2]*x*x +[1]*x +[0]");
	RejectionGraph ->Fit(f);
	RejectionGraph ->SetTitle("Rejection Across Efficiency");
	RejectionGraph ->GetXaxis()->SetTitle("Efficiency on #beta");
	RejectionGraph ->GetYaxis()->SetTitle("Rejection of #alpha");
	RejectionGraph ->Draw("a*");
	c1->Print("plots/RejectionVsEfficiency_alpha.png");

	std::vector<double> percent;
	std::vector<double> eff_Error;
	std::cout << "size of eff= "<<eff.size() << std::endl;
	std::cout << "size of rej= "<<rej.size()<< std::endl;
	std::cout << "size of rejection_Error = "<< rejection_Error.size()<< std::endl;
	std::vector<double> percent_Error;

	for (int i = 0; i < eff.size(); i++) {
		percent.push_back(1/rej[i]);	
		percent_Error.push_back(percent[i]*rejection_Error[i]/rej[i]);	
		// percent_Error.push_back(10e-7);	
		eff_Error.push_back(0);
		std::cout <<eff[i]<<"+/-"<<eff_Error[i] <<"\t"<< percent[i]<<"+/-"<<percent_Error[i] << std::endl;
	}

	TCanvas* c_new =new TCanvas();
	TGraphErrors* PercentageGraph = new TGraphErrors(percent.size(),&eff[0],&percent[0],&eff_Error[0],&percent_Error[0]);
	TF1 *f_percent_alpha = new TF1("f", "[3]*x*x*x +[2]*x*x +[1]*x +[0]");
	PercentageGraph ->Fit(f_percent_alpha);
	PercentageGraph ->SetTitle("Percentage Of #alpha's Remaining Across #beta Efficiency");
	PercentageGraph ->GetXaxis()->SetTitle("Efficiency on #beta");
	PercentageGraph ->GetYaxis()->SetTitle("Percentage Remaining");
	PercentageGraph ->Draw("ap");
	// c_new->SetLogy();
	c_new->SetGrid();
	c_new->Print("plots/PercentageVsEfficiency_alpha.tex");
	c_new->Print("plots/PercentageVsEfficiency_alpha.png");

	// TCanvas * c2= new TCanvas();
	// TGraph* RejectionGraph_rev = new TGraph(eff_rev.size(),&eff_rev[0],&rej_rev[0]);
	// TF1 *f2 = new TF1("f2", "[3]*x*x*x +[2]*x*x +[1]*x +[0]");
	// RejectionGraph_rev ->Fit(f2);
	// RejectionGraph_rev ->SetTitle("Rejection Across Efficiency");
	// RejectionGraph_rev ->GetXaxis()->SetTitle("Efficiency on #alpha");
	// RejectionGraph_rev ->GetYaxis()->SetTitle("Rejection of #beta");
	// RejectionGraph_rev ->Draw("a*");
	// c2->Print("plots/RejectionVsEfficiency_beta.png");

	// TCanvas * c3=new TCanvas();
	// c3->cd();
  //
	// TGraph* cutGraph = new TGraph(cuts_rev.size(),&energyValues_rev[0],&cuts_rev[0]);
  //       TF1 *f3 = new TF1("f3", "[3]*x*x*x +[2]*x*x +[1]*x +[0]");
  //       cutGraph->Fit(f3);
	// compareBi210->Draw();	
	// comparePo210->Draw("same");	
	// #<{(| for(int j=0;j<graphs.size();j++){ |)}>#
	// #<{(| 	graphs[j]->Draw("a* same"); |)}>#
	// #<{(| } |)}>#
	// cutGraph->Draw("a* same");
  //
	// c3->Print("plots/CutValue_rev.png");

	TFile fileout("RejectionsAndEfficiency.root","RECREATE");
	fileout.cd();
// 	-----Energy-----
	RejectionGraph->Write();
	// RejectionGraph_rev->Write();
// 	-----Radial-----
	/* compareEle_radial->Write(); */
	/* compareAlpha_radial->Write(); */
	/* RejectionRadialGraph->Write(); */
	/* RejectionRadial_can->Write(); */
	fileout.Close();

	return 0;

}

