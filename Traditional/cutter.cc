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


void FillHist(TFile* file,TH1D * hist,TH2D * hist2,TH1D *radial,TH2D *compareRadial,ofstream& outputfile){

	TTree* Tree = (TTree*) file->Get("output");
	Double_t para, mcEdepQuenched,posr;
	Bool_t  Qfit;
	Int_t pdg1, pdg2;
        Int_t parentpdg1,parentpdg2;

        Tree->SetBranchAddress("pdg1",&pdg1);
        Tree->SetBranchAddress("pdg2",&pdg2);
        Tree->SetBranchAddress("parentpdg1",&parentpdg1);
        Tree->SetBranchAddress("parentpdg2",&parentpdg2);
	
	Tree->SetBranchAddress("berkeleyAlphaBeta",&para);
	Tree->SetBranchAddress("mcEdepQuenched",&mcEdepQuenched);
	Tree->SetBranchAddress("posr",&posr);
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
			hist->Fill(para);
			hist2->Fill(mcEdepQuenched,para);
			radial->Fill(posr);
			compareRadial->Fill(posr,para);
			outputfile<< mcEdepQuenched<<","<<posr<<","<< para << endl;
		}
	}

}


void findCutsEnergy(TH2D *compareBi210,TH2D *comparePo210,vector<double>& energyValues, vector<double>& cutValues,vector<double>& mistagged_frac){
	/* This function should take a 2d histrogram split it into energy 
	 * strips and then find a cut value which retain ~99% of the signal.
	 */
	//for(double energy =0; energy<3.5;energy+=0.1){

	double startingCut=-1000;
	double step=0.1;
	double cutValue=startingCut;
	double accept=0;
	double mistagged=0;
	double mistagging_fraction=0;
	double sliceWidth=0.1;
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

	for(double energy =0; energy<1.7;energy+=0.1){
		cout<<"Energy = "<< energy<<endl;

		//TH1D* slice_bi=compareBi210->ProjectionY("slice_bi",bi_x->FindBin(0.5),bi_x->FindBin(1.0));
		TH1D* slice_bi=compareBi210->ProjectionY(SSTR(energy).c_str(),bi_x->FindBin(energy),bi_x->FindBin(energy+sliceWidth));
		TH1D* slice_po=comparePo210->ProjectionY(("Po"+SSTR(energy)).c_str(),po_x->FindBin(energy),po_x->FindBin(energy+sliceWidth));

		double bi_numEV=slice_bi->GetEntries();
		double po_numEV=slice_po->GetEntries();
		cout<<"Number of entries in bi slice = "<< bi_numEV<<endl;
		cout<<"Number of entries in po slice = "<< po_numEV<<endl;
		TAxis* bi_slice_x=slice_bi->GetXaxis();
		TAxis* po_slice_x=slice_po->GetXaxis();

		while(accept<threshold){

			double bi_remain=slice_bi->Integral(bi_slice_x->FindBin(startingCut),bi_slice_x->FindBin(cutValue+=step));
			//double po_remain=slice_po->Integral(po_slice_x->FindBin(startingCut),po_slice_x->FindBin(startingCut+step));
			cout<<"Number of entries remaining in bi slice = "<< bi_remain<<endl;
			//cout<<"Number of entries remaining in po slice = "<< po_remain<<endl;
			accept=bi_remain/bi_numEV;	
			cout<<"Accepted fraction = "<<accept<<endl;
		}
		mistagged=slice_po->Integral(po_slice_x->FindBin(startingCut),po_slice_x->FindBin(cutValue));

		mistagging_fraction= mistagged/po_numEV;
		cout<<mistagged<<" events were mistagged. "<<mistagging_fraction<<" as a fraction."<<endl;
		mistagged_frac.push_back(mistagging_fraction);
		cutValues.push_back(cutValue);
		energyValues.push_back(energy+0.05);

		TCanvas* c1 = new TCanvas();
		c1->cd();
		//TLine* cutLine = new TLine( cutValue, slice_bi->GetYaxis()->GetXmin(), cutValue, slice_bi->GetYaxis()->GetXmax() );
		TLine* cutLine = new TLine( cutValue, slice_bi->GetMinimum(), cutValue, slice_bi->GetMaximum() );
		cutLine->SetLineWidth(2);cutLine->SetLineStyle(4);
		slice_bi->Draw();
		slice_po->Draw("same");
		cutLine->Draw("same");

		cutValue=startingCut;
		accept=0;
	}

}

void findCutsRadial(TH2D *compareBi210,TH2D *comparePo210,vector<double>& energyValues, vector<double>& cutValues,vector<double>& mistagged_frac){
	/* This function should take a 2d histrogram split it into energy 
	 * strips and then find a cut value which retain ~99% of the signal.
	 */
	//for(double energy =0; energy<3.5;energy+=0.1){

	double startingCut=-1000;
	double step=0.1;
	double step_radial=500;
	double cutValue=startingCut;
	double accept=0;
	double mistagged=0;
	double mistagging_fraction=0;
	double sliceWidth=100;
	double threshold=0.9;
	TAxis* bi_x=compareBi210->GetXaxis();
	TAxis* bi_y=compareBi210->GetYaxis();
	TAxis* po_x=comparePo210->GetXaxis();
	TAxis* po_y=comparePo210->GetYaxis();

	TH1D* complete_bi=compareBi210->ProjectionY("complete_bi",bi_x->FindBin(0.),bi_x->FindBin(8000.));
	TH1D* complete_po=comparePo210->ProjectionY("complete_po",po_x->FindBin(0.),po_x->FindBin(8000.));
	double bi_numEV_complete= compareBi210->GetEntries();
	double po_numEV_complete= comparePo210->GetEntries();
	cout<<"Number of entries in complete bi = "<< bi_numEV_complete<<endl;
	cout<<"Number of entries in complete po = "<< po_numEV_complete<<endl;

	for(double energy =0; energy<6000;energy+=step_radial){
		cout<<"radius = "<< energy<<endl;

		//TH1D* slice_bi=compareBi210->ProjectionY("slice_bi",bi_x->FindBin(0.5),bi_x->FindBin(1.0));
		TH1D* slice_bi=compareBi210->ProjectionY(SSTR(energy).c_str(),bi_x->FindBin(energy),bi_x->FindBin(energy+sliceWidth));
		TH1D* slice_po=comparePo210->ProjectionY(("Po"+SSTR(energy)).c_str(),po_x->FindBin(energy),po_x->FindBin(energy+sliceWidth));

		double bi_numEV=slice_bi->GetEntries();
		double po_numEV=slice_po->GetEntries();
		cout<<"Number of entries in bi slice = "<< bi_numEV<<endl;
		cout<<"Number of entries in po slice = "<< po_numEV<<endl;
		TAxis* bi_slice_x=slice_bi->GetXaxis();
		TAxis* po_slice_x=slice_po->GetXaxis();

		while(accept<threshold){

			double bi_remain=slice_bi->Integral(bi_slice_x->FindBin(startingCut),bi_slice_x->FindBin(cutValue+=step));
			//double po_remain=slice_po->Integral(po_slice_x->FindBin(startingCut),po_slice_x->FindBin(startingCut+step));
			cout<<"Number of entries remaining in bi slice = "<< bi_remain<<endl;
			//cout<<"Number of entries remaining in po slice = "<< po_remain<<endl;
			accept=bi_remain/bi_numEV;	
			cout<<"Accepted fraction = "<<accept<<endl;
		}
		mistagged=slice_po->Integral(po_slice_x->FindBin(startingCut),po_slice_x->FindBin(cutValue));

		mistagging_fraction= mistagged/po_numEV;
		cout<<mistagged<<" events were mistagged. "<<mistagging_fraction<<" as a fraction."<<endl;
		mistagged_frac.push_back(mistagging_fraction);
		cutValues.push_back(cutValue);
		energyValues.push_back(energy+step_radial/2);

		TCanvas* c1 = new TCanvas();
		c1->cd();
		//TLine* cutLine = new TLine( cutValue, slice_bi->GetYaxis()->GetXmin(), cutValue, slice_bi->GetYaxis()->GetXmax() );
		TLine* cutLine = new TLine( cutValue, slice_bi->GetMinimum(), cutValue, slice_bi->GetMaximum() );
		cutLine->SetLineWidth(2);cutLine->SetLineStyle(4);
		slice_bi->Draw();
		slice_po->Draw("same");
		cutLine->Draw("same");
		cutValue=startingCut;
		accept=0;
	}

}

int main(){
//====================================================================	
	gStyle->SetOptStat(0);
	TH1D *hBi210   = new TH1D("hBi210","berkeleyAlphaBeta",100,-100,100);
	hBi210->SetLineColor(kBlue);hBi210->SetLineWidth(3);
	hBi210->SetFillColor(kBlue);
	TH1D *hPo210   = new TH1D("hPo210","berkeleyAlphaBeta",100,-100,100);
	hPo210->SetLineColor(kRed);hPo210->SetLineWidth(3);
	TH1D *hBlank   = new TH1D("hBlank","berkeleyAlphaBeta",100,-100,100);

	TH2D* compareBi210   = new TH2D("compareBi210","berkeleyAlphaBeta",23,0,2.3,260,-160,100);
	compareBi210->SetLineColor(kBlue);compareBi210->SetLineWidth(3);
	compareBi210->SetMarkerColor(kBlue);
	compareBi210->SetFillColor(kBlue);
	TH2D* comparePo210   = new TH2D("comparePo210","berkeleyAlphaBeta",23,0,2.3,260,-160,100);
	comparePo210->SetLineColor(kRed);comparePo210->SetLineWidth(3);
	comparePo210->SetMarkerColor(kRed);
	comparePo210->SetFillColor(kRed);
//==========================Radial Histrograms ========================
	TH1D *hEle_radial   = new TH1D("hEle_radial","berkeleyAlphaBeta",6000,0,6000);
	hEle_radial->SetLineColor(kBlue);hEle_radial->SetLineWidth(3);
	hEle_radial->SetFillColor(kBlue);
	TH1D *hAlpha_radial   = new TH1D("hAlpha_radial","berkeleyAlphaBeta",6000,0,6000);
	hAlpha_radial->SetLineColor(kRed);hAlpha_radial->SetLineWidth(3);

	TH2D* compareEle_radial  = new TH2D("compareEle_radial","berkeleyAlphaBeta",6000,0,6000,260,-160,200);
	compareEle_radial->SetLineColor(kBlue);compareEle_radial->SetLineWidth(3);
	compareEle_radial->SetMarkerColor(kBlue);
	compareEle_radial->SetFillColor(kBlue);
	TH2D* compareAlpha_radial   = new TH2D("compareAlpha_radial","berkeleyAlphaBeta",6000,0,6000,260,-160,100);
	compareAlpha_radial->SetLineColor(kRed);compareAlpha_radial->SetLineWidth(3);
	compareAlpha_radial->SetMarkerColor(kRed);
	compareAlpha_radial->SetFillColor(kRed);
	//TH2D *hBlank2   = new TH2D("hBlank2","berkeleyAlphaBeta",100,-100,100);
	
	TLegend* t1 = new TLegend( 0.6, 0.7, 0.89, 0.88 );
	TLegend* t2 = new TLegend( 0.11, 0.11, 0.21, 0.21 );
	TLegend* t3 = new TLegend( 0.11, 0.11, 0.21, 0.21 );
	
	ofstream File_bi;
	ofstream File_po;
	
	File_bi.open("Classifier_data_bi.dat");
	File_po.open("Classifier_data_po.dat");

	File_bi << "mcEdepQuenched,"<<"posr," << "BerekelyAlphaBeta" <<  endl;
	File_po << "mcEdepQuenched,"<<"posr," << "BerekelyAlphaBeta" << endl;

	vector<string> biFileList= glob("/data/snoplus/OfficialProcessing/production_5_3_0/Bi212","SolarBi212");
	vector<string> poFileList= glob("/data/snoplus/liggins/year1/fitting/fitting/alphaSims/output/ntuple","alpha");

	for( int i=0; i<biFileList.size(); i++ ){
		TFile * file= TFile::Open(biFileList[i].c_str());	
		partflag="Bi";
		FillHist( file, hBi210,compareBi210,hEle_radial,compareEle_radial,File_bi);
	  }

	for( int i=0; i<poFileList.size(); i++ ){
		TFile * file= TFile::Open(poFileList[i].c_str());	
		partflag="Po";
		FillHist( file, hPo210,comparePo210,hEle_radial,compareAlpha_radial,File_po);
	  }
	File_bi.close();
	File_po.close();

	vector<double> bi_cuts, energyValues,misstagging_frac;
	vector<double> radial_cuts, RadiusValues,misstagging_frac_radius;

	findCutsEnergy(compareBi210,comparePo210,energyValues,bi_cuts,misstagging_frac);
	findCutsRadial(compareEle_radial,compareAlpha_radial,RadiusValues,radial_cuts,misstagging_frac_radius);

//	for( int i =0; i<bi_cuts.size();i++){
//		cout<<i<<"th cut value = "<<bi_cuts[i]<<endl;
//		cout<<"The fraction of mistagged events = "<<misstagging_frac[i]<<endl;
//	}
	TGraph* cutGraph = new TGraph(bi_cuts.size(),&energyValues[0],&bi_cuts[0]);
	TF1 *f = new TF1("f", "[1]*x +[0]");
	cutGraph->Fit(f);

	TGraph* cutGraph_radial = new TGraph(radial_cuts.size(),&RadiusValues[0],&radial_cuts[0]);
	TF1 *f_Rad = new TF1("f_Rad", "[1]*x +[0]");
	cutGraph_radial->Fit(f_Rad);



//	TCanvas * c1 = new TCanvas();
//	hBlank->SetMaximum(0.3);
//	hBlank->SetTitle("Standard Classifiers From Production 5.3");
//	//hBlank->SetTitle("Classifers Recoordinated with 1 MeV #beta & 5 MeV #alpha");
//	//hBlank->SetTitle("Classifers Recoordinated with 5 MeV #beta & 50 MeV #alpha");
//	hBlank->Draw();
//	cout<<"# of entries in hBi210 = "<< hBi210->GetEntries()<<endl;
//	cout<<"# of entries in hPo210 = "<< hPo210->GetEntries()<<endl;
//	hBi210->DrawNormalized("same");
//	hPo210->DrawNormalized("same");
//
//	t1->AddEntry( hBi210, "Bi 212","f");
//	t1->AddEntry( hPo210, "Po 212","f");
//	t1->Draw();
//	//c1->Print("BerkeleyAlphaBeta_NormalWithRAT5_3.png");
//	//c1->Print("output2.png");
//	//c1->Print("production_5_3.png");
//	c1->Print("StandardClassifersFromProduction5_3_212.png");
//	//c1->Print("StandardClassifersFromProduction5_3.png");
//	//c1->Print("RecoodinatedClassifiers_1MeV_beta_5MeV_alpha.png");
//	//c1->Print("RecoodinatedClassifiers_5MeV_beta_50MeV_alpha.png");
//
	TCanvas * c2 = new TCanvas();
	cout<<"# of entries in hBi210 = "<< hBi210->GetEntries()<<endl;
	cout<<"# of entries in hPo210 = "<< hPo210->GetEntries()<<endl;
	compareBi210->SetTitle("Standard Classifiers From Production 5.3");
	compareBi210->GetXaxis()->SetTitle("mcEdepQuenched (MeV)");
	compareBi210->GetYaxis()->SetTitle("BerkeleyAlphaBeta");
	compareBi210->Draw();
	comparePo210->Draw("same");
	//cutGraph->Draw("AP*1");
	cutGraph->Draw("* same");

	t2->AddEntry( compareBi210, "Bi 210","f");
	t2->AddEntry( comparePo210, "Po 210","f");
	t2->Draw();
	c2->Print("compareBerkeleyAlphaBetaVsMCEnergy_withCuts.png");
	TCanvas* mistagged_energy = new TCanvas();
	TGraph* mistaggedEnergyGraph = new TGraph(bi_cuts.size(),&energyValues[0],&misstagging_frac[0]);
	//TF1 *f = new TF1("f", "[1]*x +[0]");
	//cutGraph->Fit(f);
	mistaggedEnergyGraph->Draw();

	TCanvas * c3 = new TCanvas();
	cout<<"# of entries in hBi210 = "<< hBi210->GetEntries()<<endl;
	cout<<"# of entries in hPo210 = "<< hPo210->GetEntries()<<endl;
	compareEle_radial->SetTitle("Standard Classifiers From Production 5.3");
	compareEle_radial->GetXaxis()->SetTitle("posr (mm)");
	compareEle_radial->GetYaxis()->SetTitle("BerkeleyAlphaBeta");
	compareEle_radial->Draw();
	compareAlpha_radial->Draw("same");
	//cutGraph->Draw("AP*1");
	cutGraph_radial->Draw("* same");

	t3->AddEntry( compareBi210, "Bi 210","f");
	t3->AddEntry( comparePo210, "Po 210","f");
	t3->Draw();
	c3->Print("compareBerkeleyAlphaBetaVsPosr_withCuts.png");

	TCanvas* mistagged_raidal = new TCanvas();
	TGraph* mistaggedRadialGraph = new TGraph(radial_cuts.size(),&RadiusValues[0],&misstagging_frac_radius[0]);
	//TF1 *f = new TF1("f", "[1]*x +[0]");
	//cutGraph->Fit(f);
	mistaggedRadialGraph->Draw();


//	TFile fileout("plot.root","RECREATE");
//	fileout.cd();
//	hBi210->Write();
//	hPo210->Write();
//	fileout.Close();
	return 0;

}

