
void plot(){

	gStyle->SetOptStat(0);
	TFile * Bi210 = TFile::Open("Bi210_100.root");
	TFile * Po210 = TFile::Open("Po210_100.root");
	TTree* Bi210Tree = (TTree*) Bi210->Get("output");
	TTree* Po210Tree = (TTree*) Po210->Get("output");
	Double_t Bi_ab, Po_ab;

	Bi210Tree->SetBranchAddress("berkeleyAlphaBeta",&Bi_ab);
	Po210Tree->SetBranchAddress("berkeleyAlphaBeta",&Po_ab);

	TH1D *hBi210   = new TH1D("hBi210","berkeleyAlphaBeta",100,-50,50);
	hBi210->SetLineColor(1);
	hBi210->SetLineWidth(3);
	TH1D *hPo210   = new TH1D("hPo210","berkeleyAlphaBeta",100,-50,50);
	hPo210->SetLineColor(4);
	hPo210->SetLineWidth(3);
	TLegend* t1 = new TLegend( 0.6, 0.7, 0.89, 0.88 );

	Int_t nBi = (Int_t)Bi210Tree->GetEntries();
	Int_t nPo = (Int_t)Po210Tree->GetEntries();

	for( Int_t i =0;i<nBi;i++){
	Bi210Tree->GetEntry(i);	
	hBi210->Fill(Bi_ab);
	}
	
	for( Int_t i =0;i<nPo;i++){
	Po210Tree->GetEntry(i);	
	hPo210->Fill(Po_ab);
	}
	
	hBi210->Draw();
	hPo210->Draw("same");

	t1->AddEntry( hBi210, "Bi 210","f");
	t1->AddEntry( hPo210, "Po 210","f");
	t1->Draw();
}


