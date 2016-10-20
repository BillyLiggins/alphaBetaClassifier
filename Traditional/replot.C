void replot(){

				TFile *file = TFile::Open("plots/scan.root");

				TGraphErrors* r4=(TGraphErrors*)file->Get("Scan_r_4000.000000;1");
				TGraphErrors* r5=(TGraphErrors*)file->Get("Scan_r_5000.000000;1");
				TGraphErrors* r6=(TGraphErrors*)file->Get("Scan_r_6000.000000;1");

				r4->SetName("rad4000");
				r5->SetName("rad5000");
				r6->SetName("rad6000");


				TCanvas* c1 = new TCanvas();
				c1->cd();
				// c1->SetLogy();
				TLegend* t1 = new TLegend( 0.1, 0.7, 0.3, 0.9 );
				gStyle->SetOptStat(kFALSE); 

				r4->SetLineColor(3);
				r5->SetLineColor(2);
				r6->SetLineColor(1);

				r4->SetLineWidth(2);
				r5->SetLineWidth(2);
				r6->SetLineWidth(2);

				r4->GetXaxis()->SetTitle("Remaining beta %");
				r4->GetYaxis()->SetTitle("Remaining alpha %");
				r4->SetTitle("Remaining #alpha's over #beta's efficiency");
				r4->GetXaxis()->CenterTitle();
				r4->GetYaxis()->CenterTitle();

				r4->SetMinimum(-1);
				r4->SetMaximum(8);
				r4->Draw("ap");
				r5->Draw("same e");
				r6->Draw("same e");

				t1->AddEntry( r4, "mcPosr<4000", "l" );
				t1->AddEntry( r5, "mcPosr<5000", "l" );
				t1->AddEntry( r6, "mcPosr<6000", "l" );

				t1->Draw();
				c1->Print("Efficiencies_vs_Rejection.png");
				c1->Print("Efficiencies_vs_Rejection.tex");

}



