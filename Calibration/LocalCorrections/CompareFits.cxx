void CompareFits(){

  TFile *f1 = TFile::Open("localEta/graphs.root");
  TFile *f2 = TFile::Open("localEta_CBfit/graphs.root");  

  string hname = "EoP_Data_19_mod4";

  f1->cd();
  TH1F* h1 = (TH1F*)f1->Get(hname.c_str());

  f2->cd();
  TH1F*  h2 = (TH1F*)f2->Get(hname.c_str());

  h1->SetMarkerStyle(20);
  h1->GetXaxis()->SetRangeUser(0.7,1.3);
  h1->Draw("e");

  for (int i =0; i<h2->GetListOfFunctions()->GetSize(); i++){
    (h2->GetListOfFunctions()->At(i))->SetBit(TF1::kNotDraw);  
  }

  h2->Draw("esame");
  gPad->Update();
  int n = h2->GetListOfFunctions()->GetSize();
  TF1* ftemp = (TF1*)(h2->GetListOfFunctions()->At(0));
  ftemp->SetLineColor(4);
  ftemp->SetLineStyle(2);
  ftemp->Draw("same");

}
