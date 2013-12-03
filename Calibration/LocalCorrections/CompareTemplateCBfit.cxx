void  ScaleGraph(TGraphErrors *g){
  
  //int midX = int(g->GetN()/2);
  //float scale = 1./g->GetY()[midX-1];
  float scale = 1./g->Eval(0.);
  for (int i=0;i<g->GetN();i++) g->GetY()[i] *= scale;
  
}

void  CheckGraphErrors(TGraphErrors *g){

  for (int i=0;i<g->GetN();i++) {
    if (g->GetErrorY(i) < 0.0005){
      float aveErr = (g->GetErrorY(0) + g->GetErrorY(g->GetN()-1))/2;
      g->SetPointError(i,0,aveErr);
    }
  }
  
}

void SkipBadFits(TGraphErrors *g){
  for (int i=0;i<g->GetN();i++) {
    if (g->GetErrorY(i) > 1){
      g->RemovePoint(i);
    }
  }
}

TGraphErrors* DiffGraph(TGraphErrors *g1, TGraphErrors *g2 ){
  TGraphErrors *g = new TGraphErrors();
  for (int i=0; i < g1->GetN(); i++){
    g->SetPoint(i,g1->GetX()[i],g1->GetY()[i]-g2->GetY()[i]);
    double err1 = g1->GetErrorY(i);
    double err2 = g2->GetErrorY(i);
    double err = sqrt(err1*err1+err2*err2);
    g->SetPointError(i,0,err);
  }
  g->SetMarkerStyle(20);
  g->SetMarkerSize(1);
  g->GetHistogram()->GetYaxis()->SetRangeUser(-0.01,0.01);
  g->GetHistogram()->GetYaxis()->SetTitle("diff");
  g->GetHistogram()->GetYaxis()->SetTitleOffset(1.4);
  return g;
}

void DrawGraphs(TGraphErrors *g1, TGraphErrors *g2, string title){
  g1   -> GetHistogram()->GetYaxis()->SetRangeUser(0.96,1.02);
  g1   -> GetXaxis()->SetTitle(title.c_str());
  g1   -> GetYaxis()->SetTitle("relative E/P scale");
  g1   -> GetYaxis()->SetTitleOffset(1.3);
  g1   -> Draw("APL");
  g2   -> Draw("PL");
}



CompareTemplateCBfit(string filename1, string filename2, bool localEta, string outdir){
  
  gROOT->SetStyle("Plain");
  gStyle->SetTextFont(42);
  gStyle->SetTextSize(0.05);
  gStyle->SetTitleFont(42,"xyz");
  gStyle->SetTitleSize(0.05);
  gStyle->SetLabelFont(42,"xyz");
  gStyle->SetLabelSize(0.05);
  gStyle->SetStatFont(42);
  gStyle->SetStatFontSize(0.05);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.1);
  gROOT->ForceStyle();

  const int Ntempl = 4;

  TFile *f[2];
  f[0] = TFile::Open(filename1.c_str());
  f[1] = TFile::Open(filename2.c_str());
  
  TGraphErrors* g_EoP_MC[Ntempl][2];
  TGraphErrors* g_EoC_MC[Ntempl][2];
  TGraphErrors* g_EoP_Data[Ntempl][2];
  TGraphErrors* g_EoC_Data[Ntempl][2];
  int style[2] = {20,24};

  for (int ifile =0; ifile <2; ifile++){
    for (int mod=0; mod<Ntempl; mod++){
      char histoName[100];
      
      sprintf(histoName, "gEoP_MC_mod%d", mod+1);
      g_EoP_MC[mod][ifile]   = (TGraphErrors*)f[ifile]->Get(histoName); 
      
      sprintf(histoName, "gEoC_MC_mod%d", mod+1);
      g_EoC_MC[mod][ifile]    = (TGraphErrors*)f[ifile] ->Get(histoName); 
      
      sprintf(histoName, "gEoP_Data_mod%d", mod+1);
      g_EoP_Data[mod][ifile]    = (TGraphErrors*)f[ifile] ->Get(histoName); 
      
      sprintf(histoName, "gEoC_Data_mod%d", mod+1);
      g_EoC_Data[mod][ifile]    = (TGraphErrors*)f[ifile] ->Get(histoName); 
      

    //scale to have middle point at 1
    ScaleGraph(g_EoP_MC[mod][ifile] );
    ScaleGraph(g_EoC_MC[mod][ifile] );
    ScaleGraph(g_EoP_Data[mod][ifile] );
    ScaleGraph(g_EoC_Data[mod][ifile] );

    CheckGraphErrors(g_EoP_MC[mod][ifile] );
    CheckGraphErrors(g_EoC_MC[mod][ifile] );
    CheckGraphErrors(g_EoP_Data[mod][ifile] );
    CheckGraphErrors(g_EoC_Data[mod][ifile] );

    SkipBadFits(g_EoP_MC[mod][ifile] );
    SkipBadFits(g_EoC_MC[mod][ifile] );
    SkipBadFits(g_EoP_Data[mod][ifile] );
    SkipBadFits(g_EoC_Data[mod][ifile] );

    // set colors for plotting
    //-- mc uncorrected
    g_EoP_MC[mod][ifile]  -> SetMarkerStyle(style[ifile]);
    g_EoP_MC[mod][ifile]  -> SetMarkerSize(1.);
    g_EoP_MC[mod][ifile]  -> SetMarkerColor(kRed); 
    g_EoP_MC[mod][ifile]  -> SetLineColor(kRed); 
    
    //-- mc corrected
    g_EoC_MC[mod][ifile]  -> SetMarkerStyle(style[ifile]);
    g_EoC_MC[mod][ifile]  -> SetMarkerSize(1.);
    g_EoC_MC[mod][ifile]  -> SetMarkerColor(kRed); 
    g_EoC_MC[mod][ifile]  -> SetLineColor(kRed); 
  
    //-- data uncorrected
    g_EoP_Data[mod][ifile] -> SetMarkerStyle(style[ifile]);
    g_EoP_Data[mod][ifile] -> SetMarkerSize(1.);
    g_EoP_Data[mod][ifile]-> SetMarkerColor(kBlue+2); 
    g_EoP_Data[mod][ifile] -> SetLineColor(kBlue+2); 

    //-- data corrected
    g_EoC_Data[mod][ifile] -> SetMarkerStyle(style[ifile]);
    g_EoC_Data[mod][ifile] -> SetMarkerSize(1.);
    g_EoC_Data[mod][ifile] -> SetMarkerColor(kBlue+2); 
    g_EoC_Data[mod][ifile] -> SetLineColor(kBlue+2); 
    }
  } 
  
  TCanvas* c_g_fit[Ntempl];
  TLegend* tl[Ntempl];
  TLegend* tlr[Ntempl];
  TLegend* tlrr[Ntempl];

  float ymin = 0.96;
  float ymax = 1.02;  
  string title = "local #eta";
  if (!localEta) title = "local #Phi";

  for(int mod=0;mod<Ntempl;mod++) {

    char padName[100];
    sprintf(padName, "g_fit_mod%d", mod+1);
    c_g_fit[mod] = new TCanvas(padName,padName,1000,500);
    c_g_fit[mod] ->Divide(4,2);

    for (int i=0; i < 8; i++){
      c_g_fit[mod] ->cd(i+1)->SetGridy();
    }
    
    for (int i=0; i < 4; i++){
      c_g_fit[mod] ->cd(i+1);
      
      // ----- Drawing--------------------
      if (i==0){
	DrawGraphs(g_EoP_MC[mod][0],g_EoP_MC[mod][1],title);
	TGraphErrors *grEoP = DiffGraph(g_EoP_MC[mod][1],g_EoP_MC[mod][0]);
	c_g_fit[mod] ->cd(i+5);
	grEoP->Draw("AP");
	c_g_fit[mod] ->cd(i+1);
	TLatex *latex2 = new TLatex(0.15,0.92,Form("module %d - MC w/o corr", mod+1));
	latex2->SetNDC();
	latex2->SetTextSize(0.05);
	latex2->Draw("same");      
      }

      if (i==1){
	DrawGraphs(g_EoC_MC[mod][0],g_EoC_MC[mod][1],title);
     	TGraphErrors *grEoC = DiffGraph(g_EoC_MC[mod][1],g_EoC_MC[mod][0]);
	c_g_fit[mod] ->cd(i+5);
	grEoC->Draw("AP");
	c_g_fit[mod] ->cd(i+1);
	TLatex *latex2 = new TLatex(0.15,0.92,Form("module %d - MC w/ regr", mod+1));
	latex2->SetNDC();
	latex2->SetTextSize(0.05);
	latex2->Draw("same");      
      }

      if (i==2){
	DrawGraphs(g_EoP_Data[mod][0],g_EoP_Data[mod][1],title);
	TGraphErrors *grEoPda = DiffGraph(g_EoP_Data[mod][1],g_EoP_Data[mod][0]);
	c_g_fit[mod] ->cd(i+5);
	grEoPda->Draw("AP");
	c_g_fit[mod] ->cd(i+1);
	TLatex *latex2 = new TLatex(0.15,0.92,Form("module %d - DATA w/o corr", mod+1));
	latex2->SetNDC();
	latex2->SetTextSize(0.05);
	latex2->Draw("same");      
      }

      if (i==3){
	DrawGraphs(g_EoC_Data[mod][0],g_EoC_Data[mod][1],title);
	TGraphErrors *grEoCda = DiffGraph(g_EoC_Data[mod][1],g_EoC_Data[mod][0]);
	c_g_fit[mod] ->cd(i+5);
	grEoCda->Draw("AP");
	c_g_fit[mod] ->cd(i+1);
	TLatex *latex2 = new TLatex(0.15,0.92,Form("module %d - DATA w/ regr", mod+1));
	latex2->SetNDC();
	latex2->SetTextSize(0.05);
	latex2->Draw("same");      
      }

      tl[mod] = new TLegend(0.15,0.71,0.65,0.90);
      tl[mod] -> SetFillStyle(0);
      tl[mod] -> SetBorderSize(0);
      tl[mod] -> AddEntry(g_EoP_MC[mod][0],"template","PL");
      tl[mod] -> AddEntry(g_EoP_MC[mod][1],"double CB","PL");
      tl[mod] -> Draw("same");

    }
   
  string cname = outdir+"/compareFits"+Form("_mod%d",mod+1);
  c_g_fit[mod]->SaveAs((cname+".pdf").c_str());
  c_g_fit[mod]->SaveAs((cname+".png").c_str());
  c_g_fit[mod]->SaveAs((cname+".C").c_str());

  }
  
}
