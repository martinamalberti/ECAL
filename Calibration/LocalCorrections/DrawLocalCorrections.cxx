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



DrawLocalCorrections(string filename, bool localEta, string outdir, string outname, bool drawFit, bool drawMC, bool drawData, string corrType){
  
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

  
  bool drawRatios = false;
  
  const int Ntempl = 4;

  TFile *f = TFile::Open(filename.c_str());
  
  TGraphErrors* g_EoP_MC[Ntempl];
  TGraphErrors* g_EoC_MC[Ntempl];
   
  TGraphErrors* g_EoP_Data[Ntempl];
  TGraphErrors* g_EoC_Data[Ntempl];

  TGraphErrors* g_ratio_MC[Ntempl];
  TGraphErrors* g_ratio_Data[Ntempl];

  TGraphErrors* g_ratio_uncorr[Ntempl];
  TGraphErrors* g_ratio_corr[Ntempl];

  for (int mod=0; mod<Ntempl; mod++){
    char histoName[100];
    
    sprintf(histoName, "gEoP_MC_mod%d", mod+1);
    g_EoP_MC[mod]   = (TGraphErrors*)f->Get(histoName); 
    
    sprintf(histoName, "gEoC_MC_mod%d", mod+1);
    g_EoC_MC[mod]   = (TGraphErrors*)f->Get(histoName); 
    
    sprintf(histoName, "gEoP_Data_mod%d", mod+1);
    g_EoP_Data[mod]   = (TGraphErrors*)f->Get(histoName); 
    
    sprintf(histoName, "gEoC_Data_mod%d", mod+1);
    g_EoC_Data[mod]   = (TGraphErrors*)f->Get(histoName); 

    sprintf(histoName, "gRatio_MC_mod%d", mod+1);
    g_ratio_MC[mod]   = (TGraphErrors*)f->Get(histoName); 
 
    sprintf(histoName, "gRatio_Data_mod%d", mod+1);
    g_ratio_Data[mod]   = (TGraphErrors*)f->Get(histoName); 
    
    sprintf(histoName, "gRatio_uncorr_mod%d", mod+1);
    g_ratio_uncorr[mod]  = (TGraphErrors*)f->Get(histoName); 
    
    sprintf(histoName, "gRatio_corr_mod%d", mod+1);
    g_ratio_corr[mod]  = (TGraphErrors*)f->Get(histoName);  

    //scale to have middle point at 1
    ScaleGraph(g_EoP_MC[mod]);
    ScaleGraph(g_EoC_MC[mod]);
    ScaleGraph(g_EoP_Data[mod]);
    ScaleGraph(g_EoC_Data[mod]);

    CheckGraphErrors(g_EoP_MC[mod]);
    CheckGraphErrors(g_EoC_MC[mod]);
    CheckGraphErrors(g_EoP_Data[mod]);
    CheckGraphErrors(g_EoC_Data[mod]);

    // set colors for plotting
    //-- mc uncorrected
    g_EoP_MC[mod] -> SetMarkerStyle(21);
    g_EoP_MC[mod] -> SetMarkerSize(1.);
    g_EoP_MC[mod] -> SetMarkerColor(kRed); 
    g_EoP_MC[mod] -> SetLineColor(kRed); 
    
    //-- mc corrected
    g_EoC_MC[mod] -> SetMarkerStyle(25);
    g_EoC_MC[mod] -> SetMarkerSize(1.);
    g_EoC_MC[mod] -> SetMarkerColor(kRed); 
    g_EoC_MC[mod] -> SetLineColor(kRed); 
  
    //-- data uncorrected
    g_EoP_Data[mod] -> SetMarkerStyle(20);
    g_EoP_Data[mod] -> SetMarkerSize(1.);
    g_EoP_Data[mod] -> SetMarkerColor(kBlue+2); 
    g_EoP_Data[mod] -> SetLineColor(kBlue+2); 

    //-- data corrected
    g_EoC_Data[mod] -> SetMarkerStyle(24);
    g_EoC_Data[mod] -> SetMarkerSize(1.);
    g_EoC_Data[mod] -> SetMarkerColor(kBlue+2); 
    g_EoC_Data[mod] -> SetLineColor(kBlue+2); 
 

    // ratio
    g_ratio_MC[mod]->SetLineColor(kMagenta+1);
    g_ratio_MC[mod]->SetMarkerColor(kMagenta+1);
    g_ratio_MC[mod]->SetMarkerStyle(21);
    g_ratio_MC[mod]->SetMarkerSize(1);
    
    g_ratio_Data[mod]->SetLineColor(kMagenta+3);
    g_ratio_Data[mod]->SetMarkerColor(kMagenta+3);
    g_ratio_Data[mod]->SetMarkerStyle(20);
    g_ratio_Data[mod]->SetMarkerSize(1);

    g_ratio_uncorr[mod]->SetLineColor(kAzure+2);
    g_ratio_uncorr[mod]->SetMarkerColor(kAzure+2);
    g_ratio_uncorr[mod]->SetMarkerStyle(20);
    g_ratio_uncorr[mod]->SetMarkerSize(1);
    
    g_ratio_corr[mod]->SetLineColor(kBlue+2);
    g_ratio_corr[mod]->SetMarkerColor(kBlue+2);
    g_ratio_corr[mod]->SetMarkerStyle(20);
    g_ratio_corr[mod]->SetMarkerSize(1);

  }
  

  TCanvas* c_g_fit[Ntempl];
  TLegend* tl[Ntempl];
  TLegend* tlr[Ntempl];
  TLegend* tlrr[Ntempl];

  TLatex *latex[Ntempl];
  
  for(int mod=0;mod<4;mod++) {
    char padName[100];
    sprintf(padName, "g_fit_mod%d", mod+1);
    latex[mod] = new TLatex(0.70,0.86,Form("module %d",mod+1));
    latex[mod]->SetNDC();
    latex[mod]->SetTextSize(0.04);


    int locmin = TMath::LocMin(g_EoP_Data[mod]->GetN(),g_EoP_Data[mod]->GetX());
    float ymin = g_EoP_Data[mod]->GetY()[locmin];
    float ymax = 1.01;

    if (!drawRatios){
      c_g_fit[mod] = new TCanvas(padName,padName,600,800);
      c_g_fit[mod] ->SetGridy();
    
      float tXoffset = 1.; 
      float tYoffset = 1.6; 
      float labSize = 0.04;
      
      TH1F *hPad;
      if (localEta){
	ymin = 0.96;
	ymax = 1.02;
	hPad = (TH1F*)gPad->DrawFrame(-0.5,ymin,0.5,ymax);
	hPad->GetXaxis()->SetTitle("local #eta");
      }
      else{
	ymin = 0.96;
	ymax = 1.02;
	hPad = (TH1F*)gPad->DrawFrame(-0.5,ymin,0.5,ymax);      
	hPad->GetXaxis()->SetTitle("local #phi");
      }

      hPad->GetYaxis()->SetTitle("Relative E/p scale");
      hPad->GetXaxis()->SetLabelSize(labSize);
      hPad->GetXaxis()->SetTitleSize(labSize);
      hPad->GetYaxis()->SetLabelSize(labSize);
      hPad->GetYaxis()->SetTitleSize(labSize);
      hPad->GetXaxis()->SetTitleOffset(tXoffset);
      hPad->GetYaxis()->SetTitleOffset(tYoffset);
      
      // ----- Drawing--------------------

      if (!drawFit){
	tl[mod] = new TLegend(0.15,0.71,0.65,0.90);
	tl[mod] -> SetFillStyle(0);
	tl[mod] -> SetBorderSize(0);

     	if (drawMC) {
	  g_EoP_MC[mod]   -> Draw("PL");
	  g_EoC_MC[mod]   -> Draw("PL");
	  tl[mod] -> AddEntry(g_EoP_MC[mod],"MC - w/o correction","PL");
	  if (corrType=="regression") tl[mod] -> AddEntry(g_EoC_MC[mod],"MC - w/ regression","PL");
	  else tl[mod] -> AddEntry(g_EoC_MC[mod],"MC - w/ correction","PL");
	}
	if (drawData){
	  g_EoP_Data[mod] -> Draw("PL");
	  g_EoC_Data[mod] -> Draw("PL");
	  tl[mod] -> AddEntry(g_EoP_Data[mod],"DATA - w/o correction","PL");
	  if (corrType=="regression") tl[mod] -> AddEntry(g_EoC_Data[mod],"DATA - w/ regression","PL");
	  else tl[mod] -> AddEntry(g_EoC_Data[mod],"DATA - w/ correction","PL");
	}
      }
      else {
	TF1 *fmc = new TF1("fmc","pol2");
	fmc->SetLineColor(kRed);	
	fmc->SetLineStyle(2);	
	g_EoP_MC[mod]   -> Fit("fmc");
	g_EoP_MC[mod]   -> Draw("PL");

	TF1 *fdata = new TF1("fdata","pol2");
	fdata->SetLineColor(kBlue+2);	
	fdata->SetLineStyle(2);	
	g_EoP_Data[mod] -> Fit("fdata");
	g_EoP_Data[mod] -> Draw("PL");
	
	tl[mod] = new TLegend(0.15,0.76,0.65,0.90);
	tl[mod] -> SetFillStyle(0);
	tl[mod] -> SetBorderSize(0);
	tl[mod] -> AddEntry(g_EoP_MC[mod],"MC - w/o correction","PL");
	tl[mod] -> AddEntry(g_EoP_Data[mod],"DATA - w/o correction","PL");
      }

      latex[mod]->Draw("same");
      tl[mod] -> Draw();

    }
    else {
      c_g_fit[mod] = new TCanvas(padName,padName,100,100,700,700);

      TPad *cLower  = new TPad("pad_0","pad_0",0.00,0.00,1.00,0.25);
      TPad *cCenter = new TPad("pad_1","pad_1",0.00,0.25,1.00,0.50);
      TPad *cUpper  = new TPad("pad_2","pad_2",0.00,0.50,1.00,1.00);
      
      cLower->SetBottomMargin(0.2); 
      cLower->SetTopMargin(0.05);
      
      cCenter->SetBottomMargin(0.1); 
      cCenter->SetTopMargin(0.05);
      
      cUpper->SetBottomMargin(0.1); 
      cUpper->SetTopMargin(0.05);
      
      cLower->Draw();
      cCenter->Draw();
      cUpper->Draw();
      
      float FontSCF = cUpper->GetHNDC()/cLower->GetHNDC(); 
      float tYoffset = 0.8; 
      float labSize = 0.05;
      
      cUpper-> cd();
      gPad->SetGrid();
      
      
      int locmin = TMath::LocMin(g_EoP_Data[mod]->GetN(),g_EoP_Data[mod]->GetX());
      //    cout<< locmin << "  " << g_EoP_Data[mod]->GetY()[locmin] << "  " << g_EoP_Data[mod]->GetX()[locmin]<<endl;
      float ymin = g_EoP_Data[mod]->GetY()[locmin];
      TH1F *hPad = (TH1F*)gPad->DrawFrame(-0.5,ymin*0.99,0.5,1.005);
      hPad->GetXaxis()->SetTitle("#eta_{SC} (deg)");
      hPad->GetYaxis()->SetTitle("Relative E/p scale");
      hPad->GetXaxis()->SetLabelSize(labSize);
      hPad->GetXaxis()->SetTitleSize(labSize);
      hPad->GetYaxis()->SetLabelSize(labSize);
      hPad->GetYaxis()->SetTitleSize(labSize);
      hPad->GetXaxis()->SetTitleOffset(tYoffset);
      hPad->GetYaxis()->SetTitleOffset(tYoffset);
      
      // ----- Drawing--------------------
      g_EoP_MC[mod]   -> Draw("PL");
      g_EoC_MC[mod]   -> Draw("PL");
      g_EoP_Data[mod] -> Draw("PL");
      g_EoC_Data[mod] -> Draw("PL");
      
      //----- legend
      tl[mod] = new TLegend(0.60,0.12,0.89,0.35);
      tl[mod] -> SetFillColor(0);
      
      tl[mod] -> AddEntry(g_EoP_MC[mod],"MC - default","PL");
      tl[mod] -> AddEntry(g_EoC_MC[mod],"MC - w/ regression","PL");
      tl[mod] -> AddEntry(g_EoP_Data[mod],"DATA - default","PL");
      tl[mod] -> AddEntry(g_EoC_Data[mod],"DATA - w/ regression","PL");
      tl[mod] -> Draw();
      
      
      //--------- RATIO PLOTS ------------------------------------------
      
      //    c_g_fit[mod]->cd(2);
      
      cCenter-> cd();
      gPad->SetGrid();
      TH1F *hPad2 = (TH1F*)gPad->DrawFrame(-0.55,0.985,0.55,1.015);
      hPad2->GetXaxis()->SetTitle("#eta_{SC} (deg)");
      hPad2->GetYaxis()->SetTitle("regr./default ratio");
      hPad2->GetXaxis()->SetLabelSize(labSize*FontSCF);
      hPad2->GetXaxis()->SetTitleSize(labSize*FontSCF);
      hPad2->GetYaxis()->SetLabelSize(labSize*FontSCF);
      hPad2->GetYaxis()->SetTitleSize(labSize*FontSCF);
      hPad2->GetYaxis()->SetTitleOffset(tYoffset/FontSCF);
      hPad2->GetYaxis()->SetNdivisions(505);
      g_ratio_MC[mod]->Draw("PL");
      g_ratio_Data[mod]->Draw("PL");
      
      tlr[mod] = new TLegend(0.60,0.15,0.89,0.35);
      tlr[mod] -> SetFillColor(0);
      tlr[mod] -> AddEntry(g_ratio_MC[mod],"MC","PL");
      tlr[mod] -> AddEntry(g_ratio_Data[mod],"DATA","PL");
      tlr[mod] -> Draw();
      
      cLower-> cd();
      gPad->SetGrid();
      TH1F *hPad3 = (TH1F*)gPad->DrawFrame(-0.55,0.985,0.55,1.015);
      hPad3->GetXaxis()->SetTitle("#eta_{SC} (deg)");
      hPad3->GetYaxis()->SetTitle("data/MC ratio");
      hPad3->GetXaxis()->SetLabelSize(labSize*FontSCF);
      hPad3->GetXaxis()->SetTitleSize(labSize*FontSCF);
      hPad3->GetYaxis()->SetLabelSize(labSize*FontSCF);
      hPad3->GetYaxis()->SetTitleSize(labSize*FontSCF);
      hPad3->GetYaxis()->SetTitleOffset(tYoffset/FontSCF);
      hPad3->GetXaxis()->SetTitleOffset(0.9);
      hPad3->GetYaxis()->SetNdivisions(505);
      g_ratio_uncorr[mod]->Draw("PL");
      g_ratio_corr[mod]->Draw("PL");
      
      
      tlrr[mod] = new TLegend(0.60,0.25,0.89,0.45);
      tlrr[mod] -> SetFillColor(0);
      tlrr[mod] -> AddEntry(g_ratio_uncorr[mod],"default","PL");
      tlrr[mod] -> AddEntry(g_ratio_corr[mod],"w/ regression","PL");
      tlrr[mod] -> Draw();
    }

    
    //--save plots
    if (mod==0) gSystem->mkdir(outdir.c_str(),true);
    //    gSystem->cd(outdir.c_str());
    string cname = outdir+"/"+outname+Form("_mod%d",mod+1);
    c_g_fit[mod]->SaveAs((cname+".pdf").c_str());
    c_g_fit[mod]->SaveAs((cname+".png").c_str());
    c_g_fit[mod]->SaveAs((cname+".C").c_str());
    
  }





}
