void plot_ALU_projections(vector<vector<vector<double>>> mean_t, 
                          vector<vector<vector<double>>> mean_t_err, 
                          vector<vector<vector<double>>> alu_t_x, 
                          vector<vector<vector<double>>> alu_t_x_err){

   TCanvas *Can = new TCanvas("Can", "Helium DVCS",0,0,750,600);
   gStyle->SetOptFit(1);
   Can->Range(0,0,1,1);
   Can->SetFillColor(0);
   Can->SetBorderMode(0);
   Can->SetBorderSize(2);
   Can->SetLeftMargin(0.16);
   Can->SetRightMargin(0.03);
   Can->SetTopMargin(0.03);
   Can->SetBottomMargin(0.15);
   Can->SetFrameBorderMode(0);

   TGraphErrors *gre;
   TGraph *graph;
   TGraph *graphCFF;
   TLine * line;
   TLegend *leg;
   Can->Draw();
   Can->cd();
   
   TH2F *HALU_He_vs_t__3 = new TH2F("HALU_He_vs_t__3","",100,0.02,0.65,100,-0.5,0.72);
   HALU_He_vs_t__3->SetStats(0);

   HALU_He_vs_t__3->GetYaxis()->SetTitle("A_{LU}^{^{4}He} (90^{#circ})");
   HALU_He_vs_t__3->GetXaxis()->SetTitle("-t [GeV^{2}]");
   HALU_He_vs_t__3->GetXaxis()->CenterTitle(true);
   HALU_He_vs_t__3->GetXaxis()->SetNdivisions(605);
   HALU_He_vs_t__3->GetXaxis()->SetLabelFont(22);
   HALU_He_vs_t__3->GetXaxis()->SetLabelSize(0.05);
   HALU_He_vs_t__3->GetXaxis()->SetTitleSize(0.07);
   HALU_He_vs_t__3->GetXaxis()->SetTitleOffset(0.9);
   HALU_He_vs_t__3->GetXaxis()->SetTitleFont(22);
   HALU_He_vs_t__3->GetYaxis()->CenterTitle(true);
   HALU_He_vs_t__3->GetYaxis()->SetNdivisions(605);
   HALU_He_vs_t__3->GetYaxis()->SetLabelFont(22);
   HALU_He_vs_t__3->GetYaxis()->SetLabelSize(0.05);
   HALU_He_vs_t__3->GetYaxis()->SetTitleSize(0.07);
   HALU_He_vs_t__3->GetYaxis()->SetTitleFont(22);
   HALU_He_vs_t__3->GetZaxis()->SetLabelFont(22);
   HALU_He_vs_t__3->GetZaxis()->SetLabelSize(0.035);
   HALU_He_vs_t__3->GetZaxis()->SetTitleSize(0.035);
   HALU_He_vs_t__3->GetZaxis()->SetTitleFont(22);
   HALU_He_vs_t__3->Draw("");
   
   Double_t HERMES_fx1036[2] = {
   0.027,
   0.105};
   Double_t HERMES_fy1036[2] = {
   0.249,
   0.28};
   Double_t HERMES_fex1036[2] = {
   0,
   0};
   Double_t HERMES_fey1036[2] = {
   0.118,
   0.146};
   gre = new TGraphErrors(2,HERMES_fx1036,HERMES_fy1036,HERMES_fex1036,HERMES_fey1036);
   gre->SetName("HERMES");
   gre->SetTitle("HERMES");
   gre->SetFillColor(1);
   gre->SetLineColor(8);
   gre->SetLineWidth(3);
   gre->SetMarkerColor(8);
   gre->SetMarkerStyle(8);
   gre->SetMarkerSize(1.5);
   
   TH1F *Graph_HERMES1036 = new TH1F("Graph_HERMES1036","HERMES",100,0.0192,0.1128);
   Graph_HERMES1036->SetMinimum(0.1015);
   Graph_HERMES1036->SetMaximum(0.4555);
   Graph_HERMES1036->SetDirectory(0);
   Graph_HERMES1036->SetStats(0);

   //ci = TColor::GetColor("#000099");
   //Graph_HERMES1036->SetLineColor(ci);
   Graph_HERMES1036->GetXaxis()->SetLabelFont(22);
   Graph_HERMES1036->GetXaxis()->SetLabelSize(0.05);
   Graph_HERMES1036->GetXaxis()->SetTitleSize(0.9);
   Graph_HERMES1036->GetXaxis()->SetTitleFont(22);
   Graph_HERMES1036->GetYaxis()->SetLabelFont(22);
   Graph_HERMES1036->GetYaxis()->SetLabelSize(0.05);
   Graph_HERMES1036->GetYaxis()->SetTitleSize(0.06);
   Graph_HERMES1036->GetYaxis()->SetTitleFont(22);
   Graph_HERMES1036->GetZaxis()->SetLabelFont(22);
   Graph_HERMES1036->GetZaxis()->SetLabelSize(0.035);
   Graph_HERMES1036->GetZaxis()->SetTitleSize(0.035);
   Graph_HERMES1036->GetZaxis()->SetTitleFont(22);
   gre->SetHistogram(Graph_HERMES1036);
   
   gre->Draw("p");
   
   Double_t Graph0_fx1037[3] = {
   0.08,
   0.094,
   0.127};
   Double_t Graph0_fy1037[3] = {
   0.376096,
   0.2452892,
   0.3182845};
   Double_t Graph0_fex1037[3] = {
   0,
   0,
   0};
   Double_t Graph0_fey1037[3] = {
   0.0416296,
   0.07200051,
   0.09482501};
   gre = new TGraphErrors(3,Graph0_fx1037,Graph0_fy1037,Graph0_fex1037,Graph0_fey1037);
   gre->SetName("Graph0");
   gre->SetTitle("This Work (<x_{B}>=0.177, <Q^{2}>= 1.492)");
   gre->SetFillColor(1);
   gre->SetLineWidth(3);
   gre->SetMarkerStyle(21);
   gre->SetMarkerSize(1.5);
   
   TH1F *Graph_Graph01037 = new TH1F("Graph_Graph01037","This Work (<x_{B}>=0.177, <Q^{2}>= 1.492)",100,0.0753,0.1317);
   Graph_Graph01037->SetMinimum(0.148845);
   Graph_Graph01037->SetMaximum(0.4421693);
   Graph_Graph01037->SetDirectory(0);
   Graph_Graph01037->SetStats(0);

   //ci = TColor::GetColor("#000099");
   //Graph_Graph01037->SetLineColor(ci);
   Graph_Graph01037->GetXaxis()->SetLabelFont(22);
   Graph_Graph01037->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph01037->GetXaxis()->SetTitleSize(0.9);
   Graph_Graph01037->GetXaxis()->SetTitleFont(22);
   Graph_Graph01037->GetYaxis()->SetLabelFont(22);
   Graph_Graph01037->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph01037->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph01037->GetYaxis()->SetTitleFont(22);
   Graph_Graph01037->GetZaxis()->SetLabelFont(22);
   Graph_Graph01037->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph01037->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph01037->GetZaxis()->SetTitleFont(22);
   gre->SetHistogram(Graph_Graph01037);
   
   gre->Draw("p");
   
   Double_t _fx25[6] = { 0.06, 0.06, 0.08, 0.094, 0.2, 0.2};
   Double_t _fy25[6] = {
   0,
   0.03269784,
   0.03269784,
   0.03091687,
   0.03456616,
   0};
   graph = new TGraph(6,_fx25,_fy25);
   graph->SetName("");
   graph->SetTitle("");
   graph->SetFillColor(40);
   
   TH1F *Graph_Graph25 = new TH1F("Graph_Graph25","",100,0.046,0.214);
   Graph_Graph25->SetMinimum(0);
   Graph_Graph25->SetMaximum(0.03802277);
   Graph_Graph25->SetDirectory(0);
   Graph_Graph25->SetStats(0);

   Graph_Graph25->GetXaxis()->SetLabelFont(22);
   Graph_Graph25->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph25->GetXaxis()->SetTitleSize(0.9);
   Graph_Graph25->GetXaxis()->SetTitleFont(22);
   Graph_Graph25->GetYaxis()->SetLabelFont(22);
   Graph_Graph25->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph25->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph25->GetYaxis()->SetTitleFont(22);
   Graph_Graph25->GetZaxis()->SetLabelFont(22);
   Graph_Graph25->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph25->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph25->GetZaxis()->SetTitleFont(22);
   graph->SetHistogram(Graph_Graph25);
   
   graph->Draw("f");
   
   Double_t Graph0_fx26[16] = {
   0.0026,
   0.0139,
   0.0336,
   0.0611,
   0.0955,
   0.1355,
   0.1796,
   0.2262,
   0.2738,
   0.3204,
   0.3645,
   0.4045,
   0.4389,
   0.4664,
   0.4861,
   0.4974};
   Double_t Graph0_fy26[16] = {
   0.0997298,
   0.100293,
   0.135753,
   0.198727,
   0.240483,
   0.266514,
   0.281026,
   0.287221,
   0.287302,
   0.282924,
   0.275682,
   0.267646,
   0.259709,
   0.252901,
   0.248061,
   0.245101};
   graph = new TGraph(16,Graph0_fx26,Graph0_fy26);
   graph->SetName("Graph0");
   graph->SetTitle("Graph");
   graph->SetFillColor(1);
   graph->SetLineColor(kBlack);
   graph->SetLineWidth(3);
   
   TH1F *Graph_Graph26 = new TH1F("Graph_Graph26","Graph",100,0,0.54688);
   Graph_Graph26->SetMinimum(0.08097258);
   Graph_Graph26->SetMaximum(0.3060592);
   Graph_Graph26->SetDirectory(0);
   Graph_Graph26->SetStats(0);

   Graph_Graph26->SetLineColor(kBlack);
   Graph_Graph26->GetXaxis()->SetLabelFont(22);
   Graph_Graph26->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph26->GetXaxis()->SetTitleSize(0.9);
   Graph_Graph26->GetXaxis()->SetTitleFont(22);
   Graph_Graph26->GetYaxis()->SetLabelFont(22);
   Graph_Graph26->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph26->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph26->GetYaxis()->SetTitleFont(22);
   Graph_Graph26->GetZaxis()->SetLabelFont(22);
   Graph_Graph26->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph26->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph26->GetZaxis()->SetTitleFont(22);
   graph->SetHistogram(Graph_Graph26);
   
   graph->Draw("csame");
   line = new TLine(0,0,1.0,0);
   line->SetLineStyle(7);
   line->Draw();
   
   leg = new TLegend(0.4,0.8,0.97,0.97);
   leg->SetBorderSize(1);
   leg->SetTextSize(0.035);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   leg->AddEntry("Graph0","CLAS-EG6: x_{B}= 0.17, Q^{2}= 1.5 GeV^{2}","P");
   leg->AddEntry("HERMES","HERMES: x_{B}= 0.07, Q^{2}= 1.8 GeV^{2}","P");
   leg->AddEntry("graph","Off-shell: x_{B}= 0.14, Q^{2}= 1.2 GeV^{2}","L");
   leg->Draw();


  // add the BSA projections here 
 double xB_lims[4] = {0.05, 0.17, 0.23, 0.5};

 for(int i=0; i<1; i++){
   for(int j=0; j<3; j++){
   
      Double_t tt[7];
      Double_t tt_err[7];
      Double_t alu[7];
      Double_t alu_err[7];

      for(int k=0; k<7; k++){
         
         tt[k]     = mean_t[i][j][k];
         tt_err[k] = mean_t_err[i][j][k];
         alu[k]    = alu_t_x[i][j][k];
         alu_err[k]= alu_t_x_err[i][j][k];
      }

      gre = new TGraphErrors(7,tt,alu,tt_err,alu_err);
      gre->SetMarkerStyle(24);
      gre->SetMarkerSize(1.3);
      gre->SetLineWidth(5);
      if(j==0) { gre->SetMarkerColor(kBlack); gre->SetLineColor(kBlack); }
      if(j==1) { gre->SetMarkerColor(kBlue);   gre->SetLineColor(kBlue);   }
      if(j==2) { gre->SetMarkerColor(kRed);  gre->SetLineColor(kRed);  }
      gre->Draw("p");
          
      TLatex *l2= new TLatex(0.023, 0.1 *(2-j)-0.32 ,Form("%.2f< x_{B} <%.2f",xB_lims[j], xB_lims[j+1]));
              l2->SetTextSize(0.04);
              if(j==0) { l2->SetTextColor(kBlack);  }
              if(j==1) { l2->SetTextColor(kBlue);   }
              if(j==2) { l2->SetTextColor(kRed);    }
 
              l2->Draw("same");
 
     }
   }

   Can->SetLogx();
   Can->Modified();
   Can->cd();
   Can->Modified();
   Can->cd();
   Can->SetSelected(Can);

   Can->Print("fig/png/BSA_Coh_90_proj.png");
   Can->Print("fig/pdf/BSA_Coh_90_proj.pdf");
}


