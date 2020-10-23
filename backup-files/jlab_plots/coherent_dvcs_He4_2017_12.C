void coherent_dvcs_He4_2017_12()
{
//=========Macro generated from canvas: c1_n22/c1_n22
//=========  (Tue May 16 11:06:04 2017) by ROOT version6.07/07
   TCanvas *c1_n22 = new TCanvas("c1_n22", "c1_n22",0,0,700,500);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   c1_n22->Range(-0.7352941,-0.2352941,0.7352941,1.235294);
   c1_n22->SetFillColor(0);
   c1_n22->SetBorderMode(0);
   c1_n22->SetBorderSize(0);
   c1_n22->SetLeftMargin(0.16);
   c1_n22->SetRightMargin(0.16);
   c1_n22->SetTopMargin(0.16);
   c1_n22->SetBottomMargin(0.16);
   c1_n22->SetFrameBorderMode(0);
   c1_n22->SetFrameBorderMode(0);
   
   THStack *deltaQ2hists = new THStack();
   deltaQ2hists->SetName("deltaQ2hists");
   deltaQ2hists->SetTitle("#Delta Q^{2}");
   
   TH1F *deltaQ2hists_stack_3 = new TH1F("deltaQ2hists_stack_3","#Delta Q^{2}",100,-0.5,0.5);
   deltaQ2hists_stack_3->SetMinimum(0);
   deltaQ2hists_stack_3->SetMaximum(0);
   deltaQ2hists_stack_3->SetDirectory(0);
   deltaQ2hists_stack_3->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   deltaQ2hists_stack_3->SetLineColor(ci);
   deltaQ2hists_stack_3->GetXaxis()->SetTitle("(Q^{2}_{thrown}-Q^{2}_{recon})/Q^{2}_{thrown}");
   deltaQ2hists_stack_3->GetXaxis()->CenterTitle(true);
   deltaQ2hists_stack_3->GetXaxis()->SetLabelFont(22);
   deltaQ2hists_stack_3->GetXaxis()->SetLabelSize(0.05);
   deltaQ2hists_stack_3->GetXaxis()->SetTitleSize(0.06);
   deltaQ2hists_stack_3->GetXaxis()->SetTitleFont(22);
   deltaQ2hists_stack_3->GetYaxis()->CenterTitle(true);
   deltaQ2hists_stack_3->GetYaxis()->SetLabelFont(22);
   deltaQ2hists_stack_3->GetYaxis()->SetLabelSize(0.05);
   deltaQ2hists_stack_3->GetYaxis()->SetTitleSize(0.06);
   deltaQ2hists_stack_3->GetYaxis()->SetTitleOffset(1.2);
   deltaQ2hists_stack_3->GetYaxis()->SetTitleFont(22);
   deltaQ2hists_stack_3->GetZaxis()->SetLabelFont(22);
   deltaQ2hists_stack_3->GetZaxis()->SetLabelSize(0.03);
   deltaQ2hists_stack_3->GetZaxis()->SetTitleSize(0.035);
   deltaQ2hists_stack_3->GetZaxis()->SetTitleFont(22);
   deltaQ2hists->SetHistogram(deltaQ2hists_stack_3);
   
   
   TH1F *hQ25__6 = new TH1F("hQ25__6","Q^{2}",100,-0.5,0.5);
   hQ25__6->SetStats(0);

   ci = TColor::GetColor("#0000ff");
   hQ25__6->SetLineColor(ci);
   hQ25__6->GetXaxis()->SetLabelFont(22);
   hQ25__6->GetXaxis()->SetLabelSize(0.05);
   hQ25__6->GetXaxis()->SetTitleSize(0.06);
   hQ25__6->GetXaxis()->SetTitleFont(22);
   hQ25__6->GetYaxis()->SetLabelFont(22);
   hQ25__6->GetYaxis()->SetLabelSize(0.05);
   hQ25__6->GetYaxis()->SetTitleSize(0.06);
   hQ25__6->GetYaxis()->SetTitleOffset(1.2);
   hQ25__6->GetYaxis()->SetTitleFont(22);
   hQ25__6->GetZaxis()->SetLabelFont(22);
   hQ25__6->GetZaxis()->SetLabelSize(0.03);
   hQ25__6->GetZaxis()->SetTitleSize(0.035);
   hQ25__6->GetZaxis()->SetTitleFont(22);
   deltaQ2hists->Add(hQ25,"");
   deltaQ2hists->Draw("nostack");
   
   TPaveText *pt = new TPaveText(0.4547989,0.9339831,0.5452011,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   TText *AText = pt->AddText("#Delta Q^{2}");
   pt->Draw();
   c1_n22->Modified();
   c1_n22->cd();
   c1_n22->SetSelected(c1_n22);
}
