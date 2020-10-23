void coherent_dvcs_He4_2017_11()
{
//=========Macro generated from canvas: c1_n21/c1_n21
//=========  (Tue May 16 11:06:04 2017) by ROOT version6.07/07
   TCanvas *c1_n21 = new TCanvas("c1_n21", "c1_n21",0,0,700,500);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   c1_n21->Range(-1.470588,-0.2352941,1.470588,1.235294);
   c1_n21->SetFillColor(0);
   c1_n21->SetBorderMode(0);
   c1_n21->SetBorderSize(0);
   c1_n21->SetLeftMargin(0.16);
   c1_n21->SetRightMargin(0.16);
   c1_n21->SetTopMargin(0.16);
   c1_n21->SetBottomMargin(0.16);
   c1_n21->SetFrameBorderMode(0);
   c1_n21->SetFrameBorderMode(0);
   
   THStack *deltathists = new THStack();
   deltathists->SetName("deltathists");
   deltathists->SetTitle("#Delta t");
   
   TH1F *deltathists_stack_2 = new TH1F("deltathists_stack_2","#Delta t",100,-1,1);
   deltathists_stack_2->SetMinimum(0);
   deltathists_stack_2->SetMaximum(0);
   deltathists_stack_2->SetDirectory(0);
   deltathists_stack_2->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   deltathists_stack_2->SetLineColor(ci);
   deltathists_stack_2->GetXaxis()->SetTitle("(t_{thrown}-t_{recon})/t_{thrown}");
   deltathists_stack_2->GetXaxis()->CenterTitle(true);
   deltathists_stack_2->GetXaxis()->SetLabelFont(22);
   deltathists_stack_2->GetXaxis()->SetLabelSize(0.05);
   deltathists_stack_2->GetXaxis()->SetTitleSize(0.06);
   deltathists_stack_2->GetXaxis()->SetTitleFont(22);
   deltathists_stack_2->GetYaxis()->CenterTitle(true);
   deltathists_stack_2->GetYaxis()->SetLabelFont(22);
   deltathists_stack_2->GetYaxis()->SetLabelSize(0.05);
   deltathists_stack_2->GetYaxis()->SetTitleSize(0.06);
   deltathists_stack_2->GetYaxis()->SetTitleOffset(1.2);
   deltathists_stack_2->GetYaxis()->SetTitleFont(22);
   deltathists_stack_2->GetZaxis()->SetLabelFont(22);
   deltathists_stack_2->GetZaxis()->SetLabelSize(0.03);
   deltathists_stack_2->GetZaxis()->SetTitleSize(0.035);
   deltathists_stack_2->GetZaxis()->SetTitleFont(22);
   deltathists->SetHistogram(deltathists_stack_2);
   
   
   TH1F *ht5__5 = new TH1F("ht5__5","t",100,-1,1);
   ht5__5->SetStats(0);

   ci = TColor::GetColor("#0000ff");
   ht5__5->SetLineColor(ci);
   ht5__5->SetLineWidth(2);
   ht5__5->GetXaxis()->SetLabelFont(22);
   ht5__5->GetXaxis()->SetLabelSize(0.05);
   ht5__5->GetXaxis()->SetTitleSize(0.06);
   ht5__5->GetXaxis()->SetTitleFont(22);
   ht5__5->GetYaxis()->SetLabelFont(22);
   ht5__5->GetYaxis()->SetLabelSize(0.05);
   ht5__5->GetYaxis()->SetTitleSize(0.06);
   ht5__5->GetYaxis()->SetTitleOffset(1.2);
   ht5__5->GetYaxis()->SetTitleFont(22);
   ht5__5->GetZaxis()->SetLabelFont(22);
   ht5__5->GetZaxis()->SetLabelSize(0.03);
   ht5__5->GetZaxis()->SetTitleSize(0.035);
   ht5__5->GetZaxis()->SetTitleFont(22);
   deltathists->Add(ht5,"");
   deltathists->Draw("nostack");
   
   TPaveText *pt = new TPaveText(0.4684483,0.94,0.5315517,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   TText *AText = pt->AddText("#Delta t");
   pt->Draw();
   c1_n21->Modified();
   c1_n21->cd();
   c1_n21->SetSelected(c1_n21);
}
