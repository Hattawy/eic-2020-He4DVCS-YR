void coherent_dvcs_He4_2017_10()
{
//=========Macro generated from canvas: c1_n20/c1_n20
//=========  (Tue May 16 11:06:04 2017) by ROOT version6.07/07
   TCanvas *c1_n20 = new TCanvas("c1_n20", "c1_n20",0,0,700,500);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   c1_n20->Range(-0.7352941,-0.2352941,0.7352941,1.235294);
   c1_n20->SetFillColor(0);
   c1_n20->SetBorderMode(0);
   c1_n20->SetBorderSize(0);
   c1_n20->SetLeftMargin(0.16);
   c1_n20->SetRightMargin(0.16);
   c1_n20->SetTopMargin(0.16);
   c1_n20->SetBottomMargin(0.16);
   c1_n20->SetFrameBorderMode(0);
   c1_n20->SetFrameBorderMode(0);
   
   THStack *deltaXhists = new THStack();
   deltaXhists->SetName("deltaXhists");
   deltaXhists->SetTitle("#Delta x");
   
   TH1F *deltaXhists_stack_1 = new TH1F("deltaXhists_stack_1","#Delta x",100,-0.5,0.5);
   deltaXhists_stack_1->SetMinimum(0);
   deltaXhists_stack_1->SetMaximum(0);
   deltaXhists_stack_1->SetDirectory(0);
   deltaXhists_stack_1->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   deltaXhists_stack_1->SetLineColor(ci);
   deltaXhists_stack_1->GetXaxis()->SetTitle("(x_{thrown}-x_{recon})/x_{thrown}");
   deltaXhists_stack_1->GetXaxis()->CenterTitle(true);
   deltaXhists_stack_1->GetXaxis()->SetLabelFont(22);
   deltaXhists_stack_1->GetXaxis()->SetLabelSize(0.05);
   deltaXhists_stack_1->GetXaxis()->SetTitleSize(0.06);
   deltaXhists_stack_1->GetXaxis()->SetTitleFont(22);
   deltaXhists_stack_1->GetYaxis()->CenterTitle(true);
   deltaXhists_stack_1->GetYaxis()->SetLabelFont(22);
   deltaXhists_stack_1->GetYaxis()->SetLabelSize(0.05);
   deltaXhists_stack_1->GetYaxis()->SetTitleSize(0.06);
   deltaXhists_stack_1->GetYaxis()->SetTitleOffset(1.2);
   deltaXhists_stack_1->GetYaxis()->SetTitleFont(22);
   deltaXhists_stack_1->GetZaxis()->SetLabelFont(22);
   deltaXhists_stack_1->GetZaxis()->SetLabelSize(0.03);
   deltaXhists_stack_1->GetZaxis()->SetTitleSize(0.035);
   deltaXhists_stack_1->GetZaxis()->SetTitleFont(22);
   deltaXhists->SetHistogram(deltaXhists_stack_1);
   
   
   TH1F *hx5__4 = new TH1F("hx5__4","x",100,-0.5,0.5);
   hx5__4->SetStats(0);

   ci = TColor::GetColor("#0000ff");
   hx5__4->SetLineColor(ci);
   hx5__4->SetLineWidth(2);
   hx5__4->GetXaxis()->SetLabelFont(22);
   hx5__4->GetXaxis()->SetLabelSize(0.05);
   hx5__4->GetXaxis()->SetTitleSize(0.06);
   hx5__4->GetXaxis()->SetTitleFont(22);
   hx5__4->GetYaxis()->SetLabelFont(22);
   hx5__4->GetYaxis()->SetLabelSize(0.05);
   hx5__4->GetYaxis()->SetTitleSize(0.06);
   hx5__4->GetYaxis()->SetTitleOffset(1.2);
   hx5__4->GetYaxis()->SetTitleFont(22);
   hx5__4->GetZaxis()->SetLabelFont(22);
   hx5__4->GetZaxis()->SetLabelSize(0.03);
   hx5__4->GetZaxis()->SetTitleSize(0.035);
   hx5__4->GetZaxis()->SetTitleFont(22);
   deltaXhists->Add(hx5,"");
   deltaXhists->Draw("nostack");
   
   TPaveText *pt = new TPaveText(0.4648563,0.94,0.5351437,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   TText *AText = pt->AddText("#Delta x");
   pt->Draw();
   c1_n20->Modified();
   c1_n20->cd();
   c1_n20->SetSelected(c1_n20);
}
