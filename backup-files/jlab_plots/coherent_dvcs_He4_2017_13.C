void coherent_dvcs_He4_2017_13()
{
//=========Macro generated from canvas: c1_n23/c1_n23
//=========  (Tue May 16 11:06:04 2017) by ROOT version6.07/07
   TCanvas *c1_n23 = new TCanvas("c1_n23", "c1_n23",0,0,700,500);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   c1_n23->Range(-0.02941176,-0.2470588,0.02941176,1.297059);
   c1_n23->SetFillColor(0);
   c1_n23->SetBorderMode(0);
   c1_n23->SetBorderSize(0);
   c1_n23->SetLeftMargin(0.16);
   c1_n23->SetRightMargin(0.16);
   c1_n23->SetTopMargin(0.16);
   c1_n23->SetBottomMargin(0.16);
   c1_n23->SetFrameBorderMode(0);
   c1_n23->SetFrameBorderMode(0);
   
   TH1F *hW2__7 = new TH1F("hW2__7","W",100,-0.02,0.02);
   hW2__7->SetStats(0);
   hW2__7->SetLineWidth(2);
   hW2__7->GetXaxis()->SetTitle("(W_{thrown}-W_{recon})/W_{thrown}");
   hW2__7->GetXaxis()->CenterTitle(true);
   hW2__7->GetXaxis()->SetLabelFont(22);
   hW2__7->GetXaxis()->SetLabelSize(0.05);
   hW2__7->GetXaxis()->SetTitleSize(0.06);
   hW2__7->GetXaxis()->SetTitleFont(22);
   hW2__7->GetYaxis()->CenterTitle(true);
   hW2__7->GetYaxis()->SetLabelFont(22);
   hW2__7->GetYaxis()->SetLabelSize(0.05);
   hW2__7->GetYaxis()->SetTitleSize(0.06);
   hW2__7->GetYaxis()->SetTitleOffset(1.2);
   hW2__7->GetYaxis()->SetTitleFont(22);
   hW2__7->GetZaxis()->SetLabelFont(22);
   hW2__7->GetZaxis()->SetLabelSize(0.03);
   hW2__7->GetZaxis()->SetTitleSize(0.035);
   hW2__7->GetZaxis()->SetTitleFont(22);
   hW2__7->Draw("");
   
   TPaveText *pt = new TPaveText(0.473477,0.94,0.526523,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   TText *AText = pt->AddText("W");
   pt->Draw();
   c1_n23->Modified();
   c1_n23->cd();
   c1_n23->SetSelected(c1_n23);
}
