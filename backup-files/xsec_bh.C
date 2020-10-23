void xsec_bh()
{
//=========Macro generated from canvas: c55/
//=========  (Sat May 13 18:02:13 2017) by ROOT version6.04/02
   TCanvas *c55 = new TCanvas("c55", "",0,0,750,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   c55->Range(-0.4705882,-3702.917,2.470588,19440.32);
   c55->SetFillColor(0);
   c55->SetBorderMode(0);
   c55->SetBorderSize(0);
   c55->SetLeftMargin(0.16);
   c55->SetRightMargin(0.16);
   c55->SetTopMargin(0.16);
   c55->SetBottomMargin(0.16);
   c55->SetFrameBorderMode(0);
   c55->SetFrameBorderMode(0);
   
   TH1D *h_xsec__1 = new TH1D("h_xsec__1","BH cross section",150,0,2);
   h_xsec__1->SetBinContent(1,14988);
   h_xsec__1->SetBinContent(2,10494);
   h_xsec__1->SetBinContent(3,9086);
   h_xsec__1->SetBinContent(4,7990);
   h_xsec__1->SetBinContent(5,7476);
   h_xsec__1->SetBinContent(6,7007);
   h_xsec__1->SetBinContent(7,6399);
   h_xsec__1->SetBinContent(8,5988);
   h_xsec__1->SetBinContent(9,5817);
   h_xsec__1->SetBinContent(10,5374);
   h_xsec__1->SetBinContent(11,5298);
   h_xsec__1->SetBinContent(12,5239);
   h_xsec__1->SetBinContent(13,4907);
   h_xsec__1->SetBinContent(14,4648);
   h_xsec__1->SetBinContent(15,4549);
   h_xsec__1->SetBinContent(16,4403);
   h_xsec__1->SetBinContent(17,4248);
   h_xsec__1->SetBinContent(18,4149);
   h_xsec__1->SetBinContent(19,4036);
   h_xsec__1->SetBinContent(20,3922);
   h_xsec__1->SetBinContent(21,3803);
   h_xsec__1->SetBinContent(22,3782);
   h_xsec__1->SetBinContent(23,3528);
   h_xsec__1->SetBinContent(24,3649);
   h_xsec__1->SetBinContent(25,3381);
   h_xsec__1->SetBinContent(26,3368);
   h_xsec__1->SetBinContent(27,3342);
   h_xsec__1->SetBinContent(28,3054);
   h_xsec__1->SetBinContent(29,3156);
   h_xsec__1->SetBinContent(30,2965);
   h_xsec__1->SetBinContent(31,3023);
   h_xsec__1->SetBinContent(32,2905);
   h_xsec__1->SetBinContent(33,2877);
   h_xsec__1->SetBinContent(34,2756);
   h_xsec__1->SetBinContent(35,2787);
   h_xsec__1->SetBinContent(36,2773);
   h_xsec__1->SetBinContent(37,2688);
   h_xsec__1->SetBinContent(38,2608);
   h_xsec__1->SetBinContent(39,2604);
   h_xsec__1->SetBinContent(40,2622);
   h_xsec__1->SetBinContent(41,2520);
   h_xsec__1->SetBinContent(42,2499);
   h_xsec__1->SetBinContent(43,2458);
   h_xsec__1->SetBinContent(44,2436);
   h_xsec__1->SetBinContent(45,2398);
   h_xsec__1->SetEntries(200000);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   h_xsec__1->SetLineColor(ci);
   h_xsec__1->GetXaxis()->SetLabelFont(22);
   h_xsec__1->GetXaxis()->SetLabelSize(0.05);
   h_xsec__1->GetXaxis()->SetTitleSize(0.06);
   h_xsec__1->GetXaxis()->SetTitleFont(22);
   h_xsec__1->GetYaxis()->SetLabelFont(22);
   h_xsec__1->GetYaxis()->SetLabelSize(0.05);
   h_xsec__1->GetYaxis()->SetTitleSize(0.06);
   h_xsec__1->GetYaxis()->SetTitleFont(22);
   h_xsec__1->GetZaxis()->SetLabelFont(22);
   h_xsec__1->GetZaxis()->SetLabelSize(0.03);
   h_xsec__1->GetZaxis()->SetTitleSize(0.035);
   h_xsec__1->GetZaxis()->SetTitleFont(22);
   h_xsec__1->Draw("");
   
   TPaveText *pt = new TPaveText(0.3318231,0.94,0.6681769,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   TText *AText = pt->AddText("BH cross section");
   pt->Draw();
   c55->Modified();
   c55->cd();
   c55->SetSelected(c55);
}
