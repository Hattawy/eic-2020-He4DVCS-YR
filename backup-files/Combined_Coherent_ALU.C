#include "Riostream.h"
#include <iostream>
#include <fstream>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <math.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TFile.h>
#include <time.h>
#include <TMultiGraph.h>
#include <TVector3.h>
#include <TMath.h>
#include <TProfile.h>
#include "TCutG.h"
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include "TGraph.h"
#include <TLine.h>
#include <TLorentzVector.h>
#include "TSystem.h"
#include "TColor.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TText.h"
#include "TPaveText.h"
#include "calculate_CFF.C"

void Combined_Coherent_ALU()
{

  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);
  gStyle->SetLabelSize(0.035, "xyz"); // size of axis value font 
  gStyle->SetTitleSize(0.035, "xyz"); // size of axis title font 
  gStyle->SetTitleFont(22, "xyz"); // font option 
  gStyle->SetLabelFont(22, "xyz");
  gStyle->SetPadBottomMargin(0.18);  
  gStyle->SetPadTopMargin(0.02);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetPaperSize(20, 24);
  gStyle->SetLabelSize(0.05, "xy");
  gStyle->SetTitleSize(0.06, "xy");
  gStyle->SetTitleSize(0.9);
 
  double xfont = 0.14;
  double lfont = 0.08;
  
  // for ALU vs phi plots
   TCanvas *c12 = new TCanvas("c12",  "", 0, 0, 1000, 900);
            c12->Divide(3,3,0,0);

   // for the CFF
   TCanvas *cc = new TCanvas("cc", "",2,51,1050,850);
            cc->Divide(3,2,3.0,0);

   TCanvas *ccc = new TCanvas("ccc", "",0,0,1250,500);
            ccc->Divide(2,1);

   // For ALU at phi = 90 deg
   TCanvas *Can = new TCanvas("Can", "Helium DVCS",0,0,1200,400);   
            Can->Divide(3,1,0,0);

   // For ALU ratios at phi = 90 deg
   TCanvas *Can_r = new TCanvas("Can_r", "Helium DVCS",0,0,1200,450);   
            Can_r->Divide(3,1); 



   TCanvas *Cc = new TCanvas("Cc", "Helium DVCS",0,0,650,500);   

////////////////////////////////////////////////////////////////////////// 
//________________                                  ____________________//
//________________      bins in Q2                  ____________________//                          
//________________                                  ____________________//

  double vQ2[3] = {1.143,  1.423,  1.902};               double vQ2_er[3] = {0.0,  0.0,  0.0};
  double vxB[3] = {0.136,  0.172,  0.224};               double vxB_er[3] = {0.0,  0.0,  0.0};
  double vt[3]  = {-0.096,  -0.099,  -0.107};            double vt_er[3]  = {0.0,  0.0,  0.0};

  double c0_BH[3],  c1_BH[3],  c2_BH[3]; 
  double A0[3],  A1[3],  A2[3],  A3[3];


  calculate_CFF(vQ2[0],  vxB[0],  vt[0],  A0[0],  A1[0],  A2[0],  A3[0],  c0_BH[0],  c1_BH[0],  c2_BH[0]);  
  calculate_CFF(vQ2[1],  vxB[1],  vt[1],  A0[1],  A1[1],  A2[1],  A3[1],  c0_BH[1],  c1_BH[1],  c2_BH[1]);  
  calculate_CFF(vQ2[2],  vxB[2],  vt[2],  A0[2],  A1[2],  A2[2],  A3[2],  c0_BH[2],  c1_BH[2],  c2_BH[2]);

   TF1 *ffit[3];
   TFitResultPtr r[3];
   double ALU_fit[3];
   double ALU_fit_err[3];

   for(int i =0; i<3; i++) {  
      ffit[i] = new TF1(Form("ffit[%u]", i),  Form("%0.7f*[0]*sin(x*3.14/180.0) / (%0.7f+ %0.7f*cos(x*3.14/180.0) + %0.7f*cos(2*x*3.14/180.0) + %0.7f*([0]*[0] + [1]*[1]) + %0.7f*[1] + %0.7f*[1]*cos(x*3.14/180.0))",  A0[i],  c0_BH[i],  c1_BH[i],  c2_BH[i],  A1[i],  A2[i],  A3[i]), 0.0, 360.0);
      ffit[i]->SetLineColor(kRed);
      ffit[i]->SetParName(0,"Im(H_{A})");
      ffit[i]->SetParName(1,"Re(H_{A})");
     } 

   TLatex *mg_title;
   double xax = 15.0; 
   double yay = -0.57;
            
// ------------>Primitives in pad: c12_1

   c12->cd(1);
   TMultiGraph *multigraph = new TMultiGraph();
   mg_title = new TLatex(xax,yay,"Q^{2}= 1.143 GeV^{2}");
   mg_title->SetTextSize(0.085); 
   mg_title->SetTextFont(62); 
   mg_title->SetTextColor(kBlack);

   TGraph *gre_sys = new TGraph(13);
   gre_sys->SetFillColor(40);
   gre_sys->SetLineColor(40);
   gre_sys->SetMarkerColor(40);
   gre_sys->SetPoint(0, 0,    0.0        );
   gre_sys->SetPoint(1, 0,    0.02571563 );
   gre_sys->SetPoint(2, 20,   0.02571563 );
   gre_sys->SetPoint(3, 60,   0.01885969 );
   gre_sys->SetPoint(4, 100,  0.04038886 );
   gre_sys->SetPoint(5, 140,  0.02678967 );
   gre_sys->SetPoint(6, 180,  0.02789155 );
   gre_sys->SetPoint(7, 220,  0.02484531 );
   gre_sys->SetPoint(8, 260,  0.03891253 );
   gre_sys->SetPoint(9, 300,  0.02632836 );
   gre_sys->SetPoint(10, 340, 0.03431105 );
   gre_sys->SetPoint(11, 360, 0.03431105 );
   gre_sys->SetPoint(12, 360, 0.0        );

   TGraphErrors *gre_st = new TGraphErrors(9);
   gre_st->SetFillColor(8);
   gre_st->SetMarkerColor(kBlack);   
   gre_st->SetMarkerSize(1.3);   
   gre_st->SetLineWidth(2);
   gre_st->SetMarkerStyle(21);
   gre_st->SetPoint(0, 24.3146  ,   0.1327475  );    gre_st->SetPointError(0, 0, sqrt(pow(  0.1088386 , 2) - pow( 0.02571563, 2) ));
   gre_st->SetPoint(1, 61.16264 ,   0.3210551  );    gre_st->SetPointError(1, 0, sqrt(pow(  0.09325182, 2) - pow( 0.01885969, 2) ));
   gre_st->SetPoint(2, 98.9465  ,   0.3709448  );    gre_st->SetPointError(2, 0, sqrt(pow(  0.1026357 , 2) - pow( 0.04038886, 2) ));
   gre_st->SetPoint(3, 140.9271 ,   0.2452102  );    gre_st->SetPointError(3, 0, sqrt(pow(  0.1517177 , 2) - pow( 0.02678967, 2) ));
   gre_st->SetPoint(4, 178.4282 ,   0.02256708 );    gre_st->SetPointError(4, 0, sqrt(pow(  0.1632185 , 2) - pow( 0.02789155, 2) ));
   gre_st->SetPoint(5, 219.1438 ,   -0.0526485 );    gre_st->SetPointError(5, 0, sqrt(pow(  0.1475768 , 2) - pow( 0.02484531, 2) ));
   gre_st->SetPoint(6, 262.95   ,   -0.2642824 );    gre_st->SetPointError(6, 0, sqrt(pow(  0.1196945 , 2) - pow( 0.03891253, 2) ));
   gre_st->SetPoint(7, 301.9212 ,   -0.1761262 );    gre_st->SetPointError(7, 0, sqrt(pow(  0.09688492, 2) - pow( 0.02632836, 2) ));
   gre_st->SetPoint(8, 337.7335 ,   -0.2791955 );    gre_st->SetPointError(8, 0, sqrt(pow(  0.1047699 , 2) - pow( 0.03431105, 2) ));
   multigraph->Add(gre_st, "");

   TGraphErrors *gre = new TGraphErrors(9);
   gre->SetFillColor(kBlue);
   gre->SetMarkerColor(kBlue);
   gre->SetMarkerStyle(21);
   gre->SetPoint(0, 24.3146  ,   0.1327475  );    gre->SetPointError(0, 0,  0.1088386  );
   gre->SetPoint(1, 61.16264 ,   0.3210551  );    gre->SetPointError(1, 0,  0.09325182 );
   gre->SetPoint(2, 98.9465  ,   0.3709448  );    gre->SetPointError(2, 0,  0.1026357  );
   gre->SetPoint(3, 140.9271 ,   0.2452102  );    gre->SetPointError(3, 0,  0.1517177  );
   gre->SetPoint(4, 178.4282 ,   0.02256708 );    gre->SetPointError(4, 0,  0.1632185  );
   gre->SetPoint(5, 219.1438 ,   -0.0526485 );    gre->SetPointError(5, 0,  0.1475768  );
   gre->SetPoint(6, 262.95   ,   -0.2642824 );    gre->SetPointError(6, 0,  0.1196945  );
   gre->SetPoint(7, 301.9212 ,   -0.1761262 );    gre->SetPointError(7, 0,  0.09688492 );
   gre->SetPoint(8, 337.7335 ,   -0.2791955 );    gre->SetPointError(8, 0,  0.1047699  );

   //multigraph->Add(gre, "");
   gre->Fit(ffit[0]);
   r[0] = gre->Fit(ffit[0],"S");

   multigraph->Draw("a2"); multigraph->Draw("p");
   multigraph->GetXaxis()->SetTitle("#phi [deg.]");
   multigraph->GetYaxis()->SetTitle("A_{LU}^{^{4}He}");
   multigraph->GetXaxis()->SetLabelSize(lfont);
   multigraph->GetYaxis()->SetLabelSize(lfont);
   multigraph->GetXaxis()->SetTitleSize(xfont);
   multigraph->GetYaxis()->SetTitleSize(xfont);
   multigraph->GetYaxis()->CenterTitle(true);
   multigraph->GetXaxis()->CenterTitle(true);
   multigraph->GetYaxis()->SetTitleOffset(0.8);
   gPad->Modified();  
   multigraph->GetYaxis()->SetRangeUser(-0.65,  0.65);
   multigraph->GetXaxis()->SetLimits(-5.0,  375);
   gre_sys->Draw("fsame");
   gre_st->Draw("psame");  
   ffit[0]->Draw("same");
   mg_title->Draw("same");

// ------------>Primitives in pad: c12_2
   c12->cd(2);   
   multigraph = new TMultiGraph();

   mg_title = new TLatex(xax,yay,"Q^{2}= 1.423 GeV^{2}");
   mg_title->SetTextSize(0.085); 
   mg_title->SetTextFont(62); 
   mg_title->SetTextColor(kBlack);

   gre_sys = new TGraph(13);
   gre_sys->SetFillColor(40);
   gre_sys->SetLineColor(40);
   gre_sys->SetMarkerColor(40);
   gre_sys->SetPoint(0, 0,    0.0        );
   gre_sys->SetPoint(1, 0,    0.02673700 );
   gre_sys->SetPoint(2, 20,   0.02673700 );
   gre_sys->SetPoint(3, 60,   0.02503328 );
   gre_sys->SetPoint(4, 100,  0.04267748 );
   gre_sys->SetPoint(5, 140,  0.02500644 );
   gre_sys->SetPoint(6, 180,  0.03026008 );
   gre_sys->SetPoint(7, 220,  0.03447198 );
   gre_sys->SetPoint(8, 260,  0.03663218 );
   gre_sys->SetPoint(9, 300,  0.03835156 );
   gre_sys->SetPoint(10, 340, 0.02613428 );
   gre_sys->SetPoint(11, 360, 0.02613428 );
   gre_sys->SetPoint(12, 360, 0.0        );
 
   gre_st = new TGraphErrors(9);
   gre_st->SetFillColor(8);
   gre_st->SetMarkerColor(kBlack); 
   gre_st->SetMarkerSize(1.3);
   gre_st->SetLineWidth(2);
   gre_st->SetMarkerStyle(21);
   gre_st->SetPoint(0,   21.10117 ,   0.1917058 );   gre_st->SetPointError(0, 0, sqrt(pow(  0.0892737 , 2) - pow( 0.02673700, 2) ));
   gre_st->SetPoint(1,   56.64015 ,   0.2819707 );   gre_st->SetPointError(1, 0, sqrt(pow(  0.08665254, 2) - pow( 0.02503328, 2) ));
   gre_st->SetPoint(2,   96.9922  ,   0.485539  );   gre_st->SetPointError(2, 0, sqrt(pow(  0.1291767 , 2) - pow( 0.04267748, 2) ));
   gre_st->SetPoint(3,   139.9828 ,   0.09969585);   gre_st->SetPointError(3, 0, sqrt(pow(  0.1679711 , 2) - pow( 0.02500644, 2) ));
   gre_st->SetPoint(4,   180.3615 ,   0.1458599 );   gre_st->SetPointError(4, 0, sqrt(pow(  0.1906658 , 2) - pow( 0.03026008, 2) ));
   gre_st->SetPoint(5,   218.7404 ,   -0.1109757);   gre_st->SetPointError(5, 0, sqrt(pow(  0.1853831 , 2) - pow( 0.03447198, 2) ));
   gre_st->SetPoint(6,   263.194  ,   -0.3516181);   gre_st->SetPointError(6, 0, sqrt(pow(  0.1371968 , 2) - pow( 0.03663218, 2) ));
   gre_st->SetPoint(7,   302.2328 ,   -0.4141039);   gre_st->SetPointError(7, 0, sqrt(pow(  0.08421253, 2) - pow( 0.03835156, 2) ));
   gre_st->SetPoint(8,   338.2142 ,   -0.2786012);   gre_st->SetPointError(8, 0, sqrt(pow(  0.08365383, 2) - pow( 0.02613428, 2) ));
   multigraph->Add(gre_st, "");

   gre = new TGraphErrors(9);
   gre->SetFillColor(1);
   gre->SetMarkerColor(kBlue);
   gre->SetMarkerStyle(21);
   gre->SetPoint(0,   21.10117 ,   0.1917058 );   gre->SetPointError(0, 0,  0.0892737 );
   gre->SetPoint(1,   56.64015 ,   0.2819707 );   gre->SetPointError(1, 0,  0.08665254);
   gre->SetPoint(2,   96.9922  ,   0.485539  );   gre->SetPointError(2, 0,  0.1291767 );
   gre->SetPoint(3,   139.9828 ,   0.09969585);   gre->SetPointError(3, 0,  0.1679711 );
   gre->SetPoint(4,   180.3615 ,   0.1458599 );   gre->SetPointError(4, 0,  0.1906658 );
   gre->SetPoint(5,   218.7404 ,   -0.1109757);   gre->SetPointError(5, 0,  0.1853831 );
   gre->SetPoint(6,   263.194  ,   -0.3516181);   gre->SetPointError(6, 0,  0.1371968 );
   gre->SetPoint(7,   302.2328 ,   -0.4141039);   gre->SetPointError(7, 0,  0.08421253);
   gre->SetPoint(8,   338.2142 ,   -0.2786012);   gre->SetPointError(8, 0,  0.08365383);
 
   
   //multigraph->Add(gre, "");
   gre->Fit(ffit[1]);
   r[1] = gre->Fit(ffit[1],"S");

   multigraph->Draw("a2"); multigraph->Draw("p");
   
   multigraph->GetXaxis()->SetTitle("#phi [deg.]");
   multigraph->GetXaxis()->CenterTitle(true);
   multigraph->GetXaxis()->SetLabelFont(22);
   multigraph->GetXaxis()->SetLabelSize(xfont);
   multigraph->GetXaxis()->SetTitleSize(xfont);
   multigraph->GetXaxis()->SetTitleFont(22);
   multigraph->GetYaxis()->SetLabelFont(22);
   multigraph->GetYaxis()->SetLabelSize(xfont);
   multigraph->GetYaxis()->SetTitleSize(xfont);
   multigraph->GetYaxis()->SetTitleOffset(1.2);
   multigraph->GetYaxis()->SetTitleFont(22);
   gPad->Modified();  
   multigraph->GetYaxis()->SetRangeUser(-0.65,  0.65);
   multigraph->GetXaxis()->SetLimits(-5.0,  375);
   gre_sys->Draw("fsame");
   gre_st->Draw("psame");  ffit[1]->Draw("same");
   mg_title->Draw("same");
  
// ------------>Primitives in pad: c12_3
   c12->cd(3);
   multigraph = new TMultiGraph();

   mg_title = new TLatex(xax,yay,"Q^{2}= 1.902 GeV^{2}");
   mg_title->SetTextSize(0.085); 
   mg_title->SetTextFont(62); 
   mg_title->SetTextColor(kBlack);

   gre_sys = new TGraph(13);
   gre_sys->SetFillColor(40);
   gre_sys->SetLineColor(40);
   gre_sys->SetMarkerColor(40);
   gre_sys->SetPoint(0, 0,    0.0        );
   gre_sys->SetPoint(1, 0,    0.02345626 );
   gre_sys->SetPoint(2, 20,   0.02345626 );
   gre_sys->SetPoint(3, 60,   0.01851792 );
   gre_sys->SetPoint(4, 100,  0.0173685  );
   gre_sys->SetPoint(5, 140,  0.01702304 );
   gre_sys->SetPoint(6, 180,  0.01265815 );
   gre_sys->SetPoint(7, 220,  0.01528429 );
   gre_sys->SetPoint(8, 260,  0.02602474 );
   gre_sys->SetPoint(9, 300,  0.02071202 );
   gre_sys->SetPoint(10, 340, 0.01635832 );
   gre_sys->SetPoint(11, 360, 0.01635832 );
   gre_sys->SetPoint(12, 360, 0.0        );

   gre_st = new TGraphErrors(9);
   gre_st->SetFillColor(8);
   gre_st->SetMarkerColor(kBlack);
   gre_st->SetMarkerSize(1.3); 
   gre_st->SetLineWidth(2);
   gre_st->SetMarkerStyle(21);
   gre_st->SetPoint(0, 21.45518 , 0.1803878 );    gre_st->SetPointError(0, 0, sqrt(pow( 0.08053756, 2) - pow( 0.02345626, 2) ));
   gre_st->SetPoint(1, 57.23776 , 0.3501509 );    gre_st->SetPointError(1, 0, sqrt(pow( 0.08242233, 2) - pow( 0.01851792, 2) ));
   gre_st->SetPoint(2, 96.22869 , 0.2701766 );    gre_st->SetPointError(2, 0, sqrt(pow( 0.1230086 , 2) - pow( 0.0173685 , 2) ));
   gre_st->SetPoint(3, 139.5664 , 0.3053819 );    gre_st->SetPointError(3, 0, sqrt(pow( 0.239036  , 2) - pow( 0.01702304, 2) ));
   gre_st->SetPoint(4, 178.2038 , 0.1028176 );    gre_st->SetPointError(4, 0, sqrt(pow( 0.2670899 , 2) - pow( 0.01265815, 2) ));
   gre_st->SetPoint(5, 221.4048 , -0.2118587);    gre_st->SetPointError(5, 0, sqrt(pow( 0.214785  , 2) - pow( 0.01528429, 2) ));
   gre_st->SetPoint(6, 263.3962 , -0.3060192);    gre_st->SetPointError(6, 0, sqrt(pow( 0.1306386 , 2) - pow( 0.02602474, 2) ));
   gre_st->SetPoint(7, 303.3615 , -0.1384293);    gre_st->SetPointError(7, 0, sqrt(pow( 0.09449148, 2) - pow( 0.02071202, 2) ));
   gre_st->SetPoint(8, 338.5945 , -0.1633788);    gre_st->SetPointError(8, 0, sqrt(pow( 0.07914896, 2) - pow( 0.01635832, 2) ));
   multigraph->Add(gre_st, "");

   gre = new TGraphErrors(9);
   gre->SetFillColor(1);
   gre->SetMarkerColor(kBlue);
   gre->SetMarkerStyle(21);
   gre->SetPoint(0, 21.45518 , 0.1803878  );   gre->SetPointError(0, 0, 0.08053756);
   gre->SetPoint(1, 57.23776 , 0.3501509  );   gre->SetPointError(1, 0, 0.08242233);
   gre->SetPoint(2, 96.22869 , 0.2701766  );   gre->SetPointError(2, 0, 0.1230086 );
   gre->SetPoint(3, 139.5664 , 0.3053819  );   gre->SetPointError(3, 0, 0.239036  );
   gre->SetPoint(4, 178.2038 , 0.1028176  );   gre->SetPointError(4, 0, 0.2670899 );
   gre->SetPoint(5, 221.4048 , -0.2118587 );   gre->SetPointError(5, 0, 0.214785  );
   gre->SetPoint(6, 263.3962 , -0.3060192 );   gre->SetPointError(6, 0, 0.1306386 );
   gre->SetPoint(7, 303.3615 , -0.1384293 );   gre->SetPointError(7, 0, 0.09449148);
   gre->SetPoint(8, 338.5945 , -0.1633788 );   gre->SetPointError(8, 0, 0.07914896);
   
   //multigraph->Add(gre, "");
   gre->Fit(ffit[2]);
   r[2] = gre->Fit(ffit[2],"S");

   multigraph->Draw("a2");   multigraph->Draw("p");
   multigraph->GetXaxis()->SetTitle("#phi [deg.]");
   multigraph->GetXaxis()->CenterTitle(true);
   multigraph->GetXaxis()->SetLabelFont(22);
   multigraph->GetXaxis()->SetLabelSize(0.05);
   multigraph->GetXaxis()->SetTitleSize(0.06);
   multigraph->GetXaxis()->SetTitleFont(22);
   multigraph->GetYaxis()->SetLabelFont(22);
   multigraph->GetYaxis()->SetLabelSize(0.05);
   multigraph->GetYaxis()->SetTitleSize(0.06);
   multigraph->GetYaxis()->SetTitleOffset(1.2);
   multigraph->GetYaxis()->SetTitleFont(22);
   gPad->Modified();  
   multigraph->GetYaxis()->SetRangeUser(-0.65,  0.65);
   multigraph->GetXaxis()->SetLimits(-5.0,  375);
   gre_sys->Draw("fsame");  
   gre_st->Draw("psame");  ffit[2]->Draw("same");
   mg_title->Draw("same");

   //__________________________________________________________________________
   // plot the Imaginary and the real parts of HA CFF vs Q2 -------------------
   cc->cd(1);

   TMultiGraph *mgCFF = new TMultiGraph();

   Double_t graphCFF_fx1001[3] = { 1.143, 1.423, 1.902};
   Double_t graphCFF_fy1001[3] = { ffit[0]->GetParameter(0), ffit[1]->GetParameter(0), ffit[2]->GetParameter(0)};
   Double_t graphCFF_fex1001[3] = { 0, 0, 0};
   Double_t graphCFF_fey1001[3] = { ffit[0]->GetParError(0), ffit[1]->GetParError(0), ffit[2]->GetParError(0)};
 
   TGraphErrors *greCFF = new TGraphErrors(3,graphCFF_fx1001,graphCFF_fy1001,graphCFF_fex1001,graphCFF_fey1001);
                 greCFF->SetFillColor(1);
                 greCFF->SetMarkerColor(kBlack);
                 greCFF->SetLineColor(kBlack); 
                 greCFF->SetMarkerStyle(21);
                 greCFF->SetMarkerSize(1.5);
                 greCFF->SetLineWidth(2);

   mgCFF->Add(greCFF,"");
   mgCFF->Draw("AP");
   mgCFF->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
   mgCFF->GetYaxis()->SetTitle("Im(H_{A})");
   mgCFF->GetXaxis()->SetLabelSize(lfont-0.01);
   mgCFF->GetYaxis()->SetLabelSize(lfont-0.01);
   mgCFF->GetXaxis()->SetTitleSize(xfont-0.012);
   mgCFF->GetYaxis()->SetTitleSize(xfont-0.012);
   mgCFF->GetXaxis()->SetNdivisions(605);
   mgCFF->GetXaxis()->CenterTitle(true);
   mgCFF->GetYaxis()->CenterTitle(true);
   mgCFF->GetXaxis()->SetLimits(0.85,2.6);
   mgCFF->GetYaxis()->SetRangeUser(-2.0, 78.0);

      // from Vadim's calculations
   TGraph * graphCFF = new TGraph(21);
   graphCFF->SetFillColor(1);
   graphCFF->SetLineColor(2);
   graphCFF->SetLineStyle(1);
   graphCFF->SetLineWidth(4);
   graphCFF->SetPoint(0,  1.0000 , 17.5160 );
   graphCFF->SetPoint(1,  1.0800 , 17.7305 );
   graphCFF->SetPoint(2,  1.1600 , 17.9450 );
   graphCFF->SetPoint(3,  1.2400 , 18.1596 );
   graphCFF->SetPoint(4,  1.3200 , 18.3741 );
   graphCFF->SetPoint(5,  1.4000 , 18.5886 );
   graphCFF->SetPoint(6,  1.4800 , 18.8031 );
   graphCFF->SetPoint(7,  1.5600 , 18.4309 );
   graphCFF->SetPoint(8,  1.6400 , 17.8632 );
   graphCFF->SetPoint(9,  1.7200 , 17.2955 );
   graphCFF->SetPoint(10, 1.8000 , 16.7278 );
   graphCFF->SetPoint(11, 1.8800 , 16.1601 );
   graphCFF->SetPoint(12, 1.9600 , 15.5923 );
   graphCFF->SetPoint(13, 2.0400 , 15.3573 );
   graphCFF->SetPoint(14, 2.1200 , 15.4550 );
   graphCFF->SetPoint(15, 2.2000 , 15.5527 );
   graphCFF->SetPoint(16, 2.2800 , 15.6504 );
   graphCFF->SetPoint(17, 2.3600 , 15.7481 );
   graphCFF->SetPoint(18, 2.4400 , 15.8458 );
   graphCFF->SetPoint(19, 2.5200 , 15.9375 );
   graphCFF->SetPoint(20, 2.6000 , 16.0110 );
   graphCFF->Draw("csame");
 


   // the real part of the CFF HA versus Q2
   cc->cd(4);
   mgCFF = new TMultiGraph();
 
   Double_t graphCFF_fx1002[3] = { 1.143, 1.423, 1.902};
   Double_t graphCFF_fy1002[3] = { ffit[0]->GetParameter(1), ffit[1]->GetParameter(1), ffit[2]->GetParameter(1)};
   Double_t graphCFF_fex1002[3] = { 0, 0, 0};
   Double_t graphCFF_fey1002[3] = { ffit[0]->GetParError(1), ffit[1]->GetParError(1), ffit[2]->GetParError(1)};

   greCFF = new TGraphErrors(3,graphCFF_fx1002,graphCFF_fy1002,graphCFF_fex1002,graphCFF_fey1002);
   greCFF->SetFillColor(1);
   greCFF->SetMarkerColor(kBlack);
   greCFF->SetMarkerStyle(25);
   greCFF->SetMarkerSize(1.5);
   greCFF->SetLineColor(kBlack); 
   greCFF->SetLineWidth(2);
   
   mgCFF->Add(greCFF,"");

   mgCFF->Draw("AP");
   mgCFF->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
   mgCFF->GetYaxis()->SetTitle("Re(H_{A})");
   mgCFF->GetXaxis()->SetLabelSize(lfont-0.01);
   mgCFF->GetYaxis()->SetLabelSize(lfont-0.01);
   mgCFF->GetXaxis()->SetTitleSize(xfont-0.017);
   mgCFF->GetYaxis()->SetTitleSize(xfont-0.017);
   mgCFF->GetXaxis()->SetTitleOffset(0.8);
   mgCFF->GetXaxis()->CenterTitle(true);
   mgCFF->GetYaxis()->CenterTitle(true);
   mgCFF->GetXaxis()->SetNdivisions(605);
   mgCFF->GetYaxis()->SetRangeUser(-60.0, 65.0);
   mgCFF->GetXaxis()->SetLimits(0.85,2.6);

   // real part from vadim vs Q2
   graphCFF = new TGraph(20);
   graphCFF->SetFillColor(1);
   graphCFF->SetLineColor(2);
   graphCFF->SetLineStyle(1);
   graphCFF->SetLineWidth(4);
   graphCFF->SetPoint(0,  1.0000 ,   -5.3312 );
   graphCFF->SetPoint(1,  1.0800 ,   -5.1177 );
   graphCFF->SetPoint(2,  1.1600 ,   -4.9042 );
   graphCFF->SetPoint(3,  1.2400 ,   -4.6907 );
   graphCFF->SetPoint(4,  1.3200 ,   -4.4771 );
   graphCFF->SetPoint(5,  1.4000 ,   -4.2636 );
   graphCFF->SetPoint(6,  1.4800 ,   -4.0501 );
   graphCFF->SetPoint(7,  1.5600 ,   -3.9338 );
   graphCFF->SetPoint(8,  1.6400 ,   -3.8500 );
   graphCFF->SetPoint(9,  1.7200 ,   -3.7661 );
   graphCFF->SetPoint(10, 1.8000 ,   -3.6822 );
   graphCFF->SetPoint(11, 1.8800 ,   -3.5984 );
   graphCFF->SetPoint(12, 1.9600 ,   -3.5145 );
   graphCFF->SetPoint(13, 2.0400 ,   -3.4239 );
   graphCFF->SetPoint(14, 2.1200 ,   -3.3266 );
   graphCFF->SetPoint(15, 2.2000 ,   -3.2293 );
   graphCFF->SetPoint(16, 2.2800 ,   -3.1320 );
   graphCFF->SetPoint(17, 2.3600 ,   -3.0347 );
   graphCFF->SetPoint(18, 2.4400 ,   -2.9374 );
   graphCFF->SetPoint(19, 2.5200 ,   -2.8453 );
   graphCFF->SetPoint(20, 2.6000 ,   -2.7688 );
   graphCFF->Draw("csame");

   // plot the coherent ALU at phi = 90 deg vs Q2  -------------------------------

   Can->cd(1);

   TH1 *HALU_He_vs_Q2 = new TH2F("HALU_He_vs_Q2","", 100,0.85, 2.6, 100,-0.1, 0.72);
        HALU_He_vs_Q2->SetStats(0);
        HALU_He_vs_Q2->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
        HALU_He_vs_Q2->GetYaxis()->SetTitle("A_{LU}^{^{4}He}(90#circ)");
        HALU_He_vs_Q2->GetXaxis()->SetTitleSize(xfont-0.06);
        HALU_He_vs_Q2->GetYaxis()->SetTitleSize(xfont-0.06);
        HALU_He_vs_Q2->GetXaxis()->SetNdivisions(205);
        HALU_He_vs_Q2->GetXaxis()->CenterTitle(true);
        HALU_He_vs_Q2->GetYaxis()->CenterTitle(true);
        HALU_He_vs_Q2->GetYaxis()->SetNdivisions(414);
        //HALU_He_vs_Q2->GetXaxis()->SetTitleOffset(1.5);
        //HALU_He_vs_Q2->GetYaxis()->SetTitleOffset(1.5);
        HALU_He_vs_Q2->Draw("");
   
   for(int i=0; i<3; i++){
        double x[1] = { 90.0};  double err[1];
        r[i]->GetConfidenceIntervals(1, 1, 1, x, err, 0.683, true);
        ALU_fit[i]     = ffit[i]->Eval(x[0]);
        ALU_fit_err[i] = err[0];
        cout<< ALU_fit[i]<<"     "<<ALU_fit_err[i]<<endl;
   }


   TGraphErrors *me_gre = new TGraphErrors(3);
                 me_gre->SetName("Graph0");
                 me_gre->SetTitle("This Work (<x_{B}>= 0.177, <-t> = 0.100)");
                 me_gre->SetMarkerColor(kBlack);
                 me_gre->SetFillColor(1);
                 me_gre->SetLineWidth(3);
                 me_gre->SetMarkerSize(1.5);
                 me_gre->SetMarkerStyle(21);
                 me_gre->SetPoint(0,1.144, ALU_fit[0]);    me_gre->SetPointError(0,0, ALU_fit_err[0]);
                 me_gre->SetPoint(1,1.422, ALU_fit[1]);    me_gre->SetPointError(1,0, ALU_fit_err[1]);
                 me_gre->SetPoint(2,1.892, ALU_fit[2]);    me_gre->SetPointError(2,0, ALU_fit_err[2]); 
                 me_gre->Draw("psame"); 


   TGraph *sys_gre = new TGraph();
           sys_gre->SetFillColor(40);
           sys_gre->SetPoint(0,1.0,0);
           sys_gre->SetPoint(1,1.0,0.101*0.3743/(1+0.1775));
           sys_gre->SetPoint(2,1.143,0.101*0.3743/(1+0.1775));
           sys_gre->SetPoint(3,1.423,0.101*0.35737/(1-0.036077));
           sys_gre->SetPoint(4,1.902,0.101*0.22352/(1-0.18));
           sys_gre->SetPoint(5,3.0,0.101*0.22352/(1-0.18));
           sys_gre->SetPoint(6,3.0,0.);
           sys_gre->Draw("fsame"); 


    TLegend* leg_Q2 = new TLegend(0.4,0.9,0.98,1.0);
             leg_Q2-> SetNColumns(1);
             leg_Q2->SetTextSize(0.05);
             leg_Q2->SetFillColor(0);
             leg_Q2->AddEntry(me_gre,"x_{B}= 0.18, -t= 0.10 GeV^{2}","P");
  leg_Q2->Draw();

  TLine *l_zero = new TLine(); 
         l_zero->SetLineColor(kBlack);
         l_zero->SetLineWidth(1);
         l_zero->SetLineStyle(7);
   l_zero->DrawLine(0.85, 0, 2.6, 0);


   // plot the coherent ALU ratio at phi = 90 deg vs Q2  -------------------------

   Can_r->cd(1);
   double y_init = 0.2;
   double y_fin = 3.9;
   TH1 *HALU_RATIO_Q2 = new TH2F("HALU_RATIO_Q2","",100,0.5,3.1,100,y_init, y_fin);
        HALU_RATIO_Q2->SetStats(0);
        HALU_RATIO_Q2->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
        HALU_RATIO_Q2->GetXaxis()->SetNdivisions(205);
        HALU_RATIO_Q2->GetYaxis()->SetNdivisions(511);
        HALU_RATIO_Q2->GetYaxis()->SetTitleSize(xfont);
        HALU_RATIO_Q2->GetXaxis()->SetTitleSize(xfont);
        HALU_RATIO_Q2->GetYaxis()->SetTitle("A_{LU}^{^{4}He}/A_{LU}^{p}(90^{#circ})");
        HALU_RATIO_Q2->GetXaxis()->CenterTitle(true);
        HALU_RATIO_Q2->GetYaxis()->CenterTitle(true);
        HALU_RATIO_Q2->GetXaxis()->SetTitleOffset(1.5);
        HALU_RATIO_Q2->GetYaxis()->SetTitleOffset(1.5);
        HALU_RATIO_Q2->Draw("");

   TGraphErrors *gre_Ratio_Q2 = new TGraphErrors(3);
                 gre_Ratio_Q2->SetName("This Work (x_{B} = 0.177");
                 gre_Ratio_Q2->SetTitle("This Work (x_{B} = 0.177");
                 gre_Ratio_Q2->SetMarkerColor(kBlack);
                 gre_Ratio_Q2->SetLineColor(kBlack);
                 gre_Ratio_Q2->SetFillColor(1);
                 gre_Ratio_Q2->SetLineWidth(3);
                 gre_Ratio_Q2->SetMarkerSize(1.6);
                 gre_Ratio_Q2->SetMarkerStyle(21);
                 gre_Ratio_Q2->SetPoint(0,1.144, 1.55242);  gre_Ratio_Q2->SetPointError(0,0,  0.857011 );
                 gre_Ratio_Q2->SetPoint(1,1.422, 2.13761 );  gre_Ratio_Q2->SetPointError(1,0, 0.431313 );
                 gre_Ratio_Q2->SetPoint(2,1.892, 1.62405 );  gre_Ratio_Q2->SetPointError(2,0, 0.808316 );
                 gre_Ratio_Q2->Draw("psame");

   line = new TLine(0.5,1,3.1,1);
   line->SetLineWidth(2);
   line->SetLineStyle(4);
   line->Draw();



    TLegend* leg_R_Q2 = new TLegend(0.45,0.92,0.96,1.0);
             leg_R_Q2-> SetNColumns(1);
             leg_R_Q2->SetFillColor(0);
             leg_R_Q2->AddEntry(gre_Ratio_Q2,"CLAS-EG6: x_{B}= 0.18, -t= 0.10 [GeV^{2}]","P");
  leg_R_Q2->Draw();



////////////////////////////////////////////////////////////////////////// 
//________________                                  ____________________//
//________________      bins in xB                  ____________________//                          
//________________                                  ____________________//


   double xvQ2[3] = {1.164, 1.439, 1.844};          
   double xvxB[3] = {0.132, 0.170, 0.225};          
   double xvt[3] =  {-0.095, -0.099, -0.107};       

  calculate_CFF(xvQ2[0], xvxB[0], xvt[0], A0[0], A1[0], A2[0], A3[0], c0_BH[0], c1_BH[0], c2_BH[0]);  
  calculate_CFF(xvQ2[1], xvxB[1], xvt[1], A0[1], A1[1], A2[1], A3[1], c0_BH[1], c1_BH[1], c2_BH[1]);  
  calculate_CFF(xvQ2[2], xvxB[2], xvt[2], A0[2], A1[2], A2[2], A3[2], c0_BH[2], c1_BH[2], c2_BH[2]);
  

   for(int i =0; i<3; i++) {  
      ffit[i] = new TF1(Form("ffit[%u]",i), Form("%0.7f*[0]*sin(x*3.14/180.0) / (%0.7f+ %0.7f*cos(x*3.14/180.0) + %0.7f*cos(2*x*3.14/180.0) + %0.7f*([0]*[0] + [1]*[1]) + %0.7f*[1] + %0.7f*[1]*cos(x*3.14/180.0))", A0[i], c0_BH[i], c1_BH[i], c2_BH[i], A1[i], A2[i], A3[i]),0.0,360.0);
      ffit[i]->SetLineColor(kRed);
      ffit[i]->SetParName(0,"Im(H_{A})");
      ffit[i]->SetParName(1,"Re(H_{A})");
 
     } 

// ------------>Primitives in pad: c13_1
   c12->cd(4);
   multigraph = new TMultiGraph();

   mg_title = new TLatex(xax,yay,"x_{B}= 0.132");
   mg_title->SetTextSize(0.085); 
   mg_title->SetTextFont(62); 
   mg_title->SetTextColor(kBlack);

   gre_sys = new TGraph(13);
   gre_sys->SetFillColor(40);
   gre_sys->SetLineColor(40);
   gre_sys->SetMarkerColor(40);
   gre_sys->SetPoint(0,0,   0.0        );
   gre_sys->SetPoint(1,0,    0.02214862);
   gre_sys->SetPoint(2,20,   0.02214862);
   gre_sys->SetPoint(3,60,   0.01984295);
   gre_sys->SetPoint(4,100,  0.04084068);
   gre_sys->SetPoint(5,140,  0.03336946);
   gre_sys->SetPoint(6,180,  0.02937379);
   gre_sys->SetPoint(7,220,  0.02438166);
   gre_sys->SetPoint(8,260,  0.0372357 );
   gre_sys->SetPoint(9,300,  0.02598807);
   gre_sys->SetPoint(10,340, 0.03198552);
   gre_sys->SetPoint(11,360, 0.03198552);
   gre_sys->SetPoint(12,360,0.0        );

   gre_st = new TGraphErrors(9);
   gre_st->SetFillColor(8);
   gre_st->SetMarkerColor(kBlack); 
   gre_st->SetMarkerSize(1.3);   
   gre_st->SetLineWidth(2);
   gre_st->SetMarkerStyle(21);
   gre_st->SetPoint(0,  25.65164,   0.01743518 );     gre_st->SetPointError(0,0, sqrt(pow( 0.1440989 ,2) - pow( 0.02214862,2) ));
   gre_st->SetPoint(1,  61.64267,   0.3478866  );     gre_st->SetPointError(1,0, sqrt(pow( 0.08670367,2) - pow( 0.01984295,2) ));
   gre_st->SetPoint(2,  99.02209,   0.380854   );     gre_st->SetPointError(2,0, sqrt(pow( 0.09484749,2) - pow( 0.04084068,2) ));
   gre_st->SetPoint(3,  141.7083,   0.2943698  );     gre_st->SetPointError(3,0, sqrt(pow( 0.1379186 ,2) - pow( 0.03336946,2) ));
   gre_st->SetPoint(4,  178.1952,   0.04334662 );     gre_st->SetPointError(4,0, sqrt(pow( 0.1519096 ,2) - pow( 0.02937379,2) ));
   gre_st->SetPoint(5,  218.5245,   -0.03520965);     gre_st->SetPointError(5,0, sqrt(pow( 0.131616  ,2) - pow( 0.02438166,2) ));
   gre_st->SetPoint(6,  262.7358,   -0.2771567 );     gre_st->SetPointError(6,0, sqrt(pow( 0.1050832 ,2) - pow( 0.0372357 ,2) ));
   gre_st->SetPoint(7,  300.7771,   -0.2138351 );     gre_st->SetPointError(7,0, sqrt(pow( 0.08368387,2) - pow( 0.02598807,2) ));
   gre_st->SetPoint(8,  335.2499,   -0.2339611 );     gre_st->SetPointError(8,0, sqrt(pow( 0.1221675 ,2) - pow( 0.03198552,2) ));
   multigraph->Add(gre_st,"");

   gre = new TGraphErrors(9);
   gre->SetFillColor(1);
   gre->SetMarkerColor(kBlue);
   gre->SetMarkerStyle(21);
   gre->SetPoint(0,  25.65164,   0.01743518 );     gre->SetPointError(0,0, 0.1440989  );
   gre->SetPoint(1,  61.64267,   0.3478866  );     gre->SetPointError(1,0, 0.08670367 );
   gre->SetPoint(2,  99.02209,   0.380854   );     gre->SetPointError(2,0, 0.09484749 );
   gre->SetPoint(3,  141.7083,   0.2943698  );     gre->SetPointError(3,0, 0.1379186  );
   gre->SetPoint(4,  178.1952,   0.04334662 );     gre->SetPointError(4,0, 0.1519096  );
   gre->SetPoint(5,  218.5245,   -0.03520965);     gre->SetPointError(5,0, 0.131616   );
   gre->SetPoint(6,  262.7358,   -0.2771567 );     gre->SetPointError(6,0, 0.1050832  );
   gre->SetPoint(7,  300.7771,   -0.2138351 );     gre->SetPointError(7,0, 0.08368387 );
   gre->SetPoint(8,  335.2499,   -0.2339611 );     gre->SetPointError(8,0, 0.1221675  );
   
   //multigraph->Add(gre,"");
   gre->Fit(ffit[0]);
   r[0] = gre->Fit(ffit[0],"S");

   multigraph->Draw("a2"); multigraph->Draw("p");
   multigraph->GetXaxis()->SetTitle("#phi [deg.]");
   multigraph->GetYaxis()->SetTitle("A_{LU}^{^{4}He}");
   multigraph->GetXaxis()->CenterTitle(true);
   multigraph->GetYaxis()->CenterTitle(true);
   multigraph->GetYaxis()->SetLabelSize(lfont);
   multigraph->GetXaxis()->SetLabelSize(lfont);
   multigraph->GetYaxis()->SetTitleSize(xfont);
   multigraph->GetXaxis()->SetTitleSize(xfont);
   multigraph->GetYaxis()->SetTitleOffset(0.8);
   gPad->Modified();  
   multigraph->GetYaxis()->SetRangeUser(-0.65, 0.65);
   multigraph->GetXaxis()->SetLimits(-5.0,  375);
   gre_sys->Draw("fsame");
   gre_st->Draw("psame");  ffit[0]->Draw("same");
   mg_title->Draw("same");
  
// ------------>Primitives in pad: c13_2

   c12->cd(5);
   multigraph = new TMultiGraph();
   
   mg_title = new TLatex(xax,yay,"x_{B}= 0.170");
   mg_title->SetTextSize(0.085); 
   mg_title->SetTextFont(62); 
   mg_title->SetTextColor(kBlack);
 
   gre_sys = new TGraph(13);
   gre_sys->SetFillColor(40);
   gre_sys->SetLineColor(40);
   gre_sys->SetMarkerColor(40);
   gre_sys->SetPoint(0,0,    0.0        );
   gre_sys->SetPoint(1,0,    0.023002136);
   gre_sys->SetPoint(2,20,   0.023002136);
   gre_sys->SetPoint(3,60,   0.01959182 );
   gre_sys->SetPoint(4,100,  0.0295878  );
   gre_sys->SetPoint(5,140,  0.01876947 );
   gre_sys->SetPoint(6,180,  0.0197908  );
   gre_sys->SetPoint(7,220,  0.0270881  );
   gre_sys->SetPoint(8,260,  0.02923329 );
   gre_sys->SetPoint(9,300,  0.03312501 );
   gre_sys->SetPoint(10,340, 0.02543044 );
   gre_sys->SetPoint(11,360, 0.02543044 );
   gre_sys->SetPoint(12,360, 0.0        );

   gre_st = new TGraphErrors(9);
   gre_st->SetFillColor(8);
   gre_st->SetMarkerColor(kBlack);   
   gre_st->SetMarkerSize(1.3);   
   gre_st->SetLineWidth(2);
   gre_st->SetMarkerStyle(21);
   gre_st->SetPoint(0,  22.96056 ,  0.1579036 );  gre_st->SetPointError(0,0, sqrt(pow(  0.08537671 ,2) - pow(0.023002136 ,2) ));
   gre_st->SetPoint(1,  57.15648 ,  0.1733943 );  gre_st->SetPointError(1,0, sqrt(pow(  0.08805212 ,2) - pow(0.01959182  ,2) ));
   gre_st->SetPoint(2,  96.39904 ,  0.2259656 );  gre_st->SetPointError(2,0, sqrt(pow(  0.1331627  ,2) - pow(0.0295878   ,2) ));
   gre_st->SetPoint(3,  139.1926 ,  0.2454004 );  gre_st->SetPointError(3,0, sqrt(pow(  0.1758729  ,2) - pow(0.01876947  ,2) ));
   gre_st->SetPoint(4,  180.1678 ,  -0.1021136);  gre_st->SetPointError(4,0, sqrt(pow(  0.19183    ,2) - pow(0.0197908   ,2) ));
   gre_st->SetPoint(5,  218.9244 ,  -0.287612 );  gre_st->SetPointError(5,0, sqrt(pow(  0.1909183  ,2) - pow(0.0270881   ,2) ));
   gre_st->SetPoint(6,  263.9158 ,  -0.2943798);  gre_st->SetPointError(6,0, sqrt(pow(  0.1360344  ,2) - pow(0.02923329  ,2) ));
   gre_st->SetPoint(7,  302.5475 ,  -0.3979002);  gre_st->SetPointError(7,0, sqrt(pow(  0.09214939 ,2) - pow(0.03312501  ,2) ));
   gre_st->SetPoint(8,  338.1138 ,  -0.2692495);  gre_st->SetPointError(8,0, sqrt(pow(  0.08347819 ,2) - pow(0.02543044  ,2) ));
   multigraph->Add(gre_st,"");

   gre = new TGraphErrors(9);
   gre->SetFillColor(1);
   gre->SetMarkerColor(kBlue);
   gre->SetMarkerStyle(21);
   gre->SetPoint(0,  22.96056 ,  0.1579036  );    gre->SetPointError(0,0,  0.08537671 );
   gre->SetPoint(1,  57.15648 ,  0.1733943  );    gre->SetPointError(1,0,  0.08805212 );
   gre->SetPoint(2,  96.39904 ,  0.2259656  );    gre->SetPointError(2,0,  0.1331627  );
   gre->SetPoint(3,  139.1926 ,  0.2454004  );    gre->SetPointError(3,0,  0.1758729  );
   gre->SetPoint(4,  180.1678 ,  -0.1021136 );    gre->SetPointError(4,0,  0.19183    );
   gre->SetPoint(5,  218.9244 ,  -0.287612  );    gre->SetPointError(5,0,  0.1909183  );
   gre->SetPoint(6,  263.9158 ,  -0.2943798 );    gre->SetPointError(6,0,  0.1360344  );
   gre->SetPoint(7,  302.5475 ,  -0.3979002 );    gre->SetPointError(7,0,  0.09214939 );
   gre->SetPoint(8,  338.1138 ,  -0.2692495 );    gre->SetPointError(8,0,  0.08347819 );
   
   //multigraph->Add(gre,"");
   gre->Fit(ffit[1]);
   r[1] = gre->Fit(ffit[1],"S");

   multigraph->Draw("a2");  multigraph->Draw("p");
   multigraph->GetXaxis()->SetTitle("#phi [deg.]");
   multigraph->GetXaxis()->CenterTitle(true);
   multigraph->GetXaxis()->SetLabelFont(22);
   multigraph->GetXaxis()->SetLabelSize(xfont);
   multigraph->GetYaxis()->SetLabelSize(xfont);
   multigraph->GetYaxis()->SetTitleSize(xfont);
   multigraph->GetXaxis()->SetTitleSize(xfont);
   multigraph->GetXaxis()->SetTitleFont(22);
   multigraph->GetYaxis()->SetLabelFont(22);
   gPad->Modified();  
   multigraph->GetYaxis()->SetRangeUser(-0.65, 0.65);
   multigraph->GetXaxis()->SetLimits(-5.0,  375);
   gre_sys->Draw("fsame");
   gre_st->Draw("psame");  ffit[1]->Draw("same");
   mg_title->Draw("same");
  
// ------------>Primitives in pad: c13_3
   c12->cd(6);
   multigraph = new TMultiGraph();

   mg_title = new TLatex(xax,yay,"x_{B}= 0.225");
   mg_title->SetTextSize(0.085); 
   mg_title->SetTextFont(62); 
   mg_title->SetTextColor(kBlack);
 
   gre_sys = new TGraph(13);
   gre_sys->SetFillColor(40);
   gre_sys->SetLineColor(40);
   gre_sys->SetMarkerColor(40);
   gre_sys->SetPoint(0,0,    0.0        );
   gre_sys->SetPoint(1,0,    0.024893846);
   gre_sys->SetPoint(2,20,   0.024893846);
   gre_sys->SetPoint(3,60,   0.02240954 );
   gre_sys->SetPoint(4,100,  0.02714954 );
   gre_sys->SetPoint(5,140,  0.01668931 );
   gre_sys->SetPoint(6,180,  0.02044162 );
   gre_sys->SetPoint(7,220,  0.02706752 );
   gre_sys->SetPoint(8,260,  0.03533464 );
   gre_sys->SetPoint(9,300,  0.02630272 );
   gre_sys->SetPoint(10,340, 0.01755089 );
   gre_sys->SetPoint(11,360, 0.01755089 );
   gre_sys->SetPoint(12,360, 0.0         );

   gre_st = new TGraphErrors(9);
   gre_st->SetFillColor(8);
   gre_st->SetMarkerColor(kBlack);   
   gre_st->SetMarkerSize(1.3);   
   gre_st->SetLineWidth(2);
   gre_st->SetMarkerStyle(21);
   gre_st->SetPoint(0,  20.13681,  0.2626366  );     gre_st->SetPointError(0,0, sqrt(pow(  0.07551056 ,2) - pow(0.024893846 ,2) ));
   gre_st->SetPoint(1,  55.99932,  0.428144   );     gre_st->SetPointError(1,0, sqrt(pow(  0.08872467 ,2) - pow(0.02240954  ,2) ));
   gre_st->SetPoint(2,  96.45656,  0.4925123  );     gre_st->SetPointError(2,0, sqrt(pow(  0.139086   ,2) - pow(0.02714954  ,2) ));
   gre_st->SetPoint(3,  137.7519,  -0.2736191 );     gre_st->SetPointError(3,0, sqrt(pow(  0.2804916  ,2) - pow(0.01668931  ,2) ));
   gre_st->SetPoint(4,  179.8777,  0.8470343  );     gre_st->SetPointError(4,0, sqrt(pow(  0.2497374  ,2) - pow(0.02044162  ,2) ));
   gre_st->SetPoint(5,  224.7934,  -0.05130243);     gre_st->SetPointError(5,0, sqrt(pow(  0.2811212  ,2) - pow(0.02706752  ,2) ));
   gre_st->SetPoint(6,  262.839 ,  -0.3416892 );     gre_st->SetPointError(6,0, sqrt(pow(  0.1694496  ,2) - pow(0.03533464  ,2) ));
   gre_st->SetPoint(7,  304.5033,  -0.135909  );     gre_st->SetPointError(7,0, sqrt(pow(  0.1025757  ,2) - pow(0.02630272  ,2) ));
   gre_st->SetPoint(8,  339.8168,  -0.1658634 );     gre_st->SetPointError(8,0, sqrt(pow(  0.07701008 ,2) - pow(0.01755089  ,2) ));
   multigraph->Add(gre_st,"");


   gre = new TGraphErrors(9);
   gre->SetName("Graph");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);
   gre->SetMarkerColor(kBlue);
   gre->SetMarkerStyle(21);
   gre->SetPoint(0,  20.13681,  0.2626366   );    gre->SetPointError(0,0,  0.07551056 );
   gre->SetPoint(1,  55.99932,  0.428144    );    gre->SetPointError(1,0,  0.08872467 );
   gre->SetPoint(2,  96.45656,  0.4925123   );    gre->SetPointError(2,0,  0.139086   );
   gre->SetPoint(3,  137.7519,  -0.2736191  );    gre->SetPointError(3,0,  0.2804916  );
   gre->SetPoint(4,  179.8777,  0.8470343   );    gre->SetPointError(4,0,  0.2497374  );
   gre->SetPoint(5,  224.7934,  -0.05130243 );    gre->SetPointError(5,0,  0.2811212  );
   gre->SetPoint(6,  262.839 ,  -0.3416892  );    gre->SetPointError(6,0,  0.1694496  );
   gre->SetPoint(7,  304.5033,  -0.135909   );    gre->SetPointError(7,0,  0.1025757  );
   gre->SetPoint(8,  339.8168,  -0.1658634  );    gre->SetPointError(8,0,  0.07701008 );
   
   //multigraph->Add(gre,"");
   gre->Fit(ffit[2]);
   r[2] = gre->Fit(ffit[2],"S");

   multigraph->Draw("a2");  multigraph->Draw("p");
   multigraph->GetXaxis()->SetTitle("#phi [deg.]");
   multigraph->GetXaxis()->CenterTitle(true);
   multigraph->GetXaxis()->SetLabelFont(22);
   multigraph->GetXaxis()->SetLabelSize(xfont);
   multigraph->GetXaxis()->SetTitleSize(xfont);
   multigraph->GetYaxis()->SetLabelSize(xfont);
   multigraph->GetYaxis()->SetTitleSize(xfont);
   multigraph->GetXaxis()->SetTitleFont(22);
   multigraph->GetYaxis()->SetLabelFont(22);
   gPad->Modified();  
   multigraph->GetYaxis()->SetRangeUser(-0.65, 0.65);
   multigraph->GetXaxis()->SetLimits(-5.0,  375);
   gre_sys->Draw("fsame");
   gre_st->Draw("psame");  ffit[2]->Draw("same");
   mg_title->Draw("same");

   //____________________________________________________________________________
   // plot the Imaginary and the real parts of HA CFF vs xB ---------------------

   cc->cd(2);
   mgCFF = new TMultiGraph();
 
   Double_t graphCFF_fx101[3] = { 0.132, 0.17, 0.225};
   Double_t graphCFF_fy101[3] = {ffit[0]->GetParameter(0),
                                 ffit[1]->GetParameter(0),
                                 ffit[2]->GetParameter(0)};   
   Double_t graphCFF_fex101[3] = { 0, 0, 0};
   Double_t graphCFF_fey101[3] = {ffit[0]->GetParError(0),
                                  ffit[1]->GetParError(0),
                                  ffit[2]->GetParError(0)};

   greCFF = new TGraphErrors(3,graphCFF_fx101,graphCFF_fy101,graphCFF_fex101,graphCFF_fey101);
   greCFF->SetFillColor(1);
   greCFF->SetMarkerColor(kBlack);
   greCFF->SetLineColor(kBlack); 
   greCFF->SetMarkerStyle(21);
   greCFF->SetMarkerSize(1.5);
   greCFF->SetLineWidth(2);
   
   mgCFF->Add(greCFF,"");

   mgCFF->Draw("AP");
   mgCFF->GetXaxis()->SetTitle("x_{B}");
   mgCFF->GetYaxis()->SetTitle("Im, Re(H_{A})");
   mgCFF->GetXaxis()->SetLabelSize(lfont);
   mgCFF->GetYaxis()->SetLabelSize(lfont);
   mgCFF->GetXaxis()->SetTitleSize(xfont);
   mgCFF->GetYaxis()->SetTitleSize(xfont);
   mgCFF->GetXaxis()->CenterTitle(true);
   mgCFF->GetYaxis()->CenterTitle(true);
   mgCFF->GetXaxis()->SetLimits(0.08,0.32);
   mgCFF->GetYaxis()->SetRangeUser(-2.0, 78.0);
   
   // imaginary part as a function of x
   graphCFF = new TGraph(21);
   graphCFF->SetFillColor(1);
   graphCFF->SetLineColor(2);
   graphCFF->SetLineStyle(1);
   graphCFF->SetLineWidth(4);
   graphCFF->SetPoint(0,  0.1000  , 36.7372 );
   graphCFF->SetPoint(1,  0.1100  , 33.9679 );
   graphCFF->SetPoint(2,  0.1200  , 31.1986 );
   graphCFF->SetPoint(3,  0.1300  , 28.4292 );
   graphCFF->SetPoint(4,  0.1400  , 25.6599 );
   graphCFF->SetPoint(5,  0.1500  , 22.8906 );
   graphCFF->SetPoint(6,  0.1600  , 21.3965 );
   graphCFF->SetPoint(7,  0.1700  , 19.9025 );
   graphCFF->SetPoint(8,  0.1800  , 18.4085 );
   graphCFF->SetPoint(9,  0.1900  , 16.9145 );
   graphCFF->SetPoint(10, 0.2000  , 15.4205 );
   graphCFF->SetPoint(11, 0.2100  , 14.4785 );
   graphCFF->SetPoint(12, 0.2200  , 13.5366 );
   graphCFF->SetPoint(13, 0.2300  , 12.5946 );
   graphCFF->SetPoint(14, 0.2400  , 11.6527 );
   graphCFF->SetPoint(15, 0.2500  , 10.7108 );
   graphCFF->SetPoint(16, 0.2600  , 10.0684 );
   graphCFF->SetPoint(17, 0.2700  ,  9.4261 );
   graphCFF->SetPoint(18, 0.2800  ,  8.7838 );
   graphCFF->SetPoint(19, 0.2900  ,  8.1414 );
   graphCFF->SetPoint(20, 0.3000  ,  7.4991 );
   graphCFF->Draw("csame");

   // based on the VGG model
   TGraph *VgraphCFF_x = new TGraph(5);
           VgraphCFF_x->SetFillColor(kBlack);
           VgraphCFF_x->SetLineColor(kBlack);
           VgraphCFF_x->SetLineStyle(7);
           VgraphCFF_x->SetLineWidth(4);
           VgraphCFF_x->SetPoint(0,  0.10 , 23.4182 );
           VgraphCFF_x->SetPoint(1,  0.15 , 16.6466 );
           VgraphCFF_x->SetPoint(2,  0.20 , 13.0322 );
           VgraphCFF_x->SetPoint(3,  0.25 , 10.6674 );
           VgraphCFF_x->SetPoint(4,  0.30 , 8.93059 );
           VgraphCFF_x->Draw("csame");



   // real part of the CFF HA vs x   -------------------------
   cc->cd(5);
   mgCFF = new TMultiGraph();

   Double_t graphCFF_fx102[3] = { 0.132, 0.17, 0.225};
   Double_t graphCFF_fy102[3] = {ffit[0]->GetParameter(1),
                                 ffit[1]->GetParameter(1),
                                 ffit[2]->GetParameter(1)};
   Double_t graphCFF_fex102[3] = {0, 0, 0};
   Double_t graphCFF_fey102[3] = {ffit[0]->GetParError(1), 
                                  ffit[1]->GetParError(1),
                                  ffit[2]->GetParError(1)};

   greCFF = new TGraphErrors(3,graphCFF_fx102,graphCFF_fy102,graphCFF_fex102,graphCFF_fey102);
   greCFF->SetFillColor(1);
   greCFF->SetMarkerColor(kBlack);
   greCFF->SetLineColor(kBlack); 
   greCFF->SetMarkerStyle(25);
   greCFF->SetMarkerSize(1.5);
   greCFF->SetLineWidth(2);
   
   mgCFF->Add(greCFF,"");

   mgCFF->Draw("AP");
   mgCFF->GetXaxis()->SetTitle("x_{B}");
   mgCFF->GetXaxis()->SetLabelSize(lfont+0.01);
   mgCFF->GetXaxis()->SetTitleSize(xfont+0.03);
   mgCFF->GetXaxis()->SetTitleOffset(0.65);
   mgCFF->GetXaxis()->CenterTitle(true);
   mgCFF->GetYaxis()->CenterTitle(true);
   mgCFF->GetYaxis()->SetRangeUser(-60.0, 65.0);
   mgCFF->GetXaxis()->SetLimits(0.08,0.32);
   
   // from Vadim's routines
   graphCFF = new TGraph(21);
   graphCFF->SetName("graphCFF");
   graphCFF->SetTitle("graphCFF");
   graphCFF->SetFillColor(1);
   graphCFF->SetLineColor(2);
   graphCFF->SetLineStyle(1);
   graphCFF->SetLineWidth(4);
   graphCFF->SetPoint(0,  0.1000  ,  -3.5974 );
   graphCFF->SetPoint(1,  0.1100  ,  -3.6032 );
   graphCFF->SetPoint(2,  0.1200  ,  -3.6089 );
   graphCFF->SetPoint(3,  0.1300  ,  -3.6147 );
   graphCFF->SetPoint(4,  0.1400  ,  -3.6204 );
   graphCFF->SetPoint(5,  0.1500  ,  -3.6261 );
   graphCFF->SetPoint(6,  0.1600  ,  -3.7634 );
   graphCFF->SetPoint(7,  0.1700  ,  -3.9007 );
   graphCFF->SetPoint(8,  0.1800  ,  -4.0379 );
   graphCFF->SetPoint(9,  0.1900  ,  -4.1752 );
   graphCFF->SetPoint(10, 0.2000  ,  -4.3125 );
   graphCFF->SetPoint(11, 0.2100  ,  -4.5291 );
   graphCFF->SetPoint(12, 0.2200  ,  -4.7457 );
   graphCFF->SetPoint(13, 0.2300  ,  -4.9623 );
   graphCFF->SetPoint(14, 0.2400  ,  -5.1789 );
   graphCFF->SetPoint(15, 0.2500  ,  -5.3955 );
   graphCFF->SetPoint(16, 0.2600  ,  -5.6701 );
   graphCFF->SetPoint(17, 0.2700  ,  -5.9447 );
   graphCFF->SetPoint(18, 0.2800  ,  -6.2193 );
   graphCFF->SetPoint(19, 0.2900  ,  -6.4938 );
   graphCFF->SetPoint(20, 0.3000  ,  -6.7684 );
   graphCFF->Draw("csame");
 
   // based on the VGG model
           VgraphCFF_x = new TGraph(5);
           VgraphCFF_x->SetFillColor(kBlack);
           VgraphCFF_x->SetLineColor(kBlack);
           VgraphCFF_x->SetLineStyle(7);
           VgraphCFF_x->SetLineWidth(4);
           VgraphCFF_x->SetPoint(0,  0.10 ,  -3.7714);
           VgraphCFF_x->SetPoint(1,  0.15 ,  -2.16029);
           VgraphCFF_x->SetPoint(2,  0.20 ,  -1.0151);
           VgraphCFF_x->SetPoint(3,  0.25 ,  -0.161793);
           VgraphCFF_x->SetPoint(4,  0.30 ,   0.484994);
           VgraphCFF_x->Draw("csame");


   // plot the coherent ALU at phi = 90 deg vs xB  -------------------------------
   Can->cd(2);
   
   TH1 *HALU_He_vs_x = new TH2F("HALU_He_vs_x","",100,0.08,0.31,100,-0.1,0.72);
        HALU_He_vs_x->SetStats(0);
        HALU_He_vs_x->GetXaxis()->SetTitle("x_{B}");
        //HALU_He_vs_x->GetYaxis()->SetTitle("A_{LU}^{^{4}He}(90#circ)");
        HALU_He_vs_x->GetXaxis()->SetTitleSize(xfont-0.05);
        HALU_He_vs_x->GetYaxis()->SetTitleSize(xfont);
        HALU_He_vs_x->GetXaxis()->SetNdivisions(205);
        HALU_He_vs_x->GetXaxis()->CenterTitle(true);
        HALU_He_vs_x->GetYaxis()->CenterTitle(true);
        HALU_He_vs_x->GetYaxis()->SetNdivisions(414);
        //HALU_He_vs_x->GetXaxis()->SetTitleOffset(1.5);
        //HALU_He_vs_x->GetYaxis()->SetTitleOffset(1.5);
        HALU_He_vs_x->Draw("");
   
   TGraphErrors *gre1 = new TGraphErrors(3);
                 gre1->SetName("LT (-t=0.100)");
                 gre1->SetTitle("LT (-t= 0.100, Q^{2}= 1.2)");
                 gre1->SetFillColor(1);
                 gre1->SetLineColor(2);
                 gre1->SetLineWidth(3);
                 gre1->SetMarkerColor(2);
                 gre1->SetMarkerStyle(21);
                 gre1->SetMarkerSize(1.5);
                 gre1->SetPoint(0,0.137113,0.5);      gre1->SetPointError(0,0,0.00603854);
                 gre1->SetPoint(1,0.189332,0.5);      gre1->SetPointError(1,0,0.00728592);
                 gre1->SetPoint(2,0.261673,0.5);      gre1->SetPointError(2,0,0.00865828);

    for(int i=0; i<3; i++){
        double x[1] = { 90.0};  double err[1];
        r[i]->GetConfidenceIntervals(1, 1, 1, x, err, 0.683, true);
        ALU_fit[i]     = ffit[i]->Eval(x[0]);
        ALU_fit_err[i] = err[0];
        cout<< ALU_fit[i]<<"     "<<ALU_fit_err[i]<<endl;
   }

    me_gre = new TGraphErrors(3);
    me_gre->SetName("Graph0");
    me_gre->SetTitle("This Work (<-t>=0.100, <Q^{2}>= 1.492)");
    me_gre->SetFillColor(1);
    me_gre->SetLineWidth(3);
    me_gre->SetMarkerSize(1.5);
    me_gre->SetMarkerStyle(21);
    me_gre->SetPoint(0,0.132, ALU_fit[0] /*0.312989*/);  me_gre->SetPointError(0,0, ALU_fit_err[0]/*0.0697107*/);
    me_gre->SetPoint(1,0.17,  ALU_fit[1] /*0.237158*/);  me_gre->SetPointError(1,0, ALU_fit_err[1]/*0.0698417*/);
    me_gre->SetPoint(2,0.225, ALU_fit[2] /*0.299174*/);  me_gre->SetPointError(2,0, ALU_fit_err[2]/*0.126198*/);
    me_gre->Draw("psame");

   sys_gre = new TGraph();
   sys_gre->SetFillColor(40);
   sys_gre->SetPoint(0,0.1,0);
   sys_gre->SetPoint(1,0.1,0.101*0.31799/(1-0.1475));
   sys_gre->SetPoint(2,0.132,0.101*0.31799/(1-0.1475));
   sys_gre->SetPoint(3,0.17,0.101*0.3263/(1+0.224));
   sys_gre->SetPoint(4,0.225,0.101*0.2448/(1-0.325));
   sys_gre->SetPoint(5,0.3,0.101*0.2448/(1-0.325));
   sys_gre->SetPoint(6,0.3,0.);
   sys_gre->Draw("fsame"); 
   
   graph = new TGraph(16);
   graph->SetName("Graph");
   graph->SetTitle("Graph");
   graph->SetFillColor(1);
   graph->SetLineColor(2);
   graph->SetLineStyle(1);
   graph->SetLineWidth(3);
   graph->SetPoint(0,0.0053,0.25272);
   graph->SetPoint(1,0.0277,0.269698);
   graph->SetPoint(2,0.0672,0.265381);
   graph->SetPoint(3,0.1223,0.240483);
   graph->SetPoint(4,0.1911,0.190882);
   graph->SetPoint(5,0.271,0.0965136);
   graph->SetPoint(6,0.3592,0.0589742);
   graph->SetPoint(7,0.4525,0.0453655);
   graph->SetPoint(8,0.5475,0.035241);
   graph->SetPoint(9,0.6408,0.0294755);
   graph->SetPoint(10,0.729,0.0294434);
   graph->SetPoint(11,0.8089,0.0407613);
   graph->SetPoint(12,0.8777,0.095138);
   graph->SetPoint(13,0.9328,0.50658);
   graph->SetPoint(14,0.9723,11.1018);
   graph->SetPoint(15,0.9947,8001.42);
   graph->Draw("csame");
   
   
   graph = new TGraph(16);
   graph->SetName("Graph");
   graph->SetTitle("Graph");
   graph->SetFillColor(1);
   graph->SetLineColor(4);
   graph->SetLineStyle(7);
   graph->SetLineWidth(4);
   graph->SetPoint(0,0.0053,0.312336);
   graph->SetPoint(1,0.0277,0.327893);
   graph->SetPoint(2,0.0672,0.316997);
   graph->SetPoint(3,0.1223,0.287221);
   graph->SetPoint(4,0.1911,0.245252);
   graph->SetPoint(5,0.271,0.194805);
   graph->SetPoint(6,0.3592,0.134905);
   graph->SetPoint(7,0.4525,0.0435876);
   graph->SetPoint(8,0.5475,0.0354022);
   graph->SetPoint(9,0.6408,0.0317877);
   graph->SetPoint(10,0.729,0.0344951);
   graph->SetPoint(11,0.8089,0.0512647);
   graph->SetPoint(12,0.8777,0.125261);
   graph->SetPoint(13,0.9328,0.684747);
   graph->SetPoint(14,0.9723,15.1177);
   graph->SetPoint(15,0.9947,10920.6);
   graph->Draw("csame");
  
 

  TLine *ll_zero = new TLine(); 
         ll_zero->SetLineColor(kBlack);
         ll_zero->SetLineWidth(1);
         ll_zero->SetLineStyle(7);
   ll_zero->DrawLine(0.075, 0, 0.31, 0);


   TLegend* leg_x = new TLegend(0.1,0.82,0.98,1.0);
            leg_x->SetNColumns(1);
            leg_x->SetFillColor(0);
            leg_x->SetTextSize(0.05);
            leg_x->AddEntry(me_gre,"-t= 0.1 GeV^{2}, Q^{2}= 1.5 GeV^{2}","P");
            leg_x->AddEntry(gre1,"-t= 0.1 GeV^{2}, Q^{2}= 1.2 GeV^{2}","L");
            leg_x->AddEntry(graph,"-t= 0.2 GeV^{2}, Q^{2}= 2.1 GeV^{2}","L");

   leg_x->Draw();


   // plot the coherent ALU ratio at phi = 90 deg vs xB  -------------------------

   Can_r->cd(2);
   TH1 *HALU_RATIO_x = new TH2F("HALU_RATIO_x","",100,-0.02,0.83,100,y_init, y_fin);
        HALU_RATIO_x->SetStats(0);
        HALU_RATIO_x->GetXaxis()->SetTitle("x_{B}");
        HALU_RATIO_x->GetXaxis()->SetNdivisions(205);
        HALU_RATIO_x->GetYaxis()->SetNdivisions(511);
        HALU_RATIO_x->GetYaxis()->SetTitleSize(xfont);
        HALU_RATIO_x->GetXaxis()->SetTitleSize(xfont);
        HALU_RATIO_x->GetYaxis()->SetTitle("A_{LU}^{^{4}He}/A_{LU}^{p} (90^{#circ})");
        HALU_RATIO_x->GetXaxis()->CenterTitle(true);
        HALU_RATIO_x->GetYaxis()->CenterTitle(true);
        HALU_RATIO_x->GetXaxis()->SetTitleOffset(1.5);
        HALU_RATIO_x->GetYaxis()->SetTitleOffset(1.5);
        HALU_RATIO_x->Draw("");

    TGraphErrors *gre_Ratio_x = new TGraphErrors(3);
                  gre_Ratio_x->SetName("This Work (-t = 0.10");
                  gre_Ratio_x->SetTitle("This Work (-t = 0.10");
                  gre_Ratio_x->SetFillColor(1);
                  gre_Ratio_x->SetLineWidth(3);
                  gre_Ratio_x->SetLineColor(kBlack);
                  gre_Ratio_x->SetMarkerColor(kBlack);
                  gre_Ratio_x->SetMarkerStyle(21);
                  gre_Ratio_x->SetMarkerSize(1.6);
                  gre_Ratio_x->SetPoint(0,0.132, 1.63394 );    gre_Ratio_x->SetPointError(0,0, 0.888503 );
                  gre_Ratio_x->SetPoint(1,0.17,  1.67078 );    gre_Ratio_x->SetPointError(1,0, 0.518536 );
                  gre_Ratio_x->SetPoint(2,0.227, 1.78805 );    gre_Ratio_x->SetPointError(2,0, 1.23573  );
                  gre_Ratio_x->Draw("psame");
 

   TGraphErrors *gre_H_R_x = new TGraphErrors(1);
                 gre_H_R_x->SetName("HERMES (-t = 0.018)");
                 gre_H_R_x->SetTitle("HERMES (-t = 0.018)");
                 gre_H_R_x->SetFillColor(1);
                 gre_H_R_x->SetLineColor(8);
                 gre_H_R_x->SetLineWidth(3);
                 gre_H_R_x->SetMarkerColor(8);
                 gre_H_R_x->SetMarkerStyle(8);
                 gre_H_R_x->SetMarkerSize(1.5);
                 gre_H_R_x->SetPoint(0,0.072,0.91);   gre_H_R_x->SetPointError(0,0,0.19);
                 gre_H_R_x->Draw("psame");

  
   TGraph *gre_LT1_R_x = new TGraph(16);
           gre_LT1_R_x->SetFillColor(1);
           gre_LT1_R_x->SetLineColor(2);
           gre_LT1_R_x->SetLineStyle(1);
           gre_LT1_R_x->SetLineWidth(3);
           gre_LT1_R_x->SetPoint(0,0.0053,1.03844);
           gre_LT1_R_x->SetPoint(1,0.0277,1.11211);
           gre_LT1_R_x->SetPoint(2,0.0672,1.11531);
           gre_LT1_R_x->SetPoint(3,0.1223,1.07369);
           gre_LT1_R_x->SetPoint(4,0.1911,1.0049);
           gre_LT1_R_x->SetPoint(5,0.271,0.927646);
           gre_LT1_R_x->SetPoint(6,0.3592,0.860616);
           gre_LT1_R_x->SetPoint(7,0.4525,0.822943);
           gre_LT1_R_x->SetPoint(8,0.5475,0.839245);
           gre_LT1_R_x->SetPoint(9,0.6408,0.960647);
           gre_LT1_R_x->SetPoint(10,0.729,1.34325);
           gre_LT1_R_x->SetPoint(11,0.8089,2.60307);
           gre_LT1_R_x->SetPoint(12,0.8777,8.3076);
           gre_LT1_R_x->SetPoint(13,0.9328,57.7025);
           gre_LT1_R_x->SetPoint(14,0.9723,1542.62);
           gre_LT1_R_x->SetPoint(15,0.9947,1.24825e+06);
           gre_LT1_R_x->Draw("csame");
   
   gre_LT2_R_x = new TGraph(16);
   gre_LT2_R_x->SetFillColor(1);
   gre_LT2_R_x->SetLineColor(4);
   gre_LT2_R_x->SetLineStyle(1);
   gre_LT2_R_x->SetLineWidth(3);
   gre_LT2_R_x->SetPoint(0,0.0053,1.15609);
   gre_LT2_R_x->SetPoint(1,0.0277,1.21547);
   gre_LT2_R_x->SetPoint(2,0.0672,1.18441);
   gre_LT2_R_x->SetPoint(3,0.1223,1.09922);
   gre_LT2_R_x->SetPoint(4,0.1911,0.993413);
   gre_LT2_R_x->SetPoint(5,0.271,0.890928);
   gre_LT2_R_x->SetPoint(6,0.3592,0.81447);
   gre_LT2_R_x->SetPoint(7,0.4525,0.790692);
   gre_LT2_R_x->SetPoint(8,0.5475,0.843083);
   gre_LT2_R_x->SetPoint(9,0.6408,1.03601);
   gre_LT2_R_x->SetPoint(10,0.729,1.57371);
   gre_LT2_R_x->SetPoint(11,0.8089,3.27383);
   gre_LT2_R_x->SetPoint(12,0.8777,10.938);
   gre_LT2_R_x->SetPoint(13,0.9328,77.9969);
   gre_LT2_R_x->SetPoint(14,0.9723,2100.63);
   gre_LT2_R_x->SetPoint(15,0.9947,1.70365e+06);
   gre_LT2_R_x->Draw("csame");

   

  TGraphErrors *gre_vd1_R_x = new TGraphErrors(8);
                gre_vd1_R_x->SetName("GS (-t=0.1)");
                gre_vd1_R_x->SetTitle("GS (-t=0.1)");
                gre_vd1_R_x->SetMarkerColor(kBlack);
                gre_vd1_R_x->SetLineColor(1);
                gre_vd1_R_x->SetLineStyle(2);
                gre_vd1_R_x->SetLineWidth(3);
                gre_vd1_R_x->SetLineStyle(1);
                gre_vd1_R_x->SetMarkerSize(1.2);
                gre_vd1_R_x->SetMarkerStyle(21);
                gre_vd1_R_x->SetPoint(0,     0.1400 ,   1.7478);    gre_vd1_R_x->SetPointError(0,0,0.);
                gre_vd1_R_x->SetPoint(1,     0.1600 ,   1.6550);    gre_vd1_R_x->SetPointError(1,0,0.);
                gre_vd1_R_x->SetPoint(2,     0.1800 ,   1.5601);    gre_vd1_R_x->SetPointError(2,0,0.);
                gre_vd1_R_x->SetPoint(3,     0.2000 ,   1.5012);    gre_vd1_R_x->SetPointError(3,0,0.);
                gre_vd1_R_x->SetPoint(4,     0.2200 ,   1.4306);    gre_vd1_R_x->SetPointError(4,0,0.);
                gre_vd1_R_x->SetPoint(5,     0.2400 ,   1.3943);    gre_vd1_R_x->SetPointError(5,0,0.);
                gre_vd1_R_x->SetPoint(6,     0.2600 ,   1.3855);    gre_vd1_R_x->SetPointError(6,0,0.);
                gre_vd1_R_x->SetPoint(7,     0.2800 ,   1.4723);    gre_vd1_R_x->SetPointError(7,0,0.);
                gre_vd1_R_x->Draw("csame");

 TGraphErrors *gre_vd2_R_x = new TGraphErrors(14);
               gre_vd2_R_x->SetName("GS (-t=0.2)");
               gre_vd2_R_x->SetTitle("GS (-t=0.2)");
               gre_vd2_R_x->SetMarkerColor(kBlack);
               gre_vd2_R_x->SetLineColor(6);
               gre_vd2_R_x->SetLineWidth(3);
               gre_vd2_R_x->SetMarkerSize(1.2);
               gre_vd2_R_x->SetMarkerStyle(20);
               gre_vd2_R_x->SetPoint(0,  0.1400, 1.4312 );    gre_vd2_R_x->SetPointError(0,0,0.);
               gre_vd2_R_x->SetPoint(1,  0.1600, 1.3556 );    gre_vd2_R_x->SetPointError(1,0,0.);
               gre_vd2_R_x->SetPoint(2,  0.1800, 1.2792 );    gre_vd2_R_x->SetPointError(2,0,0.);
               gre_vd2_R_x->SetPoint(3,  0.2000, 1.2317 );    gre_vd2_R_x->SetPointError(3,0,0.);
               gre_vd2_R_x->SetPoint(4,  0.2200, 1.1778 );    gre_vd2_R_x->SetPointError(4,0,0.);
               gre_vd2_R_x->SetPoint(5,  0.2400, 1.1427 );    gre_vd2_R_x->SetPointError(5,0,0.);
               gre_vd2_R_x->SetPoint(6,  0.2600, 1.1124 );    gre_vd2_R_x->SetPointError(6,0,0.);
               gre_vd2_R_x->SetPoint(7,  0.2800, 1.0821 );    gre_vd2_R_x->SetPointError(7,0,0.);
               gre_vd2_R_x->SetPoint(8,  0.3000, 1.0598 );    gre_vd2_R_x->SetPointError(8,0,0.);
               gre_vd2_R_x->SetPoint(9,  0.3200, 1.0300 );    gre_vd2_R_x->SetPointError(9,0,0.);
               gre_vd2_R_x->SetPoint(10, 0.3400, 1.0104 );    gre_vd2_R_x->SetPointError(10,0,0.);
               gre_vd2_R_x->SetPoint(11, 0.3600, 1.0080 );    gre_vd2_R_x->SetPointError(11,0,0.);
               gre_vd2_R_x->SetPoint(12, 0.3800, 1.0648 );    gre_vd2_R_x->SetPointError(12,0,0.);
               gre_vd2_R_x->SetPoint(13, 0.4000, 1.4436 );    gre_vd2_R_x->SetPointError(13,0,0.);
               gre_vd2_R_x->Draw("csame");

   line = new TLine(-0.05,1,0.83,1);
   line->SetLineWidth(2);
   line->SetLineStyle(4);
   line->Draw();

   TLegend* leg_R_x = new TLegend(0.32,0.79,0.95,0.95);
            leg_R_x-> SetNColumns(2);
            leg_R_x->AddEntry(gre_Ratio_x,"CLAS-EG6: -t= 0.10","P");
            leg_R_x->AddEntry(gre_H_R_x,"HERMES: -t = 0.018","P"); 
            leg_R_x->AddEntry(gre_LT1_R_x,"Off-shell cal.: -t=0.10","L");
            leg_R_x->AddEntry(gre_LT2_R_x,"Off-shell cal.: -t=0.20","L");
            leg_R_x->AddEntry(gre_vd1_R_x,"On-shell cal.: -t = 0.10","L");
            leg_R_x->AddEntry(gre_vd2_R_x,"On-shell cal.: -t = 0.20","L");


   leg_R_x->Draw();


////////////////////////////////////////////////////////////////////////// 
//________________                                  ____________________//
//________________      bins in -t                  ____________________//                          
//________________                                  ____________________//

   double tvQ2[3] = {1.36, 1.507, 1.610};     
   double tvxB[3] = {0.160, 0.179, 0.193};       
   double tvt[3]  = {-0.08, -0.094, -0.127};   

  calculate_CFF(tvQ2[0], tvxB[0], tvt[0], A0[0], A1[0], A2[0], A3[0], c0_BH[0], c1_BH[0], c2_BH[0]);  
  calculate_CFF(tvQ2[1], tvxB[1], tvt[1], A0[1], A1[1], A2[1], A3[1], c0_BH[1], c1_BH[1], c2_BH[1]);  
  calculate_CFF(tvQ2[2], tvxB[2], tvt[2], A0[2], A1[2], A2[2], A3[2], c0_BH[2], c1_BH[2], c2_BH[2]);

   for(int i =0; i<3; i++) {  
      ffit[i] = new TF1(Form("ffit[%u]",i), Form("%0.7f*[0]*sin(x*3.14/180.0) / (%0.7f+ %0.7f*cos(x*3.14/180.0) + %0.7f*cos(2*x*3.14/180.0) + %0.7f*([0]*[0] + [1]*[1]) + %0.7f*[1] + %0.7f*[1]*cos(x*3.14/180.0))", A0[i], c0_BH[i], c1_BH[i], c2_BH[i], A1[i], A2[i], A3[i]),0.0,360.0);
      ffit[i]->SetLineColor(kRed);
      ffit[i]->SetParName(0,"Im(H_{A})");
      ffit[i]->SetParName(1,"Re(H_{A})");
     } 

// ------------>Primitives in pad: c14_1
   
   c12->cd(7);
   //Cc->cd();

   multigraph = new TMultiGraph();

   mg_title = new TLatex(xax,yay,"-t= 0.08 GeV^{2}");
   mg_title->SetTextSize(0.075);
   mg_title->SetTextFont(62); 
   mg_title->SetTextColor(kBlack);
 
   gre_sys = new TGraph(13);
   gre_sys->SetFillColor(40);
   gre_sys->SetLineColor(40);
   gre_sys->SetMarkerColor(40);
   gre_sys->SetPoint(0,0,    0.0        );
   gre_sys->SetPoint(1,0,    0.026305318);
   gre_sys->SetPoint(2,20,   0.026305318);
   gre_sys->SetPoint(3,60,   0.02439207 );
   gre_sys->SetPoint(4,100,  0.03947151 );
   gre_sys->SetPoint(5,140,  0.02522018 );
   gre_sys->SetPoint(6,180,  0.03717874 );
   gre_sys->SetPoint(7,220,  0.03708316 );
   gre_sys->SetPoint(8,260,  0.04038161 );
   gre_sys->SetPoint(9,300,  0.02940459 );
   gre_sys->SetPoint(10,340, 0.02810254 );
   gre_sys->SetPoint(11,360, 0.02810254 );
   gre_sys->SetPoint(12,360, 0.0        );


   gre_st = new TGraphErrors(9);
   gre_st->SetFillColor(8);
   gre_st->SetMarkerColor(kBlack);   
   gre_st->SetMarkerSize(1.3);   
   gre_st->SetLineWidth(2);
   gre_st->SetMarkerStyle(21);
   gre_st->SetPoint(0,   22.81909,  0.2378564 );     gre_st->SetPointError(0,0, sqrt(pow(  0.09283838 ,2) - pow(0.026305318 ,2) ));
   gre_st->SetPoint(1,   57.65211,  0.3014312 );     gre_st->SetPointError(1,0, sqrt(pow(  0.08745728 ,2) - pow(0.02439207  ,2) ));
   gre_st->SetPoint(2,   98.19281,  0.4896329 );     gre_st->SetPointError(2,0, sqrt(pow(  0.1118034  ,2) - pow(0.03947151  ,2) ));
   gre_st->SetPoint(3,   139.3701,  0.1965057 );     gre_st->SetPointError(3,0, sqrt(pow(  0.1595977  ,2) - pow(0.02522018  ,2) ));
   gre_st->SetPoint(4,   178.8885,  0.05771016);     gre_st->SetPointError(4,0, sqrt(pow(  0.1923073  ,2) - pow(0.03717874  ,2) ));
   gre_st->SetPoint(5,   222.6621,  -0.1651759);     gre_st->SetPointError(5,0, sqrt(pow(  0.1640617  ,2) - pow(0.03708316  ,2) ));
   gre_st->SetPoint(6,   265.6899,  -0.3467163);     gre_st->SetPointError(6,0, sqrt(pow(  0.133935   ,2) - pow(0.04038161  ,2) ));
   gre_st->SetPoint(7,   300.0582,  -0.2891963);     gre_st->SetPointError(7,0, sqrt(pow(  0.09282593 ,2) - pow(0.02940459  ,2) ));
   gre_st->SetPoint(8,   339.3341,  -0.1849013);     gre_st->SetPointError(8,0, sqrt(pow(  0.08560304 ,2) - pow(0.02810254  ,2) ));
   multigraph->Add(gre_st,"");

   gre = new TGraphErrors(9);
   gre->SetFillColor(1);
   gre->SetMarkerColor(kBlue);
   gre->SetMarkerStyle(21);
   gre->SetPoint(0,   22.81909,  0.2378564  );   gre->SetPointError(0,0,  0.09283838 );
   gre->SetPoint(1,   57.65211,  0.3014312  );   gre->SetPointError(1,0,  0.08745728 );
   gre->SetPoint(2,   98.19281,  0.4896329  );   gre->SetPointError(2,0,  0.1118034  );
   gre->SetPoint(3,   139.3701,  0.1965057  );   gre->SetPointError(3,0,  0.1595977  );
   gre->SetPoint(4,   178.8885,  0.05771016 );   gre->SetPointError(4,0,  0.1923073  );
   gre->SetPoint(5,   222.6621,  -0.1651759 );   gre->SetPointError(5,0,  0.1640617  );
   gre->SetPoint(6,   265.6899,  -0.3467163 );   gre->SetPointError(6,0,  0.133935   );
   gre->SetPoint(7,   300.0582,  -0.2891963 );   gre->SetPointError(7,0,  0.09282593 );
   gre->SetPoint(8,   339.3341,  -0.1849013 );   gre->SetPointError(8,0,  0.08560304 );

   //multigraph->Add(gre,"");
   gre->Fit(ffit[0]);
   r[0] = gre->Fit(ffit[0],"S");

   multigraph->Draw("a2");     multigraph->Draw("p");
   multigraph->GetXaxis()->SetTitle("#phi [deg.]");
   multigraph->GetYaxis()->SetTitle("A_{LU}^{^{4}He}");
   multigraph->GetXaxis()->CenterTitle(true);
   multigraph->GetYaxis()->CenterTitle(true);
   multigraph->GetYaxis()->SetTitleSize(xfont-0.019);
   multigraph->GetXaxis()->SetTitleSize(xfont-0.019);
   multigraph->GetYaxis()->SetLabelSize(lfont-0.017);
   multigraph->GetXaxis()->SetLabelSize(lfont-0.017);
   multigraph->GetXaxis()->SetTitleOffset(0.9);
   multigraph->GetYaxis()->SetTitleOffset(0.9);
   gPad->Modified();  
   multigraph->GetYaxis()->SetRangeUser(-0.65, 0.65);
   multigraph->GetXaxis()->SetLimits(-5.0,  375);
   gre_sys->Draw("fsame");
   gre_st->Draw("psame");  ffit[0]->Draw("same");
   mg_title->Draw("same");
  
   //Cc->Print("fig_Dec2016/coh_alu_t_0.png");

// ------------>Primitives in pad: c14_2
  
   c12->cd(8); 
   multigraph = new TMultiGraph();

   mg_title = new TLatex(xax,yay,"-t= 0.094 GeV^{2}");
   mg_title->SetTextSize(0.085); 
   mg_title->SetTextFont(62); 
   mg_title->SetTextColor(kBlack);
 
   gre_sys = new TGraph(13);
   gre_sys->SetFillColor(40);
   gre_sys->SetLineColor(40);
   gre_sys->SetMarkerColor(40);
   gre_sys->SetPoint(0,0,    0.0       );
   gre_sys->SetPoint(1,0,    0.02671191);
   gre_sys->SetPoint(2,20,   0.02671191);
   gre_sys->SetPoint(3,60,   0.02842517);
   gre_sys->SetPoint(4,100,  0.03076673);
   gre_sys->SetPoint(5,140,  0.02150343);
   gre_sys->SetPoint(6,180,  0.03603272);
   gre_sys->SetPoint(7,220,  0.01477838);
   gre_sys->SetPoint(8,260,  0.02830972);
   gre_sys->SetPoint(9,300,  0.02864878);
   gre_sys->SetPoint(10,340, 0.02012911);
   gre_sys->SetPoint(11,360, 0.02012911);
   gre_sys->SetPoint(12,360, 0.0       );


   gre_st = new TGraphErrors(9);
   gre_st->SetFillColor(8);
   gre_st->SetMarkerColor(kBlack);   
   gre_st->SetMarkerSize(1.3);   
   gre_st->SetLineWidth(2);
   gre_st->SetMarkerStyle(21);
   gre_st->SetPoint(0,  21.11882, 0.2481552  );      gre_st->SetPointError(0,0, sqrt(pow( 0.09318507,2) - pow( 0.02671191,2) ));
   gre_st->SetPoint(1,  56.44725, 0.3386424  );      gre_st->SetPointError(1,0, sqrt(pow( 0.08283334,2) - pow( 0.02842517,2) ));
   gre_st->SetPoint(2,  97.5705 , 0.346783   );      gre_st->SetPointError(2,0, sqrt(pow( 0.1156562 ,2) - pow( 0.03076673,2) ));
   gre_st->SetPoint(3,  141.9895, 0.1464149  );      gre_st->SetPointError(3,0, sqrt(pow( 0.1887061 ,2) - pow( 0.02150343,2) ));
   gre_st->SetPoint(4,  179.5896, -0.2814579 );      gre_st->SetPointError(4,0, sqrt(pow( 0.1864941 ,2) - pow( 0.03603272,2) ));
   gre_st->SetPoint(5,  218.5645, 0.2104932  );      gre_st->SetPointError(5,0, sqrt(pow( 0.1998783 ,2) - pow( 0.01477838,2) ));
   gre_st->SetPoint(6,  262.8719, -0.239911  );      gre_st->SetPointError(6,0, sqrt(pow( 0.128041  ,2) - pow( 0.02830972,2) ));
   gre_st->SetPoint(7,  303.5719, -0.1989961 );      gre_st->SetPointError(7,0, sqrt(pow( 0.09553188,2) - pow( 0.02864878,2) ));
   gre_st->SetPoint(8,  338.5542, -0.2098659 );      gre_st->SetPointError(8,0, sqrt(pow( 0.08752967,2) - pow( 0.02012911,2) ));
   multigraph->Add(gre_st,"");

   gre = new TGraphErrors(9);
   gre->SetFillColor(1);
   gre->SetMarkerColor(kBlue);
   gre->SetMarkerStyle(21);
   gre->SetPoint(0,  21.11882, 0.2481552   );   gre->SetPointError(0,0, 0.09318507);
   gre->SetPoint(1,  56.44725, 0.3386424   );   gre->SetPointError(1,0, 0.08283334);
   gre->SetPoint(2,  97.5705 , 0.346783    );   gre->SetPointError(2,0, 0.1156562 );
   gre->SetPoint(3,  141.9895, 0.1464149   );   gre->SetPointError(3,0, 0.1887061 );
   gre->SetPoint(4,  179.5896, -0.2814579  );   gre->SetPointError(4,0, 0.1864941 );
   gre->SetPoint(5,  218.5645, 0.2104932   );   gre->SetPointError(5,0, 0.1998783 );
   gre->SetPoint(6,  262.8719, -0.239911   );   gre->SetPointError(6,0, 0.128041  );
   gre->SetPoint(7,  303.5719, -0.1989961  );   gre->SetPointError(7,0, 0.09553188);
   gre->SetPoint(8,  338.5542, -0.2098659  );   gre->SetPointError(8,0, 0.08752967);
 
   
   //multigraph->Add(gre,"");
   gre->Fit(ffit[1]);
   r[1] = gre->Fit(ffit[1],"S");

   multigraph->Draw("a2");    multigraph->Draw("p");
   multigraph->GetXaxis()->SetTitle("#phi [deg.]");
   multigraph->GetXaxis()->CenterTitle(true);
   multigraph->GetXaxis()->SetTitleSize(xfont+0.005);
   multigraph->GetXaxis()->SetLabelSize(lfont);
   multigraph->GetXaxis()->SetTitleOffset(0.8);
   gPad->Modified();  
   multigraph->GetYaxis()->SetRangeUser(-0.65, 0.65);
   multigraph->GetXaxis()->SetLimits(-5.0,  375);
   gre_sys->Draw("fsame"); 
   gre_st->Draw("psame");  ffit[1]->Draw("same");
   mg_title->Draw("same");
  
// ------------>Primitives in pad: c14_3
  
   c12->cd(9); 
   multigraph = new TMultiGraph();
   
   mg_title = new TLatex(xax,yay,"-t= 0.127 GeV^{2}");
   mg_title->SetTextSize(0.085); 
   mg_title->SetTextFont(62); 
   mg_title->SetTextColor(kBlack);

   gre_sys = new TGraph(13);
   gre_sys->SetFillColor(40);
   gre_sys->SetLineColor(40);
   gre_sys->SetMarkerColor(40);
   gre_sys->SetPoint(0,0,    0.0       );
   gre_sys->SetPoint(1,0,    0.02140338);
   gre_sys->SetPoint(2,20,   0.02140338);
   gre_sys->SetPoint(3,60,   0.02045513);
   gre_sys->SetPoint(4,100,  0.0313346 );
   gre_sys->SetPoint(5,140,  0.02038013);
   gre_sys->SetPoint(6,180,  0.0354366 );
   gre_sys->SetPoint(7,220,  0.01898904);
   gre_sys->SetPoint(8,260,  0.0327608 );
   gre_sys->SetPoint(9,300,  0.02802882);
   gre_sys->SetPoint(10,340, 0.02615273);
   gre_sys->SetPoint(11,360, 0.02615273);
   gre_sys->SetPoint(12,360,0.0       );

   
   gre_st = new TGraphErrors(9);
   gre_st->SetFillColor(8);
   gre_st->SetMarkerColor(kBlack);   
   gre_st->SetMarkerSize(1.3);   
   gre_st->SetLineWidth(2);
   gre_st->SetMarkerStyle(21);
   gre_st->SetPoint(0, 22.32266,  0.02807933);   gre_st->SetPointError(0,0, sqrt(pow(  0.09093945 ,2) - pow( 0.02140338,2) ));
   gre_st->SetPoint(1, 60.71096,  0.3582513 );   gre_st->SetPointError(1,0, sqrt(pow(  0.09325026 ,2) - pow( 0.02045513,2) ));
   gre_st->SetPoint(2, 96.88029,  0.2559149 );   gre_st->SetPointError(2,0, sqrt(pow(  0.1269426  ,2) - pow( 0.0313346 ,2) ));
   gre_st->SetPoint(3, 139.9055,  0.2207925 );   gre_st->SetPointError(3,0, sqrt(pow(  0.1930999  ,2) - pow( 0.02038013,2) ));
   gre_st->SetPoint(4, 178.5548,  0.5138574 );   gre_st->SetPointError(4,0, sqrt(pow(  0.1832701  ,2) - pow( 0.0354366 ,2) ));
   gre_st->SetPoint(5, 217.5716,  -0.2473028);   gre_st->SetPointError(5,0, sqrt(pow(  0.1664787  ,2) - pow( 0.01898904,2) ));
   gre_st->SetPoint(6, 260.9501,  -0.2919877);   gre_st->SetPointError(6,0, sqrt(pow(  0.1295829  ,2) - pow( 0.0327608 ,2) ));
   gre_st->SetPoint(7, 303.4526,  -0.2488257);   gre_st->SetPointError(7,0, sqrt(pow(  0.08926991 ,2) - pow( 0.02802882,2) ));
   gre_st->SetPoint(8, 336.8911,  -0.2825422);   gre_st->SetPointError(8,0, sqrt(pow(  0.09019867 ,2) - pow( 0.02615273,2) ));
   multigraph->Add(gre_st,"");

   gre = new TGraphErrors(9);
   gre->SetFillColor(1);
   gre->SetMarkerColor(kBlue);
   gre->SetMarkerStyle(21);
   gre->SetPoint(0,  22.32266,  0.02807933);    gre->SetPointError(0,0,  0.09093945 );
   gre->SetPoint(1,  60.71096,  0.3582513 );    gre->SetPointError(1,0,  0.09325026 );
   gre->SetPoint(2,  96.88029,  0.2559149 );    gre->SetPointError(2,0,  0.1269426  );
   gre->SetPoint(3,  139.9055,  0.2207925 );    gre->SetPointError(3,0,  0.1930999  );
   gre->SetPoint(4,  178.5548,  0.5138574 );    gre->SetPointError(4,0,  0.1832701  );
   gre->SetPoint(5,  217.5716,  -0.2473028);    gre->SetPointError(5,0,  0.1664787  );
   gre->SetPoint(6,  260.9501,  -0.2919877);    gre->SetPointError(6,0,  0.1295829  );
   gre->SetPoint(7,  303.4526,  -0.2488257);    gre->SetPointError(7,0,  0.08926991 );
   gre->SetPoint(8,  336.8911,  -0.2825422);    gre->SetPointError(8,0,  0.09019867 );
   
   //multigraph->Add(gre,"");
   gre->Fit(ffit[2]);
   r[2] = gre->Fit(ffit[2],"S");

   multigraph->Draw("a2");   multigraph->Draw("p");
   multigraph->GetXaxis()->SetTitle("#phi [deg.]");
   multigraph->GetXaxis()->CenterTitle(true);
   multigraph->GetXaxis()->SetTitleSize(xfont+0.005);
   multigraph->GetXaxis()->SetLabelSize(lfont);
   multigraph->GetXaxis()->SetTitleOffset(0.8);
   gPad->Modified();  
   multigraph->GetYaxis()->SetRangeUser(-0.65, 0.65);
   multigraph->GetXaxis()->SetLimits(-5.0,  375);
   gre_sys->Draw("fsame");   
   gre_st->Draw("psame");  ffit[2]->Draw("same");
   mg_title->Draw("same");


   // plot the Imaginary and the real parts of HA CFF vs -t ----------------------
  ccc->cd(1);
  
  //Cc->cd();

   mgCFF = new TMultiGraph();

   Double_t graphCFF_fx11[3] = { 0.08, 0.094, 0.127};
   Double_t graphCFF_fy11[3] = { ffit[0]->GetParameter(0), ffit[1]->GetParameter(0), ffit[2]->GetParameter(0)};
   Double_t graphCFF_fex11[3] = { 0, 0, 0};
   Double_t graphCFF_fey11[3] = { ffit[0]->GetParError(0), ffit[1]->GetParError(0), ffit[2]->GetParError(0)};

   cout<<ffit[0]->GetParameter(0)<<" ,"<< ffit[1]->GetParameter(0)<<" ,"<< ffit[2]->GetParameter(0)<<endl;
   cout<<ffit[0]->GetParError(0)<<" ,"<< ffit[1]->GetParError(0)<<" ,"<< ffit[2]->GetParError(0)<<endl;

   greCFF = new TGraphErrors(3,graphCFF_fx11,graphCFF_fy11,graphCFF_fex11,graphCFF_fey11);
   greCFF->SetFillColor(1);
   greCFF->SetMarkerColor(kBlack);
   greCFF->SetLineColor(kBlack); 
   greCFF->SetMarkerStyle(21);
   greCFF->SetMarkerSize(2.0);
   greCFF->SetLineWidth(2);
 

   mgCFF->Add(greCFF,"");
   mgCFF->Draw("AP");
   mgCFF->GetXaxis()->SetTitle("-t [GeV^{2}]");
   mgCFF->GetYaxis()->SetTitle("Im(H_{A})");
   mgCFF->GetXaxis()->SetLabelSize(lfont-0.03);
   mgCFF->GetYaxis()->SetLabelSize(lfont-0.03);
   mgCFF->GetXaxis()->SetTitleSize(xfont-0.06);
   mgCFF->GetYaxis()->SetTitleSize(xfont-0.06);
   mgCFF->GetYaxis()->SetTitleOffset(0.9);
   mgCFF->GetXaxis()->CenterTitle(true);
   mgCFF->GetYaxis()->CenterTitle(true);
   //mgCFF->GetXaxis()->SetNdivisions(605);
   mgCFF->GetXaxis()->SetLimits(0.05,0.162);
   mgCFF->GetYaxis()->SetRangeUser(-60.0, 78.0);
 
           graphCFF = new TGraph(21);
           graphCFF->SetName("graphCFF");
           graphCFF->SetTitle("graphCFF");
           graphCFF->SetFillColor(1);
           graphCFF->SetLineColor(2);
           graphCFF->SetLineStyle(1);
           graphCFF->SetLineWidth(4);
           graphCFF->SetPoint(0,  0.0500 ,  36.5249 );
           graphCFF->SetPoint(1,  0.0600 ,  31.9928 );
           graphCFF->SetPoint(2,  0.0700 ,  28.0159 );
           graphCFF->SetPoint(3,  0.0800 ,  24.5235 );
           graphCFF->SetPoint(4,  0.0900 ,  21.5069 );
           graphCFF->SetPoint(5,  0.1000 ,  18.8567 );
           graphCFF->SetPoint(6,  0.1100 ,  16.5068 );
           graphCFF->SetPoint(7,  0.1200 ,  14.4442 );
           graphCFF->SetPoint(8,  0.1300 ,  12.6338 );
           graphCFF->SetPoint(9,  0.1400 ,  11.0445 );
           graphCFF->SetPoint(10, 0.1500 ,   9.6493 );
           graphCFF->SetPoint(11, 0.1600 ,   8.4366 );
           graphCFF->SetPoint(12, 0.1700 ,   7.3834 );
           graphCFF->SetPoint(13, 0.1800 ,   6.4557 );
           graphCFF->SetPoint(14, 0.1900 ,   5.6382 );
           graphCFF->SetPoint(15, 0.2000 ,   4.9175 );
           graphCFF->SetPoint(16, 0.2100 ,   4.2736 );
           graphCFF->SetPoint(17, 0.2200 ,   3.7063 );
           graphCFF->SetPoint(18, 0.2300 ,   3.2063 );
           graphCFF->SetPoint(19, 0.2400 ,   2.7655 );
           graphCFF->SetPoint(20, 0.2500 ,   2.3768 );
           graphCFF->Draw("csame");
 
   TGraph *VgraphCFF = new TGraph(6);
           VgraphCFF->SetFillColor(kBlack);
           VgraphCFF->SetLineColor(kBlack);
           VgraphCFF->SetLineStyle(7);
           VgraphCFF->SetLineWidth(4);
           VgraphCFF->SetPoint(0,  0.05 , 30.7785 );
           VgraphCFF->SetPoint(1,  0.07 , 22.7538 );
           VgraphCFF->SetPoint(2,  0.09 , 16.8293 );
           VgraphCFF->SetPoint(3,  0.11 , 12.4511 );
           VgraphCFF->SetPoint(4,  0.13 , 9.4857  );
           VgraphCFF->SetPoint(5,  0.15 , 6.81138 );
           VgraphCFF->Draw("csame");

   TGraph*SgraphCFF = new TGraph(31);
          SgraphCFF->SetName("graphCFF");
          SgraphCFF->SetTitle("graphCFF");
          SgraphCFF->SetFillColor(1);
          SgraphCFF->SetLineColor(4);
          SgraphCFF->SetLineStyle(9);
          SgraphCFF->SetLineWidth(4);
          SgraphCFF->SetPoint(0,  0.0385,  16.59 );
          SgraphCFF->SetPoint(1,  0.0435,  22.46 );
          SgraphCFF->SetPoint(2,  0.0485,  24.85 );
          SgraphCFF->SetPoint(3,  0.0535,  24.76 );
          SgraphCFF->SetPoint(4,  0.0585,  23.29 );
          SgraphCFF->SetPoint(5,  0.0635,  22.29 );
          SgraphCFF->SetPoint(6,  0.0685,  20.98 );
          SgraphCFF->SetPoint(7,  0.0735,  19.85 );
          SgraphCFF->SetPoint(8,  0.0785,  18.68 );
          SgraphCFF->SetPoint(9,  0.0835,  17.58 );
          SgraphCFF->SetPoint(10, 0.0885,  16.59 );
          SgraphCFF->SetPoint(11, 0.0935,  15.61 );
          SgraphCFF->SetPoint(12, 0.0985,  14.69 );
          SgraphCFF->SetPoint(13, 0.1035,  13.82 );
          SgraphCFF->SetPoint(14, 0.1085,  13.04 );
          SgraphCFF->SetPoint(15, 0.1135,  12.26 );
          SgraphCFF->SetPoint(16, 0.1185,  11.54 );
          SgraphCFF->SetPoint(17, 0.1235,  10.86 );
          SgraphCFF->SetPoint(18, 0.1285,  10.21 );
          SgraphCFF->SetPoint(19, 0.1335,  9.625 );
          SgraphCFF->SetPoint(20, 0.1385,  9.054 );
          SgraphCFF->SetPoint(21, 0.1435,  8.516 );
          SgraphCFF->SetPoint(22, 0.1485,  8.009 );
          SgraphCFF->SetPoint(23, 0.1535,  7.531 );
          SgraphCFF->SetPoint(24, 0.1585,  7.081 );
          SgraphCFF->SetPoint(25, 0.1635,  6.667 );
          SgraphCFF->SetPoint(26, 0.1685,  6.267 );
          SgraphCFF->SetPoint(27, 0.1735,  5.890 );
          SgraphCFF->SetPoint(28, 0.1785,  5.534 );
          SgraphCFF->SetPoint(29, 0.1835,  5.199 );
          SgraphCFF->SetPoint(30, 0.1885,  4.883 );
          SgraphCFF->Draw("csame");         

       
   TLegend *leg_cff_Im = new TLegend(0.55,0.75,0.9,0.98);
            leg_cff_Im->SetNColumns(1);  
            //leg_cff_Im->AddEntry(greCFF,"CLAS-EG6 (Q^{2} = 1.492, x_{B} = 0.177)","P");
            //leg_cff_Im->AddEntry(graphCFF,"Off-shell cal.(Q^{2} = 1.492, x_{B} = 0.177)","L");
            //leg_cff_Im->AddEntry(VgraphCFF,"Convolution model (Q^{2} = 1.492, x_{B} = 0.177)","L");
            //leg_cff_Im->AddEntry(SgraphCFF,"On-shell cal. (Q^{2} = 0.3, x_{B} = 0.177)","L");
            leg_cff_Im->SetTextSize(0.04);
            leg_cff_Im->AddEntry(greCFF,"CLAS-EG6","P");
            leg_cff_Im->AddEntry(graphCFF,"Convolution-Dual","L");
            leg_cff_Im->AddEntry(VgraphCFF,"Convolution-VGG","L");
            leg_cff_Im->AddEntry(SgraphCFF,"Off-shell model","L");
            leg_cff_Im->Draw(); 

   //Cc->Print("fig_Dec2016/coh_im_t.png");

   // real part as a function of -t
   ccc->cd(2);
   //Cc->cd();

   mgCFF = new TMultiGraph();
   
   Double_t graphCFF_fx12[3] = { 0.08, 0.094, 0.127};
   Double_t graphCFF_fy12[3] = {ffit[0]->GetParameter(1),
                                ffit[1]->GetParameter(1), 
                                ffit[2]->GetParameter(1)};   
   Double_t graphCFF_fex12[3] = { 0, 0, 0};
   Double_t graphCFF_fey12[3] = {ffit[0]->GetParError(1), 
                                 ffit[1]->GetParError(1),
                                 ffit[2]->GetParError(1)};

   cout<<ffit[0]->GetParameter(1)<<" ,"<< ffit[1]->GetParameter(1)<<" ,"<< ffit[2]->GetParameter(1)<<endl;
   cout<<ffit[0]->GetParError(1)<<" ,"<< ffit[1]->GetParError(1)<<" ,"<< ffit[2]->GetParError(1)<<endl;




   greCFF = new TGraphErrors(3,graphCFF_fx12,graphCFF_fy12,graphCFF_fex12,graphCFF_fey12);
   greCFF->SetFillColor(1);
   greCFF->SetMarkerColor(kBlack);
   greCFF->SetLineColor(kBlack); 
   greCFF->SetMarkerStyle(21);
   greCFF->SetMarkerSize(2.0);
   greCFF->SetLineWidth(2);
 
   mgCFF->Add(greCFF,"");

   mgCFF->Draw("AP");
   mgCFF->GetXaxis()->SetTitle("-t [GeV^{2}]");
   mgCFF->GetYaxis()->SetTitle("Re(H_{A})");
   mgCFF->GetXaxis()->SetLabelSize(lfont-0.03);
   mgCFF->GetYaxis()->SetLabelSize(lfont-0.03);
   mgCFF->GetXaxis()->SetTitleSize(xfont-0.06);
   mgCFF->GetYaxis()->SetTitleSize(xfont-0.06);
   mgCFF->GetYaxis()->SetTitleOffset(0.9);
   mgCFF->GetXaxis()->CenterTitle(true);
   mgCFF->GetYaxis()->CenterTitle(true);
   mgCFF->GetYaxis()->SetRangeUser(-60.0, 65.0);
   mgCFF->GetXaxis()->SetLimits(0.05,0.162);
  

   graphCFF = new TGraph(21);
   graphCFF->SetName("graphCFF");
   graphCFF->SetTitle("graphCFF");
   graphCFF->SetFillColor(1);
   graphCFF->SetLineColor(2);
   graphCFF->SetLineStyle(1);
   graphCFF->SetLineWidth(4);
   graphCFF->SetPoint(0,   0.0500,  -6.6015 );
   graphCFF->SetPoint(1,   0.0600,  -5.9884 );
   graphCFF->SetPoint(2,   0.0700,  -5.4292 );
   graphCFF->SetPoint(3,   0.0800,  -4.9116 );
   graphCFF->SetPoint(4,   0.0900,  -4.4316 );
   graphCFF->SetPoint(5,   0.1000,  -3.9968 );
   graphCFF->SetPoint(6,   0.1100,  -3.6097 );
   graphCFF->SetPoint(7,   0.1200,  -3.2584 );
   graphCFF->SetPoint(8,   0.1300,  -2.9396 );
   graphCFF->SetPoint(9,   0.1400,  -2.6504 );
   graphCFF->SetPoint(10,  0.1500,  -2.3879 );
   graphCFF->SetPoint(11,  0.1600,  -2.1512 );
   graphCFF->SetPoint(12,  0.1700,  -1.9262 );
   graphCFF->SetPoint(13,  0.1800,  -1.7230 );
   graphCFF->SetPoint(14,  0.1900,  -1.5395 );
   graphCFF->SetPoint(15,  0.2000,  -1.3737 );
   graphCFF->SetPoint(16,  0.2100,  -1.2269 );
   graphCFF->SetPoint(17,  0.2200,  -1.0936 );
   graphCFF->SetPoint(18,  0.2300,  -0.9723 );
   graphCFF->SetPoint(19,  0.2400,  -0.8620 );
   graphCFF->SetPoint(20,  0.2500,  -0.7615 );
   graphCFF->Draw("csame");

   // from the VGG vs. -t 
   ggraphCFF = new TGraph(6);
   ggraphCFF->SetFillColor(kBlack);
   ggraphCFF->SetLineColor(kBlack);
   ggraphCFF->SetLineStyle(7);
   ggraphCFF->SetLineWidth(4);
   ggraphCFF->SetPoint(0,  0.05 , -0.958994 );
   ggraphCFF->SetPoint(1,  0.07 , -1.372    );
   ggraphCFF->SetPoint(2,  0.09 , -1.50058  );
   ggraphCFF->SetPoint(3,  0.11 , -1.46658  );
   ggraphCFF->SetPoint(4,  0.13 , -1.24124  );
   ggraphCFF->SetPoint(5,  0.15 , -1.18798  );
   ggraphCFF->Draw("csame");
 

   TLegend *leg_cff_Re = new TLegend(0.55,0.78,0.9,0.98);
            leg_cff_Re->SetNColumns(1);
            leg_cff_Re->SetTextSize(0.04);  
            leg_cff_Re->AddEntry(greCFF,"CLAS-EG6","P");   
            leg_cff_Re->AddEntry(graphCFF,"Convolution-Dual","L");
            leg_cff_Re->AddEntry(ggraphCFF,"Convolution-VGG","L");
          
   leg_cff_Re->Draw();
   //Cc->Print("fig_Dec2016/coh_re_t.png");

   // plot the coherent ALU at phi = 90 deg vs -t  -------------------------------
   Can->cd(3);
   
   //Cc->cd();

   TH1 *HALU_He_vs_t = new TH2F("HALU_He_vs_t","",100, -0.005 , 0.25,100,-0.1,0.72);
        HALU_He_vs_t->SetStats(0);
        HALU_He_vs_t->GetXaxis()->SetTitle("-t [GeV^{2}]");
        //HALU_He_vs_t->GetYaxis()->SetTitle("A_{LU}^{^{4}He}(90#circ)");
        HALU_He_vs_t->GetYaxis()->SetTitleSize(xfont);
        HALU_He_vs_t->GetXaxis()->SetTitleSize(xfont-0.05);
        HALU_He_vs_t->GetXaxis()->SetNdivisions(205);
        HALU_He_vs_t->GetXaxis()->CenterTitle(true);
        HALU_He_vs_t->GetYaxis()->CenterTitle(true);
        HALU_He_vs_t->GetYaxis()->SetNdivisions(414);
        HALU_He_vs_t->GetXaxis()->SetTitleOffset(0.9);
        //HALU_He_vs_t->GetYaxis()->SetTitleOffset(1.5);
        HALU_He_vs_t->Draw("");
   
   
   TGraphErrors *hgre = new TGraphErrors(2);
                 hgre->SetName("HERMES");
                 hgre->SetTitle("HERMES");
                 hgre->SetFillColor(1);
                 hgre->SetLineColor(8);
                 hgre->SetLineWidth(3);
                 hgre->SetMarkerColor(8);
                 hgre->SetMarkerStyle(8);
                 hgre->SetMarkerSize(1.5);
                 hgre->SetPoint(0,0.027,0.249);         hgre->SetPointError(0,0,0.107+0.011);
                 hgre->SetPoint(1,0.105,0.280);         hgre->SetPointError(1,0,0.130+0.016);
                 //   hgre->SetPoint(2,0.201,0.055);       hgre->SetPointError(2,0,0.176+0.01);
                 //   hgre->SetPoint(3,0.418,-0.009);      hgre->SetPointError(3,0,0.249+0.015);
   hgre->Draw("psame");

    for(int i=0; i<3; i++){
        double x[1] = { 90.0};  double err[1];
        r[i]->GetConfidenceIntervals(1, 1, 1, x, err, 0.683, true);
        ALU_fit[i]     = ffit[i]->Eval(x[0]);
        ALU_fit_err[i] = err[0];
        cout<< ALU_fit[i]<<"     "<<ALU_fit_err[i]<<endl;
   }

   me_gre = new TGraphErrors(3);
   me_gre->SetName("Graph0");
   me_gre->SetTitle("This Work (<x_{B}>=0.177, <Q^{2}>= 1.492)");
   me_gre->SetMarkerColor(kBlack);
   me_gre->SetFillColor(1);
   me_gre->SetLineWidth(3);
   me_gre->SetMarkerSize(1.5);
   me_gre->SetMarkerStyle(21);
   me_gre->SetPoint(0,0.08, ALU_fit[0] /*0.369816*/);   me_gre->SetPointError(0,0, ALU_fit_err[0]/*0.0526035*/);
   me_gre->SetPoint(1,0.094,ALU_fit[1] /*0.205635*/);   me_gre->SetPointError(1,0, ALU_fit_err[1]/*0.0578671*/);
   me_gre->SetPoint(2,0.127,ALU_fit[2] /*0.318257*/);   me_gre->SetPointError(2,0, ALU_fit_err[2]/*0.119681*/);
   me_gre->Draw("psame"); 


   sys_gre = new TGraph();
   sys_gre->SetFillColor(40);
   sys_gre->SetPoint(0,0.06,0);
   sys_gre->SetPoint(1,0.06,0.101*0.405/(1+0.251));
   sys_gre->SetPoint(2,0.08,0.101*0.405/(1+0.251));
   sys_gre->SetPoint(3,0.094,0.101*0.2526/(1-0.1748));
   sys_gre->SetPoint(4,0.127,0.101*0.269/(1-0.214));
   sys_gre->SetPoint(4,0.2,0.101*0.269/(1-0.214));
   sys_gre->SetPoint(5,0.2,0.);
   sys_gre->Draw("fsame"); 



   graph = new TGraph(16);
   graph->SetName("Graph");
   graph->SetTitle("Graph");
   graph->SetFillColor(1);
   graph->SetLineColor(2);
   graph->SetLineStyle(1);
   graph->SetLineWidth(3);
   graph->SetPoint(0,0.0026,0.0997298);
   graph->SetPoint(1,0.0139,0.100293);
   graph->SetPoint(2,0.0336,0.135753);
   graph->SetPoint(3,0.0611,0.198727);
   graph->SetPoint(4,0.0955,0.240483);
   graph->SetPoint(5,0.1355,0.266514);
   graph->SetPoint(6,0.1796,0.281026);
   graph->SetPoint(7,0.2262,0.287221);
   graph->SetPoint(8,0.2738,0.287302);
   graph->SetPoint(9,0.3204,0.282924);
   graph->SetPoint(10,0.3645,0.275682);
   graph->SetPoint(11,0.4045,0.267646);
   graph->SetPoint(12,0.4389,0.259709);
   graph->SetPoint(13,0.4664,0.252901);
   graph->SetPoint(14,0.4861,0.248061);
   graph->SetPoint(15,0.4974,0.245101);
   graph->Draw("csame");
   
   graph = new TGraph(16);
   graph->SetName("Graph");
   graph->SetTitle("Graph");
   graph->SetFillColor(1);
   graph->SetLineColor(4);
   graph->SetLineStyle(7);
   graph->SetLineWidth(3);
   graph->SetPoint(0,0.0026,0.0623235);
   graph->SetPoint(1,0.0139,0.0614373);
   graph->SetPoint(2,0.0336,0.0604637);
   graph->SetPoint(3,0.0611,0.059663);
   graph->SetPoint(4,0.0955,0.0589742);
   graph->SetPoint(5,0.1355,0.0581058);
   graph->SetPoint(6,0.1796,0.104233);
   graph->SetPoint(7,0.2262,0.134905);
   graph->SetPoint(8,0.2738,0.146625);
   graph->SetPoint(9,0.3204,0.148695);
   graph->SetPoint(10,0.3645,0.146013);
   graph->SetPoint(11,0.4045,0.141291);
   graph->SetPoint(12,0.4389,0.135987);
   graph->SetPoint(13,0.4664,0.13141);
   graph->SetPoint(14,0.4861,0.127973);
   graph->SetPoint(15,0.4974,0.125959);
   //graph->Draw("csame");
   

   TGraphErrors *gre_LU_1 = new TGraphErrors(3);
   gre_LU_1->SetName("");
   gre_LU_1->SetTitle("LT (x_{B}=0.137, Q^{2}= 1.2)");
   gre_LU_1->SetFillColor(1);
   gre_LU_1->SetLineColor(2);
   gre_LU_1->SetLineWidth(3);
   gre_LU_1->SetMarkerColor(2);
   gre_LU_1->SetMarkerStyle(21);
   gre_LU_1->SetMarkerSize(1.5);
   gre_LU_1->SetPoint(0,0.0996421,0.7);
   gre_LU_1->SetPointError(0,0,0.00603854);
   gre_LU_1->SetPoint(1,0.124067,0.7);
   gre_LU_1->SetPointError(1,0,0.00765358);
   gre_LU_1->SetPoint(2,0.206561,0.7);
   gre_LU_1->SetPointError(2,0,0.0121847);
   //gre_LU_1->Draw("psame");


   TGraphErrors *gre_LU_2 = new TGraphErrors(3);
   gre_LU_2->SetName("LT (x_{B}=0.307)");
   gre_LU_2->SetTitle("LT (x_{B}=0.307, Q^{2}= 2.1)");
   gre_LU_2->SetLineStyle(7);
   gre_LU_2->SetFillColor(1);
   gre_LU_2->SetLineColor(4);
   gre_LU_2->SetLineWidth(3);
   gre_LU_2->SetMarkerColor(4);
   gre_LU_2->SetMarkerStyle(22);
   gre_LU_2->SetMarkerSize(1.5);
   gre_LU_2->SetPoint(0,0.0999352,0.6);
   gre_LU_2->SetPointError(0,0,0.00865828);
   gre_LU_2->SetPoint(1,0.125002,0.6);
   gre_LU_2->SetPointError(1,0,0.00790532);
   gre_LU_2->SetPoint(2,0.208572,0.6);
   gre_LU_2->SetPointError(2,0,0.00821632);
   //gre_LU_2->Draw("psame");
   
   gre = new TGraphErrors(4);
   gre->SetName("HERMES");
   gre->SetTitle("HERMES");
   gre->SetFillColor(1);
   gre->SetLineColor(8);
   gre->SetLineWidth(3);
   gre->SetMarkerColor(8);
   gre->SetMarkerStyle(8);
   gre->SetMarkerSize(1.5);
   gre->SetPoint(0,0.02,0.21);
   gre->SetPointError(0,0,0.12);
   gre->SetPoint(1,0.1,0.15);
   gre->SetPointError(1,0,0.15);
   gre->SetPoint(2,0.2,0.3);
   gre->SetPointError(2,0,0.2);
   gre->SetPoint(3,0.42,0.1);
   gre->SetPointError(3,0,0.26);
   //gre->Draw("psame");

   
  TLine *L_zero = new TLine(); 
         L_zero->SetLineColor(kBlack);
         L_zero->SetLineWidth(1);
         L_zero->SetLineStyle(7);
   L_zero->DrawLine(0, 0, 0.25, 0);

   TLegend* leg_t = new TLegend(0.2,0.8,0.98,1.0);
            leg_t-> SetNColumns(1);
            leg_t->SetFillColor(0);
            leg_t->SetTextSize(0.05);
            leg_t->AddEntry(me_gre,"x_{B}= 0.17, Q^{2}= 1.5 GeV^{2}","P");
            leg_t->AddEntry(hgre,"x_{B}= 0.07, Q^{2}= 1.8 GeV^{2}","P");
            leg_t->AddEntry(gre_LU_1,"x_{B}= 0.14, Q^{2}= 1.2 GeV^{2}","L");
            //leg_t->AddEntry(gre_LU_2,"Off-shell cal.: x_{B}= 0.307, Q^{2}= 2.1","L");
            //leg_t->AddEntry(alpha_tt,"CLAS-E1DVCS1 free p");

   leg_t->Draw();

   //Cc->Print("fig_Dec2016/coh_alu_t_90.png");

   // plot the coherent ALU ratio at phi = 90 deg vs -t  -------------------------
   Can_r->cd(3);
   //Cc->cd();

   TH1 *HALU_RATIO_t = new TH2F("HALU_RATIO_t","",100,-0.02,0.51,100,y_init, y_fin);
        HALU_RATIO_t->SetStats(0);
        HALU_RATIO_t->GetXaxis()->SetTitle("-t [GeV^{2}]");
        HALU_RATIO_t->GetXaxis()->SetNdivisions(205);
        HALU_RATIO_t->GetYaxis()->SetNdivisions(511);
        HALU_RATIO_t->GetYaxis()->SetTitleSize(xfont);
        HALU_RATIO_t->GetXaxis()->SetTitleSize(xfont);
        HALU_RATIO_t->GetYaxis()->SetTitle("A_{LU}^{^{4}He}/A_{LU}^{p} (90^{#circ})");
        HALU_RATIO_t->GetXaxis()->CenterTitle(true);
        HALU_RATIO_t->GetYaxis()->CenterTitle(true);
        HALU_RATIO_t->GetYaxis()->SetTitleOffset(1.5);
        HALU_RATIO_t->GetXaxis()->SetTitleOffset(1.5);
        HALU_RATIO_t->Draw("");

   // CLAS-EG6 result
   TGraphErrors *gre_R_t = new TGraphErrors(1);
                 gre_R_t->SetName("This Work (x_{B} = 0.177");
                 gre_R_t->SetTitle("This Work (x_{B} = 0.177");
                 gre_R_t->SetFillColor(1);
                 gre_R_t->SetLineWidth(3);
                 gre_R_t->SetLineColor(kBlack);
                 gre_R_t->SetMarkerColor(kBlack);
                 gre_R_t->SetMarkerStyle(21);
                 gre_R_t->SetMarkerSize(1.6);
                 gre_R_t->SetPoint(0,0.1,  2.214);   gre_R_t->SetPointError(0,0,  0.515789);
                       
   // HERMES measuremenRATIO_tt
   TGraphErrors *gre_H_R_t = new TGraphErrors(1);
                 gre_H_R_t->SetName("HERMES (x_{B} = 0.072)");
                 gre_H_R_t->SetTitle("HERMES (x_{B} = 0.072)");
                 gre_H_R_t->SetFillColor(1);
                 gre_H_R_t->SetLineColor(8);
                 gre_H_R_t->SetLineWidth(3);
                 gre_H_R_t->SetMarkerColor(8);
                 gre_H_R_t->SetMarkerStyle(8);
                 gre_H_R_t->SetMarkerSize(1.5);
                 gre_H_R_t->SetPoint(0,0.018,0.91);         gre_H_R_t->SetPointError(0,0,0.19);
                 gre_H_R_t->Draw("psame");

  // calculations from vadim
  TGraph *gre_vd1_R_t = new TGraph(25);
          gre_vd1_R_t->SetName("Guzey (x_{B}=0.137)");
          gre_vd1_R_t->SetTitle("Guzey (x_{B}=0.137)");
          gre_vd1_R_t->SetMarkerColor(kBlack);
          gre_vd1_R_t->SetLineColor(6);
          gre_vd1_R_t->SetLineWidth(3);
          gre_vd1_R_t->SetMarkerSize(1.2);
          gre_vd1_R_t->SetMarkerStyle(21);
          gre_vd1_R_t->SetPoint(0,     0.0770,    1.6778   );
          gre_vd1_R_t->SetPoint(1,     0.0940,    1.5977   );
          gre_vd1_R_t->SetPoint(2,     0.1110,    1.5296   );
          gre_vd1_R_t->SetPoint(3,     0.1280,    1.4705   );
          gre_vd1_R_t->SetPoint(4,     0.1450,    1.4184   );
          gre_vd1_R_t->SetPoint(5,     0.1620,    1.3729   );
          gre_vd1_R_t->SetPoint(6,     0.1790,    1.3324   );
          gre_vd1_R_t->SetPoint(7,     0.1960,    1.2967   );
          gre_vd1_R_t->SetPoint(8,     0.2130,    1.2649   );
          gre_vd1_R_t->SetPoint(9,     0.2300,    1.2364   );
          gre_vd1_R_t->SetPoint(10,    0.2470,    1.2108   );
          gre_vd1_R_t->SetPoint(11,    0.2640,    1.1893   );
          gre_vd1_R_t->SetPoint(12,    0.2810,    1.1702   );
          gre_vd1_R_t->SetPoint(13,    0.2980,    1.1534   );
          gre_vd1_R_t->SetPoint(14,    0.3150,    1.1377   );
          gre_vd1_R_t->SetPoint(15,    0.3320,    1.1235   );
          gre_vd1_R_t->SetPoint(16,    0.3490,    1.1104   );
          gre_vd1_R_t->SetPoint(17,    0.3660,    1.1018   );
          gre_vd1_R_t->SetPoint(18,    0.3830,    1.0931   );
          gre_vd1_R_t->SetPoint(19,    0.4000,    1.0856   );
          gre_vd1_R_t->SetPoint(20,    0.4170,    1.0779   );
          gre_vd1_R_t->SetPoint(21,    0.4340,    1.0707   );
          gre_vd1_R_t->SetPoint(22,    0.4510,    1.0652   );
          gre_vd1_R_t->SetPoint(23,    0.4680,    1.0613   );
          gre_vd1_R_t->SetPoint(24,    0.4850,    1.0581   );
          gre_vd1_R_t->Draw("csame");

 
  TGraph *gre_vd2_R_t = new TGraph(25);
          gre_vd2_R_t->SetName("Guzey (x_{B}=0.177)");
          gre_vd2_R_t->SetTitle("Guzey (x_{B}=0.177)");
          gre_vd2_R_t->SetMarkerColor(kBlack);
          gre_vd2_R_t->SetFillColor(28);
          gre_vd2_R_t->SetLineColor(28);
          gre_vd2_R_t->SetLineWidth(3);
          gre_vd2_R_t->SetMarkerSize(1.2);
          gre_vd2_R_t->SetMarkerStyle(22);
          gre_vd2_R_t->SetPoint(0,     0.0770,    1.8461  ); 
          gre_vd2_R_t->SetPoint(1,     0.0940,    1.7800  ); 
          gre_vd2_R_t->SetPoint(2,     0.1110,    1.7159  ); 
          gre_vd2_R_t->SetPoint(3,     0.1280,    1.6540  ); 
          gre_vd2_R_t->SetPoint(4,     0.1450,    1.5949  ); 
          gre_vd2_R_t->SetPoint(5,     0.1620,    1.5415  ); 
          gre_vd2_R_t->SetPoint(6,     0.1790,    1.4916  ); 
          gre_vd2_R_t->SetPoint(7,     0.1960,    1.4448  ); 
          gre_vd2_R_t->SetPoint(8,     0.2130,    1.3998  ); 
          gre_vd2_R_t->SetPoint(9,     0.2300,    1.3572  ); 
          gre_vd2_R_t->SetPoint(10,    0.2470,    1.3173  ); 
          gre_vd2_R_t->SetPoint(11,    0.2640,    1.2844  ); 
          gre_vd2_R_t->SetPoint(12,    0.2810,    1.2527  ); 
          gre_vd2_R_t->SetPoint(13,    0.2980,    1.2235  ); 
          gre_vd2_R_t->SetPoint(14,    0.3150,    1.1950  ); 
          gre_vd2_R_t->SetPoint(15,    0.3320,    1.1681  ); 
          gre_vd2_R_t->SetPoint(16,    0.3490,    1.1429  ); 
          gre_vd2_R_t->SetPoint(17,    0.3660,    1.1237  ); 
          gre_vd2_R_t->SetPoint(18,    0.3830,    1.1048  ); 
          gre_vd2_R_t->SetPoint(19,    0.4000,    1.0872  ); 
          gre_vd2_R_t->SetPoint(20,    0.4170,    1.0694  ); 
          gre_vd2_R_t->SetPoint(21,    0.4340,    1.0524  ); 
          gre_vd2_R_t->SetPoint(22,    0.4510,    1.0382  ); 
          gre_vd2_R_t->SetPoint(23,    0.4680,    1.0260  ); 
          gre_vd2_R_t->SetPoint(24,    0.4850,    1.0147  ); 
          gre_vd2_R_t->Draw("csame");
 
 
 
  TGraph *gre_vd3_R_t = new TGraph(25);
          gre_vd3_R_t->SetName("Guzey (x_{B}=0.177)");
          gre_vd3_R_t->SetTitle("Guzey (x_{B}=0.177)");
          gre_vd3_R_t->SetMarkerColor(kBlack);
          gre_vd3_R_t->SetLineColor(1);
          gre_vd3_R_t->SetLineWidth(3);
          gre_vd3_R_t->SetMarkerSize(1.2);
          gre_vd3_R_t->SetMarkerStyle(20);
          gre_vd3_R_t->SetPoint(0,     0.0770,   1.7699   );
          gre_vd3_R_t->SetPoint(1,     0.0940,   1.4384   );
          gre_vd3_R_t->SetPoint(2,     0.1110,   1.3341   );
          gre_vd3_R_t->SetPoint(3,     0.1280,   1.2712   );
          gre_vd3_R_t->SetPoint(4,     0.1450,   1.2253   );
          gre_vd3_R_t->SetPoint(5,     0.1620,   1.1892   );
          gre_vd3_R_t->SetPoint(6,     0.1790,   1.1600   );
          gre_vd3_R_t->SetPoint(7,     0.1960,   1.1363   );
          gre_vd3_R_t->SetPoint(8,     0.2130,   1.1166   );
          gre_vd3_R_t->SetPoint(9,     0.2300,   1.1001   );
          gre_vd3_R_t->SetPoint(10,    0.2470,   1.0861   );
          gre_vd3_R_t->SetPoint(11,    0.2640,   1.0753   );
          gre_vd3_R_t->SetPoint(12,    0.2810,   1.0663   );
          gre_vd3_R_t->SetPoint(13,    0.2980,   1.0589   );
          gre_vd3_R_t->SetPoint(14,    0.3150,   1.0524   );
          gre_vd3_R_t->SetPoint(15,    0.3320,   1.0469   );
          gre_vd3_R_t->SetPoint(16,    0.3490,   1.0423   );
          gre_vd3_R_t->SetPoint(17,    0.3660,   1.0403   );
          gre_vd3_R_t->SetPoint(18,    0.3830,   1.0387   );
          gre_vd3_R_t->SetPoint(19,    0.4000,   1.0377   );
          gre_vd3_R_t->SetPoint(20,    0.4170,   1.0364   );
          gre_vd3_R_t->SetPoint(21,    0.4340,   1.0355   );
          gre_vd3_R_t->SetPoint(22,    0.4510,   1.0357   );
          gre_vd3_R_t->SetPoint(23,    0.4680,   1.0370   );
          gre_vd3_R_t->SetPoint(24,    0.4850,   1.0387   );
          gre_vd3_R_t->Draw("csame");
               
   // calculations from simonetta 
   TGraph *gre_LT1_R_t = new TGraph(16);
           gre_LT1_R_t->SetName("LT (x_{B}=0.137)");
           gre_LT1_R_t->SetTitle("LT (x_{B}=0.137)");
           gre_LT1_R_t->SetFillColor(1);
           gre_LT1_R_t->SetLineColor(2);
           gre_LT1_R_t->SetLineStyle(1);
           gre_LT1_R_t->SetLineWidth(3);
           gre_LT1_R_t->SetPoint(0,0.0026,1.06064);
           gre_LT1_R_t->SetPoint(1,0.0139,1.06662);
           gre_LT1_R_t->SetPoint(2,0.0336,1.07201);
           gre_LT1_R_t->SetPoint(3,0.0611,1.07401);
           gre_LT1_R_t->SetPoint(4,0.0955,1.07369);
           gre_LT1_R_t->SetPoint(5,0.1355,1.07567);
           gre_LT1_R_t->SetPoint(6,0.1796,1.08403);
           gre_LT1_R_t->SetPoint(7,0.2262,1.09922);
           gre_LT1_R_t->SetPoint(8,0.2738,1.11897);
           gre_LT1_R_t->SetPoint(9,0.3204,1.13984);
           gre_LT1_R_t->SetPoint(10,0.3645,1.15933);
           gre_LT1_R_t->SetPoint(11,0.4045,1.17837);
           gre_LT1_R_t->SetPoint(12,0.4389,1.19453);
           gre_LT1_R_t->SetPoint(13,0.4664,1.20748);
           gre_LT1_R_t->SetPoint(14,0.4861,1.21789);
           gre_LT1_R_t->SetPoint(15,0.4974,1.22326);
           gre_LT1_R_t->Draw("csame");

           
           gre_LT2_R_t = new TGraph(16);
           gre_LT2_R_t->SetName("LT (x_{B}=0.307)");
           gre_LT2_R_t->SetTitle("LT (x_{B}=0.307)");
           gre_LT2_R_t->SetFillColor(1);
           gre_LT2_R_t->SetLineColor(4);
           gre_LT2_R_t->SetLineStyle(1);
           gre_LT2_R_t->SetLineWidth(3);
           gre_LT2_R_t->SetPoint(0,0.0026,0.909493);
           gre_LT2_R_t->SetPoint(1,0.0139,0.896561);
           gre_LT2_R_t->SetPoint(2,0.0336,0.882353);
           gre_LT2_R_t->SetPoint(3,0.0611,0.870668);
           gre_LT2_R_t->SetPoint(4,0.0955,0.860616);
           gre_LT2_R_t->SetPoint(5,0.1355,0.847944);
           gre_LT2_R_t->SetPoint(6,0.1796,0.833209);
           gre_LT2_R_t->SetPoint(7,0.2262,0.81447);
           gre_LT2_R_t->SetPoint(8,0.2738,0.798337);
           gre_LT2_R_t->SetPoint(9,0.3204,0.782492);
           gre_LT2_R_t->SetPoint(10,0.3645,0.768573);
           gre_LT2_R_t->SetPoint(11,0.4045,0.756931);
           gre_LT2_R_t->SetPoint(12,0.4389,0.746667);
           gre_LT2_R_t->SetPoint(13,0.4664,0.739485);
           gre_LT2_R_t->SetPoint(14,0.4861,0.734579);
           gre_LT2_R_t->SetPoint(15,0.4974,0.731844);
           gre_LT2_R_t->Draw("csame");
   
   gre_R_t->Draw("psame");

   line = new TLine(-0.02,1,0.5,1);
   line->SetLineWidth(2);
   line->SetLineStyle(4);
   line->Draw();

   TLegend* leg_R_t = new TLegend(0.26,0.75,0.95,0.95);
            leg_R_t-> SetNColumns(2);
            leg_R_t->AddEntry(gre_R_t,"CLAS-EG6: <x_{B}> = 0.177","P");
            leg_R_t->AddEntry(gre_H_R_t,"HERMES: <x_{B}> = 0.072","P");
            leg_R_t->AddEntry(gre_LT1_R_t,"Off-shell cal.: x_{B}=0.137","L");
            leg_R_t->AddEntry(gre_LT2_R_t,"Off-shell cal.: x_{B}=0.307","L");
            leg_R_t->AddEntry(gre_vd1_R_t,"On-shell cal.: x_{B} = 0.137","L");
            leg_R_t->AddEntry(gre_vd2_R_t,"On-shell cal.: x_{B} = 0.177","L");
            leg_R_t->AddEntry(gre_vd3_R_t,"On-shell cal.: x_{B} = 0.25","L");

   leg_R_t->Draw();
   //Cc->Print("fig_Dec2016/coh_alu_t_90_ratio.png");

  ccc   ->SaveAs("fig/Coherent_CFF.pdf");
  ccc   ->SaveAs("fig/Coherent_CFF.png");
  ccc   ->SaveAs("fig/Coherent_CFF.C");




}
