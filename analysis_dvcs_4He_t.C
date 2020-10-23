#define analysis_dvcs_4He_t_cxx
#include "analysis_dvcs_4He_t.h"
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
#include <TGraphErrors.h>
#include <TLatex.h>
#include "TGraph.h"
#include<TLine.h>
#include <TLorentzVector.h>
#include "TSystem.h"
#include "TLegend.h"
#include<algorithm>
#include<vector>
#include<vector>
#include "plot_ALU_projections.C"
#include "plot_CFF_projections.C"
#include "FunProfile.C"

using namespace std;

void Find_CFF( double xB, double t, double &Im, double &Re);
void Calculate_CFF(double Q2, double xB, double t, double &A0, 
                   double &A1, double &A2, double &A3, 
                   double &c0_BH, double &c1_BH, double &c2_BH);

void analysis_dvcs_4He_t::Loop() 
{
      cout<<"starting the analysis code"<<endl;
      gStyle->SetOptStat(0);
      gStyle->SetOptFit(1);
      gStyle->SetLabelSize(0.03,"xyz"); // size of axis value font
      gStyle->SetTitleSize(0.035,"xyz"); // size of axis title font
      gStyle->SetTitleFont(22,"xyz"); // font option
      gStyle->SetLabelFont(22,"xyz");
      //gStyle->SetTitleOffSet(1.2,"y");
      gStyle->SetCanvasBorderMode(0);
      gStyle->SetCanvasBorderSize(0);
      gStyle->SetPadBottomMargin(0.16); //margins...
      gStyle->SetPadTopMargin(0.16);
      gStyle->SetPadLeftMargin(0.16);
      gStyle->SetPadRightMargin(0.16);
      gStyle->SetFrameBorderMode(0);
      gStyle->SetPaperSize(20,24);
      gStyle->SetLabelSize(0.05,"xy");
      gStyle->SetTitleSize(0.06,"xy");
   
      const int n_con = 1;
      const int n_xB = 3;
      const int n_t = 7;
    
      double xB_lims[4] = {-0.05, 0.17, 0.23, 0.5};
      double t_lims[n_t+1] = {-0.05, 0.064, 0.075, 0.086, 0.10,  0.125,  0.17, 0.95};

  // 20 days at 1.5e34
      TH1F * h_dvcs_N_p[n_con][n_xB][n_t];
      TH1F * h_dvcs_N_m[n_con][n_xB][n_t];

      TH1D *h_t_xB_Coh[n_con][n_xB][n_t];
      TH1D *h_t_Q2_Coh[n_con][n_xB][n_t];
      TH1D *h_t_t_Coh[n_con][n_xB][n_t];
 
      TH2D *h_Q2_xB_Coh[n_con];
      TH2D *h_Q2_tt_Coh[n_con];
      TH2D *h_tt_xB_Coh[n_con];

      TH1D *hh_xB[n_xB];
      TH1D *hh_Q2[n_xB];

      for(int i=0; i<n_con; i++)
      {
         hh_xB[i] = new TH1D(Form("hh_xB[%d]",i),"", 150, 0, 0.11);
         hh_Q2[i] = new TH1D(Form("hh_Q2[%d]",i),"", 150, 3.5, 32.0);

         h_Q2_tt_Coh[i]  = new TH2D(Form("h_Q2_tt_Coh[%d]",i),"",150, 3.5, 32, 150, 0.005, 0.1);
         h_Q2_xB_Coh[i]  = new TH2D(Form("h_Q2_xB_Coh[%d]",i),"",150, 0.0, 0.11, 150, 3.5, 32);
         h_tt_xB_Coh[i]  = new TH2D(Form("h_tt_xB_Coh[%d]",i)," ",150, 0.005, 1, 150, 0.0, 0.11);
         for(int j=0; j<n_xB; j++)
         {
            for(int k=0; k<n_t; k++)
            {
            
               h_dvcs_N_p[i][j][k] = new TH1F(Form("h_dvcs_N_p[%d][%d][%d]",i,j,k)," ", 12,0,360);
               h_dvcs_N_m[i][j][k] = new TH1F(Form("h_dvcs_N_m[%d][%d][%d]",i,j,k)," ", 12,0,360);
 
               h_t_Q2_Coh[i][j][k] = new TH1D(Form("h_t_Q2_Coh[%d][%d][%d]",i,j,k),"",150, 0.5, 32);
               h_t_xB_Coh[i][j][k] = new TH1D(Form("h_t_xB_Coh[%d][%d][%d]",i,j,k),"",150, 0.0, 0.1);
               h_t_t_Coh[i][j][k]  = new TH1D(Form("h_t_t_Coh[%d][%d][%d]",i,j,k),"",150, 0.0, 0.7);
          }
        }
      }

   TH1D *h_xsec = new TH1D("h_xsec","BH cross section",150, 0, 2.0);
   
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   
for (Long64_t jentry=0; jentry<nentries;jentry++)
  {

    if (jentry% 1000 == 0) printf("still running %d \n",(int)jentry);
     // if (jentry== 200000) break;

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
         int helicity = 1;
         //if(hel == 1) helicity = 1;
         //else if (hel == 0) helicity = -1;

         int which_Q2 = 0;
         int which_xB = -1;
         int which_t  = -1;
       
         for(int ii=0; ii<n_xB; ii++){
            if(xB_lims[ii]<=Xbj && Xbj <xB_lims[ii+1]) which_xB = ii ; }
       
         for(int ii=0; ii<n_t; ii++){
            if(t_lims[ii]<=t && t <t_lims[ii+1]) which_t = ii ; }
       
              h_tt_xB_Coh[which_Q2]         ->Fill(t, Xbj);
              h_Q2_xB_Coh[which_Q2]         ->Fill(Xbj, Q2);
              h_Q2_tt_Coh[which_Q2]         ->Fill( Q2,t);

              if (which_xB > -1 && which_t > -1 ) {
                   if (helicity == 1 )      h_dvcs_N_p[which_Q2][which_xB][which_t]->Fill(phih,1.0); 
                   else if (helicity == -1) h_dvcs_N_m[which_Q2][which_xB][which_t]->Fill(phih,1.0);

                   hh_xB[which_xB]        ->Fill(Xbj);
                   hh_Q2[which_xB]        ->Fill(Q2);

                   h_t_Q2_Coh[which_Q2][which_xB][which_t]->Fill(Q2);
                   h_t_xB_Coh[which_Q2][which_xB][which_t]->Fill(Xbj);
                   h_t_t_Coh[which_Q2][which_xB][which_t]->Fill(t);
                 }

   }


    TCanvas *c66 = new TCanvas("c66","",650,550 );
     for(int ii=0; ii<n_con; ii++){
        
        c66->cd();
           h_Q2_xB_Coh[ii]   ->Draw("colz");
           h_Q2_xB_Coh[ii]   ->GetYaxis()->CenterTitle(true);
           h_Q2_xB_Coh[ii]   ->SetYTitle("Q^{2} [GeV^{2}]");
           h_Q2_xB_Coh[ii]   ->GetYaxis()->SetTitleSize(0.07);
           h_Q2_xB_Coh[ii]   ->GetXaxis()->CenterTitle(true);
           h_Q2_xB_Coh[ii]   ->GetXaxis()->SetTitleSize(0.07);
           h_Q2_xB_Coh[ii]   ->SetXTitle("x_{B}");
          c66->SetLogx();
          c66->Print(Form("figs/png/coh_Q2_xB%d.png",ii));
          c66->Print(Form("figs/pdf/coh_Q2_xB%d.pdf",ii));
        
        c66->cd();
           h_Q2_tt_Coh[ii]   ->Draw("colz");
           h_Q2_tt_Coh[ii]   ->GetXaxis()->CenterTitle(true);
           h_Q2_tt_Coh[ii]   ->SetXTitle("Q^{2} [GeV^{2}]");
           h_Q2_tt_Coh[ii]   ->GetXaxis()->SetTitleSize(0.07);
           h_Q2_tt_Coh[ii]   ->GetYaxis()->CenterTitle(true);
           h_Q2_tt_Coh[ii]   ->GetYaxis()->SetTitleSize(0.07);
           h_Q2_tt_Coh[ii]   ->SetYTitle("-t [GeV^{2}]");
          c66->SetLogx();
          c66->Print(Form("figs/png/coh_Q2_tt%d.png",ii));
          c66->Print(Form("figs/pdf/coh_Q2_tt%d.pdf",ii));

        c66->cd();
           hh_xB[ii]   ->Draw();                        
           hh_xB[ii]   ->GetXaxis()->CenterTitle(true);
           hh_xB[ii]   ->SetXTitle("x_{B}");
           hh_xB[ii]   ->GetXaxis()->SetTitleSize(0.07);
           hh_xB[ii]   ->SetLineWidth(2);
          c66->SetLogy();
          c66->Print(Form("figs/png/coh_xB_%d.png",ii));
          c66->Print(Form("figs/pdf/coh_xB_%d.pdf",ii));

        c66->cd();
           hh_Q2[ii]   ->Draw();                        
           hh_Q2[ii]   ->GetXaxis()->CenterTitle(true);
           hh_Q2[ii]   ->SetXTitle("Q^{2} [GeV^{2}]");
           hh_Q2[ii]   ->GetXaxis()->SetTitleSize(0.07);
           hh_Q2[ii]   ->SetLineWidth(2);
          c66->SetLogy();
          c66->Print(Form("figs/png/coh_Q2_%d.png",ii));
          c66->Print(Form("figs/pdf/coh_Q2_%d.pdf",ii));

         }
/* 
    TCanvas *c55 = new TCanvas("c55","",750,600 ); 
             c55->cd();  
*/    
    TLine *l = new TLine();
           l->SetLineWidth(2);
           l->SetLineColor(kBlack);
  
  
    TCanvas *c6 = new TCanvas("c6","",650,550 );
     for(int ii=0; ii<n_con; ii++){
        
        c6->cd();
           h_tt_xB_Coh[ii]   ->Draw("colz");
           h_tt_xB_Coh[ii]   ->GetXaxis()->CenterTitle(true);
           h_tt_xB_Coh[ii]   ->SetXTitle("-t [GeV^{2}]");
           h_tt_xB_Coh[ii]   ->GetXaxis()->SetTitleSize(0.07);
           h_tt_xB_Coh[ii]   ->GetYaxis()->CenterTitle(true);
           h_tt_xB_Coh[ii]   ->GetYaxis()->SetTitleSize(0.07);
           h_tt_xB_Coh[ii]   ->SetYTitle("x_{B}");
              c6->SetLogz();
              c6->SetLogx();
                 for(int i=0; i<n_xB+1;i++)
                     l->DrawLine(0.05,  xB_lims[i], t_lims[n_t], xB_lims[i]);
                 for(int i=0; i<n_t+1;i++)
                        l->DrawLine(t_lims[i], 0.05, t_lims[i], 0.5);
          c6->Print(Form("figs/png/coh_t_xB%d.png",ii));
          c6->Print(Form("figs/pdf/coh_t_xB%d.pdf",ii));
         }
 

    vector<vector<vector<double>>> alu_t_x     ;
    vector<vector<vector<double>>> alu_t_x_err ;
    vector<vector<vector<double>>> Im_t_x     ;
    vector<vector<vector<double>>> Im_t_x_err ;
    vector<vector<vector<double>>> Re_t_x     ;
    vector<vector<vector<double>>> Re_t_x_err ;

    vector<vector<vector<double>>> mean_t      ;
    vector<vector<vector<double>>> mean_x      ;
    vector<vector<vector<double>>> mean_Q2     ;
    vector<vector<vector<double>>> mean_t_err  ;
    vector<vector<vector<double>>> mean_x_err  ;
    vector<vector<vector<double>>> mean_Q2_err ;

    alu_t_x    .resize(n_con); 
    alu_t_x_err.resize(n_con);
    Im_t_x     .resize(n_con); 
    Im_t_x_err .resize(n_con);
    Re_t_x     .resize(n_con); 
    Re_t_x_err .resize(n_con);

    mean_t     .resize(n_con);
    mean_x     .resize(n_con);
    mean_Q2    .resize(n_con);
    mean_t_err .resize(n_con);
    mean_x_err .resize(n_con);
    mean_Q2_err.resize(n_con);

    for(int ii=0; ii<n_con; ii++){

       alu_t_x[ii]    .resize(n_xB); 
       alu_t_x_err[ii].resize(n_xB);
       Im_t_x[ii]     .resize(n_xB); 
       Im_t_x_err[ii] .resize(n_xB);
       Re_t_x[ii]     .resize(n_xB); 
       Re_t_x_err[ii].resize(n_xB);

       mean_t[ii]     .resize(n_xB);
       mean_x[ii]     .resize(n_xB);
       mean_Q2[ii]    .resize(n_xB);
       mean_t_err[ii] .resize(n_xB);
       mean_x_err[ii] .resize(n_xB);
       mean_Q2_err[ii].resize(n_xB);
   
       for(int jj=0; jj<n_xB; jj++){

          alu_t_x[ii][jj]     .resize(n_t); 
          alu_t_x_err[ii][jj] .resize(n_t);
          Im_t_x[ii][jj]      .resize(n_t); 
          Im_t_x_err[ii][jj]  .resize(n_t);
          Re_t_x[ii][jj]      .resize(n_t); 
          Re_t_x_err[ii][jj]  .resize(n_t);

          mean_t[ii][jj]      .resize(n_t);
          mean_x[ii][jj]      .resize(n_t);
          mean_Q2[ii][jj]     .resize(n_t);
          mean_t_err[ii][jj]  .resize(n_t);
          mean_x_err[ii][jj]  .resize(n_t);
          mean_Q2_err[ii][jj] .resize(n_t);
 
       }
    }

    for(int ii=0; ii<1; ii++){
       for(int jj=0; jj<n_xB; jj++){
          for(int kk=0; kk<n_t; kk++){
             c66->cd();          
             mean_Q2[ii][jj][kk] = h_t_Q2_Coh[ii][jj][kk]->GetMean();
             mean_x[ii][jj][kk]  = h_t_xB_Coh[ii][jj][kk]->GetMean();
             mean_t[ii][jj][kk]  = h_t_t_Coh[ii][jj][kk]->GetMean();

             mean_Q2_err[ii][jj][kk] = 0.0;
             mean_x_err[ii][jj][kk] = 0.0;
             mean_t_err[ii][jj][kk] = 0.0;//(t_lims[kk+1]-t_lims[kk])/2.4;
          // h_t_Q2_Coh[ii][jj][kk]->Draw();
          // c66->Print(Form("figs/slices/coh_Q2%d_%d_%d.png",ii,jj,kk));

          // h_t_xB_Coh[ii][jj][kk]->Draw();
          // c66->Print(Form("figs/slices/coh_xB%d_%d_%d.png",ii,jj,kk));

          // h_t_t_Coh[ii][jj][kk]->Draw();
          // c66->Print(Form("figs/slices/coh_t%d_%d_%d.png",ii,jj,kk));



          }
       }
    }


   TCanvas *c3 = new TCanvas("c3","",1000,1300 );
            c3->Divide(3,7,-0.00005,-0.00005); 
            c3->SetGrid();
/*
   // as a function of -t in xB bins
   for(int ii=0; ii<1; ii++){
    for(int jj=0; jj<n_xB; jj++){
       for(int kk=0; kk<n_t; kk++){
            c3->SetGrid();
         if (jj == 0) c3->cd(3*kk +1);
         if (jj == 1) c3->cd(3*kk +2);
         if (jj == 2) c3->cd(3*kk +3);
         
         h_dvcs_N_p[ii][jj][kk]->Sumw2(); 
         h_dvcs_N_m[ii][jj][kk]->Sumw2();
         
         TH1* hsum=(TH1*)h_dvcs_N_p[ii][jj][kk]->Clone("hsum");
         TH1* hdif=(TH1*)h_dvcs_N_p[ii][jj][kk]->Clone("hdif");
         hsum->Add(h_dvcs_N_m[ii][jj][kk]);
         hdif->Add(h_dvcs_N_m[ii][jj][kk],-1);
         
         TH1* hasy=(TH1*)hdif->Clone("hasy");
              hasy->Divide(hsum);
              //hasy->SetTitle("Coherent A_{LU}");
              hasy->GetYaxis()->SetTitle("A_{LU}");
              hasy->GetYaxis()->SetTitleSize(0.07);
              hasy->GetYaxis()->CenterTitle(true);
              hasy->GetXaxis()->SetTitle("#phi_{h} [deg.]");
              hasy->GetXaxis()->CenterTitle(true);
              hasy->GetXaxis()->SetTitleSize(0.07);
              hasy->SetMarkerStyle(21);
              if(ii==0) { hasy->SetMarkerColor(kBlack); hasy->SetLineColor(kBlack);} 
              if(ii==1) { hasy->SetMarkerColor(kBlue);  hasy->SetLineColor(kBlue); }
              if(ii==2) { hasy->SetMarkerColor(kRed);   hasy->SetLineColor(kRed);  }
              hasy->GetXaxis()->SetLimits(-5.0,365.0);
              hasy->GetYaxis()->SetRangeUser(-0.7,0.7);
              if(ii==0) hasy->Draw();
              if(ii==1) hasy->Draw("same");
              if(ii==2) hasy->Draw("same");

          TLatex *l2= new TLatex(10.0,-0.25,Form("%.2f< x_{B} <%.2f",xB_lims[jj], xB_lims[jj+1]));
                  l2->Draw("same");
          TLatex *l3= new TLatex(10.0,-0.45,Form("%.2f< -t <%.2f",t_lims[kk], t_lims[kk+1]));
                  l3->Draw("same");

          TF1 *myfit;
          if(kk<7) myfit = new TF1("myfit","[0]*sin(x*3.1416/180.0)/(1 + [1]*cos(x*3.1416/180.0))",0.0,360.0);
          else if(kk==7) myfit = new TF1("myfit","[0]*sin(x*3.1416/180.0)",0.0,360.0);
          myfit->SetLineColor(kRed);
          myfit->SetLineWidth(2);

          double  A0, A1, A2, A3, c0_BH, c1_BH, c2_BH; 
          Calculate_CFF(mean_Q2[ii][jj][kk], mean_x[ii][jj][kk], -1.0*mean_t[ii][jj][kk], 
                        A0, A1, A2, A3, c0_BH, c1_BH, c2_BH);

          TF1 *ffit;
          TFitResultPtr r;
          ffit = new TF1("ffit",  Form("%0.7f*[0]*sin(x*3.14/180.0) / (%0.7f+ %0.7f*cos(x*3.14/180.0) + %0.7f*cos(2*x*3.14/180.0) + %0.7f*([0]*[0] + [1]*[1]) + %0.7f*[1] + %0.7f*[1]*cos(x*3.14/180.0))",  
                         A0, c0_BH, c1_BH, c2_BH, A1, A2, A3), 0.0, 360.0);
           ffit->SetLineColor(kRed);
           ffit->SetParName(0,"Im(H_{A})");
           ffit->SetParName(1,"Re(H_{A})");



          // fit with a sin and cosin function ------------------------------------
        //hasy->Fit("myfit");
        double IM_CFF; ;
        double RE_CFF;;
        Find_CFF(mean_x[ii][jj][kk], abs(mean_t[ii][jj][kk]), IM_CFF, RE_CFF);
          
        alu_t_x[ii][jj][kk] = 0.1 *(2-jj)-0.3;
        Im_t_x[ii][jj][kk] = IM_CFF;
        Re_t_x[ii][jj][kk] = RE_CFF;

        //
        // alu_t_x_err[ii][jj][kk] =  myfit->GetParError(0);
        //Im_t_x_err[ii][jj][kk]  = (myfit->GetParError(0)/myfit->GetParameter(0))* abs(IM_CFF);
        //Re_t_x_err[ii][jj][kk]  =  (myfit->GetParError(1)/myfit->GetParameter(1)) *abs(Re_CFF);


          // fit with the full form of the CFF

          ffit->SetParameter(0,27.0);
          ffit->SetParameter(1,-4.0);
          hasy->Fit("ffit");
          r = hasy->Fit(ffit,"S");
          double x[1] = { 90.0};  double err[1];
          r->GetConfidenceIntervals(1, 1, 1, x, err, 0.683, true);
          double ALU_fit     = ffit->Eval(x[0]);
          double ALU_fit_err = err[0];

          alu_t_x_err[ii][jj][kk] =  ALU_fit_err;
          Im_t_x_err[ii][jj][kk]  = ffit->GetParError(0); 
          Re_t_x_err[ii][jj][kk]  = ffit->GetParError(1);          

          // print a high and a low statistics Alu bins ------------------------------
          if( jj==1 && kk==1 ) 
          {
             TCanvas *c66 = new TCanvas("c6","",650,550 );
                      c66->cd();

             hasy->Draw();
             //myfit->Draw("same");
             l2->Draw("same");
             l3->Draw("same");
             c66->Print("figs/png/BSA_Coherent_Phi_high_stat.png");
             c66->Print("figs/pdf/BSA_Coherent_Phi_high_stat.pdf");
          }

          if( jj==1 && kk==5 ) 
          {
             TCanvas *c66 = new TCanvas("c6","",650,550 );
                      c66->cd();

             hasy->Draw();
             //myfit->Draw("same");
             l2->Draw("same");
             l3->Draw("same");
             c66->Print("figs/png/BSA_Coherent_Phi_low_stat.png");
             c66->Print("figs/pdf/BSA_Coherent_Phi_low_stat.pdf");
          }

       }
    }
   }


    c3->Print("figs/png/BSA_Coherent_Phi_t.png");
    c3->Print("figs/pdf/BSA_Coherent_Phi_t.pdf");
 */

  vector<double> M_XB;   M_XB.resize(n_xB);
  vector<double> M_Q2;   M_Q2.resize(n_xB);

  for(int ii=0; ii<n_xB; ii++){

    cout<<hh_xB[ii]   ->GetMean()<<endl;
    cout<<hh_Q2[ii]   ->GetMean()<<endl;
    M_XB[ii]= (hh_xB[ii]->GetMean());
    M_Q2[ii]= (hh_Q2[ii]->GetMean());
   }

//    plot_ALU_projections(mean_t, mean_t_err, alu_t_x, alu_t_x_err);
//    plot_CFF_projections(mean_t, mean_t_err, Im_t_x_err, Re_t_x_err);
//    FunProfile(mean_t, mean_t_err, M_XB, M_Q2, Im_t_x, Im_t_x_err);

}




void Find_CFF( double xB, double t, double &Im, double &Re){
   
        if( 0.10< xB && xB<0.18 ){
             Im  = (60.77 - pow(5.08 *t, 6)) * exp( -13.2*t );
             Re  = (-8.62 + 15.86*t - pow( -3.59*t, 3)) * exp( -7.27*t);
             }
       else if( 0.18< xB && xB<0.222 ){
                Im  = (39.0 - pow( 5.45*t, 5)) * exp( -12.83*t );
                Re  = (-12.53 + 33.88*t - pow( -0.11*t, 3)) * exp( -7.4*t);
                }
          else if( 0.22< xB ){
                   Im  = (26.06 - pow( 4.4*t, 6)) * exp( -12.8*t );
                   Re  = (-18.0 + 35.18*t - pow( -4.56*t, 3)) * exp( -9.5*t);
                   }

/*
 if( 0.1< xB && xB<0.145 ){
          Im  = (99.0 - pow(5.5*t,6)) * exp(-13.47 *t );
          Re  = (-5.6 + 14.7*t - pow(-1.71*t,4)) * exp(-3.9 *t );
          }
    else if( 0.145< xB && xB<0.195 ){
             Im  = (60.77 - pow(5.08 *t, 6)) * exp( -13.2*t );
             Re  = (-8.62 + 15.86*t - pow( -3.59*t, 3)) * exp( -7.27*t);
             }
       else if( 0.195< xB && xB<0.245 ){
                Im  = (39.0 - pow( 5.45*t, 5)) * exp( -12.83*t );
                Re  = (-12.53 + 33.88*t - pow( -0.11*t, 3)) * exp( -7.4*t);
                }
          else if( 0.245< xB && xB<0.295 ){
                   Im  = (26.06 - pow( 4.4*t, 6)) * exp( -12.8*t );
                   Re  = (-18.0 + 35.18*t - pow( -4.56*t, 3)) * exp( -9.5*t);
                   }
             else if( 0.295< xB && xB<0.345 ){
                      Im  = (17.01 - pow( 4.08*t, 6)) * exp( -12.5*t );
                      Re  = (-23.9 + 45.1*t - pow( -5.05*t, 3)) * exp( -10.1*t);
                      }
                else if( 0.345< xB ){
                         Im  = (11.1 - pow( 4.08*t, 6)) * exp( -12.1*t );
                         Re  = (-30.2  -154.7*t - pow( 0.037*t, 3)) * exp( -16.8*t);
                         }
 */
 }




void Calculate_CFF(double Q2, double xB, double t, 
                   double &A0, double &A1, double &A2, double &A3, 
                   double &c0_BH, double &c1_BH, double &c2_BH)
  {

  double PI = 3.1416;
  double phi = 90.0*PI/180.0;

  double MPROT = 0.93827;
  double MALPH = 3.7274;
  double EBEAM = 11.0;
  /*
  double e = 2*xA*MALPH/sqrt(Q2); // epsilon
  double e2 = e*e;
  double Tmin = -Q2 * (2*xA1*(1-sqrt(1+e2))+e2) / (4*xA*xA1 + e2);
 
  // kinematical factors
  double J = (1-y-y*e2/2) * (1+t/Q2)  - (1-xA)*(2-y)*t/Q2;
  double dt = (t - Tmin)/Q2;
  double K = -1.0*dt * xA1 * (1 -y - y*y*e2/4) * (sqrt(1+e2) + ((4*xA*xA1+e2)*dt/(4*xA1)) ); //sqrt(1 - y + e2*y*y/4)* (K_hat)/(sqrt(Q2));
  double K2 = K*K; 


  // Helium form factor
  double a=0.316;
  double b=0.681;
  double FF4He = (1-pow(a*a*Q2,6))*exp(-b*b*Q2); // 1e-15;
 
  // BH propagators
  double P1_phi = -( J + 2*K*cos(phi)) / (y * (1+e2));
  double P2_phi =  1 + t/Q2 + (1/(y*(1+e*e))) * (J + 2*K*cos(phi));


  // BH fourier coefficients
   c0_BH = ( (pow(2-y,2) + pow(y*(1+e2),2)) * (e2*Q2/t+4*(1-xA)+(4*xA+e2)*t/Q2)
                  + 2*e2*(4*(1-y)*(3+2*e2)+y*y*(2-e2*e2))
                  - 4*xA*xA*pow(2-y,2)*(2+e2)*t/Q2
                  + 8*K2*e2*Q2/t ) * pow(FF4He,2);
   c1_BH = -8*(2-y)* K * (2*xA+e*e*(1-Q2/t)) * pow(FF4He,2);
   c2_BH = 8*K2*e2*Q2*pow(FF4He,2)/t;
   */

  double xA = xB*MPROT/MALPH;
  double xA1= 1 - xA;
  double y = Q2/2/MALPH/xA/EBEAM;
  double e = 2*xA*MALPH/sqrt(Q2); // epsilon
  double e2 = e*e;
  double T0 = -Q2 * (2*xA1*(1-sqrt(1+e2))+e2) / (4*xA*xA1 + e2);

  // kinematical factors
  double J = (1-y-y*e2/2) * (1+t/Q2) - (1-xA)*(2-y)*t/Q2;
  double dt = (t - T0)/Q2;
  double K_hat = sqrt(T0 - t) * sqrt(xA1*sqrt(1+e2) + (T0 - t)*(e2 + 4*xA1*xA)/(4*Q2) );
  double K2 = -1.0*dt * xA1 * (1 -y - y*y*e2/4) * (sqrt(1+e2) + ((4*xA*xA1+e2)*dt/(4*xA1)) );
  //double K = sqrt(1 - y + e2*y*y/4)* (K_hat)/(sqrt(Q2));
  //double K2 = K*K;
  double K = sqrt(K2);

  // BH propagators
  double P1_phi = -1.0*(J + 2*K*cos(PI-phi)) / (y * (1+e2));
  double P2_phi =  1 + t/Q2 + (1/(y*(1+e*e))) * (J + 2*K*cos(PI-phi));

  
  // Helium form factor
  double a=0.316;
  double b=0.681;
  double FF4He = (1-pow(a*a*abs(t)/pow(0.197327,2),6))*exp(-b*b*abs(t)/pow(0.197327,2));//  1e-15;
  //double FF4He = (1-pow(a*a*Q2,6))*exp(-b*b*Q2);
 
  // BH fourier coefficients
   c0_BH = ( (pow(2-y,2)+pow(y * (1+e2),2)) * (e2*Q2/t + 4*xA1 + (4*xA+e2)*t/Q2) 
                   + 2*e2*(4*(1-y)*(3+2*e2) + y*y*(2-e2*e2))
                   - 4*xA*xA*pow(2-y,2)*(2+e2)*t/Q2
                   + 8*K2*e2*Q2/t) * pow(FF4He,2);
   c1_BH = -8*(2-y)* K * (2*xA + e2 - e2*Q2/t) * pow(FF4He,2);
   c2_BH =  8*K2*e2*Q2*pow(FF4He,2)/t;


  // redefine the fourier hamoinic to fit for Re and Im of HA 
  double C_DVCS_0 = 2*((2-2*y+y*y + 0.5*e2*y*y)/(1+e2)) ;

  double C_INT_plus_plus_0 = ( -4*(2-y)*(1 + sqrt(1+e2))/(pow(1+e2 ,2))) * ( (pow(K_hat*(2-y) ,2))/(Q2*sqrt(1+e2)) +  (t/Q2)*(1-y-y*y*e2/4)*(2-xA)*( 1 + (2*xA*(t/Q2)*(2-xA + ((sqrt(1+e2) -1)/(2)) + ((e2/(2*xA))) ) + e2 )/((2-xA)*(1+sqrt(1+e2)))  )) * FF4He;  
  
  double C_INT_plus_plus_1 = ((-16*K*(1 - y - e2*y*y/4))/pow(1+e2 ,5/2)) * ( ( 1 + (1-xA)*((sqrt(1+e2) -1) /(2*xA)) + e2/(4*xA))* (xA*t/Q2)  - 3.0*e2/4.0) - 4.0*K*(2-2*y+y*y + e2*y*y/2) * ( (1+sqrt(1+e2) -e2)/pow(1+e2 ,5/2) ) *(1 - (1-3.0*xA)*t/Q2 + (1-sqrt(1+e2)+3*e2)/(1+sqrt(1+e2) -e2)*(xA*t/Q2)) * FF4He;  

  double S_INT_plus_plus_1 = (8*K*(2-y)*y /(1+e2)) * ( 1 + ((1-xA+0.5*(sqrt(1+e2)-1))/(1+e2))*dt ) * FF4He;  



   A0 = xA* pow(1+e2, 2) * S_INT_plus_plus_1/y;
   A1 = 2*xA*xA*t*((1+e2)/Q2) * (2-2*y+y*y + e2*y*y/2) * P1_phi * P2_phi * C_DVCS_0; 
   A2 = xA*pow(1+e2,2) * C_INT_plus_plus_0 /y;  
   A3 = xA* pow(1+e2,2) * C_INT_plus_plus_1/y; 


  // ALU2 =  A0* Im* sin(phi)/ (c0_BH + c1_BH*cos(phi) + c2_BH*cos(2*phi) + A1*(Re*Re + Im*Im) + A2*Re + A3*Re * cos(phi));

}


