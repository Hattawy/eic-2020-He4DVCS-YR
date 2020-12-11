#define analysis_dvcs_4He_t_cxx
#include "analysis_dvcs_4He_t.h"
#include <iostream>
#include <fstream>
#include <TH1.h>
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
#include "FunProfile_0.C"
#include "FunProfile_1.C"
#include "FunProfile_2.C"

using namespace std;

int Find_Run_Number(const char*);
void Find_CFF( double xB, double t, double &Im, double &Re, int n_xB);
void Calculate_CFF(double Q2, double xB, double t, double y, 
                   double &A0, double &A1, double &A2, double &A3, 
                   double &c0_BH, double &c1_BH, double &c2_BH);
void Calculate_ALU(double Q2, double xB, double t, double y, double phi , double Im, double Re, double &ALU);

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
  
      const double xfac = 1.0/sqrt(50000.0/22.0);
      const double sys_fac = 0.01;
      const double fac_0 =0.012;
      const double fac_1 =0.0145;
      const double fac_2 =0.02175;

      const int n_con = 1;
      const int n_xB = 9;
      const int n_t_0 = 14;
      const int n_t_1 = 8;
      const int n_t_2 = 5;
    
      double xB_lims[n_xB+1] = {0.0023, 0.0028, 0.0033, 0.004, 0.0048, 0.0057, 0.0068, 0.01, 0.015, 0.08};
      double t_lims_0[n_t_0+1] = {0.01, 0.0125, 0.016, 0.02, 0.025, 0.031, 0.04, 0.05, 0.067, 0.090, 0.12, 0.155, 0.2, 0.26, 0.36};
      double t_lims_1[n_t_1+1] = {0.04, 0.05, 0.067, 0.090, 0.12, 0.155, 0.2, 0.26, 0.36};
      double t_lims_2[n_t_2+1] = {0.090, 0.12, 0.155, 0.2, 0.26, 0.36};

  // n_t_0
      TH1F * h_dvcs_N_p_0[n_con][n_xB][n_t_0];
      TH1F * h_dvcs_N_m_0[n_con][n_xB][n_t_0];

      TH1D *h_t_xB_Coh_0[n_con][n_xB][n_t_0];
      TH1D *h_t_Q2_Coh_0[n_con][n_xB][n_t_0];
      TH1D *h_t_t_Coh_0[n_con][n_xB][n_t_0];
      TH1D *h_t_y_Coh_0[n_con][n_xB][n_t_0];
 
      TH2D *h_Q2_xB_Coh_0[n_con];
      TH2D *h_Q2_tt_Coh_0[n_con];
      TH2D *h_tt_xB_Coh_0[n_con];

      TH1D *hh_xB_0[n_xB];
      TH1D *hh_Q2_0[n_xB];

      for(int i=0; i<n_con; i++) {
         h_Q2_tt_Coh_0[i]  = new TH2D(Form("h_Q2_tt_Coh_0[%d]",i),"",500, 3.5, 32, 500, 0.006, 1);
         h_Q2_xB_Coh_0[i]  = new TH2D(Form("h_Q2_xB_Coh_0[%d]",i),"",500, 0.001, 0.11, 500, 3.5, 32);
         h_tt_xB_Coh_0[i]  = new TH2D(Form("h_tt_xB_Coh_0[%d]",i),"",500, 0.006, 1, 500, 0.0017, 0.11);
        
         for(int j=0; j<n_xB; j++) {
           hh_xB_0[j] = new TH1D(Form("hh_xB_0[%d]",j),"", 300, 0, 0.11);
           hh_Q2_0[j] = new TH1D(Form("hh_Q2_0[%d]",j),"", 300, 3.5, 32.0);

          for(int k=0; k<n_t_0; k++) {
            
               h_dvcs_N_p_0[i][j][k] = new TH1F(Form("h_dvcs_N_p_0[%d][%d][%d]",i,j,k)," ", 12,0,360);
               h_dvcs_N_m_0[i][j][k] = new TH1F(Form("h_dvcs_N_m_0[%d][%d][%d]",i,j,k)," ", 12,0,360);
 
               h_t_Q2_Coh_0[i][j][k] = new TH1D(Form("h_t_Q2_Coh_0[%d][%d][%d]",i,j,k),"",300, 0.5, 32);
               h_t_xB_Coh_0[i][j][k] = new TH1D(Form("h_t_xB_Coh_0[%d][%d][%d]",i,j,k),"",300, 0.0, 0.1);
               h_t_t_Coh_0[i][j][k]  = new TH1D(Form("h_t_t_Coh_0[%d][%d][%d]",i,j,k),"",300, 0.0, 0.7);
               h_t_y_Coh_0[i][j][k]  = new TH1D(Form("h_t_y_Coh_0[%d][%d][%d]",i,j,k),"",300, 0.0, 1.0);
          }}}



  // n_t_1
      TH1F * h_dvcs_N_p_1[n_con][n_xB][n_t_1];
      TH1F * h_dvcs_N_m_1[n_con][n_xB][n_t_1];

      TH1D *h_t_xB_Coh_1[n_con][n_xB][n_t_1];
      TH1D *h_t_Q2_Coh_1[n_con][n_xB][n_t_1];
      TH1D *h_t_t_Coh_1[n_con][n_xB][n_t_1];
      TH1D *h_t_y_Coh_1[n_con][n_xB][n_t_1];
 
      TH2D *h_Q2_xB_Coh_1[n_con];
      TH2D *h_Q2_tt_Coh_1[n_con];
      TH2D *h_tt_xB_Coh_1[n_con];

      TH1D *hh_xB_1[n_xB];
      TH1D *hh_Q2_1[n_xB];

      for(int i=0; i<n_con; i++) {
         h_Q2_tt_Coh_1[i]  = new TH2D(Form("h_Q2_tt_Coh_1[%d]",i),"",500, 3.5, 32, 500, 0.006, 1);
         h_Q2_xB_Coh_1[i]  = new TH2D(Form("h_Q2_xB_Coh_1[%d]",i),"",500, 0.001, 0.11, 500, 3.5, 32);
         h_tt_xB_Coh_1[i]  = new TH2D(Form("h_tt_xB_Coh_1[%d]",i),"",500, 0.006, 1, 500, 0.0017, 0.11);
        
         for(int j=0; j<n_xB; j++) {
           hh_xB_1[j] = new TH1D(Form("hh_xB_1[%d]",j),"", 300, 0, 0.11);
           hh_Q2_1[j] = new TH1D(Form("hh_Q2_1[%d]",j),"", 300, 3.5, 32.0);

          for(int k=0; k<n_t_1; k++) {
            
               h_dvcs_N_p_1[i][j][k] = new TH1F(Form("h_dvcs_N_p_1[%d][%d][%d]",i,j,k)," ", 12,0,360);
               h_dvcs_N_m_1[i][j][k] = new TH1F(Form("h_dvcs_N_m_1[%d][%d][%d]",i,j,k)," ", 12,0,360);
 
               h_t_Q2_Coh_1[i][j][k] = new TH1D(Form("h_t_Q2_Coh_1[%d][%d][%d]",i,j,k),"",300, 0.5, 32);
               h_t_xB_Coh_1[i][j][k] = new TH1D(Form("h_t_xB_Coh_1[%d][%d][%d]",i,j,k),"",300, 0.0, 0.1);
               h_t_t_Coh_1[i][j][k]  = new TH1D(Form("h_t_t_Coh_1[%d][%d][%d]",i,j,k),"",300, 0.0, 0.7);
               h_t_y_Coh_1[i][j][k]  = new TH1D(Form("h_t_y_Coh_1[%d][%d][%d]",i,j,k),"",300, 0.0, 1.0);
          }}}

 

  // n_t_2
      TH1F * h_dvcs_N_p_2[n_con][n_xB][n_t_2];
      TH1F * h_dvcs_N_m_2[n_con][n_xB][n_t_2];

      TH1D *h_t_xB_Coh_2[n_con][n_xB][n_t_2];
      TH1D *h_t_Q2_Coh_2[n_con][n_xB][n_t_2];
      TH1D *h_t_t_Coh_2[n_con][n_xB][n_t_2];
      TH1D *h_t_y_Coh_2[n_con][n_xB][n_t_2];
 
      TH2D *h_Q2_xB_Coh_2[n_con];
      TH2D *h_Q2_tt_Coh_2[n_con];
      TH2D *h_tt_xB_Coh_2[n_con];

      TH1D *hh_xB_2[n_xB];
      TH1D *hh_Q2_2[n_xB];

      for(int i=0; i<n_con; i++) {
         h_Q2_tt_Coh_2[i]  = new TH2D(Form("h_Q2_tt_Coh_2[%d]",i),"",500, 3.5, 32, 500, 0.006, 1);
         h_Q2_xB_Coh_2[i]  = new TH2D(Form("h_Q2_xB_Coh_2[%d]",i),"",500, 0.001, 0.11, 500, 3.5, 32);
         h_tt_xB_Coh_2[i]  = new TH2D(Form("h_tt_xB_Coh_2[%d]",i),"",500, 0.006, 1, 500, 0.0017, 0.11);
        
         for(int j=0; j<n_xB; j++) {
           hh_xB_2[j] = new TH1D(Form("hh_xB_2[%d]",j),"", 300, 0, 0.11);
           hh_Q2_2[j] = new TH1D(Form("hh_Q2_2[%d]",j),"", 300, 3.5, 32.0);

          for(int k=0; k<n_t_2; k++) {
            
               h_dvcs_N_p_2[i][j][k] = new TH1F(Form("h_dvcs_N_p_2[%d][%d][%d]",i,j,k)," ", 12,0,360);
               h_dvcs_N_m_2[i][j][k] = new TH1F(Form("h_dvcs_N_m_2[%d][%d][%d]",i,j,k)," ", 12,0,360);
 
               h_t_Q2_Coh_2[i][j][k] = new TH1D(Form("h_t_Q2_Coh_2[%d][%d][%d]",i,j,k),"",300, 0.5, 32);
               h_t_xB_Coh_2[i][j][k] = new TH1D(Form("h_t_xB_Coh_2[%d][%d][%d]",i,j,k),"",300, 0.0, 0.1);
               h_t_t_Coh_2[i][j][k]  = new TH1D(Form("h_t_t_Coh_2[%d][%d][%d]",i,j,k),"",300, 0.0, 0.7);
               h_t_y_Coh_2[i][j][k]  = new TH1D(Form("h_t_y_Coh_2[%d][%d][%d]",i,j,k),"",300, 0.0, 1.0);
          }}}

 




   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   const int n_runs = 22;
   bool run_tag[n_runs];
   for (int i=0; i<n_runs; i++){ run_tag[i]  = false; }
   
for (Long64_t jentry=0; jentry<nentries;jentry++)
  {

    if (jentry% 1000000 == 0) printf("still running %d \n",(int)jentry);
     // if (jentry== 200000) break;

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
         
      int RunNumber = -999;
      const char *input_file_name;
      input_file_name = fChain->GetCurrentFile()->GetName();
      RunNumber = Find_Run_Number(input_file_name);

      if (run_tag[RunNumber-1] == false ){
         cout<<input_file_name<<" run number "<< RunNumber<<endl;
         run_tag[RunNumber-1] = true;
      }
      
      int helicity = 0;
      if(RunNumber% 2 == 1) helicity = 1;
      else if (RunNumber% 2 == 0) helicity = -1;

      int which_Q2 = 0;
      int which_xB = -1;
      int which_t_0  = -1;
      int which_t_1  = -1;
      int which_t_2  = -1;
      
      for(int ii=0; ii<n_xB; ii++){ if(xB_lims[ii]<=Xbj && Xbj <xB_lims[ii+1]) which_xB = ii ; }
      for(int ii=0; ii<n_t_0; ii++){ if(t_lims_0[ii]<=t && t <t_lims_0[ii+1]) which_t_0 = ii ; }
      for(int ii=0; ii<n_t_1; ii++){ if(t_lims_1[ii]<=t && t <t_lims_1[ii+1]) which_t_1 = ii ; }
      for(int ii=0; ii<n_t_2; ii++){ if(t_lims_2[ii]<=t && t <t_lims_2[ii+1]) which_t_2 = ii ; }
      
      if(t>= t_lims_0[0]) {     
         h_tt_xB_Coh_0[which_Q2]  ->Fill(t, Xbj);
         h_Q2_xB_Coh_0[which_Q2]  ->Fill(Xbj, Q2);
         h_Q2_tt_Coh_0[which_Q2]  ->Fill( Q2,t);
      
         if (which_xB > -1 && which_t_0 > -1 ) {
              if (helicity == 1 )      h_dvcs_N_p_0[which_Q2][which_xB][which_t_0]->Fill(phih,1.0); 
              else if (helicity == -1) h_dvcs_N_m_0[which_Q2][which_xB][which_t_0]->Fill(phih,1.0);
      
              h_t_Q2_Coh_0[which_Q2][which_xB][which_t_0]->Fill(Q2);
              h_t_xB_Coh_0[which_Q2][which_xB][which_t_0]->Fill(Xbj);
              h_t_t_Coh_0[which_Q2][which_xB][which_t_0]->Fill(t);
              h_t_y_Coh_0[which_Q2][which_xB][which_t_0]->Fill(y);
     
              hh_xB_0[which_xB]        ->Fill(Xbj);
              hh_Q2_0[which_xB]        ->Fill(Q2); 
         } }


      if(t>= t_lims_1[0]) {     
         h_tt_xB_Coh_1[which_Q2]  ->Fill(t, Xbj);
         h_Q2_xB_Coh_1[which_Q2]  ->Fill(Xbj, Q2);
         h_Q2_tt_Coh_1[which_Q2]  ->Fill( Q2,t);
      
         if (which_xB > -1 && which_t_1 > -1 ) {
              if (helicity == 1 )      h_dvcs_N_p_1[which_Q2][which_xB][which_t_1]->Fill(phih,1.0); 
              else if (helicity == -1) h_dvcs_N_m_1[which_Q2][which_xB][which_t_1]->Fill(phih,1.0);
      
              h_t_Q2_Coh_1[which_Q2][which_xB][which_t_1]->Fill(Q2);
              h_t_xB_Coh_1[which_Q2][which_xB][which_t_1]->Fill(Xbj);
              h_t_t_Coh_1[which_Q2][which_xB][which_t_1]->Fill(t);
              h_t_y_Coh_1[which_Q2][which_xB][which_t_1]->Fill(y);
     
              hh_xB_1[which_xB]        ->Fill(Xbj);
              hh_Q2_1[which_xB]        ->Fill(Q2); 
         } }


      if(t>= t_lims_2[0]) {     
         h_tt_xB_Coh_2[which_Q2]  ->Fill(t, Xbj);
         h_Q2_xB_Coh_2[which_Q2]  ->Fill(Xbj, Q2);
         h_Q2_tt_Coh_2[which_Q2]  ->Fill( Q2,t);
      
         if (which_xB > -1 && which_t_2 > -1 ) {
              if (helicity == 1 )      h_dvcs_N_p_2[which_Q2][which_xB][which_t_2]->Fill(phih,1.0); 
              else if (helicity == -1) h_dvcs_N_m_2[which_Q2][which_xB][which_t_2]->Fill(phih,1.0);
      
              h_t_Q2_Coh_2[which_Q2][which_xB][which_t_2]->Fill(Q2);
              h_t_xB_Coh_2[which_Q2][which_xB][which_t_2]->Fill(Xbj);
              h_t_t_Coh_2[which_Q2][which_xB][which_t_2]->Fill(t);
              h_t_y_Coh_2[which_Q2][which_xB][which_t_2]->Fill(y);
     
              hh_xB_2[which_xB]        ->Fill(Xbj);
              hh_Q2_2[which_xB]        ->Fill(Q2); 
         } }

   }  // end the loop over the events 

////// n__0 ///////////////////////////////////////////////////////
    TCanvas *c66 = new TCanvas("c66","",650,550 );
     for(int ii=0; ii<n_con; ii++){
      
        c66->cd();
           h_Q2_xB_Coh_0[ii]   ->Draw("colz");
           h_Q2_xB_Coh_0[ii]   ->GetYaxis()->CenterTitle(true);
           h_Q2_xB_Coh_0[ii]   ->SetYTitle("Q^{2} [GeV^{2}]");
           h_Q2_xB_Coh_0[ii]   ->GetYaxis()->SetTitleSize(0.07);
           h_Q2_xB_Coh_0[ii]   ->GetXaxis()->CenterTitle(true);
           h_Q2_xB_Coh_0[ii]   ->GetXaxis()->SetTitleSize(0.07);
           h_Q2_xB_Coh_0[ii]   ->SetXTitle("x_{B}");
          c66->SetLogx();
          c66->Print("figs/png/coh_Q2_xB_0.png");
          c66->Print("figs/pdf/coh_Q2_xB_0.pdf");
        
        c66->cd();
           h_Q2_tt_Coh_0[ii]   ->Draw("colz");
           h_Q2_tt_Coh_0[ii]   ->GetXaxis()->CenterTitle(true);
           h_Q2_tt_Coh_0[ii]   ->SetXTitle("Q^{2} [GeV^{2}]");
           h_Q2_tt_Coh_0[ii]   ->GetXaxis()->SetTitleSize(0.07);
           h_Q2_tt_Coh_0[ii]   ->GetYaxis()->CenterTitle(true);
           h_Q2_tt_Coh_0[ii]   ->GetYaxis()->SetTitleSize(0.07);
           h_Q2_tt_Coh_0[ii]   ->SetYTitle("-t [GeV^{2}]");
          c66->SetLogy();
          c66->SetLogx();
          c66->Print("figs/png/coh_Q2_tt_0.png");
          c66->Print("figs/pdf/coh_Q2_tt_0.pdf");

         }



    TLine *l = new TLine();
           l->SetLineWidth(2);
           l->SetLineColor(kBlack);
  
  
    TCanvas *c6 = new TCanvas("c6","",650,550 );
     for(int ii=0; ii<n_con; ii++){
        
        c6->cd();
           h_tt_xB_Coh_0[ii]   ->Draw("colz");
           h_tt_xB_Coh_0[ii]   ->GetXaxis()->CenterTitle(true);
           h_tt_xB_Coh_0[ii]   ->SetXTitle("-t [GeV^{2}]");
           h_tt_xB_Coh_0[ii]   ->GetXaxis()->SetTitleSize(0.07);
           h_tt_xB_Coh_0[ii]   ->GetYaxis()->CenterTitle(true);
           h_tt_xB_Coh_0[ii]   ->GetYaxis()->SetTitleSize(0.07);
           h_tt_xB_Coh_0[ii]   ->SetYTitle("x_{B}");
              c6->SetLogz();
              c6->SetLogx();
              c6->SetLogy();
                 for(int i=0; i<n_xB+1;i++)
                     l->DrawLine(t_lims_0[0],  xB_lims[i], t_lims_0[n_t_0], xB_lims[i]);
                 for(int i=0; i<n_t_0+1;i++)
                        l->DrawLine(t_lims_0[i], xB_lims[0], t_lims_0[i], xB_lims[n_xB]);
          c6->Print("figs/png/coh_t_xB_0.png");
          c6->Print("figs/pdf/coh_t_xB_0.pdf");
         }
 

    vector<vector<vector<double>>> alu_t_x     ;
    vector<vector<vector<double>>> alu_t_x_err ;
    vector<vector<vector<double>>> Im_t_x     ;
    vector<vector<vector<double>>> Im_t_x_err ;
    vector<vector<vector<double>>> Re_t_x     ;
    vector<vector<vector<double>>> Re_t_x_err ;
    vector<vector<vector<double>>> modle_Im ;
    vector<vector<vector<double>>> modle_Re ;

    vector<vector<vector<double>>> mean_y      ;
    vector<vector<vector<double>>> mean_t_0      ;
    vector<vector<vector<double>>> mean_x      ;
    vector<vector<vector<double>>> mean_Q2     ;
    vector<vector<vector<double>>> mean_y_err  ;
    vector<vector<vector<double>>> mean_t_0_err  ;
    vector<vector<vector<double>>> mean_x_err  ;
    vector<vector<vector<double>>> mean_Q2_err ;

    alu_t_x    .resize(n_con); 
    alu_t_x_err.resize(n_con);
    Im_t_x     .resize(n_con); 
    Im_t_x_err .resize(n_con);
    Re_t_x     .resize(n_con); 
    Re_t_x_err .resize(n_con);
    modle_Im   .resize(n_con);
    modle_Re   .resize(n_con);

    mean_y     .resize(n_con);
    mean_t_0     .resize(n_con);
    mean_x     .resize(n_con);
    mean_Q2    .resize(n_con);
    mean_y_err .resize(n_con);
    mean_t_0_err .resize(n_con);
    mean_x_err .resize(n_con);
    mean_Q2_err.resize(n_con);

    for(int ii=0; ii<n_con; ii++){

       alu_t_x[ii]    .resize(n_xB); 
       alu_t_x_err[ii].resize(n_xB);
       Im_t_x[ii]     .resize(n_xB); 
       Im_t_x_err[ii] .resize(n_xB);
       Re_t_x[ii]     .resize(n_xB); 
       Re_t_x_err[ii].resize(n_xB);
       modle_Im[ii]   .resize(n_xB);
       modle_Re[ii]   .resize(n_xB);

       mean_y[ii]     .resize(n_xB);
       mean_t_0[ii]     .resize(n_xB);
       mean_x[ii]     .resize(n_xB);
       mean_Q2[ii]    .resize(n_xB);
       mean_y_err[ii] .resize(n_xB);
       mean_t_0_err[ii] .resize(n_xB);
       mean_x_err[ii] .resize(n_xB);
       mean_Q2_err[ii].resize(n_xB);
   
       for(int jj=0; jj<n_xB; jj++){

          alu_t_x[ii][jj]     .resize(n_t_0); 
          alu_t_x_err[ii][jj] .resize(n_t_0);
          Im_t_x[ii][jj]      .resize(n_t_0); 
          Im_t_x_err[ii][jj]  .resize(n_t_0);
          Re_t_x[ii][jj]      .resize(n_t_0); 
          Re_t_x_err[ii][jj]  .resize(n_t_0);
          modle_Im[ii][jj]    .resize(n_t_0);
          modle_Re[ii][jj]    .resize(n_t_0);

          mean_y[ii][jj]      .resize(n_t_0);
          mean_t_0[ii][jj]      .resize(n_t_0);
          mean_x[ii][jj]      .resize(n_t_0);
          mean_Q2[ii][jj]     .resize(n_t_0);
          mean_y_err[ii][jj]  .resize(n_t_0);
          mean_t_0_err[ii][jj]  .resize(n_t_0);
          mean_x_err[ii][jj]  .resize(n_t_0);
          mean_Q2_err[ii][jj] .resize(n_t_0);
 
       }
    }


    ifstream infile;
    infile.open("imcff_recff_moh_4he.dat");
    int n_row = 126;
    int n_col = 6;
    double parameters[n_row][n_col]; 
    for(int i =0; i<n_row; i++){
      for(int j =0; j<n_col; j++){
          infile>>parameters[i][j];
      }}
    infile.close();
    
    for(int i =0; i<n_row; i++){
    cout<< parameters[i][0]<<"   "<<parameters[i][1]<<"   "<<parameters[i][2]<<"   "<<parameters[i][3]<<"   "<<parameters[i][4]<<"   "<<parameters[i][5]<<endl;
    }

    int ncc=1;
    ofstream outfile;
    outfile.open ("Q2_xB_t_mean_values.txt");
    outfile << "bin #"<<"  "<<"<Q2>"<<"   "<<"<xB>"<<"   "<<"<-t>"<<"\n";

    for(int ii=0; ii<1; ii++){
       for(int jj=0; jj<n_xB; jj++){
          for(int kk=0; kk<n_t_0; kk++){
             c66->cd();          
             mean_Q2[ii][jj][kk] = h_t_Q2_Coh_0[ii][jj][kk]->GetMean();
             mean_x[ii][jj][kk]  = h_t_xB_Coh_0[ii][jj][kk]->GetMean();
             mean_t_0[ii][jj][kk]  = h_t_t_Coh_0[ii][jj][kk]->GetMean();
             mean_y[ii][jj][kk]  = h_t_y_Coh_0[ii][jj][kk]->GetMean();

             mean_Q2_err[ii][jj][kk] = 0.0;
             mean_x_err[ii][jj][kk] = 0.0;
             mean_t_0_err[ii][jj][kk] = 0.0;
             mean_y_err[ii][jj][kk] = 0.0;
    
             outfile <<ncc<<"  "<<h_t_Q2_Coh_0[ii][jj][kk]->GetMean()<<"   "<<h_t_xB_Coh_0[ii][jj][kk]->GetMean()<<"   "<<h_t_t_Coh_0[ii][jj][kk]->GetMean()<<"\n";    
             modle_Im[ii][jj][kk] = parameters[ncc-1][4]; 
             modle_Re[ii][jj][kk] = parameters[ncc-1][5];
             ncc++;
          }
       }
    }

outfile.close();


   TCanvas *c3 = new TCanvas("c3","",1300,1000 );
            c3->Divide(n_t_0,n_xB,-0.00005,-0.00005); 
            c3->SetGrid();

   // as a function of -t in xB bins
   int counter =1;
   for(int ii=0; ii<1; ii++){
    for(int jj=n_xB-1; jj>=0; jj--){
    //for(int jj=0; jj<n_xB; jj++){
       for(int kk=0; kk<n_t_0; kk++){
            
          cout<< modle_Im[ii][jj][kk]<<"    "<<modle_Re[ii][jj][kk]<<endl; 
            
          c3->SetGrid();
          c3->cd(counter);
          counter++;
         
          h_dvcs_N_p_0[ii][jj][kk]->Sumw2(); 
          h_dvcs_N_m_0[ii][jj][kk]->Sumw2();
         
          TH1* hsum=(TH1*)h_dvcs_N_p_0[ii][jj][kk]->Clone("hsum");
          TH1* hdif=(TH1*)h_dvcs_N_p_0[ii][jj][kk]->Clone("hdif");
          hsum->Add(h_dvcs_N_m_0[ii][jj][kk]);
          hdif->Add(h_dvcs_N_m_0[ii][jj][kk],-1);
         
          TH1* hasy=(TH1*)hdif->Clone("hasy");
               hasy->Divide(hsum);
              
              cout<<"n bins in phi = "<< hasy->GetNbinsX()<<endl;
              for(int mm=1; mm<hasy->GetNbinsX()+1; mm++){
                  double phi_hh = hasy->GetBinCenter(mm);
                  double ALU_phi = -0.1;
                  double Im, Re;
                  //Find_CFF( mean_x[ii][jj][kk], mean_t[ii][jj][kk], Im, Re, jj);
                  Im = modle_Im[ii][jj][kk];
                  Re = modle_Re[ii][jj][kk];
                  //Re = 0.0;

                  Calculate_ALU(mean_Q2[ii][jj][kk], mean_x[ii][jj][kk], -1.0*mean_t_0[ii][jj][kk], mean_y[ii][jj][kk],  phi_hh, Im, Re, ALU_phi); 

                  hasy->SetBinContent(mm, ALU_phi);
                  double bin_error_st = hasy->GetBinError(mm);
                  double bin_error_sys = sys_fac*ALU_phi; 
                  double bin_error_all = bin_error_st + bin_error_sys; 
                  hasy->SetBinError(mm, bin_error_all); 
                  cout<<hasy->GetBinCenter(mm)<<"    "<< hasy->GetBinContent(mm)<<"     "<<bin_error_st<<"     "<<hasy->GetBinError(mm)<<endl;
              }



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
              hasy->GetYaxis()->SetRangeUser(-0.6,0.6);
              if(ii==0) hasy->Draw();
              if(ii==1) hasy->Draw("same");
              if(ii==2) hasy->Draw("same");

          TLatex *l2= new TLatex(10.0,-0.25,Form("%.4f< x_{B} <%.4f",xB_lims[jj], xB_lims[jj+1]));
                  
                  l2->Draw("same");
          TLatex *l3= new TLatex(10.0,-0.4,Form("%.4f< -t <%.4f",t_lims_0[kk], t_lims_0[kk+1]));
                  l3->Draw("same");

          // TF1 *myfit;
          // if(kk<7) myfit = new TF1("myfit","[0]*sin(x*3.1416/180.0)/(1 + [1]*cos(x*3.1416/180.0))",0.0,360.0);
          // else if(kk==7) myfit = new TF1("myfit","[0]*sin(x*3.1416/180.0)",0.0,360.0);
          // myfit->SetLineColor(kRed);
          // myfit->SetLineWidth(2);

          double  A0, A1, A2, A3, c0_BH, c1_BH, c2_BH; 
          Calculate_CFF(mean_Q2[ii][jj][kk], mean_x[ii][jj][kk], -1.0*mean_t_0[ii][jj][kk], mean_y[ii][jj][kk], 
                        A0, A1, A2, A3, c0_BH, c1_BH, c2_BH);

          
          TF1 *ffit;
          TFitResultPtr r;
          ffit = new TF1("ffit",  Form("%0.7f*[0]*sin(x*3.14/180.0) / (%0.7f+ %0.7f*cos(x*3.14/180.0) + %0.7f*cos(2*x*3.14/180.0) + %0.7f*([0]*[0] + [1]*[1]) + %0.7f*[1] + %0.7f*[1]*cos(x*3.14/180.0))",  
                         A0, c0_BH, c1_BH, c2_BH, A1, A2, A3), 0.0, 360.0);
           ffit->SetLineColor(kRed);
           ffit->SetParName(0,"Im(H_{A})");
           ffit->SetParName(1,"Re(H_{A})");


          // fit with a sin and cosin function ------------------------------------
       // hasy->Fit("myfit");
     
        double IM_CFF; ;
        double RE_CFF;;
        Find_CFF(mean_x[ii][jj][kk], abs(mean_t_0[ii][jj][kk]), IM_CFF, RE_CFF, jj);
          
        alu_t_x[ii][jj][kk] = 0.1 *(2-jj)-0.3;
        Im_t_x[ii][jj][kk] = fac_0*modle_Im[ii][jj][kk];
        Re_t_x[ii][jj][kk] = modle_Re[ii][jj][kk];
        //Im_t_x[ii][jj][kk] = IM_CFF;
        //Re_t_x[ii][jj][kk] = RE_CFF;

        //
        // alu_t_x_err[ii][jj][kk] =  myfit->GetParError(0);
        //Im_t_x_err[ii][jj][kk]  = (myfit->GetParError(0)/myfit->GetParameter(0))* abs(IM_CFF);
        //Re_t_x_err[ii][jj][kk]  =  (myfit->GetParError(1)/myfit->GetParameter(1)) *abs(Re_CFF);
        

          // fit with the full form of the CFF

          //ffit->SetParameter(0,27.0);
          //ffit->SetParameter(1,-4.0);
          hasy->Fit("ffit");
          r = hasy->Fit(ffit,"S");
          double x[1] = { 90.0};  double err[1];
          r->GetConfidenceIntervals(1, 1, 1, x, err, 0.683, true);
          double ALU_fit     = ffit->Eval(x[0]);
          double ALU_fit_err = err[0];

          alu_t_x_err[ii][jj][kk] =  ALU_fit_err;
          Im_t_x_err[ii][jj][kk]  = xfac * ffit->GetParError(0); 
          Re_t_x_err[ii][jj][kk]  = xfac * ffit->GetParError(1);          

          }
    }
   }

    c3->Print("figs/png/BSA_Coherent_xB_t_phi_0.png");
    c3->Print("figs/pdf/BSA_Coherent_xB_t_phi_0.pdf");


  vector<double> M_XB;   M_XB.resize(n_xB);
  vector<double> M_Q2;   M_Q2.resize(n_xB);

  for(int ii=0; ii<n_xB; ii++){

    M_XB[ii]= (hh_xB_0[ii]->GetMean());
    M_Q2[ii]= (hh_Q2_0[ii]->GetMean());
   }

//    plot_ALU_projections(mean_t_0, mean_t_0_err, alu_t_x, alu_t_x_err);
//    plot_CFF_projections(mean_t_0, mean_t_0_err, Im_t_x_err, Re_t_x_err);
    FunProfile_0(mean_t_0, mean_t_0_err, M_XB, M_Q2, Im_t_x, Im_t_x_err);




////////   n_t_1 //////////////////////////////////////////////////////////////////
    TCanvas *c661 = new TCanvas("c661","",650,550 );
     for(int ii=0; ii<n_con; ii++){
      
        c661->cd();
           h_Q2_xB_Coh_1[ii]   ->Draw("colz");
           h_Q2_xB_Coh_1[ii]   ->GetYaxis()->CenterTitle(true);
           h_Q2_xB_Coh_1[ii]   ->SetYTitle("Q^{2} [GeV^{2}]");
           h_Q2_xB_Coh_1[ii]   ->GetYaxis()->SetTitleSize(0.07);
           h_Q2_xB_Coh_1[ii]   ->GetXaxis()->CenterTitle(true);
           h_Q2_xB_Coh_1[ii]   ->GetXaxis()->SetTitleSize(0.07);
           h_Q2_xB_Coh_1[ii]   ->SetXTitle("x_{B}");
          c661->SetLogx();
          c661->Print("figs/png/coh_Q2_xB_1.png");
          c661->Print("figs/pdf/coh_Q2_xB_1.pdf");
        
        c661->cd();
           h_Q2_tt_Coh_1[ii]   ->Draw("colz");
           h_Q2_tt_Coh_1[ii]   ->GetXaxis()->CenterTitle(true);
           h_Q2_tt_Coh_1[ii]   ->SetXTitle("Q^{2} [GeV^{2}]");
           h_Q2_tt_Coh_1[ii]   ->GetXaxis()->SetTitleSize(0.07);
           h_Q2_tt_Coh_1[ii]   ->GetYaxis()->CenterTitle(true);
           h_Q2_tt_Coh_1[ii]   ->GetYaxis()->SetTitleSize(0.07);
           h_Q2_tt_Coh_1[ii]   ->SetYTitle("-t [GeV^{2}]");
          c661->SetLogy();
          c661->SetLogx();
          c661->Print("figs/png/coh_Q2_tt_1.png");
          c661->Print("figs/pdf/coh_Q2_tt_1.pdf");

         }



    TCanvas *c61 = new TCanvas("c61","",650,550 );
     for(int ii=0; ii<n_con; ii++){
        
        c61->cd();
           h_tt_xB_Coh_1[ii]   ->Draw("colz");
           h_tt_xB_Coh_1[ii]   ->GetXaxis()->CenterTitle(true);
           h_tt_xB_Coh_1[ii]   ->SetXTitle("-t [GeV^{2}]");
           h_tt_xB_Coh_1[ii]   ->GetXaxis()->SetTitleSize(0.07);
           h_tt_xB_Coh_1[ii]   ->GetYaxis()->CenterTitle(true);
           h_tt_xB_Coh_1[ii]   ->GetYaxis()->SetTitleSize(0.07);
           h_tt_xB_Coh_1[ii]   ->SetYTitle("x_{B}");
              c61->SetLogz();
              c61->SetLogx();
              c61->SetLogy();
                 for(int i=0; i<n_xB+1;i++)
                     l->DrawLine(t_lims_1[0],  xB_lims[i], t_lims_1[n_t_1], xB_lims[i]);
                 for(int i=0; i<n_t_1+1;i++)
                        l->DrawLine(t_lims_1[i], xB_lims[0], t_lims_1[i], xB_lims[n_xB]);
          c61->Print("figs/png/coh_t_xB_1.png");
          c61->Print("figs/pdf/coh_t_xB_1.pdf");
         }
 

    vector<vector<vector<double>>> alu_t_x_1     ;
    vector<vector<vector<double>>> Im_t_x_1     ;
    vector<vector<vector<double>>> Re_t_x_1     ;
    vector<vector<vector<double>>> alu_t_x_1_err ;
    vector<vector<vector<double>>> Im_t_x_1_err ;
    vector<vector<vector<double>>> Re_t_x_1_err ;
    vector<vector<vector<double>>> modle_Im_1 ;
    vector<vector<vector<double>>> modle_Re_1 ;

    vector<vector<vector<double>>> mean_y_1      ;
    vector<vector<vector<double>>> mean_t_1      ;
    vector<vector<vector<double>>> mean_x_1      ;
    vector<vector<vector<double>>> mean_Q2_1     ;
    vector<vector<vector<double>>> mean_t_1_err  ;

    alu_t_x_1    .resize(n_con); 
    Im_t_x_1     .resize(n_con); 
    Re_t_x_1     .resize(n_con); 
    alu_t_x_1_err.resize(n_con);
    Im_t_x_1_err .resize(n_con);
    Re_t_x_1_err .resize(n_con);
    modle_Im_1   .resize(n_con);
    modle_Re_1   .resize(n_con);

    mean_y_1     .resize(n_con);
    mean_t_1     .resize(n_con);
    mean_x_1     .resize(n_con);
    mean_Q2_1    .resize(n_con);
    mean_t_1_err .resize(n_con);

    for(int ii=0; ii<n_con; ii++){

       alu_t_x_1[ii]    .resize(n_xB); 
       Im_t_x_1[ii]     .resize(n_xB); 
       Re_t_x_1[ii]     .resize(n_xB); 
       alu_t_x_1_err[ii].resize(n_xB);
       Im_t_x_1_err[ii] .resize(n_xB);
       Re_t_x_1_err[ii].resize(n_xB);
       modle_Im_1[ii]   .resize(n_xB);
       modle_Re_1[ii]   .resize(n_xB);

       mean_y_1[ii]     .resize(n_xB);
       mean_t_1[ii]     .resize(n_xB);
       mean_x_1[ii]     .resize(n_xB);
       mean_Q2_1[ii]    .resize(n_xB);
       mean_t_1_err[ii] .resize(n_xB);
   
       for(int jj=0; jj<n_xB; jj++){

          alu_t_x_1[ii][jj]     .resize(n_t_1); 
          Im_t_x_1[ii][jj]      .resize(n_t_1); 
          Re_t_x_1[ii][jj]      .resize(n_t_1); 
          alu_t_x_1_err[ii][jj] .resize(n_t_1);
          Im_t_x_1_err[ii][jj]  .resize(n_t_1);
          Re_t_x_1_err[ii][jj]  .resize(n_t_1);
          modle_Im_1[ii][jj]    .resize(n_t_1);
          modle_Re_1[ii][jj]    .resize(n_t_1);

          mean_y_1[ii][jj]      .resize(n_t_1);
          mean_t_1[ii][jj]      .resize(n_t_1);
          mean_x_1[ii][jj]      .resize(n_t_1);
          mean_Q2_1[ii][jj]     .resize(n_t_1);
          mean_t_1_err[ii][jj]  .resize(n_t_1);
 
       }
    }


    ifstream infile_1;
    infile_1.open("imcff_recff_moh_4he_1.dat");
    int n_row_1 = 72;
    int n_col_1 = 6;
    double parameters_1[n_row_1][n_col_1]; 
    for(int i =0; i<n_row_1; i++){
      for(int j =0; j<n_col_1; j++){
          infile_1>>parameters_1[i][j];
      }}
    infile_1.close();
    
    for(int i =0; i<n_row_1; i++){
    cout<< parameters_1[i][0]<<"   "<<parameters_1[i][1]<<"   "<<parameters_1[i][2]<<"   "<<parameters_1[i][3]<<"   "<<parameters_1[i][4]<<"   "<<parameters_1[i][5]<<endl;
    }

    int ncc_1=1;
    for(int ii=0; ii<1; ii++){
       for(int jj=0; jj<n_xB; jj++){
          for(int kk=0; kk<n_t_1; kk++){
             c66->cd();          
             mean_Q2_1[ii][jj][kk] = h_t_Q2_Coh_1[ii][jj][kk]->GetMean();
             mean_x_1[ii][jj][kk]  = h_t_xB_Coh_1[ii][jj][kk]->GetMean();
             mean_t_1[ii][jj][kk]  = h_t_t_Coh_1[ii][jj][kk]->GetMean();
             mean_y_1[ii][jj][kk]  = h_t_y_Coh_1[ii][jj][kk]->GetMean();

             mean_t_1_err[ii][jj][kk] = 0.0;
    
             modle_Im_1[ii][jj][kk] = parameters_1[ncc_1-1][4]; 
             modle_Re_1[ii][jj][kk] = parameters_1[ncc_1-1][5];
             ncc_1++;
          }
       }
    }



   TCanvas *c31 = new TCanvas("c31","",1300,1000 );
            c31->Divide(n_t_1,n_xB,-0.00005,-0.00005); 
            c31->SetGrid();

   // as a function of -t in xB bins
   int counter_1 =1;
   for(int ii=0; ii<1; ii++){
    for(int jj=n_xB-1; jj>=0; jj--){
    //for(int jj=0; jj<n_xB; jj++){
       for(int kk=0; kk<n_t_1; kk++){
            
          cout<< modle_Im_1[ii][jj][kk]<<"    "<<modle_Re_1[ii][jj][kk]<<endl; 
            
          c31->SetGrid();
          c31->cd(counter_1);
          counter_1++;
         
          h_dvcs_N_p_1[ii][jj][kk]->Sumw2(); 
          h_dvcs_N_m_1[ii][jj][kk]->Sumw2();
         
          TH1* hsum=(TH1*)h_dvcs_N_p_1[ii][jj][kk]->Clone("hsum");
          TH1* hdif=(TH1*)h_dvcs_N_p_1[ii][jj][kk]->Clone("hdif");
          hsum->Add(h_dvcs_N_m_1[ii][jj][kk]);
          hdif->Add(h_dvcs_N_m_1[ii][jj][kk],-1);
         
          TH1* hasy=(TH1*)hdif->Clone("hasy");
               hasy->Divide(hsum);
              
              cout<<"n bins in phi = "<< hasy->GetNbinsX()<<endl;
              for(int mm=1; mm<hasy->GetNbinsX()+1; mm++){
                  double phi_hh = hasy->GetBinCenter(mm);
                  double ALU_phi = -0.1;
                  double Im, Re;
                  //Find_CFF( mean_x[ii][jj][kk], mean_t[ii][jj][kk], Im, Re, jj);
                  Im = modle_Im_1[ii][jj][kk];
                  Re = modle_Re_1[ii][jj][kk];
                  //Re = 0.0;

                  Calculate_ALU(mean_Q2_1[ii][jj][kk], mean_x_1[ii][jj][kk], -1.0*mean_t_1[ii][jj][kk], mean_y_1[ii][jj][kk],  phi_hh, Im, Re, ALU_phi); 

                  hasy->SetBinContent(mm, ALU_phi);
                  double bin_error_st = hasy->GetBinError(mm);
                  double bin_error_sys = sys_fac*ALU_phi; 
                  double bin_error_all = bin_error_st + bin_error_sys; 
                  hasy->SetBinError(mm, bin_error_all); 
                  cout<<hasy->GetBinCenter(mm)<<"    "<< hasy->GetBinContent(mm)<<"     "<<bin_error_st<<"     "<<hasy->GetBinError(mm)<<endl;
              
              }



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
              hasy->GetYaxis()->SetRangeUser(-0.6,0.6);
              if(ii==0) hasy->Draw();
              if(ii==1) hasy->Draw("same");
              if(ii==2) hasy->Draw("same");

          TLatex *l2= new TLatex(10.0,-0.2,Form("%.4f< x_{B} <%.4f",xB_lims[jj], xB_lims[jj+1]));
                  l2->SetTextSize(0.07);
                  l2->Draw("same");
          TLatex *l3= new TLatex(10.0,-0.4,Form("%.4f< -t <%.4f",t_lims_1[kk], t_lims_1[kk+1]));
                  l3->SetTextSize(0.07);
                  l3->Draw("same");

          double  A0, A1, A2, A3, c0_BH, c1_BH, c2_BH; 
          Calculate_CFF(mean_Q2_1[ii][jj][kk], mean_x_1[ii][jj][kk], -1.0*mean_t_1[ii][jj][kk], mean_y_1[ii][jj][kk], 
                        A0, A1, A2, A3, c0_BH, c1_BH, c2_BH);

          
          TF1 *ffit;
          TFitResultPtr r;
          ffit = new TF1("ffit",  Form("%0.7f*[0]*sin(x*3.14/180.0) / (%0.7f+ %0.7f*cos(x*3.14/180.0) + %0.7f*cos(2*x*3.14/180.0) + %0.7f*([0]*[0] + [1]*[1]) + %0.7f*[1] + %0.7f*[1]*cos(x*3.14/180.0))",  
                         A0, c0_BH, c1_BH, c2_BH, A1, A2, A3), 0.0, 360.0);
           ffit->SetLineColor(kRed);
           ffit->SetParName(0,"Im(H_{A})");
           ffit->SetParName(1,"Re(H_{A})");


          // fit with a sin and cosin function ------------------------------------
     
        double IM_CFF; ;
        double RE_CFF;;
        Find_CFF(mean_x[ii][jj][kk], abs(mean_t_1[ii][jj][kk]), IM_CFF, RE_CFF, jj);
          
        alu_t_x_1[ii][jj][kk] = 0.1 *(2-jj)-0.3;
        Im_t_x_1[ii][jj][kk] = fac_1*modle_Im_1[ii][jj][kk];
        Re_t_x_1[ii][jj][kk] = modle_Re_1[ii][jj][kk];
        //Im_t_x[ii][jj][kk] = IM_CFF;
        //Re_t_x[ii][jj][kk] = RE_CFF;

        //
        // alu_t_x_err[ii][jj][kk] =  myfit->GetParError(0);
        //Im_t_x_err[ii][jj][kk]  = (myfit->GetParError(0)/myfit->GetParameter(0))* abs(IM_CFF);
        //Re_t_x_err[ii][jj][kk]  =  (myfit->GetParError(1)/myfit->GetParameter(1)) *abs(Re_CFF);
        

          // fit with the full form of the CFF

          //ffit->SetParameter(0,27.0);
          //ffit->SetParameter(1,-4.0);
          hasy->Fit("ffit");
          r = hasy->Fit(ffit,"S");
          double x[1] = { 90.0};  double err[1];
          r->GetConfidenceIntervals(1, 1, 1, x, err, 0.683, true);
          double ALU_fit     = ffit->Eval(x[0]);
          double ALU_fit_err = err[0];

          alu_t_x_1_err[ii][jj][kk] =  ALU_fit_err;
          Im_t_x_1_err[ii][jj][kk]  = xfac * ffit->GetParError(0); 
          Re_t_x_1_err[ii][jj][kk]  = xfac * ffit->GetParError(1);          

          }
    }
   }

    c31->Print("figs/png/BSA_Coherent_xB_t_phi_1.png");
    c31->Print("figs/pdf/BSA_Coherent_xB_t_phi_1.pdf");


  vector<double> M_XB_1;   M_XB_1.resize(n_xB);
  vector<double> M_Q2_1;   M_Q2_1.resize(n_xB);

  for(int ii=0; ii<n_xB; ii++){

    M_XB_1[ii]= (hh_xB_1[ii]->GetMean());
    M_Q2_1[ii]= (hh_Q2_1[ii]->GetMean());
   }

    FunProfile_1(mean_t_1, mean_t_1_err, M_XB_1, M_Q2_1, Im_t_x_1, Im_t_x_1_err);




///////////////   n_t_2  ///////////////////////////////////

    TCanvas *c662 = new TCanvas("c662","",650,550 );
     for(int ii=0; ii<n_con; ii++){
      
        c662->cd();
           h_Q2_xB_Coh_2[ii]   ->Draw("colz");
           h_Q2_xB_Coh_2[ii]   ->GetYaxis()->CenterTitle(true);
           h_Q2_xB_Coh_2[ii]   ->SetYTitle("Q^{2} [GeV^{2}]");
           h_Q2_xB_Coh_2[ii]   ->GetYaxis()->SetTitleSize(0.07);
           h_Q2_xB_Coh_2[ii]   ->GetXaxis()->CenterTitle(true);
           h_Q2_xB_Coh_2[ii]   ->GetXaxis()->SetTitleSize(0.07);
           h_Q2_xB_Coh_2[ii]   ->SetXTitle("x_{B}");
          c662->SetLogx();
          c662->Print("figs/png/coh_Q2_xB_2.png");
          c662->Print("figs/pdf/coh_Q2_xB_2.pdf");
        
        c662->cd();
           h_Q2_tt_Coh_2[ii]   ->Draw("colz");
           h_Q2_tt_Coh_2[ii]   ->GetXaxis()->CenterTitle(true);
           h_Q2_tt_Coh_2[ii]   ->SetXTitle("Q^{2} [GeV^{2}]");
           h_Q2_tt_Coh_2[ii]   ->GetXaxis()->SetTitleSize(0.07);
           h_Q2_tt_Coh_2[ii]   ->GetYaxis()->CenterTitle(true);
           h_Q2_tt_Coh_2[ii]   ->GetYaxis()->SetTitleSize(0.07);
           h_Q2_tt_Coh_2[ii]   ->SetYTitle("-t [GeV^{2}]");
          c662->SetLogy();
          c662->SetLogx();
          c662->Print("figs/png/coh_Q2_tt_2.png");
          c662->Print("figs/pdf/coh_Q2_tt_2.pdf");

         }


    TCanvas *c62 = new TCanvas("c62","",650,550 );
     for(int ii=0; ii<n_con; ii++){
        
        c62->cd();
           h_tt_xB_Coh_2[ii]   ->Draw("colz");
           h_tt_xB_Coh_2[ii]   ->GetXaxis()->CenterTitle(true);
           h_tt_xB_Coh_2[ii]   ->SetXTitle("-t [GeV^{2}]");
           h_tt_xB_Coh_2[ii]   ->GetXaxis()->SetTitleSize(0.07);
           h_tt_xB_Coh_2[ii]   ->GetYaxis()->CenterTitle(true);
           h_tt_xB_Coh_2[ii]   ->GetYaxis()->SetTitleSize(0.07);
           h_tt_xB_Coh_2[ii]   ->SetYTitle("x_{B}");
              c62->SetLogz();
              c62->SetLogx();
              c62->SetLogy();
                 for(int i=0; i<n_xB+1;i++)
                     l->DrawLine(t_lims_2[0],  xB_lims[i], t_lims_2[n_t_2], xB_lims[i]);
                 for(int i=0; i<n_t_2+1;i++)
                        l->DrawLine(t_lims_2[i], xB_lims[0], t_lims_2[i], xB_lims[n_xB]);
          c62->Print("figs/png/coh_t_xB_2.png");
          c62->Print("figs/pdf/coh_t_xB_2.pdf");
         }
 
   vector<vector<vector<double>>> alu_t_x_2     ;
    vector<vector<vector<double>>> Im_t_x_2     ;
    vector<vector<vector<double>>> Re_t_x_2     ;
    vector<vector<vector<double>>> alu_t_x_2_err ;
    vector<vector<vector<double>>> Im_t_x_2_err ;
    vector<vector<vector<double>>> Re_t_x_2_err ;
    vector<vector<vector<double>>> modle_Im_2 ;
    vector<vector<vector<double>>> modle_Re_2 ;

    vector<vector<vector<double>>> mean_y_2      ;
    vector<vector<vector<double>>> mean_t_2      ;
    vector<vector<vector<double>>> mean_x_2      ;
    vector<vector<vector<double>>> mean_Q2_2     ;
    vector<vector<vector<double>>> mean_t_2_err  ;

    alu_t_x_2    .resize(n_con); 
    Im_t_x_2     .resize(n_con); 
    Re_t_x_2     .resize(n_con); 
    alu_t_x_2_err.resize(n_con);
    Im_t_x_2_err .resize(n_con);
    Re_t_x_2_err .resize(n_con);
    modle_Im_2   .resize(n_con);
    modle_Re_2   .resize(n_con);

    mean_y_2     .resize(n_con);
    mean_t_2     .resize(n_con);
    mean_x_2     .resize(n_con);
    mean_Q2_2    .resize(n_con);
    mean_t_2_err .resize(n_con);

    for(int ii=0; ii<n_con; ii++){

       alu_t_x_2[ii]    .resize(n_xB); 
       Im_t_x_2[ii]     .resize(n_xB); 
       Re_t_x_2[ii]     .resize(n_xB); 
       alu_t_x_2_err[ii].resize(n_xB);
       Im_t_x_2_err[ii] .resize(n_xB);
       Re_t_x_2_err[ii].resize(n_xB);
       modle_Im_2[ii]   .resize(n_xB);
       modle_Re_2[ii]   .resize(n_xB);

       mean_y_2[ii]     .resize(n_xB);
       mean_t_2[ii]     .resize(n_xB);
       mean_x_2[ii]     .resize(n_xB);
       mean_Q2_2[ii]    .resize(n_xB);
       mean_t_2_err[ii] .resize(n_xB);
   
       for(int jj=0; jj<n_xB; jj++){

          alu_t_x_2[ii][jj]     .resize(n_t_2); 
          Im_t_x_2[ii][jj]      .resize(n_t_2); 
          Re_t_x_2[ii][jj]      .resize(n_t_2); 
          alu_t_x_2_err[ii][jj] .resize(n_t_2);
          Im_t_x_2_err[ii][jj]  .resize(n_t_2);
          Re_t_x_2_err[ii][jj]  .resize(n_t_2);
          modle_Im_2[ii][jj]    .resize(n_t_2);
          modle_Re_2[ii][jj]    .resize(n_t_2);

          mean_y_2[ii][jj]      .resize(n_t_2);
          mean_t_2[ii][jj]      .resize(n_t_2);
          mean_x_2[ii][jj]      .resize(n_t_2);
          mean_Q2_2[ii][jj]     .resize(n_t_2);
          mean_t_2_err[ii][jj]  .resize(n_t_2);
 
       }
    }


    ifstream infile_2;
    infile_2.open("imcff_recff_moh_4he_2.dat");
    int n_row_2 = 45;
    int n_col_2 = 6;
    double parameters_2[n_row_2][n_col_2]; 
    for(int i =0; i<n_row_2; i++){
      for(int j =0; j<n_col_2; j++){
          infile_2>>parameters_2[i][j];
      }}
    infile_2.close();
    
    for(int i =0; i<n_row_2; i++){
    cout<< parameters_2[i][0]<<"   "<<parameters_2[i][1]<<"   "<<parameters_2[i][2]<<"   "<<parameters_2[i][3]<<"   "<<parameters_2[i][4]<<"   "<<parameters_2[i][5]<<endl;
    }

    int ncc_2=1;
    for(int ii=0; ii<1; ii++){
       for(int jj=0; jj<n_xB; jj++){
          for(int kk=0; kk<n_t_2; kk++){
             c66->cd();          
             mean_Q2_2[ii][jj][kk] = h_t_Q2_Coh_2[ii][jj][kk]->GetMean();
             mean_x_2[ii][jj][kk]  = h_t_xB_Coh_2[ii][jj][kk]->GetMean();
             mean_t_2[ii][jj][kk]  = h_t_t_Coh_2[ii][jj][kk]->GetMean();
             mean_y_2[ii][jj][kk]  = h_t_y_Coh_2[ii][jj][kk]->GetMean();

             mean_t_2_err[ii][jj][kk] = 0.0;
    
             modle_Im_2[ii][jj][kk] = parameters_2[ncc_2-1][4]; 
             modle_Re_2[ii][jj][kk] = parameters_2[ncc_2-1][5];
             ncc_2++;
          }
       }
    }



   TCanvas *c32 = new TCanvas("c32","",1300,1000 );
            c32->Divide(n_t_2,n_xB,-0.00005,-0.00005); 
            c32->SetGrid();

   // as a function of -t in xB bins
   int counter_2 =1;
   for(int ii=0; ii<1; ii++){
    for(int jj=n_xB-1; jj>=0; jj--){
    //for(int jj=0; jj<n_xB; jj++){
       for(int kk=0; kk<n_t_2; kk++){
            
          cout<< modle_Im_2[ii][jj][kk]<<"    "<<modle_Re_2[ii][jj][kk]<<endl; 
            
          c32->SetGrid();
          c32->cd(counter_2);
          counter_2++;
         
          h_dvcs_N_p_2[ii][jj][kk]->Sumw2(); 
          h_dvcs_N_m_2[ii][jj][kk]->Sumw2();
         
          TH1* hsum=(TH1*)h_dvcs_N_p_2[ii][jj][kk]->Clone("hsum");
          TH1* hdif=(TH1*)h_dvcs_N_p_2[ii][jj][kk]->Clone("hdif");
          hsum->Add(h_dvcs_N_m_2[ii][jj][kk]);
          hdif->Add(h_dvcs_N_m_2[ii][jj][kk],-1);
         
          TH1* hasy=(TH1*)hdif->Clone("hasy");
               hasy->Divide(hsum);
              
              cout<<"n bins in phi = "<< hasy->GetNbinsX()<<endl;
              for(int mm=1; mm<hasy->GetNbinsX()+1; mm++){
                  double phi_hh = hasy->GetBinCenter(mm);
                  double ALU_phi = -0.1;
                  double Im, Re;
                  //Find_CFF( mean_x[ii][jj][kk], mean_t[ii][jj][kk], Im, Re, jj);
                  Im = modle_Im_2[ii][jj][kk];
                  Re = modle_Re_2[ii][jj][kk];
                  //Re = 0.0;

                  Calculate_ALU(mean_Q2_2[ii][jj][kk], mean_x_2[ii][jj][kk], -1.0*mean_t_2[ii][jj][kk], mean_y_2[ii][jj][kk],  phi_hh, Im, Re, ALU_phi); 

                  hasy->SetBinContent(mm, ALU_phi);
                  double bin_error_st = hasy->GetBinError(mm);
                  double bin_error_sys = sys_fac*ALU_phi; 
                  double bin_error_all = bin_error_st + bin_error_sys; 
                  hasy->SetBinError(mm, bin_error_all); 
                  cout<<hasy->GetBinCenter(mm)<<"    "<< hasy->GetBinContent(mm)<<"     "<<bin_error_st<<"     "<<hasy->GetBinError(mm)<<endl;
 
              }



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
              hasy->GetYaxis()->SetRangeUser(-0.6,0.6);
              if(ii==0) hasy->Draw();
              if(ii==1) hasy->Draw("same");
              if(ii==2) hasy->Draw("same");

          TLatex *l2= new TLatex(10.0,-0.2,Form("%.4f< x_{B} <%.4f",xB_lims[jj], xB_lims[jj+1]));
                  l2->SetTextSize(0.09);
                  l2->Draw("same");
          TLatex *l3= new TLatex(10.0,-0.4,Form("%.4f< -t <%.4f",t_lims_2[kk], t_lims_2[kk+1]));
                  l3->SetTextSize(0.09);
                  l3->Draw("same");

          double  A0, A1, A2, A3, c0_BH, c1_BH, c2_BH; 
          Calculate_CFF(mean_Q2_2[ii][jj][kk], mean_x_2[ii][jj][kk], -1.0*mean_t_2[ii][jj][kk], mean_y_2[ii][jj][kk], 
                        A0, A1, A2, A3, c0_BH, c1_BH, c2_BH);

          
          TF1 *ffit;
          TFitResultPtr r;
          ffit = new TF1("ffit",  Form("%0.7f*[0]*sin(x*3.14/180.0) / (%0.7f+ %0.7f*cos(x*3.14/180.0) + %0.7f*cos(2*x*3.14/180.0) + %0.7f*([0]*[0] + [1]*[1]) + %0.7f*[1] + %0.7f*[1]*cos(x*3.14/180.0))",  
                         A0, c0_BH, c1_BH, c2_BH, A1, A2, A3), 0.0, 360.0);
           ffit->SetLineColor(kRed);
           ffit->SetParName(0,"Im(H_{A})");
           ffit->SetParName(1,"Re(H_{A})");


          // fit with a sin and cosin function ------------------------------------
     
        double IM_CFF; ;
        double RE_CFF;;
        Find_CFF(mean_x[ii][jj][kk], abs(mean_t_0[ii][jj][kk]), IM_CFF, RE_CFF, jj);
          
        alu_t_x_2[ii][jj][kk] = 0.1 *(2-jj)-0.3;
        Im_t_x_2[ii][jj][kk] = fac_2*modle_Im_2[ii][jj][kk];
        Re_t_x_2[ii][jj][kk] = modle_Re_2[ii][jj][kk];
        //Im_t_x[ii][jj][kk] = IM_CFF;
        //Re_t_x[ii][jj][kk] = RE_CFF;

        //
        // alu_t_x_err[ii][jj][kk] =  myfit->GetParError(0);
        //Im_t_x_err[ii][jj][kk]  = (myfit->GetParError(0)/myfit->GetParameter(0))* abs(IM_CFF);
        //Re_t_x_err[ii][jj][kk]  =  (myfit->GetParError(1)/myfit->GetParameter(1)) *abs(Re_CFF);
        

          // fit with the full form of the CFF

          //ffit->SetParameter(0,27.0);
          //ffit->SetParameter(1,-4.0);
          hasy->Fit("ffit");
          r = hasy->Fit(ffit,"S");
          double x[1] = { 90.0};  double err[1];
          r->GetConfidenceIntervals(1, 1, 1, x, err, 0.683, true);
          double ALU_fit     = ffit->Eval(x[0]);
          double ALU_fit_err = err[0];

          alu_t_x_2_err[ii][jj][kk] =  ALU_fit_err;
          Im_t_x_2_err[ii][jj][kk]  = xfac * ffit->GetParError(0); 
          Re_t_x_2_err[ii][jj][kk]  = xfac * ffit->GetParError(1);          

          }
    }
   }

    c32->Print("figs/png/BSA_Coherent_xB_t_phi_2.png");
    c32->Print("figs/pdf/BSA_Coherent_xB_t_phi_2.pdf");


  vector<double> M_XB_2;   M_XB_2.resize(n_xB);
  vector<double> M_Q2_2;   M_Q2_2.resize(n_xB);

  for(int ii=0; ii<n_xB; ii++){

    M_XB_2[ii]= (hh_xB_2[ii]->GetMean());
    M_Q2_2[ii]= (hh_Q2_2[ii]->GetMean());
   }

    FunProfile_2(mean_t_2, mean_t_2_err, M_XB_2, M_Q2_2, Im_t_x_2, Im_t_x_2_err);







}




void Find_CFF( double xB, double t, double &Im, double &Re, int n_xB){
             Im  = (1 + 0.1/(1+n_xB))*(50.77 - pow(5.08 *t, 6)) * exp( -13.2*t );
             Re  = (-8.62 + 15.86*t - pow( -3.59*t, 3)) * exp( -7.27*t);
  
}




void Calculate_CFF(double Q2, double xB, double t, double y, 
                   double &A0, double &A1, double &A2, double &A3, 
                   double &c0_BH, double &c1_BH, double &c2_BH)
  {

  double PI = 3.1416;
  double phi = 90.0*PI/180.0;

  double MPROT = 0.93827;
  double MALPH = 3.7274;
  double EBEAM = 18.0;

  double xA = xB*MPROT/MALPH;
  double xA1= 1 - xA;
  //double y = Q2/2/MALPH/xA/EBEAM;
  double e = 2*xA*MALPH/sqrt(Q2); // epsilon
  double e2 = e*e;
  double T0 = -Q2 * (2*xA1*(1-sqrt(1+e2))+e2) / (4*xA*xA1 + e2);

  // kinematical factors
  double J = (1-y-y*e2/2) * (1+t/Q2) - (1-xA)*(2-y)*t/Q2;
  double dt = (t - T0)/Q2;
  double K_hat = sqrt(T0 - t) * sqrt(xA1*sqrt(1+e2) + (T0 - t)*(e2 + 4*xA1*xA)/(4*Q2) );
  double K2 = -1.0*dt * xA1 * (1 -y - y*y*e2/4) * (sqrt(1+e2) + ((4*xA*xA1+e2)*dt/(4*xA1)) );
  double K = sqrt(1 - y + e2*y*y/4)* (K_hat)/(sqrt(Q2));
  //double K2 = K*K;
  //double K = sqrt(-1.0*K2);

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
   A1 = xA*xA*t*(pow((1+e2),2)/Q2) * P1_phi * P2_phi * C_DVCS_0; 
   A2 = xA*pow(1+e2,2) * C_INT_plus_plus_0 /y;  
   A3 = xA* pow(1+e2,2) * C_INT_plus_plus_1/y; 

  // ALU2 =  A0* Im* sin(phi)/ (c0_BH + c1_BH*cos(phi) + c2_BH*cos(2*phi) + A1*(Re*Re + Im*Im) + A2*Re + A3*Re * cos(phi));

}

int Find_Run_Number(const char* input_file_name){

int RunNumber= 0;
int fullNameLength = strlen(input_file_name);
string StringRunNumber;

   int  first_numberPosition = fullNameLength - 8;

   for (int i = first_numberPosition; i<first_numberPosition+3; i++)
     StringRunNumber += input_file_name[i];

     std::istringstream iss(StringRunNumber);
     iss >> RunNumber;
return RunNumber;
}


void Calculate_ALU(double Q2, double xB, double t, double y, double phi , double Im, double Re, double &ALU)
  {

  double PI = 3.1416;
   phi = phi*PI/180.0;

  double MPROT = 0.93827;
  double MALPH = 3.7274;
  double EBEAM = 18.0;
  double EHad = 220.0; 

  double xA = xB*MPROT/MALPH;
  double xA1= 1 - xA;
  //double y = Q2/2/MALPH/xA/EBEAM;

  double s = MALPH*MALPH +2.0*EBEAM*(EHad+sqrt(EHad*EHad-MALPH*MALPH));
  double yy = Q2/(s-MALPH*MALPH)/xA;
  //cout<< y<<"     "<<yy<<endl;
  double e = 2*xA*MALPH/sqrt(Q2); // epsilon
  double e2 = e*e;
  double T0 = -Q2 * (2*xA1*(1-sqrt(1+e2))+e2) / (4*xA*xA1 + e2);

  // kinematical factors
  double J = (1-y-y*e2/2) * (1+t/Q2) - (1-xA)*(2-y)*t/Q2;
  double dt = (t - T0)/Q2;
  double K_hat = sqrt(T0 - t) * sqrt(xA1*sqrt(1+e2) + (T0 - t)*(e2 + 4*xA1*xA)/(4*Q2) );
  double K2 = -1.0*dt * xA1 * (1 -y - y*y*e2/4) * (sqrt(1+e2) + ((4*xA*xA1+e2)*dt/(4*xA1)) );
  double K = sqrt(1 - y + e2*y*y/4)* (K_hat)/(sqrt(Q2));
  //double K2 = K*K;
  //double K = sqrt(-1.0*K2);

  // BH propagators
  double P1_phi = -1.0*(J + 2*K*cos(PI-phi)) / (y * (1+e2));
  double P2_phi =  1 + t/Q2 + (1/(y*(1+e*e))) * (J + 2*K*cos(PI-phi));

  
  // Helium form factor
  double a=0.316;
  double b=0.681;
  double FF4He = (1-pow(a*a*abs(t)/pow(0.197327,2),6))*exp(-b*b*abs(t)/pow(0.197327,2)); //  1e-15;
  //double FF4He = (1-pow(a*a*Q2,6))*exp(-b*b*Q2);
 
  // BH fourier coefficients
  double c0_BH = ( (pow(2-y,2)+pow(y * (1+e2),2)) * (e2*Q2/t + 4*xA1 + (4*xA+e2)*t/Q2) 
                   + 2*e2*(4*(1-y)*(3+2*e2) + y*y*(2-e2*e2))
                   - 4*xA*xA*pow(2-y,2)*(2+e2)*t/Q2
                   + 8*K2*e2*Q2/t) * pow(FF4He,2);
  double c1_BH = -8*(2-y)* K * (2*xA + e2 - e2*Q2/t) * pow(FF4He,2);
  double c2_BH =  8*K2*e2*Q2*pow(FF4He,2)/t;


  // redefine the fourier hamoinic to fit for Re and Im of HA 
  double C_DVCS_0 = 2*((2-2*y+y*y + 0.5*e2*y*y)/(1+e2)) ;

  double C_INT_plus_plus_0 = ( -4*(2-y)*(1 + sqrt(1+e2))/(pow(1+e2 ,2))) * ( (pow(K_hat*(2-y) ,2))/(Q2*sqrt(1+e2)) +  (t/Q2)*(1-y-y*y*e2/4)*(2-xA)*( 1 + (2*xA*(t/Q2)*(2-xA + ((sqrt(1+e2) -1)/(2)) + ((e2/(2*xA))) ) + e2 )/((2-xA)*(1+sqrt(1+e2)))  )) * FF4He;  
  
  double C_INT_plus_plus_1 = ((-16*K*(1 - y - e2*y*y/4))/pow(1+e2 ,5/2)) * ( ( 1 + (1-xA)*((sqrt(1+e2) -1) /(2*xA)) + e2/(4*xA))* (xA*t/Q2)  - 3.0*e2/4.0) - 4.0*K*(2-2*y+y*y + e2*y*y/2) * ( (1+sqrt(1+e2) -e2)/pow(1+e2 ,5/2) ) *(1 - (1-3.0*xA)*t/Q2 + (1-sqrt(1+e2)+3*e2)/(1+sqrt(1+e2) -e2)*(xA*t/Q2)) * FF4He;  


  double S_INT_plus_plus_1 = (8*K*(2-y)*y /(1+e2)) * ( 1 + ((1-xA+0.5*(sqrt(1+e2)-1))/(1+e2))*dt ) * FF4He;  


  double A0 = xA* pow(1+e2, 2) * S_INT_plus_plus_1/y;
  double A1 = xA*xA*t*(pow((1+e2),2)/Q2) * P1_phi * P2_phi * C_DVCS_0; 
  double A2 = xA*pow(1+e2,2) * C_INT_plus_plus_0 /y;  
  double A3 = xA* pow(1+e2,2) * C_INT_plus_plus_1/y; 

   ALU =  A0* Im* sin(phi)/ (c0_BH + c1_BH*cos(phi) + c2_BH*cos(2*phi) + A1*(Re*Re + Im*Im) + A2*Re + A3*Re * cos(phi));
}



