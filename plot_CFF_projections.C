
void plot_CFF_projections( vector<vector<vector<double>>> mean_t, 
                           vector<vector<vector<double>>> mean_t_err,
                           vector<vector<vector<double>>> Im_t_x_err,
                           vector<vector<vector<double>>> Re_t_x_err  )
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

  TCanvas *ccc = new TCanvas("ccc", "",0,0,750,600);
  TCanvas *cc = new TCanvas("cc", "",0,0,750,600);

   ccc->cd();
   ccc->SetLogx();
  
   double xfont = 0.14;
   double lfont = 0.08;
   TMultiGraph *mgCFF = new TMultiGraph();
   TGraphErrors *greCFF;;
   TGraph *graphCFF;
   TGraph *ggraphCFF;
   mgCFF = new TMultiGraph();

   double fx11[3] = { 0.08, 0.094, 0.127};
   double fy11[3] = { 35.1955 ,15.1153 ,18.2116};
   double fex11[3] = { 0, 0, 0};
   double fey11[3] = { 11.8244 ,4.29386 ,7.54076};

   greCFF = new TGraphErrors(3, fx11, fy11, fex11, fey11);
   greCFF->SetFillColor(1);
   greCFF->SetMarkerColor(kBlack);
   greCFF->SetLineColor(kBlack); 
   greCFF->SetMarkerStyle(21);
   greCFF->SetMarkerSize(2.0);
   greCFF->SetLineWidth(4);


cout<<"hello   "<<endl;
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
   mgCFF->GetXaxis()->SetNdivisions(105);
   mgCFF->GetXaxis()->SetLimits(0.04,1.0);
   mgCFF->GetYaxis()->SetRangeUser(-25.0, 50.0);

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

   TGraph*SgraphCFF = new TGraph(29);
          SgraphCFF->SetName("graphCFF");
          SgraphCFF->SetTitle("graphCFF");
          SgraphCFF->SetFillColor(1);
          SgraphCFF->SetLineColor(4);
          SgraphCFF->SetLineStyle(9);
          SgraphCFF->SetLineWidth(4);
        //SgraphCFF->SetPoint(,  0.0385,  16.59 );
        //SgraphCFF->SetPoint(,  0.0435,  22.46 );
          SgraphCFF->SetPoint(0 ,  0.0485,  24.85 );
          SgraphCFF->SetPoint(1 ,  0.0535,  24.76 );
          SgraphCFF->SetPoint(2 ,  0.0585,  23.29 );
          SgraphCFF->SetPoint(3 ,  0.0635,  22.29 );
          SgraphCFF->SetPoint(4 ,  0.0685,  20.98 );
          SgraphCFF->SetPoint(5 ,  0.0735,  19.85 );
          SgraphCFF->SetPoint(6 ,  0.0785,  18.68 );
          SgraphCFF->SetPoint(7 ,  0.0835,  17.58 );
          SgraphCFF->SetPoint(8 , 0.0885,  16.59 );
          SgraphCFF->SetPoint(9 , 0.0935,  15.61 );
          SgraphCFF->SetPoint(10, 0.0985,  14.69 );
          SgraphCFF->SetPoint(11, 0.1035,  13.82 );
          SgraphCFF->SetPoint(12, 0.1085,  13.04 );
          SgraphCFF->SetPoint(13, 0.1135,  12.26 );
          SgraphCFF->SetPoint(14, 0.1185,  11.54 );
          SgraphCFF->SetPoint(15, 0.1235,  10.86 );
          SgraphCFF->SetPoint(16, 0.1285,  10.21 );
          SgraphCFF->SetPoint(17, 0.1335,  9.625 );
          SgraphCFF->SetPoint(18, 0.1385,  9.054 );
          SgraphCFF->SetPoint(19, 0.1435,  8.516 );
          SgraphCFF->SetPoint(20, 0.1485,  8.009 );
          SgraphCFF->SetPoint(21, 0.1535,  7.531 );
          SgraphCFF->SetPoint(22, 0.1585,  7.081 );
          SgraphCFF->SetPoint(23, 0.1635,  6.667 );
          SgraphCFF->SetPoint(24, 0.1685,  6.267 );
          SgraphCFF->SetPoint(25, 0.1735,  5.890 );
          SgraphCFF->SetPoint(26, 0.1785,  5.534 );
          SgraphCFF->SetPoint(27, 0.1835,  5.199 );
          SgraphCFF->SetPoint(28, 0.1885,  4.883 );
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

    double xB_lims[4] = {0.05, 0.17, 0.23, 0.5};
    const int n_tt = 7;
    for(int i=0; i<1; i++){
      for(int j=0; j<3; j++){
      
       Double_t tt[n_tt];
       Double_t tt_err[n_tt];
       Double_t IM[n_tt];
       Double_t IM_err[n_tt];
   
       for(int k=0; k<n_tt; k++){
         
         tt[k]     = mean_t[i][j][k];
         if(j==1) tt[k]     = mean_t[i][j][k] -0.005;
         tt_err[k] = mean_t_err[i][j][k];
         IM[k]    = 7 *(2-j) - 15;
         IM_err[k]= Im_t_x_err[i][j][k];
       }
   
         TGraphErrors * gre = new TGraphErrors(8,tt,IM,tt_err,IM_err);
         gre->SetMarkerStyle(24);
         gre->SetMarkerSize(1.3);
         gre->SetLineWidth(5);
         if(j==0) { gre->SetMarkerColor(kBlack);  gre->SetLineColor(kBlack); }
         if(j==1) { gre->SetMarkerColor(kBlue);   gre->SetLineColor(kBlue);   }
         if(j==2) { gre->SetMarkerColor(kRed);    gre->SetLineColor(kRed);  }
         gre->Draw("psame");
             
         TLatex *l2= new TLatex(0.38, 7 *(2-j) - 16,Form("%.2f< x_{B} <%.2f",xB_lims[j], xB_lims[j+1]));
                 l2->SetTextSize(0.045);
                 if(j==0) { l2->SetTextColor(kBlack);  }
                 if(j==1) { l2->SetTextColor(kBlue);   }
                 if(j==2) { l2->SetTextColor(kRed);    }
                 l2->Draw("same");
        }
      }





// // real part as a function of -t
   cc->cd();
   cc->SetLogx();
 
   mgCFF = new TMultiGraph();
   
   Double_t graphCFF_fx12[3] = { 0.08, 0.094, 0.127};
   Double_t graphCFF_fy12[3] = {-4.88354 ,19.6659 ,-2.28985};   
   Double_t graphCFF_fex12[3] = { 0, 0, 0};
   Double_t graphCFF_fey12[3] = {22.4969 ,11.1467 ,16.0775};

   greCFF = new TGraphErrors(3,graphCFF_fx12,graphCFF_fy12,graphCFF_fex12,graphCFF_fey12);
   greCFF->SetFillColor(1);
   greCFF->SetMarkerColor(kBlack);
   greCFF->SetLineColor(kBlack); 
   greCFF->SetMarkerStyle(21);
   greCFF->SetMarkerSize(2.0);
   greCFF->SetLineWidth(4);

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
   mgCFF->GetYaxis()->SetRangeUser(-60.0, 40.0);
   mgCFF->GetXaxis()->SetLimits(0.04,1.0);
  

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

   ccc->SetLogx(); 

   TLegend *leg_cff_Re = new TLegend(0.55,0.78,0.9,0.98);
            leg_cff_Re->SetNColumns(1);
            leg_cff_Re->SetTextSize(0.04);  
            leg_cff_Re->AddEntry(greCFF,"CLAS-EG6","P");   
            leg_cff_Re->AddEntry(graphCFF,"Convolution-Dual","L");
            leg_cff_Re->AddEntry(ggraphCFF,"Convolution-VGG","L");
          
   leg_cff_Re->Draw();

  for(int i=0; i<1; i++){
   for(int j=0; j<3; j++){
   
      Double_t tt[8];
      Double_t tt_err[8];
      Double_t RE[8];
      Double_t RE_err[8];

      for(int k=0; k<8; k++){
         
         tt[k]     = mean_t[i][j][k];
         if(j==1) tt[k]     = mean_t[i][j][k] -0.005;
         tt_err[k] = mean_t_err[i][j][k];
         RE[k]    = 8 *(2-j) - 50;
         RE_err[k]= Re_t_x_err[i][j][k];
      }

      TGraphErrors *gre = new TGraphErrors(8, tt, RE, tt_err, RE_err);
      gre->SetMarkerStyle(24);
      gre->SetMarkerSize(1.3);
      gre->SetLineWidth(5);
      if(j==0) { gre->SetMarkerColor(kBlack); gre->SetLineColor(kBlack); }
      if(j==1) { gre->SetMarkerColor(kBlue);   gre->SetLineColor(kBlue);   }
      if(j==2) { gre->SetMarkerColor(kRed);  gre->SetLineColor(kRed);  }
      gre->Draw("psame");
          
      TLatex *l2= new TLatex(0.38, 8 *(2-j)-52 ,Form("%.2f< x_{B} <%.2f",xB_lims[j], xB_lims[j+1]));
              l2->SetTextSize(0.045);
              if(j==0) { l2->SetTextColor(kBlack);  }
              if(j==1) { l2->SetTextColor(kBlue);   }
              if(j==2) { l2->SetTextColor(kRed);    }
              l2->Draw("same");

     }
   }

   ccc->Print("fig/png/Coherent_CFF_Im.png");
   ccc->Print("fig/pdf/Coherent_CFF_Im.pdf");

   cc->Print("fig/png/Coherent_CFF_Re.png");
   cc->Print("fig/pdf/Coherent_CFF_Re.pdf");

}
