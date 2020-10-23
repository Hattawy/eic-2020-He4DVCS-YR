#include "TRandom3.h"
#include <TH1D.h>
TRandom3 RANDOM;

void combine_one_bin_function_integration(){

         gStyle->SetLabelSize(0.03,"xyz"); // size of axis value font
         gStyle->SetTitleSize(0.035,"xyz"); // size of axis title font
         gStyle->SetTitleFont(22,"xyz"); // font option
         gStyle->SetLabelFont(22,"xyz");
         gStyle->SetTitleOffset(1.2,"y");
         gStyle->SetCanvasBorderMode(0);
         gStyle->SetCanvasBorderSize(0);
         gStyle->SetPadBottomMargin(0.16); //margins...
         gStyle->SetPadTopMargin(0.16);
         gStyle->SetPadLeftMargin(0.16);
         gStyle->SetPadRightMargin(0.08);
         gStyle->SetFrameBorderMode(0);
         gStyle->SetPaperSize(20,24);
         gStyle->SetLabelSize(0.05,"xy");
         gStyle->SetTitleSize(0.06,"xy");


   const int np = 50;
    TH1D *hh[np]; 
    for( int j=0; j<np; j++){
        hh[j] = new TH1D(Form("hh[%d]",j),Form("point b_{#perp} = %f ", 0.1*j+0.05) ,150, -0.1, 0.1);
       }

   double b_perp[np];
   double b_perp_err[np];
   double rho_val[np];
   double rho_val_err[np];

   TLegend* leg = new TLegend(0.65,0.65,0.84,0.84);
   leg-> SetNColumns(1);

   TCanvas *c3 = new TCanvas("c3","",750,600 ); 
            c3->cd(); 
            //c3->SetGrid();

   TMultiGraph * mg = new TMultiGraph();
                 mg->SetTitle("Projected charge profile precisions");

   auto f = [](double *x, double *p) {
          double t = x[0]/0.197;
          double b = p[0]/0.197;
          double p0 = p[1];
          double p1 = p[2];
          double Delta = TMath::Sqrt(t);
          double ImH = p0*(1-TMath::Power(pow(0.316/0.197,2)*t,6))*TMath::Exp(p1*t);
          //double integrand = ImH*TMath::BesselJ0(b*Delta)*Delta/(2 *(TMath::Pi()));
          double integrand = ImH*TMath::BesselJ0(b*Delta)/(4.0 *(TMath::Pi()));

          return integrand;
       };
 
    
    TGraph * gr0;
    
    for (int m =0; m<3; m++)
     {   
      for (int i =0; i<3000; i++)
       {
          if (i% 300 == 0) printf("still running %d \n",(int)i);
       double xpp[3][3] = { 1.0, 130.08, -14.127, 
                            1.0, 54.77, -15.67,
                            1.0, 17.11, -11.34  };

       double xpp_sig[3][3] = { 0.0, 25.5, 3.5, 
                                0.0, 8.92, 1.98,
                                0.0, 3.13, 1.65  };
 

       double pp[3][3] = {1.0, RANDOM.Gaus(xpp[0][1],xpp_sig[0][1]),  RANDOM.Gaus(xpp[0][2],xpp_sig[0][2]),
                          1.0, RANDOM.Gaus(xpp[1][1],xpp_sig[1][1]),  RANDOM.Gaus(xpp[1][2],xpp_sig[1][2]),
                          1.0, RANDOM.Gaus(xpp[2][1],xpp_sig[2][1]),  RANDOM.Gaus(xpp[2][2],xpp_sig[2][2])
                          };

          //cout<<pp[m][1]<<"     "<<xpp[m][1]<<endl; 
    
           TF1 * func = new TF1("func", f, 0, 10.0, 3);
           func->SetParameter(0,pp[m][0]);
           func->SetParameter(1,pp[m][1]);
           func->SetParameter(2,pp[m][2]);
        
           auto rho = [&](double *x, double *p){
              double b = x[0];
              func->SetParameter(0,b);
              double res = func->Integral(0,10);
              return res; };
    
           TF1 * rho_func = new TF1("rho_func", rho, 0, 5.0, 0);
           TGraph * gr = new TGraph(rho_func->DrawCopy("goff")->GetHistogram());
                    gr->SetLineColor(kGreen);
    
           // find the origional rho to find the difference
           TF1 * func_0 = new TF1("func_0", f, 0, 10.0, 3);
           func_0->SetParameter(0,xpp[m][0]);
           func_0->SetParameter(1,xpp[m][1]);
           func_0->SetParameter(2,xpp[m][2]);
        
           auto rho_0 = [&](double *x, double *p){
              double b = x[0];
              func_0->SetParameter(0,b);
              double res_0 = func_0->Integral(0,100);
              return res_0; };
           
    
           TF1 * rho_func_0 = new TF1("rho_func_0", rho_0, 0, 5.0, 0);
                  gr0 = new TGraph(rho_func_0->DrawCopy("goff")->GetHistogram()); ;
    
           for( int j=0; j< np; j++){
              hh[j]->Fill(rho_func_0->Eval(0.1*j+0.05) - rho_func->Eval(0.1*j+0.05));
              b_perp[j]     =  0.1*j+0.05;
              b_perp_err[j] = 0.0;
              rho_val[j]   = rho_func_0->Eval(0.1*j+0.05);
              }
       }
    
       for( int j=0; j< np; j++){
          hh[j]->Draw();
          hh[j]->SetLineColor(kBlue);
          hh[j]->SetLineWidth(2);
          hh[j]->Fit("gaus","Q","",-0.1, 0.1);
          TF1 *myfit= (TF1*)hh[j]->GetFunction("gaus");
          c3->Print(Form("fig_density/point%d.png", j));
          rho_val_err[j] = myfit->GetParameter(2);
         }
      
       TGraphErrors *gk;
              gk  = new TGraphErrors(np, b_perp, rho_val, b_perp_err, rho_val_err);
              gk->SetMarkerStyle(20);
              if(m == 0) {gk->SetMarkerColor(kBlack);  gk->SetLineColor(kBlack);    leg->AddEntry(gk,"0.1<x_{B}< 0.18 (<x_{B}> = 0.15)","P");} 
              if(m == 1) {gk->SetMarkerColor(kBlue);   gk->SetLineColor(kBlue);     leg->AddEntry(gk,"0.18<x_{B}< 0.24 (<x_{B}> = 0.21)","P");} 
              if(m == 2) {gk->SetMarkerColor(kRed);    gk->SetLineColor(kRed);      leg->AddEntry(gk,"0.24<x_{B}< 0.65 (<x_{B}> = 0.28)","P");}  
              gk->SetLineWidth(2);
    
      //mg->Add(gr0,"l");
      mg->Add(gk,"p");

   }


    // impulse approximation calculations --------------------------------------------------

   auto cal_f = [](double *x, double *p) {
          double t = x[0]/0.197;
          double b = p[0]/0.197;
          double p0 = p[1];
          double p1 = p[2];
          double Delta = TMath::Sqrt(t);
          double ImH = p0*(1-TMath::Power(pow(0.316/0.197,2)*t,6))*TMath::Exp(p1*t);
          //double integrand = ImH*TMath::BesselJ0(b*Delta)*Delta/(2 *(TMath::Pi()));
          double integrand = ImH*TMath::BesselJ0(b*Delta)/(4 *(TMath::Pi()));
          return integrand;
   };

   for (int i =0; i<3; i++)
   {
     double cal_pp[3][3] = {
     1.0,  104.62, -13.34,
     1.0,  48.78, -13.12,
     1.0,  23.99, -12.75
     }; 


     //double cal_pp[3][3] = {
     //1.0,  14.579, -13.6168,
     //1.0,  8.25877, -13.3859,
     //1.0,  4.09494, -13.0487
     //}; 
 

       TF1 * cal_func = new TF1("cal_func", cal_f, 0, 10.0, 3);
       cal_func->SetParameter(0,cal_pp[i][0]);
       cal_func->SetParameter(1,cal_pp[i][1]);
       cal_func->SetParameter(2,cal_pp[i][2]);
    
       auto cal_rho = [&](double *x, double *p){
          double b = x[0];
          cal_func->SetParameter(0,b);
          double res = cal_func->Integral(0,10);
          return res;
       };
    
       TF1 * cal_rho_func = new TF1("cal_rho_func", cal_rho, 0, 5.0, 0);
    
       TGraph * cal_gr = new TGraph(cal_rho_func->DrawCopy("goff")->GetHistogram());
                if (i == 0) { cal_gr->SetLineColor(kBlack);     leg->AddEntry(cal_gr,"IA (x_{B} = 0.15)","L");}      
                if (i == 1) { cal_gr->SetLineColor(kBlue);      leg->AddEntry(cal_gr,"IA (x_{B} = 0.21)","L");}
                if (i == 2) { cal_gr->SetLineColor(kRed);       leg->AddEntry(cal_gr,"IA (x_{B} = 0.28)","L");}

       mg->Add(cal_gr,"l");

   }



   mg->Draw("a");
    leg->Draw("same");
   mg->GetYaxis()->SetTitleOffset(1.4);
   mg->GetYaxis()->CenterTitle(true);
   mg->GetYaxis()->SetTitleSize(0.05);
   mg->GetXaxis()->CenterTitle(true);
   mg->GetXaxis()->SetTitleSize(0.05);
   mg->GetYaxis()->SetTitle("#rho [fm^{-2}])");
   mg->GetXaxis()->SetTitle("b_{#perp}   [fm]");
   gPad->Modified(); 
   mg->GetXaxis()->SetLimits(-0.05, 5.0);
   //mg->GetYaxis()->SetRangeUser(-0.005,0.04);
   //mg->Draw("a");
    c3->SaveAs("fig_density/projected_density_profile.pdf");
    c3->SaveAs("fig_density/projected_density_profile.png");
    c3->SaveAs("fig_density/projected_density_profile.C");
   
}
