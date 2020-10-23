void FunProfile(vector<vector<vector<double>>> mean_t, 
                vector<vector<vector<double>>> mean_t_err,
                vector<double> M_XB,
                vector<double> M_Q2,
                vector<vector<vector<double>>> Im_t_x,
                vector<vector<vector<double>>> Im_t_x_err
){

const int N=3; // number of bins in xB
int npt[N] = {7,7,7}; // number of -t bind in each bin in xB

double X[N] = {M_XB[0], M_XB[1], M_XB[2]};
double Q[N] = {M_Q2[0], M_Q2[1], M_Q2[2]};
double T[N][7];
double Him[N][7]; 
double dHim[N][7]; 
 for(int j=0; j<3; j++){
  for(int k=0; k<7; k++){
   T[j][k] = mean_t[0][j][k];
   Him[j][k] = Im_t_x[0][j][k];
   dHim[j][k] = Im_t_x_err[0][j][k];

   }
 }

const int Npx = 100;


 Float_t Bnorm[N], dBnorm[N];
 Float_t Bslopes[N], dBslopes[N];

 // plot the imaginary part of H -------------------------------
 TCanvas *CANfits = new TCanvas("CANfits","CANfits",3*400,1*400);
          CANfits->Divide(3,1);

 for(int ifit=0;ifit<N;ifit++){
  CANfits->cd(ifit+1);
  TGraphErrors *g_fit_slope = new TGraphErrors(npt[ifit],T[ifit],Him[ifit],0,dHim[ifit]);
                g_fit_slope->SetName(Form("t_dist_%d",ifit+1));
                g_fit_slope->SetTitle(Form("t_dist_%d",ifit+1));
                g_fit_slope->SetMarkerStyle(29);
  TF1 *fitfun = new TF1(Form("fit_exp_%d",ifit+1),"[0]*TMath::Exp(-[1]*x)",0,0.35);
       fitfun->SetParameters(Him[ifit][0],2.5);

  g_fit_slope->Fit(Form("fit_exp_%d",ifit+1));
  g_fit_slope->Draw("AP");

  Bnorm[ifit]    = fitfun->GetParameter(0);
  dBnorm[ifit]   = fitfun->GetParError(0);
  Bslopes[ifit]  = fitfun->GetParameter(1);
  dBslopes[ifit] = fitfun->GetParError(1);
 }
 CANfits->Print("fig/png/clas12_Him_slopes.png");
 CANfits->Print("fig/pdf/clas12_Him_slopes.pdf");


 // plot the profile densities ------------------------------------
 TCanvas *CANprofs = new TCanvas("CANprofs","CANprofs",3*400,1*400);
          CANprofs->Divide(3,1);

 for(int ifit=0;ifit<3;ifit++){
    CANprofs->cd(ifit+1);
    TF1 *fun_prof = new TF1(Form("prof_%d",ifit+1),"[0]*TMath::Exp(-x*x/(4*[1]))",0,0.35);
         fun_prof->SetParameters(Bnorm[ifit],0.197*0.197*Bslopes[ifit]/2.);
         fun_prof->Draw();

    TF1 *fun_profL = new TF1(Form("profL_%d",ifit+1),"[0]*TMath::Exp(-x*x/(4*[1]))",0,0.35);
         fun_profL->SetParameters(Bnorm[ifit]-dBnorm[ifit],0.197*0.197*(Bslopes[ifit]-dBslopes[ifit])/2.);
         fun_profL->Draw("same");

    TF1 *fun_profH = new TF1(Form("profH_%d",ifit+1),"[0]*TMath::Exp(-x*x/(4*[1]))",0,0.35);
         fun_profH->SetParameters(Bnorm[ifit]+dBnorm[ifit],0.197*0.197*(Bslopes[ifit]+dBslopes[ifit])/2.);
         fun_profH->Draw("same");
 }
 CANprofs->Print("fig/png/clas12_profiles.png");
 CANprofs->Print("fig/pdf/clas12_profiles.pdf");

 // plot the 3D profiles ------------------------------------------  
 gStyle->SetOptTitle(0);
 gStyle->SetOptStat(0);
 Float_t obs_x =  2.0;
 Float_t obs_y = -0.25;
 Float_t obs_z =  3.0;

 TCanvas *c0 = new TCanvas("c1","3D",1200,800);
 TH2F *Hemb0 = new TH2F("emb0","emb0",200,-2.0,-0.7,200,-3,-1.3);
       Hemb0->GetXaxis()->SetLabelSize(0);
       Hemb0->GetYaxis()->SetLabelSize(0);
       Hemb0->GetXaxis()->SetTickLength(0);
       Hemb0->GetYaxis()->SetTickLength(0);
       Hemb0->GetXaxis()->SetAxisColor(0);
       Hemb0->GetYaxis()->SetAxisColor(0);
       Hemb0->Draw();

 c0->Update();
 Float_t perspective = 0.95;
 c0->GetFrame()->SetLineColor(0);

 Float_t X0 = (0.-obs_x)*(X[0]-obs_y)/(X[0]-obs_y);
 Float_t Z0 = (0.-obs_z)*(X[0]-obs_y)/(X[0]-obs_y);
 Float_t X1 = (1.-obs_x)*(X[0]-obs_y)/(X[0]-obs_y);
 Float_t Z1 = (0.-obs_z)*(X[0]-obs_y)/(X[0]-obs_y);

 // b axis
 TGaxis *baxis = new TGaxis(X0,Z0,X1,Z1,0,3.999,605);
         baxis->SetLineWidth(3);
         baxis->Draw("same");

 X0 = (0.-obs_x)*(X[0]-obs_y)/(X[0]-obs_y);
 Z0 = (0.-obs_z)*(X[0]-obs_y)/(X[0]-obs_y);
 X1 = (0.-obs_x)*(X[0]-obs_y)/(X[0]-obs_y);
 Z1 = (1.-obs_z)*(X[0]-obs_y)/(X[0]-obs_y);
 TGaxis *zaxis = new TGaxis(X0,Z0,X1,Z1,0,0.20,605);
         zaxis->SetLineWidth(3);
         zaxis->Draw("same");

 // out z axis
 X0 = (1.-obs_x)*(X[0]-obs_y)/(X[0]-obs_y);
 Z0 = (0.-obs_z)*(X[0]-obs_y)/(X[0]-obs_y);
 X1 = (1.-obs_x)*(X[0]-obs_y)/(X[2]-obs_y);
 Z1 = (0.-obs_z)*(X[0]-obs_y)/(X[2]-obs_y);
  TF1 *funaxis = new TF1("funaxis",Form("(0.-%f)*(%f-%f)/(x-%f)",obs_x,X[0],obs_y,obs_y),X[0],X[2]);
  TGaxis *xaxis = new TGaxis(X1,Z1,X0,Z0,"funaxis",605);
          xaxis->SetLabelOffset(0.02);
          xaxis->SetLineWidth(3);
          xaxis->Draw("same");

 //middle z axis
 TLine *ob_ax;
        X0 = (0.-obs_x)*(X[0]-obs_y)/(X[0]-obs_y);
        Z0 = (0.-obs_z)*(X[0]-obs_y)/(X[0]-obs_y);
        X1 = (0.-obs_x)*(X[0]-obs_y)/(X[2]-obs_y);
        Z1 = (0.-obs_z)*(X[0]-obs_y)/(X[2]-obs_y);
        ob_ax = new TLine(X0,Z0,X1,Z1);
        ob_ax->SetLineWidth(3); 
        ob_ax->Draw();

 // last x-axis 
 TLine *vert_ax;
        X0 = (0.-obs_x)*(X[0]-obs_y)/(X[2]-obs_y);
        Z0 = (0.-obs_z)*(X[0]-obs_y)/(X[2]-obs_y);
        X1 = (0.-obs_x)*(X[0]-obs_y)/(X[2]-obs_y);
        Z1 = (1.-obs_z)*(X[0]-obs_y)/(X[2]-obs_y);
        vert_ax = new TLine(X0,Z0,X1,Z1);
        //vert_ax->Draw();
        
 TLine *hor_ax;
 // bt vs. xB division lines
                for(int idiv=1;idiv<5;idiv++)
                 {
                    X0 = (0.2*idiv-obs_x)*(X[0]-obs_y)/(X[0]-obs_y);
                    Z0 = (0.-obs_z)*(X[0]-obs_y)/(X[0]-obs_y);
                    X1 = (0.2*idiv-obs_x)*(X[0]-obs_y)/(X[2]-obs_y);
                    Z1 = (0.-obs_z)*(X[0]-obs_y)/(X[2]-obs_y);
                    vert_ax = new TLine(X0,Z0,X1,Z1);
                    //vert_ax->Draw();
                  }

  // xB division lines 
 for(int idiv=1;idiv<3;idiv++)
  {
     X0 = (0.-obs_x)*(X[0]-obs_y)/(X[idiv]-obs_y);
     Z0 = (0.-obs_z)*(X[0]-obs_y)/(X[idiv]-obs_y);
     X1 = (0.-obs_x)*(X[0]-obs_y)/(X[idiv]-obs_y);
     Z1 = (1.-obs_z)*(X[0]-obs_y)/(X[idiv]-obs_y);
     vert_ax = new TLine(X0,Z0,X1,Z1);
     vert_ax->Draw();

     zaxis = new TGaxis(X0,Z0,X1,Z1,0,0.20,605);
     zaxis->SetLineWidth(3);
     zaxis->Draw("same");
  }
 
  // rho division lines 
               for(int idiv=1;idiv<6;idiv++)
               {
                  X0 = (0.-obs_x)*(X[0]-obs_y)/(X[0]-obs_y);
                  Z0 = (0.2*idiv-obs_z)*(X[0]-obs_y)/(X[0]-obs_y);
                  X1 = (0.-obs_x)*(X[0]-obs_y)/(X[2]-obs_y);
                  Z1 = (0.2*idiv-obs_z)*(X[0]-obs_y)/(X[2]-obs_y);
                  ob_ax = new TLine(X0,Z0,X1,Z1);
                  //ob_ax->Draw();
               }

 for(int idiv=1;idiv<3;idiv++)
  {
     X0 = (0.-obs_x)*(X[0]-obs_y)/(X[idiv]-obs_y);
     Z0 = (0.-obs_z)*(X[0]-obs_y)/(X[idiv]-obs_y);
     X1 = (1.-obs_x)*(X[0]-obs_y)/(X[idiv]-obs_y);
     Z1 = (0.-obs_z)*(X[0]-obs_y)/(X[idiv]-obs_y);
     hor_ax = new TLine(X0,Z0,X1,Z1);
     // hor_ax->Draw();
     baxis = new TGaxis(X0,Z0,X1,Z1,0,3.999,605);
     baxis->SetLineWidth(3);
     baxis->Draw("same");
  }


 TF1 *profile_fun2d[7], *Hprofile_fun2d[7];
 TGraphErrors *profile_G2D[7];
 Float_t gb[Npx], gr[Npx], gdr[Npx];
 Float_t midX[Npx], midY[Npx];

 for(int ifit=0;ifit<3;ifit++)
   {
     Float_t persp_fact;
     //if(ifit ==0 || ifit==2)  persp_fact = (X[ifit]-obs_y)/(X[2]-obs_y);
     //else if(ifit==1)  persp_fact = (X[0]-obs_y)/(X[ifit]-obs_y);
      persp_fact = (X[0]-obs_y)/(X[0]-obs_y);
     
     X0 = (0.-obs_x)*persp_fact;
     Z0 = (0.-obs_z)*persp_fact;
     X1 = (1.0-obs_x)*persp_fact;
     Z1 = (1.0-obs_z)*persp_fact;

     Float_t xi = X[ifit];
     profile_fun2d[ifit] = new TF1(Form("prof_fun%d",ifit+1),"[2]+[0]*TMath::Exp(-TMath::Power((x-[3])/[4],2.)/(4*[1]))/90",X0,X1);
     profile_fun2d[ifit]->SetParameters(persp_fact*Bnorm[ifit],0.197*0.197*Bslopes[ifit]/(2.*(15.0-xi)),
                                        Z0,X0,persp_fact);
     Hprofile_fun2d[ifit] = new TF1(Form("Hprof_fun%d",ifit+1),"[2]+[0]*TMath::Exp(-TMath::Power((x-[3])/[4],2.)/(4*[1]))/90",X0,X1);
     Hprofile_fun2d[ifit]->SetParameters(persp_fact*(Bnorm[ifit]-dBnorm[ifit]),0.197*0.197*(Bslopes[ifit]-dBslopes[ifit])/(2.*(15.-xi)),Z0,X0,persp_fact);
     midX[ifit] = X0 + persp_fact* 2.*TMath::Sqrt(0.197*0.197*(Bslopes[ifit]-dBslopes[ifit])/(2.*(15.0-xi)));
     midY[ifit] = profile_fun2d[ifit]->Eval(midX[ifit]);
  
   for(int point=0;point<Npx;point++)
    {
      gb[point] = X0 + (point+0.5)*(X1-X0)/Npx;
      gr[point] = profile_fun2d[ifit]->Eval(gb[point]);
      gdr[point] = TMath::Abs(Hprofile_fun2d[ifit]->Eval(gb[point])-profile_fun2d[ifit]->Eval(gb[point]));
     }

    profile_G2D[ifit] = new TGraphErrors(Npx,gb,gr,0,gdr);
    profile_G2D[ifit]->SetMarkerStyle(20);
    profile_G2D[ifit]->SetMarkerSize(0.3);
    profile_G2D[ifit]->SetMarkerColor(4);
    profile_G2D[ifit]->SetFillColor(40); profile_G2D[ifit]->SetLineColor(40);
    
    profile_G2D[ifit]->SetLineWidth(2);
    profile_G2D[ifit]->Draw("LPsameE3");
    profile_fun2d[ifit]->SetLineColor(2);
    profile_fun2d[ifit]->Draw("same");
  
 }

 TLatex *comments = new TLatex();
         comments->SetTextColor(1);
         comments->SetTextSize(0.06);
         comments->DrawLatex(-0.7,-2.8,"x_{B}");
         comments->DrawLatex(-1.5,-3.22,"b_{#perp}  (fm)");
         comments->DrawLatex(-2.27,-2.42,"#rho(fm^{-2})");
         comments->SetTextAngle(28);
         comments->SetTextSize(0.055);
         comments->DrawLatex(-2.15,-1.96,"^{4}He quark density profiles");

   TLegend* leg = new TLegend(0.58,0.78,0.99,0.95);
   leg-> SetNColumns(1);
   leg->SetTextSize(0.04);
   leg->AddEntry(profile_fun2d[0]," Impulse Appr. Model","L");
   leg->AddEntry(profile_G2D[0]," Projected uncertainties","f");
   leg->Draw("same");


c0->Print("fig/nucl_profile.png");
c0->Print("fig/nucl_profile.pdf");
c0->Print("fig/nucl_profile2.png");
c0->Print("fig/nucl_profile.ps");


 //gSystem->Exit(0);
}
