/***********************************************************************
*root -l 
*.x read_draw.cxx
*.q
*
*
************************************************************************/
#include<iostream>
#include<fstream>
#include<cstdlib>
#include"TH1F.h"
#include"TCanvas.h"
#include"TMath.h"
using namespace std;
void hqrh()
{
   TCanvas *c1=new TCanvas("hhraa","a Canvas",10,10,600,500);
   c1->Divide(1,1);

   TH1D *h1=new TH1D("h1","",10,0.0,100.0);

   h1->SetLineColor(1);
   h1->SetLineWidth(2);
   h1->SetLineStyle(2);

   ifstream fort;
   Double_t x0[10],y0[10];
   Double_t x1[9],y11[9],y12[9];
   Double_t x2[9],y21[9],y22[9];
   Double_t x3[9],y31[9],y32[9];
   Double_t x4[9],y41[9],y42[9];
   Double_t x5[100],y5[100];
   Double_t x6[100],y6[100];
   Double_t x7[100],y7[100];
   Double_t x8[100],y8[100];
   Double_t x9[100],y9[100];
   Double_t xe[100],ye[100];
   Double_t ce,te,eaverE;

   Int_t d=0;
   Int_t nbin;
   fort.open("referline3.dat");
   if(fort.fail()){
      cout<<"sorry,couldn't open file"<<endl;
      exit(1);
}
   while(!fort.eof())
   {
      fort>>x0[d]>>y0[d];
      d++;
}
   fort.close();

   d=0;
   fort.open("energy-aver-y.dat");
   if(fort.fail()){
      cout<<"sorry,couldn't open file"<<endl;
      exit(1);
}
   while(!fort.eof())
   {
      fort>>x1[d]>>y11[d]>>y12[d];
      d++;
}
   fort.close();

    d=0;
   fort.open("energy-aver-s.dat");
   if(fort.fail()){
      cout<<"sorry,couldn't open file"<<endl;
      exit(1);
}
   while(!fort.eof())
   {
      fort>>x2[d]>>y21[d]>>y22[d];
      d++;
}
   fort.close();  
   
   d=0;
   fort.open("energy-aver-sy.dat");
   if(fort.fail()){
      cout<<"sorry,couldn't open file"<<endl;
      exit(1);
}
   while(!fort.eof())
   {
      fort>>x3[d]>>y31[d]>>y32[d];
      d++;
}
   fort.close();
   
   d=0;
   fort.open("energy-aver-d4.dat");
   if(fort.fail()){
      cout<<"sorry,couldn't open file"<<endl;
      exit(1);
}
   while(!fort.eof())
   {
      fort>>x4[d]>>y41[d]>>y42[d];
      d++;
}
   fort.close();
   
   
   
   
   
   
 
   for(Int_t i=0;i<10;i++)
   {
      nbin=(x0[i]+5.1)/10.0;
      h1->SetBinContent(nbin,y0[i]);

}


//   ce=1.0;
//   te=0.5; 
//   for(Int_t i=0;i<100;i++)
//   {
//     xe[i]=x8[i];
//     aver=aver+
//     ye[i]=ce*exp(-xe[i]/te);
//}




     c1->cd(1);
     h1->Draw();
     //gPad->SetLogx(1);
    // gPad->SetLogy(1);
     h1->GetYaxis()->SetRangeUser(0.0,1500);
     h1->GetXaxis()->SetRangeUser(1.0,50.0);
     h1->GetYaxis()->SetTitle("T (MeV)");
     h1->GetXaxis()->SetTitle("t (fm/c)");
     h1->GetXaxis()->SetLabelSize(0.05);
     h1->GetYaxis()->SetLabelSize(0.05);
     h1->GetYaxis()->SetNdivisions(010);
     h1->GetXaxis()->SetNdivisions(010);
     h1->GetYaxis()->CenterTitle();
     h1->GetYaxis()->SetTitleOffset(1.1);
     h1->GetYaxis()->SetTitleSize(0.06);
     h1->GetXaxis()->CenterTitle();
     h1->GetXaxis()->SetTitleOffset(0.85);
     h1->GetXaxis()->SetTitleSize(0.07);

     TGraph *gr1=new TGraph(9,x1,y12);
     gr1->SetLineColor(8);
     gr1->SetMarkerStyle(20);
     gr1->SetMarkerSize(1.5);
     gr1->SetLineWidth(5);
     gr1->SetLineStyle(9);
     gr1->Draw("l");
     gr1->SetMarkerColor(8);
    
     TGraph *gr2=new TGraph(9,x2,y22);
     gr2->SetLineColor(4);
     gr2->SetMarkerStyle(20);
     gr2->SetMarkerSize(1.5);
     gr2->SetLineWidth(5);
     gr2->SetLineStyle(6);
     gr2->Draw("l");
     gr2->SetMarkerColor(4);

     TGraph *gr3=new TGraph(9,x3,y32);
     gr3->SetLineColor(2);
     gr3->SetMarkerStyle(20);
     gr3->SetMarkerSize(1.5);
     gr3->SetLineWidth(5);
     gr3->SetLineStyle(1);
     gr3->Draw("l");
     gr3->SetMarkerColor(2);

     TGraph *gr4=new TGraph(9,x4,y42);
     gr4->SetLineColor(7);
     gr4->SetMarkerStyle(20);
     gr4->SetMarkerSize(1.5);
     gr4->SetLineWidth(5);
     gr4->SetLineStyle(5);
     gr4->Draw("l");
     gr4->SetMarkerColor(7);




 
     TLegend *l1 = new TLegend(0.48,0.61,0.87,0.83);
     l1->SetTextSize(0.05);   
     //l1->AddEntry(gr1, "From the energy spectrum",""); 
     l1->AddEntry(gr1,"Yukawa","l");
     l1->AddEntry(gr2,"string","l");
     l1->AddEntry(gr3,"Yukawa + string","l"); 
     l1->AddEntry(gr4,"D_{s}(2#piT) = 4","l");      
     l1->SetFillColor(0);
     l1->SetLineColor(0);
     l1->Draw();




}
