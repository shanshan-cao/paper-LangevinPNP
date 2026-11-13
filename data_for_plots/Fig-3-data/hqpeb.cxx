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
void hqpeb()
{
   TCanvas *c1=new TCanvas("hhraa","a Canvas",10,10,600,500);
   c1->Divide(1,1);

   TH1D *h1=new TH1D("h1","",100,0.0,10.0);

   h1->SetLineColor(0);
   h1->SetLineWidth(1);
   h1->SetLineStyle(1);

   ifstream fort;
   Double_t x0[100],y0[100];
   Double_t x1[87],y1[87];
   Double_t x2[87],y2[87];
   Double_t x3[87],y3[87];
   Double_t x4[87],y4[87];
   Double_t x5[87],y5[87];
   Double_t x6[87],y6[87];
   Double_t x7[87],y7[87];
   Double_t x8[87],y8[87];
   Double_t x9[87],y9[87];
   Double_t xe[87],ye[87];
   Double_t ce,te,eaverE;

   Int_t d=0;
   Int_t nbin;
   fort.open("referline1.dat");
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
   fort.open("hqpet2fm.dat");
   if(fort.fail()){
      cout<<"sorry,couldn't open file"<<endl;
      exit(1);
}
   while(!fort.eof())
   {
      fort>>x1[d]>>y1[d];
      d++;
}
   fort.close();

   d=0;
   fort.open("hqpet4fm.dat");
   if(fort.fail()){
      cout<<"sorry,couldn't open file"<<endl;
      exit(1);
}
   while(!fort.eof())
   {
      fort>>x2[d]>>y2[d];
      d++;
}
   fort.close();

   d=0;
   fort.open("hqpet6fm.dat");
   if(fort.fail()){
      cout<<"sorry,couldn't open file"<<endl;
      exit(1);
}
   while(!fort.eof())
   {
      fort>>x3[d]>>y3[d];
      d++;
}
   fort.close();

   d=0;
   fort.open("hqpet8fm.dat");
   if(fort.fail()){
      cout<<"sorry,couldn't open file"<<endl;
      exit(1);
}
   while(!fort.eof())
   {
      fort>>x4[d]>>y4[d];
      d++;
}
   fort.close();

   d=0;
   fort.open("hqpet10fm.dat");
   if(fort.fail()){
      cout<<"sorry,couldn't open file"<<endl;
      exit(1);
}
   while(!fort.eof())
   {
      fort>>x5[d]>>y5[d];
      d++;
}
   fort.close();

   d=0;
   fort.open("hqpet15fm.dat");
   if(fort.fail()){
      cout<<"sorry,couldn't open file"<<endl;
      exit(1);
}
   while(!fort.eof())
   {
      fort>>x6[d]>>y6[d];
      d++;
}
   fort.close();


   d=0;
   fort.open("hqpet20fm.dat");
   if(fort.fail()){
      cout<<"sorry,couldn't open file"<<endl;
      exit(1);
}
   while(!fort.eof())
   {
      fort>>x7[d]>>y7[d];
      d++;
}
   fort.close();


   d=0;
   fort.open("hqpet30fm.dat");
   if(fort.fail()){
      cout<<"sorry,couldn't open file"<<endl;
      exit(1);
}
   while(!fort.eof())
   {
      fort>>x8[d]>>y8[d];
      d++;
}
   fort.close();
   
   d=0;
   fort.open("hqpet40fm.dat");
   if(fort.fail()){
      cout<<"sorry,couldn't open file"<<endl;
      exit(1);
}
   while(!fort.eof())
   {
      fort>>x9[d]>>y9[d];
      d++;
}
   fort.close();
 
   for(Int_t i=0;i<100;i++)
   {
      nbin=(x0[i]+0.051)/0.1;
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
     gPad->SetLogy(1);
     h1->GetYaxis()->SetRangeUser(1.0e-6,2.5e1);
     h1->GetXaxis()->SetRangeUser(1.0,8.0);
     h1->GetYaxis()->SetTitle("dN/pEdE");
     h1->GetXaxis()->SetTitle("E (GeV) ");
     h1->GetXaxis()->SetLabelSize(0.05);
     h1->GetYaxis()->SetLabelSize(0.05);
     h1->GetYaxis()->SetNdivisions(010);
     h1->GetXaxis()->SetNdivisions(010);
     h1->GetYaxis()->CenterTitle();
     h1->GetYaxis()->SetTitleOffset(1.1);
     h1->GetYaxis()->SetTitleSize(0.06);
     h1->GetXaxis()->CenterTitle();
     h1->GetXaxis()->SetTitleOffset(0.85);
     h1->GetXaxis()->SetTitleSize(0.06);

     TGraph *gr1=new TGraph(87,x1,y1);
     gr1->SetLineColor(8);
     gr1->SetMarkerStyle(20);
     gr1->SetMarkerSize(1.0);
     gr1->SetLineWidth(5);
     gr1->SetLineStyle(9);
    // gr1->Draw("l");
     gr1->SetMarkerColor(8);
    
     TGraph *gr2=new TGraph(87,x2,y2);
     gr2->SetLineColor(4);
     gr2->SetMarkerStyle(21);
     gr2->SetMarkerSize(1.0);
     gr2->SetLineWidth(5);
     gr2->SetLineStyle(6);
    // gr2->Draw("l");
     gr2->SetMarkerColor(4);

     TGraph *gr3=new TGraph(87,x3,y3);
     gr3->SetLineColor(7);
     gr3->SetMarkerStyle(22);
     gr3->SetMarkerSize(1.0);
     gr3->SetLineWidth(5);
     gr3->SetLineStyle(2);
    // gr3->Draw("l");
     gr3->SetMarkerColor(7);

     TGraph *gr4=new TGraph(87,x4,y4);
     gr4->SetLineColor(6);
     gr4->SetMarkerStyle(23);
     gr4->SetMarkerSize(1.0);
     gr4->SetLineWidth(5);
     gr4->SetLineStyle(5);
     //gr4->Draw("l");
     gr4->SetMarkerColor(6);

     TGraph *gr5=new TGraph(87,x5,y5);
     gr5->SetLineColor(8);
     gr5->SetMarkerStyle(24);
     gr5->SetMarkerSize(1.0);
     gr5->SetLineWidth(5);
     gr5->SetLineStyle(9);
     gr5->Draw("l");
     gr5->SetMarkerColor(8);

     TGraph *gr6=new TGraph(87,x6,y6);
     gr6->SetLineColor(4);
     gr6->SetMarkerStyle(33);
     gr6->SetMarkerSize(1.5);
     gr6->SetLineWidth(5);
     gr6->SetLineStyle(6);
     gr6->Draw("l");
     gr6->SetMarkerColor(7);

     TGraph *gr7=new TGraph(87,x7,y7);
     gr7->SetLineColor(7);
     gr7->SetMarkerStyle(23);
     gr7->SetMarkerSize(1.0);
     gr7->SetLineWidth(5);
     gr7->SetLineStyle(2);
     gr7->Draw("l");
     gr7->SetMarkerColor(7);

     TGraph *gr8=new TGraph(87,x8,y8);
     gr8->SetLineColor(6);
     gr8->SetMarkerStyle(35);
     gr8->SetMarkerSize(1.0);
     gr8->SetLineWidth(5);
     gr8->SetLineStyle(5);
     gr8->Draw("l");
     gr8->SetMarkerColor(6);
     
     TGraph *gr9=new TGraph(87,x9,y9);
     gr9->SetLineColor(2);
     gr9->SetMarkerStyle(35);
     gr9->SetMarkerSize(1.0);
     gr9->SetLineWidth(5);
     gr9->SetLineStyle(1);
     gr9->Draw("l");
     gr9->SetMarkerColor(2);
     
     
     TLegend *l = new TLegend(0.48,0.61,0.87,0.83);   
     l->SetTextSize(0.05);
     l->AddEntry(gr1,"Yukawa + string"," ");
     l->SetFillColor(0);
     l->SetLineColor(0);
     l->Draw();   
     
     

     TLegend *l1 = new TLegend(0.48,0.61,0.87,0.83);   
     l1->SetTextSize(0.05);
     l1->AddEntry(gr5, "t = 10 fm/c","l");
     l1->AddEntry(gr6, "t = 15 fm/c","l");
     l1->AddEntry(gr7, "t = 20 fm/c","l");
     l1->AddEntry(gr8, "t = 30 fm/c","l");
     l1->AddEntry(gr9, "t = 40 fm/c","l");
     l1->SetFillColor(0);
     l1->SetLineColor(0);
     l1->Draw();




}
