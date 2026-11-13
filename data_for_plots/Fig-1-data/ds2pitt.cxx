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
void ds2pitt()
{
   TCanvas *c1=new TCanvas("hhraa","a Canvas",10,10,600,500);
   c1->Divide(1,1);

   TH1D *h1=new TH1D("h1","",10,0.0,5.0);
   TH1D *h2=new TH1D("h2","",10,0.0,5.0);

   h1->SetLineColor(0);
   h1->SetLineWidth(1);
   h1->SetLineStyle(1);

   h2->SetLineColor(0);
   h2->SetLineWidth(2);
   h2->SetLineStyle(2);


   ifstream fort;

   Double_t x[10],y1[10],y2[10];
   Double_t x2[30],y21[30];
   Double_t x3[30],y31[30];
   Double_t x4[30],y41[30];
   Double_t x5[30],y51[30];
   Double_t x6[30],y61[30];
   Double_t x7[30],y71[30];
   Double_t x8[30],y81[30];


   Double_t temsy[1500],pesy[1500],kpsy[1500],kpdsy[1500];
   Double_t tems[1500],pes[1500],kps[1500],kpds[1500];
   Double_t temy[1500],pey[1500],kpy[1500],kpdy[1500];

   Double_t dtpsy[1500],dtps[1500],dtpy[1500];

   Double_t tem1[30],pe1[30],dtp1[30];
   Double_t tem2[30],pe2[30],dtp2[30];
   Double_t tem3[30],pe3[30],dtp3[30];

   Double_t tem1a[50],pe1a[50],dtp1a[50];
   Double_t tem2a[50],pe2a[50],dtp2a[50];
   Double_t tem3a[30],pe3a[30],dtp3a[50];




   Int_t d=0;
   Int_t nbin;
   fort.open("referline3.dat");
   if(fort.fail()){
      cout<<"sorry,couldn't open file"<<endl;
      exit(1);
}
   while(!fort.eof())
   {
      fort>>x[d]>>y1[d]>>y2[d];
      d++;
}
   fort.close();

   d=0;
   fort.open("Ds_tp-sy-c-p5.dat");
   if(fort.fail()){
      cout<<"sorry,couldn't open file"<<endl;
      exit(1);
}
   while(!fort.eof())
   {
      fort>>x2[d]>>y21[d];
      d++;
}
   fort.close();

   d=0;
   fort.open("Ds_tp-sy-c-p10.dat");
   if(fort.fail()){
      cout<<"sorry,couldn't open file"<<endl;
      exit(1);
}
   while(!fort.eof())
   {
      fort>>x3[d]>>y31[d];
      d++;
}
   fort.close();

   d=0;
   fort.open("Ds_tp-sy-c-p50.dat");
   if(fort.fail()){
      cout<<"sorry,couldn't open file"<<endl;
      exit(1);
}
   while(!fort.eof())
   {
      fort>>x4[d]>>y41[d];
      d++;
}
   fort.close();

   d=0;
   fort.open("Ds_tp-sy-c-p5.dat");
   if(fort.fail()){
      cout<<"sorry,couldn't open file"<<endl;
      exit(1);
}
   while(!fort.eof())
   {
      fort>>x5[d]>>y51[d];
      d++;
}
   fort.close();

   d=0;
   fort.open("Ds_tp-sy-c-p10.dat");
   if(fort.fail()){
      cout<<"sorry,couldn't open file"<<endl;
      exit(1);
}
   while(!fort.eof())
   {
      fort>>x6[d]>>y61[d];
      d++;
}
   fort.close();
  
   d=0;
   fort.open("Ds_tp-sy-c-p50.dat");
   if(fort.fail()){
      cout<<"sorry,couldn't open file"<<endl;
      exit(1);
}
   while(!fort.eof())
   {
      fort>>x7[d]>>y71[d];
      d++;
}
   fort.close();




   for(Int_t i=0;i<10;i++)
   {
      nbin=(x[i]+0.25)/0.5;
      h1->SetBinContent(nbin,y1[i]);
//      h2->SetBinContent(nbin,y2[i]);

}



     c1->cd(1);
     h1->Draw();
     //h2->Draw("same");
     //gPad->SetLogx(1);
     h1->GetYaxis()->SetRangeUser(0.0,25.0);
     h1->GetXaxis()->SetRangeUser(0.5,4.5);
     h1->GetYaxis()->SetTitle("D_{s}(2#piT)");
     h1->GetXaxis()->SetTitle("T/T_{c}");
     h1->GetXaxis()->SetLabelSize(0.05);
     h1->GetYaxis()->SetLabelSize(0.05);
     h1->GetYaxis()->SetNdivisions(010);
     h1->GetXaxis()->SetNdivisions(010);
     h1->GetYaxis()->CenterTitle();
     h1->GetYaxis()->SetTitleOffset(0.8);
     h1->GetYaxis()->SetTitleSize(0.06);
     h1->GetXaxis()->CenterTitle();
     h1->GetXaxis()->SetTitleOffset(0.9);
     h1->GetXaxis()->SetTitleSize(0.07);


    
     TGraph *gr11=new TGraph(30,x2,y21);
     gr11->SetLineColor(4);
     gr11->SetMarkerStyle(20);
     gr11->SetMarkerSize(0.7);
     gr11->SetLineWidth(4);
     gr11->SetLineStyle(6);
     //gr11->Draw("l");
     gr11->SetMarkerColor(4);

     TGraph *gr12=new TGraph(30,x3,y31);
     gr12->SetLineColor(2);
     gr12->SetMarkerStyle(21);
     gr12->SetMarkerSize(0.7);
     gr12->SetLineWidth(4);
     gr12->SetLineStyle(1);
     //gr12->Draw("l");
     gr12->SetMarkerColor(4);

    
     TGraph *gr14=new TGraph(30,x4,y41);
     gr14->SetLineColor(8);
     gr14->SetMarkerStyle(20);
     gr14->SetMarkerSize(0.7);
     gr14->SetLineWidth(4);
     gr14->SetLineStyle(9);
     //gr14->Draw("l");
     gr14->SetMarkerColor(1);

     TGraph *gr15=new TGraph(30,x5,y51);
     gr15->SetLineColor(8);
     gr15->SetMarkerStyle(21);
     gr15->SetMarkerSize(0.7);
     gr15->SetLineWidth(5);
     gr15->SetLineStyle(9);
     gr15->Draw("l");
     gr15->SetMarkerColor(2);

     TGraph *gr16=new TGraph(30,x6,y61);
     gr16->SetLineColor(4);
     gr16->SetMarkerStyle(22);
     gr16->SetMarkerSize(0.7);
     gr16->SetLineWidth(5);
     gr16->SetLineStyle(6);
     gr16->Draw("l");
     gr16->SetMarkerColor(3);

     TGraph *gr17=new TGraph(30,x7,y71);
     gr17->SetLineColor(2);
     gr17->SetMarkerStyle(22);
     gr17->SetMarkerSize(0.7);
     gr17->SetLineWidth(5);
     gr17->SetLineStyle(1);
     gr17->Draw("l");
     gr17->SetMarkerColor(3);


     TLegend *l1 = new TLegend(0.15,0.75,0.85,0.85);
     l1->SetTextSize(0.05);      
    // l1->AddEntry(gr11, "Bottom"," ");
    // l1->AddEntry(gr11, "p = 5 GeV","l");
    // l1->AddEntry(gr12, "p = 10 GeV","l");
    // l1->AddEntry(gr14, "p = 50 GeV","l");
     l1->AddEntry(gr11, "charm"," ");
     l1->AddEntry(gr15, "p = 5 GeV","l");
     l1->AddEntry(gr16, "p = 10 GeV","l");
     l1->AddEntry(gr17, "p = 50 GeV","l");
     l1->SetFillColor(0);
     l1->SetLineColor(0);
     l1->Draw();




}
