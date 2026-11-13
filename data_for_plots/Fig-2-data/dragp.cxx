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
void dragp()
{
   TCanvas *c1=new TCanvas("hhraa","a Canvas",10,10,600,500);
   c1->Divide(1,1);

   TH1D *h1=new TH1D("h1","",10,0.0,50.0);
   TH1D *h2=new TH1D("h2","",10,0.0,50.0);

   h1->SetLineColor(1);
   h1->SetLineWidth(1);
   h1->SetLineStyle(1);

   h2->SetLineColor(0);
   h2->SetLineWidth(2);
   h2->SetLineStyle(2);


   ifstream fort;

   Double_t x[100],y1[100],y2[100];
   Double_t x2[30],y21[30];
   Double_t x3[30],y31[30];
   Double_t x4[30],y41[30];
   Double_t x5[30],y51[30];
   Double_t x6[30],y61[30];
   Double_t x7[30],y71[30];
   Double_t x8[30],y81[30];


   Double_t tem1[1500],pe1[1500],drag1[1500];
   Double_t tem2[1500],pe2[1500],drag2[1500];
   Double_t tem3[1500],pe3[1500],drag3[1500];

   Double_t tem1b[50],pe1b[50],drag1b[50];
   Double_t tem2b[50],pe2b[50],drag2b[50];
   Double_t tem3b[50],pe3b[50],drag3b[50];
   Double_t tem1c[50],pe1c[50],drag1c[50];
   Double_t tem2c[50],pe2c[50],drag2c[50];
   Double_t tem3c[50],pe3c[50],drag3c[50];

   Int_t d=0;
   Int_t nbin;
   fort.open("referline4.dat");
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
   fort.open("drag-sy-c.dat");
   if(fort.fail()){
      cout<<"sorry,couldn't open file"<<endl;
      exit(1);
}
   while(!fort.eof())
   {
      fort>>tem1[d]>>pe1[d]>>drag1[d];
      d++;
}
   fort.close();

   d=0;
   fort.open("drag-sy-c.dat");
   if(fort.fail()){
      cout<<"sorry,couldn't open file"<<endl;
      exit(1);
}
   while(!fort.eof())
   {
      fort>>tem2[d]>>pe2[d]>>drag2[d];
      d++;
}
   fort.close();



   for(Int_t i=0;i<100;i++)
   {
      nbin=(x[i]+0.25)/0.5;
      h1->SetBinContent(nbin,y1[i]);
//      h2->SetBinContent(nbin,y2[i]);

}

 for(Int_t i=200;i<250;i++)
   {
      
      pe1b[i-200]=pe1[i];
      drag1b[i-200]=drag1[i];
      pe1c[i-200]=pe2[i];
      drag1c[i-200]=drag2[i];
//      cout<<temsy[i]<<endl;
}

   for(Int_t i=450;i<500;i++)
   {

      pe2b[i-450]=pe1[i];
      drag2b[i-450]=drag1[i];
      pe2c[i-450]=pe2[i];
      drag2c[i-450]=drag2[i];
//      cout<<temsy[i]<<endl;
}

   for(Int_t i=1200;i<1250;i++)
   {
      
      pe3b[i-1200]=pe1[i];
      drag3b[i-1200]=drag1[i];
      pe3c[i-1200]=pe2[i];
      drag3c[i-1200]=drag2[i];
//      cout<<temsy[i]<<endl;
}




     c1->cd(1);
     h1->Draw();
     //h2->Draw("same");
     //gPad->SetLogx(1);
     h1->GetYaxis()->SetRangeUser(0.0,0.14);
     h1->GetXaxis()->SetRangeUser(0.0,50);
     h1->GetYaxis()->SetTitle("#eta_{D} (GeV)");
     h1->GetXaxis()->SetTitle("p (GeV)");
     h1->GetXaxis()->SetLabelSize(0.05);
     h1->GetYaxis()->SetLabelSize(0.05);
     h1->GetYaxis()->SetNdivisions(010);
     h1->GetXaxis()->SetNdivisions(010);
     h1->GetYaxis()->CenterTitle();
     h1->GetYaxis()->SetTitleOffset(1.2);
     h1->GetYaxis()->SetTitleSize(0.06);
     h1->GetXaxis()->CenterTitle();
     h1->GetXaxis()->SetTitleOffset(0.9);
     h1->GetXaxis()->SetTitleSize(0.07);


    
     TGraph *gr11=new TGraph(50,pe1b,drag1b);
     gr11->SetLineColor(4);
     gr11->SetMarkerStyle(20);
     gr11->SetMarkerSize(0.7);
     gr11->SetLineWidth(4);
     gr11->SetLineStyle(6);
    // gr11->Draw("l");
     gr11->SetMarkerColor(4);

     TGraph *gr12=new TGraph(50,pe2b,drag2b);
     gr12->SetLineColor(2);
     gr12->SetMarkerStyle(21);
     gr12->SetMarkerSize(0.7);
     gr12->SetLineWidth(4);
     gr12->SetLineStyle(1);
   //  gr12->Draw("l");
     gr12->SetMarkerColor(4);

    
     TGraph *gr14=new TGraph(50,pe3b,drag3b);
     gr14->SetLineColor(8);
     gr14->SetMarkerStyle(20);
     gr14->SetMarkerSize(0.7);
     gr14->SetLineWidth(4);
     gr14->SetLineStyle(9);
   //  gr14->Draw("l");
     gr14->SetMarkerColor(1);

    
     TGraph *gr15=new TGraph(50,pe1c,drag1c);
     gr15->SetLineColor(8);
     gr15->SetMarkerStyle(20);
     gr15->SetMarkerSize(0.7);
     gr15->SetLineWidth(5);
     gr15->SetLineStyle(9);
     gr15->Draw("l");
     gr15->SetMarkerColor(4);

     TGraph *gr16=new TGraph(50,pe2c,drag2c);
     gr16->SetLineColor(4);
     gr16->SetMarkerStyle(21);
     gr16->SetMarkerSize(0.7);
     gr16->SetLineWidth(5);
     gr16->SetLineStyle(6);
     gr16->Draw("l");
     gr16->SetMarkerColor(4);

    
     TGraph *gr17=new TGraph(50,pe3c,drag3c);
     gr17->SetLineColor(2);
     gr17->SetMarkerStyle(20);
     gr17->SetMarkerSize(0.7);
     gr17->SetLineWidth(5);
     gr17->SetLineStyle(1);
     gr17->Draw("l");
     gr17->SetMarkerColor(2);






     TLegend *l1 = new TLegend(0.15,0.75,0.85,0.85);
     l1->SetTextSize(0.05); 
     //l1->AddEntry(gr11, "Bottom"," ");
     //l1->AddEntry(gr11, "T = 200 MeV","l");
    // l1->AddEntry(gr12, "T = 300 MeV","l");
    // l1->AddEntry(gr14, "T = 600 MeV","l");
     l1->AddEntry(gr11, "charm"," ");
     l1->AddEntry(gr15, "T = 200 MeV","l");
     l1->AddEntry(gr16, "T = 300 MeV","l");
     l1->AddEntry(gr17, "T = 600 MeV","l");
     l1->SetFillColor(0);
     l1->SetLineColor(0);
     l1->Draw();




}
