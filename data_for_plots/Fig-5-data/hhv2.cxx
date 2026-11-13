/***********************************************************************
 * root -l
 * .x hhv2.cxx
 * .q
 ***********************************************************************/
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "TH1F.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TMath.h"
#include "TStyle.h"
using namespace std;

void hhv2()
{
   // ================= 画布 =================
   TCanvas *c1 = new TCanvas("hhpt", "a Canvas", 10, 10, 600, 500);
   c1->Divide(1, 1);

   // 去掉上、右边框的刻度线
   c1->SetTickx(0);
   c1->SetTicky(0);

   // 保证四个边框线粗细一样
   gStyle->SetFrameLineWidth(1);
   gStyle->SetFrameLineColor(1);
   gStyle->SetPadTickX(0);
   gStyle->SetPadTickY(0);

   // ================= 直方图（仅用于做坐标框架）=================
   TH1D *h1 = new TH1D("h1", "", 100, 0.0, 50.0);
   TH1D *h2 = new TH1D("h2", "", 100, 0.0, 50.0);
   TH1D *h3 = new TH1D("h3", "", 100, 0.0, 50.0);

   h1->SetLineColor(1);
   h1->SetLineWidth(1);
   h1->SetLineStyle(1);

   h2->SetLineColor(2);
   h2->SetLineWidth(1);
   h2->SetLineStyle(1);

   h3->SetLineColor(4);
   h3->SetLineWidth(1);
   h3->SetLineStyle(2);

   // ================= 数据读取 =================
   ifstream fort;

   Double_t x1[100], y1[100];
   Double_t x2[40],  y2[40];
   Double_t x3[40],  y3[40];
   Double_t x4[40],  y4[40];
   Double_t x5[100], y5[100];
   Double_t xdat1[8], ydat1[8], eyl1[8], eyh1[8];
   Double_t exl1[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
   Double_t exh1[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

   Int_t d = 0;
   Int_t nbin;

   fort.open("referline1.dat");
   if (fort.fail()) {
      cout << "sorry, couldn't open file" << endl;
      exit(1);
   }
   d = 0;
   while (fort >> x1[d] >> y1[d]) d++;
   fort.close();

   fort.open("D0v2-auau200-10-40y.dat");
   if (fort.fail()) {
      cout << "sorry, couldn't open file" << endl;
      exit(1);
   }
   d = 0;
   while (fort >> x2[d] >> y2[d]) d++;
   fort.close();

   fort.open("D0v2-auau200-10-40s.dat");
   if (fort.fail()) {
      cout << "sorry, couldn't open file" << endl;
      exit(1);
   }
   d = 0;
   while (fort >> x3[d] >> y3[d]) d++;
   fort.close();

   fort.open("D0v2-auau200-10-40sy.dat");
   if (fort.fail()) {
      cout << "sorry, couldn't open file" << endl;
      exit(1);
   }
   d = 0;
   while (fort >> x4[d] >> y4[d]) d++;
   fort.close();

   fort.open("D0v2-10-40-d3.dat");
   if (fort.fail()) {
      cout << "sorry, couldn't open file" << endl;
      exit(1);
   }
   d = 0;
   while (fort >> x5[d] >> y5[d]) d++;
   fort.close();

   fort.open("D0v2(10-40).txt");
   if (fort.fail()) {
      cout << "sorry, couldn't open file" << endl;
      exit(1);
   }
   d = 0;
   while (fort >> xdat1[d] >> ydat1[d] >> eyl1[d] >> eyh1[d]) d++;
   fort.close();

   // ================= 填充直方图（仅用于 frame）=================
   for (Int_t i = 0; i < 100; i++) {
      nbin = (x1[i] + 0.251) / 0.5;   // 原逻辑
      if (nbin < 1) nbin = 1;         // 防止 0 号 bin 越界
      if (nbin > 100) nbin = 100;
      h1->SetBinContent(nbin, y1[i]);
      h2->SetBinContent(nbin, y2[i]);
      h3->SetBinContent(nbin, y3[i]);
   }

   // ================= 绘制 frame =================
   c1->cd(1);
   h1->Draw();

   // 四边框线粗细一致（关键）
   h1->SetLineWidth(1);

   h1->GetYaxis()->SetRangeUser(0.0, 0.4);
   h1->GetXaxis()->SetRangeUser(0.0, 10.0);
   h1->GetYaxis()->SetTitle("v_{2}");
   h1->GetXaxis()->SetTitle("p_{T} (GeV) ");
   h1->GetXaxis()->SetLabelSize(0.05);
   h1->GetYaxis()->SetLabelSize(0.05);
   h1->GetYaxis()->SetNdivisions(005);
   h1->GetXaxis()->SetNdivisions(010);
   h1->GetYaxis()->CenterTitle();
   h1->GetYaxis()->SetTitleOffset(0.75);
   h1->GetYaxis()->SetTitleSize(0.07);
   h1->GetXaxis()->CenterTitle();
   h1->GetXaxis()->SetTitleOffset(0.9);
   h1->GetXaxis()->SetTitleSize(0.06);

   // ================= 理论曲线 =================
   TGraph *gr1 = new TGraph(40, x2, y2);
   gr1->SetLineColor(8);
   gr1->SetMarkerStyle(20);
   gr1->SetMarkerSize(1.0);
   gr1->SetLineWidth(5);
   gr1->SetLineStyle(9);
   gr1->SetMarkerColor(8);
   gr1->Draw("l same");

   TGraph *gr2 = new TGraph(40, x3, y3);
   gr2->SetLineColor(4);
   gr2->SetMarkerStyle(22);
   gr2->SetMarkerSize(1.0);
   gr2->SetLineWidth(5);
   gr2->SetLineStyle(6);
   gr2->SetMarkerColor(4);
   gr2->Draw("l same");

   TGraph *gr4 = new TGraph(40, x4, y4);
   gr4->SetLineColor(2);
   gr4->SetMarkerStyle(20);
   gr4->SetMarkerSize(1.0);
   gr4->SetLineWidth(5);
   gr4->SetLineStyle(1);
   gr4->SetMarkerColor(1);
   gr4->Draw("l same");

   TGraph *gr5 = new TGraph(100, x5, y5);
   gr5->SetLineColor(7);
   gr5->SetMarkerStyle(20);
   gr5->SetMarkerSize(1.0);
   gr5->SetLineWidth(5);
   gr5->SetLineStyle(5);
   gr5->SetMarkerColor(1);
   gr5->Draw("l same");

   // ================= 实验数据点 =================
   TGraphAsymmErrors *gr13 = new TGraphAsymmErrors(
      8, xdat1, ydat1, exl1, exh1, eyl1, eyh1);
   gr13->SetLineColor(1);
   gr13->SetMarkerStyle(20);
   gr13->SetMarkerSize(1.8);
   gr13->SetLineWidth(2);
   gr13->SetMarkerColor(1);
   gr13->Draw("P same");

   // ================= 图例 =================
   
   TLegend *leg1 = new TLegend(0.55, 0.80, 0.65, 0.90);
   leg1->SetFillColor(0);
   leg1->SetLineColor(0);
   leg1->SetTextSize(0.05);
   leg1->AddEntry(gr2,"Au-Au #sqrt{s_{NN}} = 200 GeV","");
   leg1->AddEntry(gr1, "10-40%  D^{0} ", "");
   leg1->Draw();
   
   
   TLegend *l = new TLegend(0.4, 0.55, 0.88, 0.88);
   l->SetTextSize(0.05);
   l->AddEntry(gr13, "STAR",          "p");
   l->AddEntry(gr1,  "Yukawa",        "l");
   l->AddEntry(gr2,  "string",        "l");
   l->AddEntry(gr4,  "Yukawa + string", "l");
   l->AddEntry(gr5,  "D_{s}(2#piT) = 3",    "l");
   l->SetFillColor(0);
   l->SetLineColor(0);
   l->Draw();
}
