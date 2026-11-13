/***********************************************************************
 * root -l
 * .x hhraa.cxx
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

void hhraa()
{
   // ================= 画布 =================
   TCanvas *c1 = new TCanvas("hhraa", "a Canvas", 10, 10, 600, 500);
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
   TH1D *h1 = new TH1D("h1", "", 100, 0.0, 200.0);
   TH1D *h2 = new TH1D("h2", "", 100, 0.0, 200.0);

   h1->SetLineColor(0);
   h1->SetLineWidth(1);
   h1->SetLineStyle(2);

   h2->SetLineColor(1);
   h2->SetLineWidth(1);
   h2->SetLineStyle(1);

   // ================= 数据读取 =================
   ifstream fort;
   Double_t x1[100], y1[100];
   Double_t x2[100], y2[100];
   Double_t x3[100], y3[100];
   Double_t x4[100], y4[100];
   Double_t x5[100], y5[100];
   Double_t x6[100], y6[100];
   Double_t x7[100], y7[100];
   Double_t x8[100], y8[100];

   Double_t xdat1[14], ydat1[14], eyl1[14], eyh1[14];
   Double_t xdat2[14], ydat2[14], eyl2[14], eyh2[14];
   Double_t xdat3[14], ydat3[14], eyl3[14], eyh3[14];
   Double_t exl1[14] = {0}, exh1[14] = {0};
   Double_t exl2[14] = {0}, exh2[14] = {0};
   Double_t exl3[14] = {0}, exh3[14] = {0};

   Int_t d = 0;
   Int_t nbin;

   fort.open("referline1.dat");
   if (fort.fail()) { cout << "sorry, couldn't open file" << endl; exit(1); }
   d = 0;
   while (fort >> x1[d] >> y1[d]) d++;
   fort.close();

   fort.open("D0raa-pbpb5020-0-10y.dat");
   if (fort.fail()) { cout << "sorry, couldn't open file" << endl; exit(1); }
   d = 0;
   while (fort >> x2[d] >> y2[d]) d++;
   fort.close();

   fort.open("D0raa-pbpb5020-0-10s.dat");
   if (fort.fail()) { cout << "sorry, couldn't open file" << endl; exit(1); }
   d = 0;
   while (fort >> x3[d] >> y3[d]) d++;
   fort.close();

   fort.open("D0raa-pbpb5020-0-10sy.dat");
   if (fort.fail()) { cout << "sorry, couldn't open file" << endl; exit(1); }
   d = 0;
   while (fort >> x4[d] >> y4[d]) d++;
   fort.close();

   fort.open("D0raafnaver0-10.dat");
   if (fort.fail()) { cout << "sorry, couldn't open file" << endl; exit(1); }
   d = 0;
   while (fort >> x5[d] >> y5[d]) d++;
   fort.close();

   fort.open("D0raa_alice0-10.txt");
   if (fort.fail()) { cout << "sorry, couldn't open file" << endl; exit(1); }
   d = 0;
   while (fort >> xdat1[d] >> ydat1[d] >> eyl1[d] >> eyh1[d]) d++;
   fort.close();

   fort.open("D0raa_cms0-10stat.txt");
   if (fort.fail()) { cout << "sorry, couldn't open file" << endl; exit(1); }
   d = 0;
   while (fort >> xdat2[d] >> ydat2[d] >> eyl2[d] >> eyh2[d]) d++;
   fort.close();

   fort.open("D0raa_cms0-10sys.txt");
   if (fort.fail()) { cout << "sorry, couldn't open file" << endl; exit(1); }
   d = 0;
   while (fort >> xdat3[d] >> ydat3[d] >> eyl3[d] >> eyh3[d]) d++;
   fort.close();

   // ================= 填充直方图（仅用于 frame）=================
   for (Int_t i = 0; i < 100; i++) {
      nbin = (x1[i] + 1.0) / 2.0;   // 原逻辑
      if (nbin < 1) nbin = 1;       // 防止 0 号 bin 越界
      if (nbin > 100) nbin = 100;
      h1->SetBinContent(nbin, y1[i]);
      // h2->SetBinContent(nbin, y5[i]);
   }

   // ================= 绘制 frame =================
   c1->cd(1);
   h1->Draw();

   // 四边框线粗细一致（关键）
   h1->SetLineWidth(1);

   // gPad->SetLogx(1);
   h1->GetYaxis()->SetRangeUser(0.0, 2.0);
   h1->GetXaxis()->SetRangeUser(0.0, 100.0);
   h1->GetYaxis()->SetTitle("R_{AA}");
   h1->GetXaxis()->SetTitle("p_{T} (GeV) ");
   h1->GetXaxis()->SetLabelSize(0.05);
   h1->GetYaxis()->SetLabelSize(0.05);
   h1->GetYaxis()->SetNdivisions(010);
   h1->GetXaxis()->SetNdivisions(010);
   h1->GetYaxis()->CenterTitle();
   h1->GetYaxis()->SetTitleOffset(0.7);
   h1->GetYaxis()->SetTitleSize(0.06);
   h1->GetXaxis()->CenterTitle();
   h1->GetXaxis()->SetTitleOffset(0.9);
   h1->GetXaxis()->SetTitleSize(0.06);

   // ================= 理论曲线 =================
   TGraph *gr2 = new TGraph(100, x2, y2);
   gr2->SetLineColor(8);
   gr2->SetMarkerStyle(20);
   gr2->SetMarkerSize(1.0);
   gr2->SetLineWidth(5);
   gr2->SetLineStyle(9);
   gr2->SetMarkerColor(2);
   gr2->Draw("l same");

   TGraph *gr3 = new TGraph(100, x3, y3);
   gr3->SetLineColor(4);
   gr3->SetMarkerStyle(20);
   gr3->SetMarkerSize(1.0);
   gr3->SetLineWidth(5);
   gr3->SetLineStyle(6);
   gr3->SetMarkerColor(4);
   gr3->Draw("l same");

   TGraph *gr4 = new TGraph(100, x4, y4);
   gr4->SetLineColor(2);
   gr4->SetMarkerStyle(20);
   gr4->SetMarkerSize(1.0);
   gr4->SetLineWidth(5);
   gr4->SetLineStyle(1);
   gr4->SetMarkerColor(2);
   gr4->Draw("l same");

   TGraph *gr5 = new TGraph(100, x5, y5);
   gr5->SetLineColor(7);
   gr5->SetMarkerStyle(20);
   gr5->SetMarkerSize(1.0);
   gr5->SetLineWidth(5);
   gr5->SetLineStyle(5);
   gr5->SetMarkerColor(8);
   gr5->Draw("l same");

   // ================= 实验数据点 =================
   TGraphAsymmErrors *gr11 = new TGraphAsymmErrors(
      14, xdat1, ydat1, exl1, exh1, eyl1, eyh1);
   gr11->SetLineColor(1);
   gr11->SetMarkerStyle(20);
   gr11->SetMarkerSize(1.8);
   gr11->SetLineWidth(2);
   gr11->SetMarkerColor(1);
   gr11->Draw("P same");

   TGraphAsymmErrors *gr12 = new TGraphAsymmErrors(
      14, xdat2, ydat2, exl2, exh2, eyl2, eyh2);
   gr12->SetLineColor(1);
   gr12->SetMarkerStyle(24);
   gr12->SetMarkerSize(1.8);
   gr12->SetLineWidth(2);
   gr12->SetMarkerColor(1);
   gr12->Draw("P same");

   TGraphAsymmErrors *gr13 = new TGraphAsymmErrors(
      14, xdat3, ydat3, exl3, exh3, eyl3, eyh3);
   gr13->SetLineColor(1);
   gr13->SetMarkerStyle(24);
   gr13->SetMarkerSize(1.8);
   gr13->SetLineWidth(2);
   gr13->SetMarkerColor(1);
   gr13->Draw("P same");

   // ================= 图例 =================
   TLegend *l1 = new TLegend(0.15,0.75,0.85,0.85);
   l1->SetTextSize(0.05);
   l1->AddEntry(gr2,"Pb-Pb #sqrt{s_{NN}} = 5.02 TeV","");
   l1->AddEntry(gr3,"0-10%  D^{0} ","");
   l1->SetFillColor(0);
   l1->SetLineColor(0);
   l1->Draw();
   
   TLegend *l2 = new TLegend(0.15, 0.55, 0.88, 0.88);
   l2->SetTextSize(0.05);
   l2->AddEntry(gr11, "ALICE",          "p");   // PLB 782 (2018) 474-496, |y|<1.0
   l2->AddEntry(gr13, "CMS",            "p");
   l2->AddEntry(gr2,  "Yukawa",         "l");
   l2->AddEntry(gr3,  "string",         "l");
   l2->AddEntry(gr4,  "Yukawa + string",  "l");
   l2->AddEntry(gr5,  "D_{s}(2#piT) = 4",     "l");
   l2->SetFillColor(0);
   l2->SetLineColor(0);
   l2->Draw();
}
