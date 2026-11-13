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

   h1->SetLineColor(1);
   h1->SetLineWidth(1);   // 由 1 改为 2，保证和全局 frame 一致
   h1->SetLineStyle(1);

   h2->SetLineColor(2);
   h2->SetLineWidth(2);
   h2->SetLineStyle(1);

   // ================= 数据读取 =================
   ifstream fort;
   Double_t x1[100], y1[100];
   Double_t x2[100], y2[100];
   Double_t x3[100], y3[100];
   Double_t x4[100], y4[100];
   Double_t x8[100], y8[100];

   Double_t xdat1[10], ydat1[10], eyl1[10], eyh1[10];
   Double_t exl1[10] = {0}, exh1[10] = {0};

   Double_t xdat2[11], ydat2[11], eyl2[11], eyh2[11];
   Double_t exl2[11] = {0}, exh2[11] = {0};

   Int_t d = 0;
   Int_t nbin;

   fort.open("referline1.dat");
   if (fort.fail()) { cout << "sorry, couldn't open file" << endl; exit(1); }
   d = 0;
   while (fort >> x8[d] >> y8[d]) d++;
   fort.close();

   fort.open("D0v2-pbpb5020-30-50y.dat");
   if (fort.fail()) { cout << "sorry, couldn't open file" << endl; exit(1); }
   d = 0;
   while (fort >> x1[d] >> y1[d]) d++;
   fort.close();

   fort.open("D0v2-pbpb5020-30-50s.dat");
   if (fort.fail()) { cout << "sorry, couldn't open file" << endl; exit(1); }
   d = 0;
   while (fort >> x2[d] >> y2[d]) d++;
   fort.close();

   fort.open("D0v2-pbpb5020-30-50sy.dat");
   if (fort.fail()) { cout << "sorry, couldn't open file" << endl; exit(1); }
   d = 0;
   while (fort >> x3[d] >> y3[d]) d++;
   fort.close();

   fort.open("D0v2fnaver30-50.dat");
   if (fort.fail()) { cout << "sorry, couldn't open file" << endl; exit(1); }
   d = 0;
   while (fort >> x4[d] >> y4[d]) d++;
   fort.close();

   fort.open("D0v2_cms30-50.txt");
   if (fort.fail()) { cout << "sorry, couldn't open file" << endl; exit(1); }
   d = 0;
   while (fort >> xdat1[d] >> ydat1[d] >> eyl1[d] >> eyh1[d]) d++;
   fort.close();

   fort.open("D0v2_alice30-50.txt");
   if (fort.fail()) { cout << "sorry, couldn't open file" << endl; exit(1); }
   d = 0;
   while (fort >> xdat2[d] >> ydat2[d] >> eyl2[d] >> eyh2[d]) d++;
   fort.close();

   // ================= 填充直方图（仅用于 frame）=================
   for (Int_t i = 0; i < 100; i++) {
      nbin = (x8[i] + 1.0) / 2.0;   // 原逻辑
      if (nbin < 1) nbin = 1;       // 防止 0 号 bin 越界
      if (nbin > 100) nbin = 100;
      h1->SetBinContent(nbin, y8[i]);
   }

   // ================= 绘制 frame =================
   c1->cd(1);
   h1->Draw();

   // 四边框线粗细一致（关键）
   h1->SetLineWidth(1);

   // gPad->SetLogx(1);
   h1->GetYaxis()->SetRangeUser(-0.1, 0.4);
   h1->GetXaxis()->SetRangeUser(0.0, 40.0);
   h1->GetYaxis()->SetTitle("v_{2}");
   h1->GetXaxis()->SetTitle("p_{T} (GeV) ");
   h1->GetXaxis()->SetLabelSize(0.05);
   h1->GetYaxis()->SetLabelSize(0.05);
   h1->GetYaxis()->SetNdivisions(010);
   h1->GetXaxis()->SetNdivisions(010);
   h1->GetYaxis()->CenterTitle();
   h1->GetYaxis()->SetTitleOffset(0.7);
   h1->GetYaxis()->SetTitleSize(0.07);
   h1->GetXaxis()->CenterTitle();
   h1->GetXaxis()->SetTitleOffset(0.9);
   h1->GetXaxis()->SetTitleSize(0.06);

   // ================= 理论曲线 =================
   TGraph *gr1 = new TGraph(100, x1, y1);
   gr1->SetLineColor(8);
   gr1->SetMarkerStyle(20);
   gr1->SetMarkerSize(1.0);
   gr1->SetLineWidth(5);
   gr1->SetLineStyle(9);
   gr1->SetMarkerColor(2);
   gr1->Draw("l same");

   TGraph *gr2 = new TGraph(100, x2, y2);
   gr2->SetLineColor(4);
   gr2->SetMarkerStyle(20);
   gr2->SetMarkerSize(1.0);
   gr2->SetLineWidth(5);
   gr2->SetLineStyle(6);
   gr2->SetMarkerColor(4);
   gr2->Draw("l same");

   TGraph *gr3 = new TGraph(100, x3, y3);
   gr3->SetLineColor(2);
   gr3->SetMarkerStyle(20);
   gr3->SetMarkerSize(1.0);
   gr3->SetLineWidth(5);
   gr3->SetLineStyle(1);
   gr3->SetMarkerColor(7);
   gr3->Draw("l same");

   TGraph *gr4 = new TGraph(100, x4, y4);
   gr4->SetLineColor(7);
   gr4->SetMarkerStyle(20);
   gr4->SetMarkerSize(1.0);
   gr4->SetLineWidth(5);
   gr4->SetLineStyle(5);
   gr4->SetMarkerColor(1);
   gr4->Draw("l same");

   // ================= 实验数据点 =================
   TGraphAsymmErrors *gr11 = new TGraphAsymmErrors(
      10, xdat1, ydat1, exl1, exh1, eyl1, eyh1);
   gr11->SetLineColor(1);
   gr11->SetMarkerStyle(24);
   gr11->SetMarkerSize(1.8);
   gr11->SetLineWidth(2);
   gr11->SetMarkerColor(1);
   gr11->Draw("P same");

   TGraphAsymmErrors *gr12 = new TGraphAsymmErrors(
      11, xdat2, ydat2, exl2, exh2, eyl2, eyh2);
   gr12->SetLineColor(1);
   gr12->SetMarkerStyle(20);
   gr12->SetMarkerSize(1.8);
   gr12->SetLineWidth(2);
   gr12->SetMarkerColor(1);
   gr12->Draw("P same");

   // ================= 图例 =================
   // 标题性图例（左上角空白处）
   TLegend *leg1 = new TLegend(0.15, 0.78, 0.45, 0.88);
   leg1->SetFillColor(0);
   leg1->SetLineColor(0);
   leg1->SetTextSize(0.05);
   leg1->AddEntry(gr1, "Pb-Pb #sqrt{s_{NN}} = 5.02 TeV", "");
   leg1->AddEntry(gr1, "30-50%  D^{0} ", "");
   leg1->Draw();

   // 实验数据 + 理论曲线图例
   TLegend *leg2 = new TLegend(0.50, 0.55, 0.90, 0.80);
   leg2->SetFillColor(0);
   leg2->SetLineColor(0);
   leg2->SetTextSize(0.05);
   leg2->AddEntry(gr12, "ALICE",            "p");
   leg2->AddEntry(gr11, "CMS",          "p");
   leg2->AddEntry(gr1,  "Yukawa",         "l");
   leg2->AddEntry(gr2,  "string",         "l");
   leg2->AddEntry(gr3,  "Yukawa + string",  "l");
   leg2->AddEntry(gr4,  "D_{s}(2#piT) = 4",     "l");
   leg2->Draw();
}
