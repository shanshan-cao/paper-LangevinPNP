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

// ---------- 读取两列数据文件，返回实际点数 ----------
Int_t Read2Col(const char *fname, Double_t *x, Double_t *y)
{
   ifstream fort(fname);
   if (fort.fail()) {
      cout << "sorry, couldn't open file: " << fname << endl;
      exit(1);
   }
   Int_t d = 0;
   while (fort >> x[d] >> y[d]) d++;
   fort.close();
   return d;
}

// ---------- 读取四列数据文件 (x, y, eyl, eyh) ----------
Int_t Read4Col(const char *fname, Double_t *x, Double_t *y,
               Double_t *eyl, Double_t *eyh)
{
   ifstream fort(fname);
   if (fort.fail()) {
      cout << "sorry, couldn't open file: " << fname << endl;
      exit(1);
   }
   Int_t d = 0;
   while (fort >> x[d] >> y[d] >> eyl[d] >> eyh[d]) d++;
   fort.close();
   return d;
}

void hhraa()
{
   // ================= 画布 =================
   TCanvas *c1 = new TCanvas("hhpt", "a Canvas", 10, 10, 600, 500);

   // 去掉上、右边框的刻度线
   c1->SetTickx(0);
   c1->SetTicky(0);

   // 保证四个边框线粗细一样
   gStyle->SetFrameLineWidth(1);
   gStyle->SetFrameLineColor(1);
   gStyle->SetPadTickX(0);
   gStyle->SetPadTickY(0);

   // ================= 数据读取 =================
   Double_t x1[200], y1[200];   // referline1
   Double_t x2[200], y2[200];   // Yukawa
   Double_t x3[200], y3[200];   // String
   Double_t x4[200], y4[200];   // Yukawa + String
   Double_t x5[200], y5[200];   // D(2piT)=3

   Double_t xdat1[50], ydat1[50], eyl1[50], eyh1[50];  // STAR 2014
   Double_t xdat2[50], ydat2[50], eyl2[50], eyh2[50];  // STAR 2010/11

   Int_t n1 = Read2Col("referline1.dat",            x1, y1);
   Int_t n2 = Read2Col("D0raa-auau200-0-10y.dat",   x2, y2);
   Int_t n3 = Read2Col("D0raa-auau200-0-10s.dat",   x3, y3);
   Int_t n4 = Read2Col("D0raa-auau200-0-10sy.dat",  x4, y4);
   Int_t n5 = Read2Col("D0raa-0-10-d3.dat",         x5, y5);

   Int_t nd1 = Read4Col("D0raa(0-10)2014.1812.txt", xdat1, ydat1, eyl1, eyh1);
   Int_t nd2 = Read4Col("D0raa(0-10)2010.1812.txt", xdat2, ydat2, eyl2, eyh2);

   // 对称误差的 x 误差设为 0
   Double_t exl1[50] = {0}, exh1[50] = {0};
   Double_t exl2[50] = {0}, exh2[50] = {0};

   // ================= 参考线（用直方图做坐标框架）=================
   TH1D *hframe = new TH1D("hframe", "", 100, 0.0, 50.0);
   hframe->SetLineColor(0);
   hframe->GetYaxis()->SetRangeUser(0.0, 2.5);
   hframe->GetXaxis()->SetRangeUser(0.0, 10.0);
   hframe->GetYaxis()->SetTitle("R_{AA}");
   hframe->GetXaxis()->SetTitle("p_{T} (GeV)");
   hframe->GetYaxis()->SetTitleOffset(0.7);
   hframe->GetXaxis()->SetTitleOffset(0.9);
   hframe->GetYaxis()->SetTitleSize(0.06);
   hframe->GetXaxis()->SetTitleSize(0.06);
   hframe->GetYaxis()->SetLabelSize(0.05);
   hframe->GetXaxis()->SetLabelSize(0.05);
   hframe->GetYaxis()->CenterTitle();
   hframe->GetXaxis()->CenterTitle();
   hframe->GetYaxis()->SetNdivisions(010);
   hframe->GetXaxis()->SetNdivisions(010);
   hframe->SetStats(0);

   // 四边框线粗细一致
   hframe->SetLineWidth(1);

   c1->cd();
   hframe->Draw();

   // ================= 理论曲线 =================
   // 颜色 + 线型都区分开

   // (1) Yukawa：绿色，点划线
   TGraph *gr1 = new TGraph(n2, x2, y2);
   gr1->SetLineColor(8);
   gr1->SetLineWidth(5);
   gr1->SetLineStyle(9);      // 虚线
   gr1->Draw("l same");

   // (2) String：蓝色，点线
   TGraph *gr2 = new TGraph(n3, x3, y3);
   gr2->SetLineColor(4);
   gr2->SetLineWidth(5);
   gr2->SetLineStyle(6);      // 实线
   gr2->Draw("l same");

   // (3) Yukawa+String：红色，实线
   TGraph *gr3 = new TGraph(n4, x4, y4);
   gr3->SetLineColor(2);
   gr3->SetLineWidth(5);
   gr3->SetLineStyle(1);      // 点划线
   gr3->Draw("l same");

   // (4) D(2piT)=3：橙色，长虚线
   TGraph *gr4 = new TGraph(n5, x5, y5);
   gr4->SetLineColor(7);
   gr4->SetLineWidth(5);
   gr4->SetLineStyle(5);     // 长虚线
   gr4->Draw("l same");

   // ================= 实验数据点：统一黑色，区分实心/空心 =================

   // STAR 2014：黑色实心圆 + 误差棒
   TGraphAsymmErrors *gr11 = new TGraphAsymmErrors(
      nd1, xdat1, ydat1, exl1, exh1, eyl1, eyh1);
   gr11->SetLineColor(1);
   gr11->SetMarkerColor(1);
   gr11->SetMarkerStyle(20);   // 实心圆
   gr11->SetMarkerSize(1.8);
   gr11->SetLineWidth(2);
   gr11->Draw("P same");

   // STAR 2010/11：黑色空心方块 + 误差棒
   TGraphAsymmErrors *gr12 = new TGraphAsymmErrors(
      nd2, xdat2, ydat2, exl2, exh2, eyl2, eyh2);
   gr12->SetLineColor(1);
   gr12->SetMarkerColor(1);
   gr12->SetMarkerStyle(24);   // 空心方圆
   gr12->SetMarkerSize(1.8);
   gr12->SetLineWidth(2);
   gr12->Draw("P same");

   // ================= 图例 =================
   // 理论曲线图例
   TLegend *leg1 = new TLegend(0.55, 0.80, 0.65, 0.90);
   leg1->SetFillColor(0);
   leg1->SetLineColor(0);
   leg1->SetTextSize(0.05);
   leg1->AddEntry(gr2,"Au-Au #sqrt{s_{NN}} = 200 GeV","");
   leg1->AddEntry(gr1, "0-10%  D^{0} ", "");
   leg1->Draw();

   // 实验数据图例
   TLegend *leg2 = new TLegend(0.5, 0.55, 0.90, 0.80);
   leg2->SetFillColor(0);
   leg2->SetLineColor(0);
   leg2->SetTextSize(0.05);
   leg2->AddEntry(gr11, "STAR 2014",    "p");
   leg2->AddEntry(gr12, "STAR 2010/11", "p");
   leg2->AddEntry(gr1, "Yukawa",          "l");
   leg2->AddEntry(gr2, "string",          "l");
   leg2->AddEntry(gr3, "Yukawa + string", "l");
   leg2->AddEntry(gr4, "D_{s}(2#piT) = 3",    "l");
   leg2->Draw();
}
