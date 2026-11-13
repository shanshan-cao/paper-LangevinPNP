/***********************************************************************
 * root -l
 * .x read_draw_opt.C
 * .q
 ***********************************************************************/
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "TH1F.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TPad.h"
#include "TLine.h"
#include "TMath.h"
using namespace std;

// ---------- 读取两列数据文件 ----------
// 返回实际读入的点数
Int_t ReadXY(const char *fname, Double_t *x, Double_t *y, Int_t maxN)
{
    ifstream fin(fname);
    if (fin.fail()) {
        cout << "sorry, couldn't open file: " << fname << endl;
        exit(1);
    }
    Int_t n = 0;
    while (n < maxN && (fin >> x[n] >> y[n])) n++;
    fin.close();
    return n;
}

void hqptratios()
{
    // ================= 1. 读取数据 =================
    const Int_t MAXN = 100;
    Double_t x1[MAXN], y1[MAXN];   // charm
    Double_t x2[MAXN], y2[MAXN];   // D0
    Double_t x3[MAXN], y3[MAXN];   // Ds
    Double_t x4[MAXN], y4[MAXN];   // Lamc
    Double_t x5[MAXN], y5[MAXN];   // Xi0c
    Double_t x6[MAXN], y6[MAXN];   // Omega0c
    Double_t x8[MAXN], y8[MAXN];   // referline

    Int_t n1 = ReadXY("charm-pt0-10s.dat", x1, y1, MAXN);
    Int_t n2 = ReadXY("D0-pt0-10s.dat",    x2, y2, MAXN);
    Int_t n3 = ReadXY("Ds-pt0-10s.dat",    x3, y3, MAXN);
    Int_t n4 = ReadXY("Lamc-pt0-10s.dat",  x4, y4, MAXN);
    Int_t n5 = ReadXY("cas0-pt0-10s.dat",  x5, y5, MAXN);
    Int_t n6 = ReadXY("Ome-pt0-10s.dat",   x6, y6, MAXN);
    Int_t n8 = ReadXY("referline1.dat",    x8, y8, MAXN);

    // ================= 2. Canvas 与 Pad（3:1）=================
    TCanvas *c1 = new TCanvas("hhraa", "a Canvas", 10, 10, 600, 500);

    // 上图占 75%，下图占 25%，即 3:1
    TPad *padTop = new TPad("padTop", "", 0.0, 0.4, 1.0, 1.0);
    padTop->SetBottomMargin(0.00);
    padTop->SetTopMargin(0.08);
    padTop->SetLeftMargin(0.14);
    padTop->SetRightMargin(0.05);
    padTop->Draw();

    TPad *padBot = new TPad("padBot", "", 0.0, 0.0, 1.0, 0.4);
    padBot->SetTopMargin(0.00);
    padBot->SetBottomMargin(0.30);
    padBot->SetLeftMargin(0.14);
    padBot->SetRightMargin(0.05);
    padBot->Draw();

    // ================= 3. 上图：p_T 谱 =================
    padTop->cd();

    // 参考线直方图作为框架
    TH1D *hFrame = new TH1D("hFrame", "", 100, 0.0, 200.0);
    for (Int_t i = 0; i < n8; i++) {
        Int_t nbin = (x8[i] + 1.0) / 2.0;
        hFrame->SetBinContent(nbin, y8[i]);
    }
    hFrame->SetLineColor(0);
    hFrame->SetLineWidth(1);
    hFrame->SetLineStyle(1);
    hFrame->SetStats(0);
    hFrame->Draw("HIST");
    gPad->SetLogy(1);
    hFrame->GetYaxis()->SetRangeUser(1.0e-8, 1.0e2);
    hFrame->GetXaxis()->SetRangeUser(0.0, 20.0);
    hFrame->GetYaxis()->SetTitle("dN^{2}/dp_{T}dy (GeV^{-1})");
    hFrame->GetYaxis()->SetTitleSize(0.09);
    hFrame->GetYaxis()->SetTitleOffset(0.6);
    hFrame->GetYaxis()->SetLabelSize(0.06);
    hFrame->GetYaxis()->SetNdivisions(505);
    hFrame->GetYaxis()->CenterTitle();
    hFrame->GetXaxis()->SetLabelSize(0);   // 上图不显示 X 标签
    hFrame->GetXaxis()->SetNdivisions(010);

    // 各粒子曲线
    TGraph *gr1 = new TGraph(n1, x1, y1);
    gr1->SetLineColor(2);          gr1->SetLineWidth(5);  gr1->SetLineStyle(1);
    gr1->Draw("l");

    TGraph *gr2 = new TGraph(n2, x2, y2);
    gr2->SetLineColor(4);          gr2->SetLineWidth(5);  gr2->SetLineStyle(9);
    gr2->Draw("l");

    TGraph *gr3 = new TGraph(n3, x3, y3);
    gr3->SetLineColor(6);          gr3->SetLineWidth(5);  gr3->SetLineStyle(6);
    gr3->Draw("l");

    TGraph *gr4 = new TGraph(n4, x4, y4);
    gr4->SetLineColor(7);          gr4->SetLineWidth(5);  gr4->SetLineStyle(5);
    gr4->Draw("l");

    TGraph *gr5 = new TGraph(n5, x5, y5);
    gr5->SetLineColor(kGreen - 4); gr5->SetLineWidth(5);  gr5->SetLineStyle(2);
    gr5->Draw("l");

    TGraph *gr6 = new TGraph(n6, x6, y6);
    gr6->SetLineColor(kOrange + 1);gr6->SetLineWidth(5);  gr6->SetLineStyle(7);
    gr6->Draw("l");

    // 图例
    TLegend *leg1 = new TLegend(0.20, 0.62, 0.58, 0.90);
    leg1->AddEntry(gr1, "Pb-Pb #sqrt{s_{NN}} = 5.02 TeV   string", "");
    leg1->AddEntry(gr1, "0-10%", "");    
    leg1->SetFillColor(0);
    leg1->SetLineColor(0);
    leg1->Draw();
    
    TLegend *leg = new TLegend(0.20, 0.62, 0.58, 0.90);
    leg->AddEntry(gr1, "charm",          "l");
    leg->AddEntry(gr2, "D^{0}",          "l");
    leg->AddEntry(gr3, "D^{+}_{s}",      "l");
    leg->SetFillColor(0);
    leg->SetLineColor(0);
    leg->Draw();    
    
    TLegend *leg2 = new TLegend(0.20, 0.62, 0.58, 0.90);
    leg2->AddEntry(gr4, "#Lambda^{+}_{c}","l");
    leg2->AddEntry(gr5, "#Xi^{0}_{c}",    "l");
    leg2->AddEntry(gr6, "#Omega^{0}_{c}", "l");
    leg2->SetFillColor(0);
    leg2->SetLineColor(0);
    leg2->Draw();  
    
    
    
    
    

    // ================= 4. 计算比值 =================
    Double_t r2[MAXN], r3[MAXN], r4[MAXN], r5[MAXN], r6[MAXN];
    for (Int_t i = 0; i < n2; i++) r2[i] = (y1[i] != 0) ? y2[i] / y1[i] : 0;
    for (Int_t i = 0; i < n3; i++) r3[i] = (y1[i] != 0) ? y3[i] / y1[i] : 0;
    for (Int_t i = 0; i < n4; i++) r4[i] = (y1[i] != 0) ? y4[i] / y1[i] : 0;
    for (Int_t i = 0; i < n5; i++) r5[i] = (y1[i] != 0) ? y5[i] / y1[i] : 0;
    for (Int_t i = 0; i < n6; i++) r6[i] = (y1[i] != 0) ? y6[i] / y1[i] : 0;

    // ================= 5. 下图：Ratio =================
    padBot->cd();

    TH1D *hRatioFrame = new TH1D("hRatioFrame", "", 100, 0.0, 200.0);
    hRatioFrame->SetStats(0);
    hRatioFrame->Draw("HIST");
    hRatioFrame->GetYaxis()->SetRangeUser(0.0, 0.5);
    hRatioFrame->GetXaxis()->SetRangeUser(0.0, 20.0);
    hRatioFrame->GetYaxis()->SetTitle("h/c");
    hRatioFrame->GetYaxis()->SetTitleSize(0.16);
    hRatioFrame->GetYaxis()->SetTitleOffset(0.35);
    hRatioFrame->GetYaxis()->SetLabelSize(0.11);
    hRatioFrame->GetYaxis()->SetNdivisions(010);
    hRatioFrame->GetYaxis()->CenterTitle();
    hRatioFrame->GetXaxis()->SetTitle("p_{T} (GeV)");
    hRatioFrame->GetXaxis()->SetTitleSize(0.16);
    hRatioFrame->GetXaxis()->SetTitleOffset(0.9);
    hRatioFrame->GetXaxis()->SetLabelSize(0.10);
    hRatioFrame->GetXaxis()->SetNdivisions(505);
    hRatioFrame->GetXaxis()->CenterTitle();

    TGraph *gr22 = new TGraph(n2, x2, r2);
    gr22->SetLineColor(4);          gr22->SetLineWidth(5);  gr22->SetLineStyle(9);
    gr22->Draw("l");

    TGraph *gr23 = new TGraph(n3, x3, r3);
    gr23->SetLineColor(6);          gr23->SetLineWidth(5);  gr23->SetLineStyle(6);
    gr23->Draw("l");

    TGraph *gr24 = new TGraph(n4, x4, r4);
    gr24->SetLineColor(7);          gr24->SetLineWidth(5);  gr24->SetLineStyle(5);
    gr24->Draw("l");

    TGraph *gr25 = new TGraph(n5, x5, r5);
    gr25->SetLineColor(kGreen - 4); gr25->SetLineWidth(5);  gr25->SetLineStyle(2);
    gr25->Draw("l");

    TGraph *gr26 = new TGraph(n6, x6, r6);
    gr26->SetLineColor(kOrange + 1);gr26->SetLineWidth(5);  gr26->SetLineStyle(7);
    gr26->Draw("l");

    // 参考线 y = 1
    TLine *line = new TLine(0, 1, 20, 1);
    line->SetLineStyle(2);
    line->SetLineColor(kRed);
    line->Draw();

    c1->cd();
    c1->Update();
}
