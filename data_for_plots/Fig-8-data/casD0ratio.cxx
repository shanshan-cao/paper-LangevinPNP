/***********************************************************************
 * LamD0ratio.cxx
 *   Λc+ / D0 ratio vs pT  (0-10% Pb-Pb 5.02 TeV)
 *   Compare models (Yukawa / String / Yukawa+String) with ALICE data
 *
 *   Usage:
 *     root -l -b -q 'LamD0ratio.cxx+'
 ***********************************************************************/
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>

#include "TCanvas.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TError.h"

using namespace std;

// ---------------------------------------------------------------------
// 读入两列 (x, y) 数据文件
// ---------------------------------------------------------------------
static bool ReadXY(const string& fname,
                   vector<double>& x,
                   vector<double>& y)
{
    ifstream fin(fname);
    if (!fin) {
        Error("ReadXY", "cannot open %s", fname.c_str());
        return false;
    }
    x.clear(); y.clear();
    double a, b;
    while (fin >> a >> b) {
        x.push_back(a);
        y.push_back(b);
    }
    return true;
}

// ---------------------------------------------------------------------
// 读入四列 (x, y, eyl, eyh) ALICE 数据文件
// ---------------------------------------------------------------------
static bool ReadXYErr(const string& fname,
                      vector<double>& x,
                      vector<double>& y,
                      vector<double>& eyl,
                      vector<double>& eyh)
{
    ifstream fin(fname);
    if (!fin) {
        Error("ReadXYErr", "cannot open %s", fname.c_str());
        return false;
    }
    x.clear(); y.clear(); eyl.clear(); eyh.clear();
    double a, b, c, d;
    while (fin >> a >> b >> c >> d) {
        x.push_back(a);
        y.push_back(b);
        eyl.push_back(c);
        eyh.push_back(d);
    }
    return true;
}

// ---------------------------------------------------------------------
// 统一坐标轴样式
// ---------------------------------------------------------------------
static void StyleAxis(TH1* h)
{
    h->GetXaxis()->SetTitle("p_{T} (GeV)");
    h->GetYaxis()->SetTitle("#Xi^{0}_{c}/D^{0}");
    h->GetXaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetTitleSize(0.06);
    h->GetXaxis()->SetTitleOffset(0.9);
    h->GetYaxis()->SetTitleOffset(0.9);
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetLabelSize(0.05);
    h->GetXaxis()->SetNdivisions(010);
    h->GetYaxis()->SetNdivisions(010);
    h->GetXaxis()->CenterTitle();
    h->GetYaxis()->CenterTitle();
}

// ---------------------------------------------------------------------
// 主函数
// ---------------------------------------------------------------------
void casD0ratio()
{
    // ---------------- 全局样式：去掉右边 & 上边刻度线 ----------------
    gStyle->SetPadTickX(0);   // 顶边不画 x 刻度线
    gStyle->SetPadTickY(0);   // 右边不画 y 刻度线
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetFrameLineWidth(1);
    

    // 如需彻底隐藏上/右边框（不要框），再加：
    // gStyle->SetFrameBorderMode(0);

    // gROOT->SetBatch(kTRUE);
    


    // ---------- 1. 读数据 ----------
    vector<double> x1,y1, x2,y2, x3,y3, x4,y4, x5,y5, x6,y6, x8,y8;
    vector<double> xdat,ydat,eyl,eyh;
    Double_t exl[4]={0.5,1.0,1.0,2.0}; 
    Double_t exh[4]={0.5,1.0,1.0,2.0}; 

    if (!ReadXY   ("cas0-pt0-10y.dat",  x1, y1)) return;
    if (!ReadXY   ("D0-pt0-10y.dat",    x2, y2)) return;
    if (!ReadXY   ("cas0-pt0-10s.dat",  x3, y3)) return;
    if (!ReadXY   ("D0-pt0-10s.dat",    x4, y4)) return;
    if (!ReadXY   ("cas0-pt0-10sy.dat", x5, y5)) return;
    if (!ReadXY   ("D0-pt0-10sy.dat",   x6, y6)) return;
    if (!ReadXY   ("referline1.dat",    x8, y8)) return;
    if (!ReadXYErr("ALICE-cas0D0-ratio0-10.txt",
                   xdat, ydat, eyl, eyh)) return;

    // ---------- 2. 计算比值 ----------
    const size_t N = min({x1.size(), x2.size(), x3.size(),
                          x4.size(), x5.size(), x6.size()});
    vector<double> rY(N), rS(N), rSY(N);
    for (size_t i = 0; i < N; ++i) {
        rY [i] = (y2[i] != 0.0) ? y1[i] / y2[i] : 0.0;
        rS [i] = (y4[i] != 0.0) ? y3[i] / y4[i] : 0.0;
        rSY[i] = (y6[i] != 0.0) ? y5[i] / y6[i] : 0.0;
    }

    // ---------- 3. 画布 ----------
    TCanvas* c1 = new TCanvas("hhraa", "#Xi^{0}_{c}/D^{0}", 600, 500);
    // 注意：c1->SetTicks(1,1) 会强制四边都画刻度，
    // 与下面 gStyle 的 SetPadTickX/Y(0) 效果相反，所以这里去掉或设为 (0,0)
    c1->SetTicks(0, 0);
    c1->SetLeftMargin(0.14);
    c1->SetBottomMargin(0.13);

    TH1D* hFrame = new TH1D("hFrame", "", 100, 0.0, 200.0);
    hFrame->SetLineColor(0);
    StyleAxis(hFrame);
    hFrame->GetYaxis()->SetRangeUser(0.0, 2.5);
    hFrame->GetXaxis()->SetRangeUser(0.0, 14.0);

    // 单独对直方图再关一次顶/右刻度（防止 TCanvas::SetTicks 覆盖全局）
    hFrame->GetXaxis()->SetTickLength(0.0);
    hFrame->GetYaxis()->SetTickLength(0.0);
    // 恢复底/左刻度长度（只关顶/右）：
    hFrame->GetXaxis()->SetTickLength(0.03);
    hFrame->GetYaxis()->SetTickLength(0.03);

    hFrame->Draw();

    // ---------- 4. 模型曲线 ----------
    auto MakeGraph = [](const vector<double>& x, const vector<double>& y,
                        int color, int style) -> TGraph*
    {
        TGraph* g = new TGraph((int)x.size(), x.data(), y.data());
        g->SetLineColor(color);
        g->SetLineWidth(5);
        g->SetLineStyle(style);
        return g;
    };

    TGraph* grY  = MakeGraph(x2, rY,  8, 9);   // Yukawa
    TGraph* grS  = MakeGraph(x4, rS,  4, 6);   // String
    TGraph* grSY = MakeGraph(x6, rSY, 2, 1);   // Yukawa+String

    grY ->Draw("L same");
    grS ->Draw("L same");
    grSY->Draw("L same");

    // ---------- 5. ALICE 数据点 ----------
    TGraphAsymmErrors* grAlice =
        new TGraphAsymmErrors((int)xdat.size(),
                              xdat.data(), ydat.data(),
                               exl,exh,
                              eyl.data(), eyh.data());
    grAlice->SetLineColor(1);
    grAlice->SetMarkerColor(1);
    grAlice->SetMarkerStyle(20);
    grAlice->SetMarkerSize(1.8);
    grAlice->SetLineWidth(2);
    grAlice->Draw("P same");

    // ---------- 6. 图例 ----------
    TLegend* l1= new TLegend(0.42, 0.8, 0.90, 0.88);
    l1->SetBorderSize(0);
    l1->SetFillStyle(0);
    l1->SetTextSize(0.05);
    l1->AddEntry(grAlice,"Pb-Pb #sqrt{s_{NN}} = 5.02 TeV"," "); 
    l1->AddEntry(grAlice,"0-10%"," ");     
    l1->Draw();
    
    
    TLegend* leg = new TLegend(0.60, 0.60, 0.90, 0.79);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.05);
    leg->AddEntry(grAlice,"ALICE ","p");    //arXiv:2607.17903 |y|<0.5 
    leg->AddEntry(grY,  "Yukawa",          "l");
    leg->AddEntry(grS,  "string",          "l");
    leg->AddEntry(grSY, "Yukawa + string",   "l");
    leg->Draw();

    // ---------- 7. 可选保存 ----------
    // c1->SaveAs("LamD0ratio.pdf");
    // c1->SaveAs("LamD0ratio.png");
}
