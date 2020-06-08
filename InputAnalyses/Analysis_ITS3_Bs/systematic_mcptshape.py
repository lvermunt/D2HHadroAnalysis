import pickle
from array import array
from ROOT import TFile, TH1F, kBlue, kBlack, kOrange, kRed, \
    TCanvas, TGraph, gPad, gStyle, TLegend
from root_numpy import fill_hist

def extract_Bs_eff(filereco, filegen, mlcut, fout):
    print("Opening reco file", filereco)
    dfreco = pickle.load(open(filereco,'rb'))
    print("Opening gen file", filegen)
    dfgen = pickle.load(open(filegen,'rb'))
    dfgen.query("abs(y_cand) < 0.5 and abs(z_vtx_gen) < 10")

    h_gen_pr = TH1F("h_gen_pr", "Prompt Generated in acceptance |y|<0.5", 1600, 0, 16)
    h_presel_pr = TH1F("h_presel_pr", "Prompt Reco in acc |#eta|<0.8 and presel", 1600, 0, 16)
    h_sel_pr = TH1F("h_sel_pr", "Prompt Reco in acc |#eta|<0.8 and sel", 1600, 0, 16)

    dfgen_pr = dfgen[dfgen.ismcprompt == 1]
    dfreco_presel_pr = dfreco[dfreco.ismcprompt == 1]
    dfreco_sel_pr = dfreco_presel_pr.query("y_test_probxgboost > @mlcut")
    print("Check lengths (FD rejected):", len(dfgen), len(dfgen_pr))

    fill_hist(h_gen_pr, dfgen_pr["pt_cand"])
    fill_hist(h_presel_pr, dfreco_presel_pr["pt_cand"])
    fill_hist(h_sel_pr, dfreco_sel_pr["pt_cand"])

    heff_presel = h_presel_pr.Clone("heff_presel")
    heff_sel = h_sel_pr.Clone("heff_sel")
    heff_presel.Divide(heff_presel, h_gen_pr, 1, 1, "B")
    heff_sel.Divide(heff_sel, h_gen_pr, 1, 1, "B")

    print("Are histos filled?", h_gen_pr.GetEntries(), h_sel_pr.GetEntries())
    print("Creating file", fout)
    feffDs = TFile.Open(fout, "recreate")
    feffDs.cd()
    h_gen_pr.Write()
    h_presel_pr.Write()
    h_sel_pr.Write()
    heff_presel.Write()
    heff_sel.Write()
    feffDs.Close()

def pythia_generated(pathgen, fout, isits2):
    filegen04 = pathgen + "AnalysisResultsGen_pt_cand0_4.pkl"
    filegen48 = pathgen + "AnalysisResultsGen_pt_cand4_8.pkl"
    filegen812 = pathgen + "AnalysisResultsGen_pt_cand8_12.pkl"
    filegen1216 = pathgen + "AnalysisResultsGen_pt_cand12_16.pkl"
    filegen1624 = pathgen + "AnalysisResultsGen_pt_cand16_24.pkl"

    print("Opening gen file", filegen04)
    dfgen04 = pickle.load(open(filegen04,'rb'))
    dfgen04.query("abs(y_cand) < 0.5 and abs(z_vtx_gen) < 10 and ismcprompt == 1")
    print("Opening gen file", filegen48)
    dfgen48 = pickle.load(open(filegen48,'rb'))
    dfgen48.query("abs(y_cand) < 0.5 and abs(z_vtx_gen) < 10 and ismcprompt == 1")
    print("Opening gen file", filegen812)
    dfgen812 = pickle.load(open(filegen812,'rb'))
    dfgen812.query("abs(y_cand) < 0.5 and abs(z_vtx_gen) < 10 and ismcprompt == 1")
    print("Opening gen file", filegen1216)
    dfgen1216 = pickle.load(open(filegen1216,'rb'))
    dfgen1216.query("abs(y_cand) < 0.5 and abs(z_vtx_gen) < 10 and ismcprompt == 1")
    print("Opening gen file", filegen1624)
    dfgen1624 = pickle.load(open(filegen1624,'rb'))
    dfgen1624.query("abs(y_cand) < 0.5 and abs(z_vtx_gen) < 10 and ismcprompt == 1")

    h_gen_pr = TH1F("h_gen_pr", "Prompt Generated in acceptance |y|<0.5", 2001, 0, 100.05)

    fill_hist(h_gen_pr, dfgen04["pt_cand"])
    fill_hist(h_gen_pr, dfgen48["pt_cand"])
    fill_hist(h_gen_pr, dfgen812["pt_cand"])
    fill_hist(h_gen_pr, dfgen1216["pt_cand"])
    fill_hist(h_gen_pr, dfgen1624["pt_cand"])

    print("Are histos filled?", h_gen_pr.GetEntries())
    print("Creating file", fout)
    fileout = TFile.Open(fout, "recreate")
    fileout.cd()
    h_gen_pr.Write()
    fileout.Close()

    filereco04 = pathgen + "AnalysisResultsReco0_4_0.00.pkl"
    filereco48 = pathgen + "AnalysisResultsReco4_8_0.00.pkl"
    filereco812 = pathgen + "AnalysisResultsReco8_12_0.00.pkl"
    filereco1216 = pathgen + "AnalysisResultsReco12_16_0.00.pkl"
    filereco1624 = pathgen + "AnalysisResultsReco16_24_0.00.pkl"

    #ITS2 = {4988, 4800, 4400, 4050, 4800}
    #ITS3 = {4840, 4800, 4850, 4800, 4550}
    if isits2:
        mlcut04 = 0.5 + 0.4988
        mlcut48 = 0.5 + 0.4800
        mlcut812 = 0.5 + 0.4400
        mlcut1216 = 0.5 + 0.4050
        mlcut1624 = 0.5 + 0.4800
    else:
        mlcut04 = 0.5 + 0.4840
        mlcut48 = 0.5 + 0.4800
        mlcut812 = 0.5 + 0.4850
        mlcut1216 = 0.5 + 0.4800
        mlcut1624 = 0.5 + 0.4550

    print("Opening gen file", filereco04)
    dfreco04 = pickle.load(open(filereco04, 'rb'))
    dfreco04 = dfreco04.query("y_test_probxgboost > @mlcut04 and ismcprompt == 1")
    print("Opening gen file", filereco48)
    dfreco48 = pickle.load(open(filereco48, 'rb'))
    dfreco48 = dfreco48.query("y_test_probxgboost > @mlcut48 and ismcprompt == 1")
    print("Opening gen file", filereco812)
    dfreco812 = pickle.load(open(filereco812, 'rb'))
    dfreco812 = dfreco812.query("y_test_probxgboost > @mlcut812 and ismcprompt == 1")
    print("Opening gen file", filereco1216)
    dfreco1216 = pickle.load(open(filereco1216, 'rb'))
    dfreco1216 = dfreco1216.query("y_test_probxgboost > @mlcut1216 and ismcprompt == 1")
    print("Opening gen file", filereco1624)
    dfreco1624 = pickle.load(open(filereco1624, 'rb'))
    dfreco1624 = dfreco1624.query("y_test_probxgboost > @mlcut1624 and ismcprompt == 1")

    h_sel_pr = TH1F("h_sel_pr", "Prompt Reco in acc |#eta|<0.8 and sel", 2001, 0, 100.05)

    fill_hist(h_sel_pr, dfreco04["pt_cand"])
    fill_hist(h_sel_pr, dfreco48["pt_cand"])
    fill_hist(h_sel_pr, dfreco812["pt_cand"])
    fill_hist(h_sel_pr, dfreco1216["pt_cand"])
    fill_hist(h_sel_pr, dfreco1624["pt_cand"])

    print("Are histos filled?", h_sel_pr.GetEntries())
    print("Creating file", fout)
    fileout = TFile.Open(fout, "update")
    fileout.cd()
    h_sel_pr.Write()
    fileout.Close()

def define_weights(filepythia, its):
    ftamu = TFile("theory/input_RAA_TAMU_Bs_B.root")
    grTAMURAABs = ftamu.Get("RAA_TAMU_Bs")

    ffonllcent = TFile("fonll/DfromB_FONLLcentPythia8_FFee_yDcut_sqrts5500.root")
    hfonllB = ffonllcent.Get("hfonllB")
    hfonllB.SetLineColor(kBlue)
    hfonllB.SetLineWidth(2)
    hfonllB.Scale(1./hfonllB.Integral())

    hfonllBTAMU = hfonllB.Clone("hfonllBTAMU")
    for i in range(hfonllB.GetNbinsX()):
        hfonllBTAMU.SetBinContent(i, hfonllB.GetBinContent(i) * grTAMURAABs.Eval(hfonllB.GetBinCenter(i)))
    hfonllBTAMU.SetLineColor(kRed+1)
    hfonllBTAMU.SetLineWidth(2)

    hflat = hfonllB.Clone("flat")
    for i in range(hflat.GetNbinsX()):
        hflat.SetBinContent(i, 0.001)
    hflat.SetLineColor(kRed+1)
    hflat.SetLineWidth(2)
    hflat.SetLineStyle(2)

    fpythia = TFile(filepythia)
    hgen = fpythia.Get("h_gen_pr")
    hpythia = hgen.Clone("hpythia")
    hpythia.SetLineColor(kBlack)
    hpythia.Scale(1./hpythia.Integral())
    hpythia.SetLineWidth(2)
    hpythia.SetTitle(";#it{p}_{T} (GeV/#it{c});d#it{N}/d#it{p}_{T} (a.u.)")

    c = TCanvas("c", "c", 450, 400)
    c.cd()
    hpythia.Draw("hist")
    gStyle.SetOptStat(0)
    hpythia.GetXaxis().SetRangeUser(0,24)
    hpythia.GetYaxis().SetRangeUser(0.00001,0.1)
    gPad.SetLogy()
    gPad.SetTickx()
    gPad.SetTicky()
    hfonllB.Draw("same hist")
    hfonllBTAMU.Draw("same hist")
    #hflat.Draw("same hist")

    leg = TLegend(.6, .69, .75, .85)
    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.034)
    leg.AddEntry(hpythia, "Pythia", "l")
    leg.AddEntry(hfonllB, "FONLL", "l")
    leg.AddEntry(hfonllBTAMU, "FONLL x TAMU", "l")
    #leg.AddEntry(hflat, "Flat", "l")
    leg.Draw()

    c.SaveAs("tempfig_" + its + "/MCpTshape_dNdpT_" + its + ".eps")

    hpythia_weight = hpythia.Clone("hpythia_weight")
    hfonllB_weight = hfonllB.Clone("hfonllB_weight")
    hfonllBTAMU_weight = hfonllBTAMU.Clone("hfonllBTAMU_weight")
    hflat_weight = hflat.Clone("hflat_weight")

    hpythia_weight.SetTitle(";#it{p}_{T} (GeV/#it{c});#it{p}_{T} weights")
    hpythia_weight.SetLineStyle(9)
    hflat_weight.SetLineStyle(2)

    hpythia_weight.Divide(hpythia_weight, hpythia, 1, 1, "B")
    hfonllB_weight.Divide(hfonllB_weight, hpythia, 1, 1, "B")
    hfonllBTAMU_weight.Divide(hfonllBTAMU_weight, hpythia, 1, 1, "B")
    hflat_weight.Divide(hflat_weight, hpythia, 1, 1, "B")

    cw = TCanvas("cw", "cw", 450, 400)
    cw.cd()
    hpythia_weight.Draw("hist")
    gStyle.SetOptStat(0)
    hpythia_weight.GetXaxis().SetRangeUser(0,23.95)
    hpythia_weight.GetYaxis().SetRangeUser(0.1,10)
    gPad.SetLogy()
    gPad.SetTickx()
    gPad.SetTicky()
    hfonllB_weight.Draw("same hist")
    hfonllBTAMU_weight.Draw("same hist")
    #hflat_weight.Draw("same hist")
    leg.Draw()

    cw.SaveAs("tempfig_" + its + "/MCpTshape_weights_" + its + ".eps")

    hreco = fpythia.Get("h_sel_pr")
    hgen = fpythia.Get("h_gen_pr")
    hreco.SetLineWidth(2)
    hgen.SetLineWidth(2)

    hgenpythia = hgen.Clone("hgenpythia")
    hgenfonll = hgen.Clone("hgenfonll")
    hgenfonlltamu = hgen.Clone("hgenfonlltamu")
    hgenflat = hgen.Clone("hgenflat")

    hrecopythia = hreco.Clone("hrecopythia")
    hrecofonll = hreco.Clone("hrecofonll")
    hrecofonlltamu = hreco.Clone("hrecofonlltamu")
    hrecoflat = hreco.Clone("hrecoflat")

    hgenpythia.SetLineColor(kBlack)
    hgenfonll.SetLineColor(kBlue)
    hgenfonlltamu.SetLineColor(kRed+1)
    hgenflat.SetLineColor(kOrange+1)
    hrecopythia.SetLineColor(kBlack)
    hrecofonll.SetLineColor(kBlue)
    hrecofonlltamu.SetLineColor(kRed+1)
    hrecoflat.SetLineColor(kOrange+1)

    hgenpythia.Multiply(hgenpythia, hpythia_weight, 1, 1, "B")
    hgenfonll.Multiply(hgenfonll, hfonllB_weight, 1, 1, "B")
    hgenfonlltamu.Multiply(hgenfonlltamu, hfonllBTAMU_weight, 1, 1, "B")
    hgenflat.Multiply(hgenflat, hflat_weight, 1, 1, "B")

    hrecopythia.Multiply(hrecopythia, hpythia_weight, 1, 1, "B")
    hrecofonll.Multiply(hrecofonll, hfonllB_weight, 1, 1, "B")
    hrecofonlltamu.Multiply(hrecofonlltamu, hfonllBTAMU_weight, 1, 1, "B")
    hrecoflat.Multiply(hrecoflat, hflat_weight, 1, 1, "B")

    arr_ptbins = [2, 4, 8, 12, 16, 24] #[0.5 + i*0.25 for i in range(20)]
    heff_cent = TH1F("heff_cent", ";#it{p}_{T} (GeV/#it{c});efficiency ratio", len(arr_ptbins) - 1, array("d", arr_ptbins))
    heff_pythia = TH1F("heff_pythia", ";#it{p}_{T} (GeV/#it{c});efficiency ratio", len(arr_ptbins) - 1, array("d", arr_ptbins))
    heff_fonll = TH1F("heff_fonll", ";#it{p}_{T} (GeV/#it{c});efficiency ratio", len(arr_ptbins) - 1, array("d", arr_ptbins))
    heff_fonlltamu = TH1F("heff_fonlltamu", ";#it{p}_{T} (GeV/#it{c});efficiency ratio", len(arr_ptbins) - 1, array("d", arr_ptbins))
    heff_flat = TH1F("heff_flat", ";#it{p}_{T} (GeV/#it{c});efficiency ratio", len(arr_ptbins) - 1, array("d", arr_ptbins))
    arr_ptbins2 = [0, 2, 4, 8, 12, 16, 24]
    href = TH1F("href", ";#it{p}_{T} (GeV/#it{c});efficiency ratio", len(arr_ptbins2) - 1, array("d", arr_ptbins2))

    for i in range(len(arr_ptbins) - 1):
        binmin = arr_ptbins[i]
        binmax = arr_ptbins[i+1]
        eff_cent = hrecopythia.Integral(hrecopythia.FindBin(binmin), hrecopythia.FindBin(binmax-0.001)) / \
                   hgenpythia.Integral(hgenpythia.FindBin(binmin), hgenpythia.FindBin(binmax-0.001))
        eff_pythia = hrecopythia.Integral(hrecopythia.FindBin(binmin), hrecopythia.FindBin(binmax-0.001)) / \
                    hgenpythia.Integral(hgenpythia.FindBin(binmin), hgenpythia.FindBin(binmax-0.001))
        eff_fonll = hrecofonll.Integral(hrecofonll.FindBin(binmin), hrecofonll.FindBin(binmax-0.001)) / \
                    hgenfonll.Integral(hgenfonll.FindBin(binmin), hgenfonll.FindBin(binmax-0.001))
        eff_fonlltamu = hrecofonlltamu.Integral(hrecofonlltamu.FindBin(binmin), hrecofonlltamu.FindBin(binmax-0.001)) / \
                        hgenfonlltamu.Integral(hgenfonlltamu.FindBin(binmin), hgenfonlltamu.FindBin(binmax-0.001))
        eff_flat = hrecoflat.Integral(hrecoflat.FindBin(binmin), hrecoflat.FindBin(binmax-0.001)) / \
                   hgenflat.Integral(hgenflat.FindBin(binmin), hgenflat.FindBin(binmax-0.001))
        heff_cent.SetBinContent(i+1, eff_cent)
        heff_pythia.SetBinContent(i + 1, eff_pythia)
        heff_fonll.SetBinContent(i + 1, eff_fonll)
        heff_fonlltamu.SetBinContent(i + 1, eff_fonlltamu)
        heff_flat.SetBinContent(i + 1, eff_flat)
        heff_cent.SetBinError(i+1, 0.0001)
        heff_pythia.SetBinError(i+1, 0.0001)
        heff_fonll.SetBinError(i+1, 0.0001)
        heff_fonlltamu.SetBinError(i+1, 0.0001)
        heff_flat.SetBinError(i+1, 0.0001)

    heff_pythia.Divide(heff_pythia, heff_cent, 1, 1, "B")
    heff_fonll.Divide(heff_fonll, heff_cent, 1, 1, "B")
    heff_fonlltamu.Divide(heff_fonlltamu, heff_cent, 1, 1, "B")
    heff_flat.Divide(heff_flat, heff_cent, 1, 1, "B")

    for i in range(len(arr_ptbins) - 1): heff_pythia.SetBinError(i, 0.0001)

    heff_pythia.SetLineColor(kBlack)
    heff_fonll.SetLineColor(kBlue)
    heff_fonlltamu.SetLineColor(kRed+1)
    heff_flat.SetLineColor(kOrange+1)
    heff_pythia.SetMarkerColor(kBlack)
    heff_fonll.SetMarkerColor(kBlue)
    heff_fonlltamu.SetMarkerColor(kRed+1)
    heff_flat.SetMarkerColor(kOrange+1)
    heff_pythia.SetLineWidth(2)
    heff_fonll.SetLineWidth(2)
    heff_fonlltamu.SetLineWidth(2)
    heff_flat.SetLineWidth(2)
    heff_flat.SetLineStyle(2)
    heff_pythia.SetMarkerStyle(20)
    heff_fonll.SetMarkerStyle(20)
    heff_fonlltamu.SetMarkerStyle(20)
    heff_flat.SetMarkerStyle(20)

    ce = TCanvas("ce", "ce", 450, 400)
    ce.cd()
    href.Draw()
    gStyle.SetOptStat(0)
    href.GetXaxis().SetRangeUser(0, 23.95)
    href.GetYaxis().SetRangeUser(0.8, 1.2)
    gPad.SetTickx()
    gPad.SetTicky()
    heff_pythia.Draw("same ep")
    heff_fonll.Draw("same ep")
    heff_fonlltamu.Draw("same ep")
    #heff_flat.Draw("same ep")
    leg.Draw()

    ce.SaveAs("tempfig_" + its + "/MCpTshape_efficiency_ratio_" + its + ".eps")

pythia_generated("/Users/lvermunt/cernbox/Analyses/ML/input/Derived/BskAnyITS2_new/vAN-20200116_ROOT6-1/ITS2_19h1b2_sequence/337_20200116-2345/skpkldecmerged_coarse/analysis/",
                 "./theory/pythia_Bs_ITS2.root", True)
define_weights("./theory/pythia_Bs_ITS2.root", "ITS2")

pythia_generated("/Users/lvermunt/cernbox/Analyses/ML/input/Derived/BskAnyITS3_new/vAN-20200116_ROOT6-1/ITS2_19h1b2_sequence/337_20200116-2345/skpkldecmerged_coarse/analysis/",
                 "./theory/pythia_Bs_ITS3.root", False)
define_weights("./theory/pythia_Bs_ITS3.root", "ITS3")

#extract_Bs_eff("/Users/lvermunt/cernbox/Analyses/ML/input/Derived/BskAnyITS2_new/vAN-20200116_ROOT6-1/ITS2_19h1b2_sequence/337_20200116-2345/skpkldecmerged_coarse/analysis/AnalysisResultsReco0_4_0.00.pkl","/Users/lvermunt/cernbox/Analyses/ML/input/Derived/BskAnyITS2_new/vAN-20200116_ROOT6-1/ITS2_19h1b2_sequence/337_20200116-2345/skpkldecmerged_coarse/analysis/AnalysisResultsGen_pt_cand0_4.pkl", 0.96, "./test.root")
