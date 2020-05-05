import pickle
import pandas
from ROOT import TFile, TH1F
from root_numpy import fill_hist

def extract_Ds_eff(filereco, filegen, fout):
    print("Opening reco file", filereco)
    dfreco = pickle.load(open(filereco,'rb'))
    print("Opening gen file", filegen)
    dfgen = pickle.load(open(filegen,'rb'))
    dfgen.query("abs(y_cand) < 0.5 and abs(z_vtx_gen) < 10")

    h_gen_pr = TH1F("h_gen_pr", "Prompt Generated in acceptance |y|<0.5", 1000, 0, 50)
    h_gen_fd = TH1F("h_gen_fd", "FD Generated in acceptance |y|<0.5", 1000, 0, 50)
    h_sel_pr = TH1F("h_presel_pr", "Prompt Reco in acc |#eta|<0.8 and sel", 1000, 0, 50)
    h_sel_fd = TH1F("h_presel_fd", "FD Reco in acc |#eta|<0.8 and sel", 1000, 0, 50)

    dfgen_pr = dfgen[dfgen.ismcprompt == 1]
    dfreco_presel_pr = dfreco[dfreco.ismcprompt == 1]
    dfgen_fd = dfgen[dfgen.ismcfd == 1]
    dfreco_presel_fd = dfreco[dfreco.ismcfd == 1]

    fill_hist(h_gen_pr, dfgen_pr["pt_cand"])
    fill_hist(h_gen_fd, dfgen_fd["pt_cand"])
    fill_hist(h_sel_pr, dfreco_presel_pr["pt_cand"])
    fill_hist(h_sel_fd, dfreco_presel_fd["pt_cand"])

    print("Are histos filled?", h_gen_pr.GetEntries(), h_gen_fd.GetEntries(), h_sel_pr.GetEntries(), h_sel_fd.GetEntries())
    print("Creating file", fout)
    feffDs = TFile.Open(fout, "recreate")
    feffDs.cd()
    h_gen_pr.Write()
    h_gen_fd.Write()
    h_sel_pr.Write()
    h_sel_fd.Write()
    feffDs.Close()

def match_Dseff_toBs_binning(path, fout):
    df04 = pickle.load(open(path + "AnalysisResultsReco_pt_cand0_4.pkl", 'rb'))
    df48 = pickle.load(open(path + "AnalysisResultsReco_pt_cand4_8.pkl", 'rb'))
    df812 = pickle.load(open(path + "AnalysisResultsReco_pt_cand8_12.pkl", 'rb'))
    df1216 = pickle.load(open(path + "AnalysisResultsReco_pt_cand12_16.pkl", 'rb'))
    df1624 = pickle.load(open(path + "AnalysisResultsReco_pt_cand16_24.pkl", 'rb'))

    pr04 = df04.query("isdsprompt == 1")["pt_Ds"]
    pr48 = df48.query("isdsprompt == 1")["pt_Ds"]
    pr812 = df812.query("isdsprompt == 1")["pt_Ds"]
    pr1216 = df1216.query("isdsprompt == 1")["pt_Ds"]
    pr1624 = df1624.query("isdsprompt == 1")["pt_Ds"]
    fd04 = df04.query("isdsprompt != 1")["pt_Ds"]
    fd48 = df48.query("isdsprompt != 1")["pt_Ds"]
    fd812 = df812.query("isdsprompt != 1")["pt_Ds"]
    fd1216 = df1216.query("isdsprompt != 1")["pt_Ds"]
    fd1624 = df1624.query("isdsprompt != 1")["pt_Ds"]

    h_match_pr = TH1F("h_match_pr", "0 < #it{p}_{T}(B_{s}) < 24 GeV/#it{c};#it{p}_{T}(pr. D_{s}) (GeV/#it{c});Counts (norm.)", 1000, 0, 50)
    h_match_pr04 = TH1F("h_match_pr04", "0 < #it{p}_{T}(B_{s}) < 4 GeV/#it{c};#it{p}_{T}(pr. D_{s}) (GeV/#it{c});Counts (norm.)", 1000, 0, 50)
    h_match_pr48 = TH1F("h_match_pr48", "4 < #it{p}_{T}(B_{s}) < 8 GeV/#it{c};#it{p}_{T}(pr. D_{s}) (GeV/#it{c});Counts (norm.)", 1000, 0, 50)
    h_match_pr812 = TH1F("h_match_pr812", "8 < #it{p}_{T}(B_{s}) < 12 GeV/#it{c};#it{p}_{T}(pr. D_{s}) (GeV/#it{c});Counts (norm.)", 1000, 0, 50)
    h_match_pr1216 = TH1F("h_match_pr1216", "12 < #it{p}_{T}(B_{s}) < 16 GeV/#it{c};#it{p}_{T}(pr. D_{s}) (GeV/#it{c});Counts (norm.)", 1000, 0, 50)
    h_match_pr1624 = TH1F("h_match_pr1624", "16 < #it{p}_{T}(B_{s}) < 24 GeV/#it{c};#it{p}_{T}(pr. D_{s}) (GeV/#it{c});Counts (norm.)", 1000, 0, 50)

    h_match_fd = TH1F("h_match_fd", "0 < #it{p}_{T}(B_{s}) < 4 GeV/#it{c};#it{p}_{T}(fd. D_{s}) (GeV/#it{c});Counts (norm.)", 1000, 0, 50)
    h_match_fd04 = TH1F("h_match_fd04", "0 < #it{p}_{T}(B_{s}) < 4 GeV/#it{c};#it{p}_{T}(fd. D_{s}) (GeV/#it{c});Counts (norm.)", 1000, 0, 50)
    h_match_fd48 = TH1F("h_match_fd48", "4 < #it{p}_{T}(B_{s}) < 8 GeV/#it{c};#it{p}_{T}(fd. D_{s}) (GeV/#it{c});Counts (norm.)", 1000, 0, 50)
    h_match_fd812 = TH1F("h_match_fd812", "8 < #it{p}_{T}(B_{s}) < 12 GeV/#it{c};#it{p}_{T}(fd. D_{s}) (GeV/#it{c});Counts (norm.)", 1000, 0, 50)
    h_match_fd1216 = TH1F("h_match_fd1216", "12 < #it{p}_{T}(B_{s}) < 16 GeV/#it{c};#it{p}_{T}(fd. D_{s}) (GeV/#it{c});Counts (norm.)", 1000, 0, 50)
    h_match_fd1624 = TH1F("h_match_fd1624", "16 < #it{p}_{T}(B_{s}) < 24 GeV/#it{c};#it{p}_{T}(fd. D_{s}) (GeV/#it{c});Counts (norm.)", 1000, 0, 50)

    fill_hist(h_match_pr, pr04)
    fill_hist(h_match_pr, pr48)
    fill_hist(h_match_pr, pr812)
    fill_hist(h_match_pr, pr1216)
    fill_hist(h_match_pr, pr1624)
    fill_hist(h_match_pr04, pr04)
    fill_hist(h_match_pr48, pr48)
    fill_hist(h_match_pr812, pr812)
    fill_hist(h_match_pr1216, pr1216)
    fill_hist(h_match_pr1624, pr1624)

    fill_hist(h_match_fd, fd04)
    fill_hist(h_match_fd, fd48)
    fill_hist(h_match_fd, fd812)
    fill_hist(h_match_fd, fd1216)
    fill_hist(h_match_fd, fd1624)
    fill_hist(h_match_fd04, fd04)
    fill_hist(h_match_fd48, fd48)
    fill_hist(h_match_fd812, fd812)
    fill_hist(h_match_fd1216, fd1216)
    fill_hist(h_match_fd1624, fd1624)

    h_match_pr04.Divide(h_match_pr)
    h_match_pr48.Divide(h_match_pr)
    h_match_pr812.Divide(h_match_pr)
    h_match_pr1216.Divide(h_match_pr)
    h_match_pr1624.Divide(h_match_pr)

    h_match_fd04.Divide(h_match_fd)
    h_match_fd48.Divide(h_match_fd)
    h_match_fd812.Divide(h_match_fd)
    h_match_fd1216.Divide(h_match_fd)
    h_match_fd1624.Divide(h_match_fd)

    print("Are histos filled?", h_match_pr04.GetEntries(), h_match_fd812.GetEntries(), h_match_pr1216.GetEntries(), h_match_fd1624.GetEntries())
    print("Creating file", fout)
    fmatcheffDs = TFile.Open(fout, "recreate")
    fmatcheffDs.cd()
    h_match_pr.Write()
    h_match_pr04.Write()
    h_match_pr48.Write()
    h_match_pr812.Write()
    h_match_pr1216.Write()
    h_match_pr1624.Write()
    h_match_fd.Write()
    h_match_fd04.Write()
    h_match_fd48.Write()
    h_match_fd812.Write()
    h_match_fd1216.Write()
    h_match_fd1624.Write()
    fmatcheffDs.Close()

match_Dseff_toBs_binning("/data/Derived/BskAnyITS2_new/vAN-20200426_ROOT6-1/ITS2_19h1a2_sequence/412_20200427-1444/pklskml_coarse/","theory/matchDspT_toBspt_ITS2_050_50MeVbins.root")
#extract_Ds_eff("/data/Derived/BskAnyITS2_new/vAN-20200426_ROOT6-1/ITS2_19h1a2_sequence/412_20200427-1444/pklDsskml_coarse/AnalysisResultsReco_pt_cand0_50.pkl", "/data/Derived/BskAnyITS2_new/vAN-20200426_ROOT6-1/ITS2_19h1a2_sequence/412_20200427-1444/pklDsskml_coarse/AnalysisResultsGen_pt_cand0_50.pkl", "theory/efficiency_Ds_050_50MeVbins.root")
