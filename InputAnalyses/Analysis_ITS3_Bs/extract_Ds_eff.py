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

    h_gen_pr = TH1F("h_gen_pr", "Prompt Generated in acceptance |y|<0.5", 200, 0, 50)
    h_gen_fd = TH1F("h_gen_fd", "FD Generated in acceptance |y|<0.5", 200, 0, 50)
    h_sel_pr = TH1F("h_presel_pr", "Prompt Reco in acc |#eta|<0.8 and sel", 200, 0, 50)
    h_sel_fd = TH1F("h_presel_fd", "FD Reco in acc |#eta|<0.8 and sel", 200, 0, 50)

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

extract_Ds_eff("/data/Derived/BskAnyITS2_new/vAN-20200426_ROOT6-1/ITS2_19h1a2_sequence/412_20200427-1444/pklDsskml_coarse/AnalysisResultsReco_pt_cand0_50.pkl", "/data/Derived/BskAnyITS2_new/vAN-20200426_ROOT6-1/ITS2_19h1a2_sequence/412_20200427-1444/pklDsskml_coarse/AnalysisResultsGen_pt_cand0_50.pkl", "theory/efficiency_Ds_050_250MeVbins.root")
