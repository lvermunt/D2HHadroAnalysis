import pickle
import pandas
from ROOT import TFile, TH1F
from root_numpy import fill_hist

def match_dspt_to_bspt(path, outpath):

    print("Opening files from:", path)
    print("Saving histograms in:", outpath)
    fout = TFile.Open(outpath, "recreate")

    df04 = pickle.load(open(path + "/skpkldecmerged_coarse/analysis/AnalysisResultsReco0_4_0.00.pkl", 'rb'))
    df48 = pickle.load(open(path + "/skpkldecmerged_coarse/analysis/AnalysisResultsReco4_8_0.00.pkl", 'rb'))
    df812 = pickle.load(open(path + "/skpkldecmerged_coarse/analysis/AnalysisResultsReco8_12_0.00.pkl", 'rb'))
    df1216 = pickle.load(open(path + "/skpkldecmerged_coarse/analysis/AnalysisResultsReco12_16_0.00.pkl", 'rb'))
    df1624 = pickle.load(open(path + "/skpkldecmerged_coarse/analysis/AnalysisResultsReco16_24_0.00.pkl", 'rb'))

    frames = [df04, df48, df812, df1216, df1624]
    dfbs = pandas.concat(frames)
    dfbs_pr = dfbs.query("isdsprompt == 1")
    dfbs_fd = dfbs.query("isdsprompt != 1")

    dfds10 = pickle.load(open(path + "/pklDs/child_1/pack_0/AnalysisResultsReco.pkl", 'rb'))
    dfds11 = pickle.load(open(path + "/pklDs/child_1/pack_1/AnalysisResultsReco.pkl", 'rb'))
    dfds20 = pickle.load(open(path + "/pklDs/child_2/pack_0/AnalysisResultsReco.pkl", 'rb'))
    dfds30 = pickle.load(open(path + "/pklDs/child_3/pack_0/AnalysisResultsReco.pkl", 'rb'))
    dfds40 = pickle.load(open(path + "/pklDs/child_4/pack_0/AnalysisResultsReco.pkl", 'rb'))
    dfds41 = pickle.load(open(path + "/pklDs/child_4/pack_1/AnalysisResultsReco.pkl", 'rb'))

    framesds = [dfds10, dfds11, dfds20, dfds30, dfds40, dfds41]
    dfds = pandas.concat(framesds)
    dfds_pr = dfds.query("ismcprompt == 1")
    dfds_fd = dfds.query("ismcfd == 1")

    nbins = 1000 #from pT=0 till pT=50
    binwidth = 0.05 #50 MeV
    for i in range(nbins):
        ptmin = i * binwidth
        ptmax = (i + 1) * binwidth
        print("Checking pTDs: ", ptmin, "-", ptmax)
        dfbs_pr_cutds = dfbs_pr.query("pt_Ds >= @ptmin and pt_Ds < @ptmax")
        dfbs_fd_cutds = dfbs_fd.query("pt_Ds >= @ptmin and pt_Ds < @ptmax")

        dfds_pr_cutds = dfds_pr.query("pt_cand >= @ptmin and pt_cand < @ptmax")
        dfds_fd_cutds = dfds_fd.query("pt_cand >= @ptmin and pt_cand < @ptmax")

        h_bs_pr = TH1F("h_bs_pr_%d" % i, "B_{s} in TTree for %.3f < p_{T}(pr. D_{s}) < %.3f" % (ptmin, ptmax), 100, 0, 50)
        h_bs_fd = TH1F("h_bs_fd_%d" % i, "B_{s} in TTree for %.3f < p_{T}(fd. D_{s}) < %.3f" % (ptmin, ptmax), 100, 0, 50)

        fill_hist(h_bs_pr, dfbs_pr_cutds["pt_cand"])
        fill_hist(h_bs_fd, dfbs_fd_cutds["pt_cand"])
        if len(dfds_pr_cutds) > 0: h_bs_pr.Scale(1. / len(dfds_pr_cutds))
        if len(dfds_fd_cutds) > 0: h_bs_fd.Scale(1. / len(dfds_fd_cutds))

        fout.cd()
        h_bs_pr.Write()
        h_bs_fd.Write()
    fout.Close()

match_dspt_to_bspt("/data/Derived/BskAnyITS2_new/vAN-20200426_ROOT6-1/ITS2_19h1a2_sequence/412_20200427-1444", "/home/lvermunt/D2HHadroAnalysis/InputAnalyses/Analysis_ITS3_Bs/BspT_for_50MeVDsbins_ITS2_Tree.root")