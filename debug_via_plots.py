import pickle
from hipe4ml import plot_utils
import matplotlib.pyplot as plt
from machine_learning_hep.utilities import openfile

fndata = "/mnt/temp/OngoingAnalysis_LcPbPb/Derived/LcPbPb/vAN-20210110_ROOT6-1/PbPb_data_010/566_20210117-1935/pklsk/child_1/pack_0/AnalysisResults"
fnmc = "/mnt/temp/OngoingAnalysis_LcPbPb/Derived/LcPbPb/vAN-20210201_ROOT6-1/PbPb_sim_010/573_20210202-0756/mlprob_applied/pklsk/child_1/pack_0/AnalysisResults"

binmin = [1, 2, 4, 6, 8, 12]
binmax = [2, 4, 6, 8, 12, 24]

var_evt = ["centrality", "z_vtx_reco", "n_vtx_contributors", "n_tracks", "is_ev_rej", "run_number", "ev_id", "ev_id_ext", "z_vtx_gen"]

var_reco_1 = ["ml_prob", "cos_t_star", "inv_mass", "pt_cand", "phi_cand", "eta_cand", "y_cand", "cand_type", "dca_K0s", "imp_par_K0s", "d_len_K0s", "armenteros_K0s", "ctau_K0s", "cos_p_K0s", "inv_mass_K0s", "pt_K0s", "signd0"]
var_reco_2 = ["nsigTPC_Pr_0", "nsigTOF_Pr_0", "nsigComb_Pr_0", "pt_prong0", "pt_prong1", "pt_prong2", "imp_par_prong0", "imp_par_prong1", "imp_par_prong2", "spdhits_prong0", "spdhits_prong1", "spdhits_prong2", "run_number", "ev_id", "ev_id_ext"]

var_gen = ["y_cand", "pt_cand", "eta_cand", "phi_cand", "cand_type", "run_number", "ev_id", "ev_id_ext"]

for ipt in range(6):
    print("opening ptbin", binmin[ipt], binmax[ipt])
    df_data = pickle.load(openfile(f"{fndata}Reco_pt_cand{binmin[ipt]}_{binmax[ipt]}.pkl.lz4", "rb"))
    df_mc = pickle.load(openfile(f"{fnmc}Reco_pt_cand{binmin[ipt]}_{binmax[ipt]}.pkl.lz4", "rb"))
    df_mcgen = pickle.load(openfile(f"{fnmc}Gen_pt_cand{binmin[ipt]}_{binmax[ipt]}.pkl.lz4", "rb"))

    leglabels_rec = ["Data", "Signal"]
    listdf_rec = [df_data, df_mc]

    plot_utils.plot_distr(listdf_rec, var_reco_1, 100, leglabels_rec, figsize=(12, 7), alpha=0.3, log=True, grid=True, density=True)
    plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
    figname = f'debug_plots/DistributionsReco1_pT_{binmin[ipt]}_{binmax[ipt]}.pdf'
    plt.savefig(figname)
    plt.close('all')

    plot_utils.plot_distr(listdf_rec, var_reco_2, 100, leglabels_rec, figsize=(12, 7), alpha=0.3, log=True, grid=True, density=True)
    plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
    figname = f'debug_plots/DistributionsReco2_pT_{binmin[ipt]}_{binmax[ipt]}.pdf'
    plt.savefig(figname)
    plt.close('all')

    leglabels_gen = ["Generated1", "Generated2"]
    listdf_gen = [df_mcgen, df_mcgen]

    plot_utils.plot_distr(listdf_gen, var_gen, 100, leglabels_gen, figsize=(12, 7), alpha=0.3, log=True, grid=True, density=True)
    plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
    figname = f'debug_plots/DistributionsGen_pT_{binmin[ipt]}_{binmax[ipt]}.pdf'
    plt.savefig(figname)
    plt.close('all')

    leglabels_evt = ["Data", "Signal", "Generated"]
    listdf_evt = [df_data, df_mc, df_mcgen]

    plot_utils.plot_distr(listdf_evt, var_evt, 100, leglabels_evt, figsize=(12, 7), alpha=0.3, log=True, grid=True, density=True)
    plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
    figname = f'debug_plots/DistributionsEvt_pT_{binmin[ipt]}_{binmax[ipt]}.pdf'
    plt.savefig(figname)
    plt.close('all')

