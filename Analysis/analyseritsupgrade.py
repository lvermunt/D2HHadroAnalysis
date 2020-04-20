#############################################################################
##  © Copyright CERN 2018. All rights not expressly granted are reserved.  ##
##                 Author: Gian.Michele.Innocenti@cern.ch                  ##
## This program is free software: you can redistribute it and/or modify it ##
##  under the terms of the GNU General Public License as published by the  ##
## Free Software Foundation, either version 3 of the License, or (at your  ##
## option) any later version. This program is distributed in the hope that ##
##  it will be useful, but WITHOUT ANY WARRANTY; without even the implied  ##
##     warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    ##
##           See the GNU General Public License for more details.          ##
##    You should have received a copy of the GNU General Public License    ##
##   along with this program. if not, see <https://www.gnu.org/licenses/>. ##
#############################################################################

"""
main script for doing final stage analysis
"""
# pylint: disable=too-many-lines
import os
# pylint: disable=unused-wildcard-import, wildcard-import
from array import array

# pylint: disable=import-error, no-name-in-module, unused-import
from ROOT import TFile, TH1F, TCanvas

from Analysis.analyser import Analyser
from machine_learning_hep.logger import get_logger
from machine_learning_hep.utilities_selection import getnormforselevt


class AnalyserITSUpgrade(Analyser): # pylint: disable=invalid-name
    species = "analyser"
    def __init__(self, datap, typean, period, run_param):
        super().__init__(datap, typean, period, run_param)
        self.logger = get_logger()

        self.run_param = run_param
        self.typean = typean
        self.p_period = datap["multi"]["data"]["period"][period] if period is not None \
                else "merged"

        # Analysis binning
        self.v_var_binning = datap["var_binning"]
        self.lpt_finbinmin = datap["analysis"][self.typean]["sel_an_binmin"]
        self.lpt_finbinmax = datap["analysis"][self.typean]["sel_an_binmax"]
        self.p_nptfinbins = len(self.lpt_finbinmin)
        self.lvar2_binmin = datap["analysis"][self.typean].get("sel_binmin2", None)

        #Input for probability scan
        self.probscan_min = datap["analysis"][self.typean]["probscan_min"]
        self.probscan_max = datap["analysis"][self.typean]["probscan_max"]
        self.n_probscan = datap["analysis"][self.typean]["n_probscan"]
        self.a_probscan = [[] for _ in range(self.p_nptfinbins)]

        #ML model variables
        self.p_modelname = datap["mlapplication"]["modelname"]
        self.lpt_probcutpre_mc = datap["mlapplication"]["probcutpresel"]["mc"]
        self.lpt_probcutpre_data = datap["mlapplication"]["probcutpresel"]["data"]
        self.bin_matching = datap["analysis"][self.typean]["binning_matching"]

        #input filenames
        #Build names for input pickle files (data, mc_reco, mc_gen)
        self.n_reco = datap["files_names"]["namefile_reco"]
        self.n_gen = datap["files_names"]["namefile_gen"]
        self.lpt_recodec_data = [self.n_reco.replace(".pkl", "%d_%d_%.2f.pkl" % \
                                 (self.lpt_finbinmin[i], self.lpt_finbinmax[i],
                                  self.lpt_probcutpre_data[i])) for i in range(self.p_nptfinbins)]
        self.lpt_recodec_mc = [self.n_reco.replace(".pkl", "%d_%d_%.2f.pkl" % \
                               (self.lpt_finbinmin[i], self.lpt_finbinmax[i],
                                self.lpt_probcutpre_mc[i])) for i in range(self.p_nptfinbins)]
        self.lpt_gensk = [self.n_gen.replace(".pkl", "_%s%d_%d.pkl" % \
                          (self.v_var_binning, self.lpt_finbinmin[i], self.lpt_finbinmax[i])) \
                          for i in range(self.p_nptfinbins)]

        self.d_pkl_decmerged_mc = datap["mlapplication"]["mc"]["pkl_skimmed_decmerged"][period]
        self.d_pkl_decmerged_data = datap["mlapplication"]["data"]["pkl_skimmed_decmerged"][period]
        self.lpt_recodecmerged_data = [join(self.d_pkl_decmerged_data,
                                       self.lpt_recodec_data[ipt]) for ipt in range(self.p_nptfinbins)]
        self.lpt_recodecmerged_mc = [join(self.d_pkl_decmerged_mc, self.lpt_recodec_mc[ipt])
                                     for ipt in range(self.p_nptfinbins)]
        self.lpt_gendecmerged = [os.path.join(self.d_pkl_decmerged_mc, self.lpt_gensk[ipt])
                                     for ipt in range(self.p_nptfinbins)]

        #output files
        self.d_resultsallpdata = datap["analysis"][typean]["data"]["results"][period] \
                if period is not None else datap["analysis"][typean]["data"]["resultsallp"]
        self.n_filemass_name = datap["files_names"]["histofilename"]
        self.n_fileeff_temp = "efficiencies.root"
        self.n_fileeff_name = datap["files_names"]["efffilename"]
        self.n_filemass_probscan = os.path.join(self.d_resultsallpdata, self.n_filemass_name)
        self.n_fileeff_probscan = os.path.join(self.d_resultsallpdata, self.n_fileeff_temp)
        self.n_fileeff_probscanfinal = os.path.join(self.d_resultsallpdata, self.n_fileeff_name)

        # Extra pre-selections
        self.s_evtsel = datap["analysis"][self.typean]["evtsel"]
        self.s_presel_gen_eff = datap["analysis"][self.typean]["presel_gen_eff"]
        self.s_trigger_mc = datap["analysis"][self.typean]["triggersel"]["mc"]
        self.s_trigger_data = datap["analysis"][self.typean]["triggersel"]["data"]

        #to get number of analysed events (actually not analysed, but used when merging)
        self.n_evt = data_param["files_names"]["namefile_evt"]
        self.n_evtorig = datap["files_names"]["namefile_evtorig"]
        self.d_pklevt_mergedallp_data = datap["multi"]["data"]["pkl_evtcounter_all"]
        self.f_evt_mergedallp = os.path.join(self.d_pklevt_mergedallp_data, self.n_evt)
        self.f_evtorig_mergedallp = os.path.join(self.d_pklevt_mergedallp_data, self.n_evtorig)

        #histogram settings
        self.p_mass_fit_lim = datap["analysis"][self.typean]['mass_fit_lim']
        self.p_bin_width = datap["analysis"][self.typean]['bin_width']
        self.p_num_bins = int(round((self.p_mass_fit_lim[1] - self.p_mass_fit_lim[0]) / \
                                    self.p_bin_width))


    def define_probscan_limits(self):
        """
        Sets the array with the probability cuts used for the scan
        """

        #Define these values once
        if len(self.a_probscan[0]) is not 0:
            return

        for ipt in range(self.p_nptfinbins):
            step = (self.probscan_max[ipt] - self.probscan_min[ipt])/(self.n_probscan)
            for iscan in range(self.n_probscan + 1):
                self.a_probscan[ipt].append(self.probscan_min[ipt] + iscan * step)

        self.logger.info("Array for probability scan set. Using the following ML selections:")
        for ipt in range(self.p_nptfinbins):
             print("pT bin:", ipt, ":", self.a_probscan[ipt])


    def probscan_mass_histo(self):
        """
        Invariant mass histogram for each variation, for each pT bin (2nd binning to be added)

        Similar as process_histomass_single(self, index) in processor.py
        """

        # Define limits first
        self.define_probscan_limits()

        myfile = TFile.Open(self.n_filemass_probscan, "recreate")

        self.logger.info("Doing mass histo for period", self.p_period)

        df_evt_all = pickle.load(openfile(self.f_evt_mergedallp, "rb"))
        nselevt = len(df_evt_all.query("is_ev_rej==0"))
        normevt = getnormforselevt(df_evt_all)
        hNorm = TH1F("hEvForNorm", ";;Normalisation", 2, 0.5, 2.5)
        hNorm.GetXaxis().SetBinLabel(1, "normsalisation factor")
        hNorm.GetXaxis().SetBinLabel(2, "selected events")
        hNorm.SetBinContent(1, normevt)
        hNorm.SetBinContent(2, nselevt)
        hNorm.Write()

        for ipt in range(self.p_nptfinbins):
            bin_id = self.bin_matching[ipt]
            df = pickle.load(openfile(self.lpt_recodecmerged_data[bin_id], "rb"))

            if self.s_evtsel is not None:
                df = df.query(self.s_evtsel)
            if self.s_trigger_data is not None:
                df = df.query(self.s_trigger_data)
            df = seldf_singlevar(df, self.v_var_binning,
                                 self.lpt_finbinmin[ipt], self.lpt_finbinmax[ipt])

            ptstring = "%s%d_%d" % (self.v_var_binning,
                                    self.lpt_finbinmin[ipt],
                                    self.lpt_finbinmax[ipt])
            h_isc_match = TH1F("h_match" + ptstring, ";scan iterator;ML prob.",
                               self.n_probscan + 2, -0.5, self.n_probscan + 1.5)

            for isc in range(len(self.a_probscan[ipt])):
                h_isc_match.Fill(isc, self.a_probscan[ipt][isc])
                selml_cv = "y_test_prob%s>%s" % (self.p_modelname, self.a_probscan[ipt][isc])
                df = df.query(selml_cv)

                if self.lvar2_binmin is None:
                    suffix = "%s%d_%d_%d" % (self.v_var_binning,
                                             self.lpt_finbinmin[ipt],
                                             self.lpt_finbinmax[ipt], isc)
                    h_invmass = TH1F("hmass" + suffix, "", self.p_num_bins,
                                     self.p_mass_fit_lim[0], self.p_mass_fit_lim[1])

                    fill_hist(h_invmass, df_bin.inv_mass)

                    myfile.cd()
                    h_invmass.Write()
                else:
                    self.logger.fatal("2nd binning to be implemented (code exists)")
            h_isc_match.Write()

        self.done_mass = True


    def probscan_eff(self):
        """
        Efficiency histogram for each variation, for each pT bin (2nd binning to be added)

        Similar as process_efficiency_single(self, index) in processor.py
        """

        # Define limits first
        self.define_probscan_limits()

        myfile = TFile.Open(self.n_fileeff_probscan, "recreate")

        h_gen_pr = []
        h_sel_pr = []
        h_gen_fd = []
        h_sel_fd = []

        self.logger.info("Doing eff histo for period", self.p_period)

        for ipt in range(self.p_nptfinbins):
            bin_id = self.bin_matching[ipt]

            df = pickle.load(openfile(self.lpt_recodecmerged_mc[bin_id], "rb"))
            if self.s_evtsel is not None:
                df = df.query(self.s_evtsel)
            if self.s_trigger_mc is not None:
                df = df.query(self.s_trigger_mc)
            df = seldf_singlevar(df, self.v_var_binning, self.lpt_finbinmin[ipt],
                                 self.lpt_finbinmax[ipt])

            df_gen = pickle.load(openfile(self.lpt_gendecmerged[bin_id], "rb"))
            df_gen = df_gen.query(self.s_presel_gen_eff)
            df_gen = seldf_singlevar(df_gen, self.v_var_binning, self.lpt_finbinmin[ipt],
                                     self.lpt_finbinmax[ipt])

            idx = 0
            for isc in range(len(self.a_probscan[ipt])):
                selml_cv = "y_test_prob%s>%s" % (self.p_modelname, self.a_probscan[ipt][isc])
                df = df.query(selml_cv)

                if self.lvar2_binmin is None:
                    suffix = "_%d" % (isc)
                    if ipt == 0:
                        n_bins = len(self.lpt_finbinmin)
                        analysis_bin_lims_temp = self.lpt_finbinmin.copy()
                        analysis_bin_lims_temp.append(self.lpt_finbinmax[n_bins-1])
                        analysis_bin_lims = array('f', analysis_bin_lims_temp)
                        h_gen_pr.append(TH1F("h_gen_pr" + suffix,
                                             "Prompt Generated in acceptance |y|<0.5",
                                             n_bins, analysis_bin_lims))
                        h_sel_pr.append(TH1F("h_sel_pr" + suffix,
                                             "Prompt Reco and sel in acc |#eta|<0.8 and sel",
                                             n_bins, analysis_bin_lims))
                        h_gen_fd.append(TH1F("h_gen_fd" + suffix,
                                             "FD Generated in acceptance |y|<0.5",
                                             n_bins, analysis_bin_lims))
                        h_sel_fd.append(TH1F("h_sel_fd" + suffix,
                                             "FD Reco and sel in acc |#eta|<0.8 and sel",
                                             n_bins, analysis_bin_lims))

                    df_sel_pr = df[df.ismcprompt == 1]
                    df_gen_pr = df_gen[df_gen.ismcprompt == 1]
                    df_sel_fd = df[df.ismcfd == 1]
                    df_gen_fd = df_gen[df_gen.ismcfd == 1]

                    h_gen_pr[idx].SetBinContent(ipt + 1, len(df_gen_pr))
                    h_gen_pr[idx].SetBinError(ipt + 1, math.sqrt(len(df_gen_pr)))
                    h_sel_pr[idx].SetBinContent(ipt + 1, len(df_sel_pr))
                    h_sel_pr[idx].SetBinError(ipt + 1, math.sqrt(len(df_sel_pr)))

                    h_gen_fd[idx].SetBinContent(ipt + 1, len(df_gen_fd))
                    h_gen_fd[idx].SetBinError(ipt + 1, math.sqrt(len(df_gen_fd)))
                    h_sel_fd[idx].SetBinContent(ipt + 1, len(df_sel_fd))
                    h_sel_fd[idx].SetBinError(ipt + 1, math.sqrt(len(df_sel_fd)))
                    idx = idx + 1

                else:
                    self.logger.fatal("2nd binning to be implemented (code exists)")

        myfile.cd()
        for i in range(idx):
            h_gen_pr[i].Write()
            h_sel_pr[i].Write()
            h_gen_fd[i].Write()
            h_sel_fd[i].Write()

        self.done_eff = True


    def probscan_eff_histo(self):
        """
        Efficiency histogram for each variation, for each pT bin (2nd binning to be added)

        Similar as efficiency(self) in analyzer.py
        """
        lfileeff = TFile.Open(self.n_fileeff_probscan, "READ")
        fileout = TFile(self.n_fileeff_probscanfinal, "RECREATE")

        for icv in range(len(self.a_probscan[ipt])):

            if self.lvar2_binmin is None:
                suffix = "_%d" % (icv)

                h_gen_pr = lfileeff.Get("h_gen_pr" + suffix)
                h_sel_pr = lfileeff.Get("h_sel_pr" + suffix)
                h_sel_pr.Divide(h_sel_pr, h_gen_pr, 1.0, 1.0, "B")

                h_gen_fd = lfileeff.Get("h_gen_fd" + suffix)
                h_sel_fd = lfileeff.Get("h_sel_fd" + suffix)
                h_sel_fd.Divide(h_sel_fd, h_gen_fd, 1.0, 1.0, "B")

                fileout.cd()
                h_sel_pr.SetName("eff_mult" + suffix)
                h_sel_fd.SetName("eff_fd_mult" + suffix)
                h_sel_pr.Write()
                h_sel_fd.Write()
            else:
                self.logger.fatal("2nd binning to be implemented (code exists)")

        fileout.Close()
