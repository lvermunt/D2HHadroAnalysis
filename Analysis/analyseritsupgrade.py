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
import pickle
import math
# pylint: disable=unused-wildcard-import, wildcard-import
from array import array
from numpy import isnan, average
from root_numpy import fill_hist  # pylint: disable=import-error, no-name-in-module

# pylint: disable=import-error, no-name-in-module, unused-import
from ROOT import TFile, TH1F, TCanvas, TF1, TGraph, gPad, TGaxis, gStyle
from ROOT import kBlack, kRed, kGreen, kBlue, kOrange, kViolet, kAzure

from Analysis.analyser import Analyser
from machine_learning_hep.logger import get_logger
from machine_learning_hep.utilities import openfile
from machine_learning_hep.utilities_selection import getnormforselevt, seldf_singlevar


class AnalyserITSUpgrade(Analyser): # pylint: disable=invalid-name
    species = "analyser"
    def __init__(self, datap, case, typean, period):
        super().__init__(datap, case, typean, period)
        self.logger = get_logger()

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
        self.probscan_prefilter = datap["analysis"][self.typean].get("probscan_prefilter", None)
        self.probscan_min = datap["analysis"][self.typean]["probscan_min"]
        self.probscan_max = datap["analysis"][self.typean]["probscan_max"]
        self.n_probscan = datap["analysis"][self.typean]["n_probscan"]
        self.a_probscan = [[] for _ in range(self.p_nptfinbins)]
        self.a_probscan_bin = [[] for _ in range(self.p_nptfinbins)]
        self.fpar1func = datap["analysis"][self.typean].get("probscan_func_par1", "pol2")
        self.n_bkgentrfits = datap["analysis"][self.typean].get("probscan_nfits", 8)
        self.v_fitsplit = datap["analysis"][self.typean].get("probscan_length_fits",
                                                             [0.005 for _ in range(self.p_nptfinbins)])
        self.minbc = datap["analysis"][self.typean].get("probscan_min_binc", 100)

        #ML model variables
        self.p_modelname = datap["mlapplication"]["modelname"]
        self.lpt_probcutpre_mc = datap["mlapplication"]["probcutpresel"]["mc"]
        self.lpt_probcutpre_data = datap["mlapplication"]["probcutpresel"]["data"]
        self.bin_matching = datap["analysis"][self.typean]["binning_matching"]
        self.lpt_skimbinmin = datap["sel_skim_binmin"]
        self.lpt_skimbinmax = datap["sel_skim_binmax"]
        self.p_nptskimbins = len(self.lpt_skimbinmin)

        #input filenames
        #Build names for input pickle files (data, mc_reco, mc_gen)
        self.n_reco = datap["files_names"]["namefile_reco"]
        self.n_gen = datap["files_names"]["namefile_gen"]
        self.lpt_recodec_data = [self.n_reco.replace(".pkl", "%d_%d_%.2f.pkl" % \
                                 (self.lpt_skimbinmin[i], self.lpt_skimbinmax[i],
                                  self.lpt_probcutpre_data[i])) for i in range(self.p_nptskimbins)]
        self.lpt_recodec_mc = [self.n_reco.replace(".pkl", "%d_%d_%.2f.pkl" % \
                               (self.lpt_skimbinmin[i], self.lpt_skimbinmax[i],
                                self.lpt_probcutpre_mc[i])) for i in range(self.p_nptskimbins)]
        self.lpt_gensk = [self.n_gen.replace(".pkl", "_%s%d_%d.pkl" % \
                          (self.v_var_binning, self.lpt_skimbinmin[i], self.lpt_skimbinmax[i])) \
                          for i in range(self.p_nptskimbins)]

        #FIXME
        self.d_pkl_decmerged_mc = datap["mlapplication"]["mc"]["pkl_skimmed_decmerged"][period] + "/analysis" \
            if period is not None else "./"
        self.d_pkl_decmerged_data = datap["mlapplication"]["data"]["pkl_skimmed_decmerged"][period] + "/analysis" \
            if period is not None else "./"
        self.lpt_recodecmerged_data = [os.path.join(self.d_pkl_decmerged_data, self.lpt_recodec_data[ipt])
                                       for ipt in range(self.p_nptskimbins)]
        self.lpt_recodecmerged_mc = [os.path.join(self.d_pkl_decmerged_mc, self.lpt_recodec_mc[ipt])
                                     for ipt in range(self.p_nptskimbins)]
        self.lpt_gendecmerged = [os.path.join(self.d_pkl_decmerged_mc, self.lpt_gensk[ipt])
                                 for ipt in range(self.p_nptskimbins)]

        #output files
        self.d_resultsallpdata = datap["analysis"][self.typean]["data"]["resultsallp"]
        self.n_filemass_name = datap["files_names"]["histofilename"]
        self.n_filemass_namefit = "backgroundparameterfit.root"
        self.n_filepars_probscanfit = "parametrised_bkgpars.root"
        self.n_fileeff_temp = "efficiencies.root"
        self.n_fileeff_name = datap["files_names"]["efffilename"]
        self.n_filemass_probscan = os.path.join(self.d_resultsallpdata, self.n_filemass_name)
        self.n_filemass_bkgshape = self.n_filemass_probscan.replace(".root", "_bkgshape.root")
        self.n_filemass_bkgshapescaled = self.n_filemass_probscan.replace(".root", "_bkgshape_scaled.root")
        self.n_filemass_probscanfit = os.path.join(self.d_resultsallpdata, self.n_filemass_namefit)
        self.n_file_params = os.path.join(self.d_resultsallpdata, self.n_filepars_probscanfit)
        self.n_fileeff_probscan = os.path.join(self.d_resultsallpdata, self.n_fileeff_temp)
        self.n_fileeff_probscanfinal = os.path.join(self.d_resultsallpdata, self.n_fileeff_name)

        # Extra pre-selections
        self.s_evtsel = datap["analysis"][self.typean]["evtsel"]
        self.s_presel_gen_eff = datap["analysis"][self.typean]["presel_gen_eff"]
        self.s_trigger_mc = datap["analysis"][self.typean]["triggersel"]["mc"]
        self.s_trigger_data = datap["analysis"][self.typean]["triggersel"]["data"]

        #Names for bitmap selections
        self.v_isstd = datap["bitmap_sel"]["var_isstd"]
        self.v_ismcsignal = datap["bitmap_sel"]["var_ismcsignal"]
        self.v_ismcprompt = datap["bitmap_sel"]["var_ismcprompt"]
        self.v_ismcfd = datap["bitmap_sel"]["var_ismcfd"]
        self.v_ismcbkg = datap["bitmap_sel"]["var_ismcbkg"]
        self.v_ismcrefl = datap["bitmap_sel"]["var_ismcrefl"]
        self.v_dsprompt = datap["bitmap_sel"].get("var_isdsprompt", None)
        self.v_dsfdbplus = datap["bitmap_sel"].get("var_isdsfdbplus", None)
        self.v_dsfdbzero = datap["bitmap_sel"].get("var_isdsfdbzero", None)
        self.v_dsfdlambdab = datap["bitmap_sel"].get("var_isdsfdlambdab", None)
        self.v_dsfdbs = datap["bitmap_sel"].get("var_isdsfdbs", None)

        #to get number of analysed events (actually not analysed, but used when merging)
        self.n_evt = datap["files_names"]["namefile_evt"]
        self.n_evtorig = datap["files_names"]["namefile_evtorig"]
        self.d_pklevt_mergedallp_data = datap["multi"]["data"]["pkl_evtcounter_all"]
        self.f_evt_mergedallp = os.path.join(self.d_pklevt_mergedallp_data, self.n_evt)
        self.f_evtorig_mergedallp = os.path.join(self.d_pklevt_mergedallp_data, self.n_evtorig)

        #histogram settings
        self.p_mass_fit_lim = datap["analysis"][self.typean]['mass_fit_lim']
        self.p_bin_width = datap["analysis"][self.typean]['bin_width']
        self.p_num_bins = int(round((self.p_mass_fit_lim[1] - self.p_mass_fit_lim[0]) / \
                                    self.p_bin_width))
        self.p_bkgfunc = datap["analysis"][self.typean]["bkgfunc"]
        self.bkg_fmap = {"kExpo": "expo", "kLin": "pol1", "Pol2": "pol2"}
        self.rebins = datap["analysis"][self.typean]["rebin"]
        self.binwidth = [None for _ in range(self.p_nptfinbins)]

        #FIXME, move to a database
        self.n_filebkgcorrection = "/Users/lvermunt/PycharmProjects/D2HHadroAnalysis/InputAnalyses/Analysis_ITS3_Bs/theory/BkgCorrFactor_Bs_1DataFile_25MCFile.root"
        self.n_filefonllinput = "/Users/lvermunt/PycharmProjects/D2HHadroAnalysis/InputAnalyses/Analysis_ITS3_Bs/theory/inputFONLL_02468121624.root"
        self.n_filetamuinput = "/Users/lvermunt/PycharmProjects/D2HHadroAnalysis/InputAnalyses/Analysis_ITS3_Bs/theory/input_RAA_TAMU_Bs_02468121624.root"
        self.v_ninterestingtrials = 500
        self.brbs = 0.00304
        self.brbs_err = 0.00023
        self.brds = 0.0227
        self.brds_err = 0.0008
        self.nevexpected = 8000000000
        self.nevanalysed = 852644
        self.taa = 0.02307
        self.gauss3sigmafac = 0.9973

    def define_probscan_limits(self):
        """
        Sets the array with the probability cuts used for the scan
        """

        #Define these values once
        if len(self.a_probscan[0]) is not 0:
            return

        for ipt in range(self.p_nptfinbins):
            step = (self.probscan_max[ipt] - self.probscan_min[ipt])/(self.n_probscan)
            self.a_probscan_bin[ipt].append(self.probscan_min[ipt] - 0.5 * step)
            for iscan in range(self.n_probscan + 1):
                self.a_probscan[ipt].append(self.probscan_min[ipt] + iscan * step)
                self.a_probscan_bin[ipt].append(self.probscan_min[ipt] + (iscan + 0.5) * step)

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

        self.logger.info("Doing mass histo for period %s", self.p_period)

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
                    h_invmass_sig = TH1F("hmass_sig" + suffix, "", self.p_num_bins,
                                         self.p_mass_fit_lim[0], self.p_mass_fit_lim[1])
                    h_invmass_bkg = TH1F("hmass_bkg" + suffix, "", self.p_num_bins,
                                         self.p_mass_fit_lim[0], self.p_mass_fit_lim[1])

                    df_sig = df[df[self.v_ismcbkg] == 0]
                    df_bkg = df[df[self.v_ismcbkg] == 1]

                    fill_hist(h_invmass, df.inv_mass)
                    fill_hist(h_invmass_sig, df_sig.inv_mass)
                    fill_hist(h_invmass_bkg, df_bkg.inv_mass)

                    myfile.cd()
                    h_invmass.Write()
                    h_invmass_sig.Write()
                    h_invmass_bkg.Write()
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

        self.logger.info("Doing eff histo for period %s", self.p_period)

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
            df_gen_pr = df_gen[df_gen.ismcprompt == 1]
            df_gen_fd = df_gen[df_gen.ismcfd == 1]

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
                    df_sel_fd = df[df.ismcfd == 1]

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

        # Define limits first
        self.define_probscan_limits()

        lfileeff = TFile.Open(self.n_fileeff_probscan, "READ")
        fileout = TFile(self.n_fileeff_probscanfinal, "RECREATE")

        self.logger.info("Doing final eff histo for period %s", self.p_period)

        for icv in range(len(self.a_probscan[0])):

            if self.lvar2_binmin is None:
                suffix = "_%d" % (icv)

                h_gen_pr = lfileeff.Get("h_gen_pr" + suffix)
                h_sel_pr = lfileeff.Get("h_sel_pr" + suffix)
                h_sel_pr.Divide(h_sel_pr, h_gen_pr, 1.0, 1.0, "B")

                h_gen_fd = lfileeff.Get("h_gen_fd" + suffix)
                h_sel_fd = lfileeff.Get("h_sel_fd" + suffix)
                h_sel_fd.Divide(h_sel_fd, h_gen_fd, 1.0, 1.0, "B")

                fileout.cd()
                h_sel_pr.SetName("eff" + suffix)
                h_sel_fd.SetName("eff_fd" + suffix)
                h_sel_pr.Write()
                h_sel_fd.Write()
            else:
                self.logger.fatal("2nd binning to be implemented (code exists)")

        fileout.Close()


    def probability_scan(self):
        """
        Function to run full probability scan in one go
        """

        if self.p_period is "merged":
            self.logger.warning("Invalid option for ITSUpgrade analyser, skipping")
            return
        self.define_probscan_limits()
        self.probscan_mass_histo()
        self.probscan_eff()
        self.probscan_eff_histo()


    def fit_invmassbkg_scan(self):
        """
        Fit background with TF1, and parametrise the background parameters
        """

        if self.p_period is "merged":
            self.logger.warning("Invalid option for ITSUpgrade analyser, skipping")
            return

        # Define limits first
        self.define_probscan_limits()
        self.logger.info("\n\nArray for probability histo binning set:")
        for ipt in range(self.p_nptfinbins):
             print("pT bin:", ipt, ":", self.a_probscan_bin[ipt])

        fmass = TFile.Open(self.n_filemass_probscan, "READ")
        fileout = TFile(self.n_filemass_probscanfit, "RECREATE")

        self.logger.info("Fitting background mass spectra for %s", self.p_period)
        hpar0 = []
        hpar1 = []
        hpar2 = []
        hentr = []
        for ipt in range(self.p_nptfinbins):

            hpar0.append(TH1F("hpar0_%d" % ipt,
                              "%d < #it{p}_{T} < %d GeV/#it{c};ML probability;Background par0" %
                              (self.lpt_finbinmin[ipt], self.lpt_finbinmax[ipt]),
                              self.n_probscan + 1, array("d", self.a_probscan_bin[ipt])))
            hpar1.append(TH1F("hpar1_%d" % ipt,
                              "%d < #it{p}_{T} < %d GeV/#it{c};ML probability;Background par1" %
                              (self.lpt_finbinmin[ipt], self.lpt_finbinmax[ipt]),
                              self.n_probscan + 1, array("d", self.a_probscan_bin[ipt])))
            hpar2.append(TH1F("hpar2_%d" % ipt,
                              "%d < #it{p}_{T} < %d GeV/#it{c};ML probability;Background par2" %
                              (self.lpt_finbinmin[ipt], self.lpt_finbinmax[ipt]),
                              self.n_probscan + 1, array("d", self.a_probscan_bin[ipt])))
            hentr.append(TH1F("hentr_%d" % ipt,
                              "%d < #it{p}_{T} < %d GeV/#it{c};ML probability;Background Entries" %
                              (self.lpt_finbinmin[ipt], self.lpt_finbinmax[ipt]),
                              self.n_probscan + 1, array("d", self.a_probscan_bin[ipt])))

            hpar0[ipt].SetStats(0)
            hpar1[ipt].SetStats(0)
            hpar2[ipt].SetStats(0)
            hentr[ipt].SetStats(0)

            for isc in range(len(self.a_probscan[ipt])):
                suffix = "%s%d_%d_%d" % (self.v_var_binning,
                                         self.lpt_finbinmin[ipt],
                                         self.lpt_finbinmax[ipt], isc)
                hbkg = fmass.Get("hmass_bkg" + suffix)
                fbkg = TF1("fbkg" + suffix, self.bkg_fmap[self.p_bkgfunc[ipt]],
                           self.p_mass_fit_lim[0], self.p_mass_fit_lim[1])

                hbkg.Rebin(self.rebins[ipt])
                self.binwidth[ipt] = hbkg.GetBinWidth(1)
                hbkg.Fit("fbkg" + suffix, "R,E,+,0")

                bkgentries = hbkg.GetEntries()
                bkgpar0 = fbkg.GetParameter(0)
                bkgparerr0 = fbkg.GetParError(0)
                bkgpar1 = fbkg.GetParameter(1)
                bkgparerr1 = fbkg.GetParError(1)
                if fbkg.GetNpar() == 3:
                    bkgpar2 = fbkg.GetParameter(2)
                    bkgparerr2 = fbkg.GetParError(2)
                if fbkg.GetNpar() > 3:
                    self.logger.warning("N parameters > 3 not supported!")

                hentr[ipt].SetBinContent(hentr[ipt].FindBin(self.a_probscan[ipt][isc]), bkgentries)
                hentr[ipt].SetBinError(hentr[ipt].FindBin(self.a_probscan[ipt][isc]), math.sqrt(bkgentries))
                if not isnan(bkgpar0) and not isnan(bkgpar1):
                    hpar0[ipt].SetBinContent(hpar0[ipt].FindBin(self.a_probscan[ipt][isc]), bkgpar0)
                    hpar0[ipt].SetBinError(hpar0[ipt].FindBin(self.a_probscan[ipt][isc]), bkgparerr0)
                    hpar1[ipt].SetBinContent(hpar1[ipt].FindBin(self.a_probscan[ipt][isc]), bkgpar1)
                    hpar1[ipt].SetBinError(hpar1[ipt].FindBin(self.a_probscan[ipt][isc]), bkgparerr1)
                    if fbkg.GetNpar() == 3:
                        hpar2[ipt].SetBinContent(hpar2[ipt].FindBin(self.a_probscan[ipt][isc]), bkgpar2)
                        hpar2[ipt].SetBinError(hpar2[ipt].FindBin(self.a_probscan[ipt][isc]), bkgparerr2)
                print("Fitting pT bin", ipt, "for iscan", isc)

        fileout.cd()
        for ipt in range(self.p_nptfinbins):
            hentr[ipt].Write()
            hpar0[ipt].Write()
            hpar1[ipt].Write()
            hpar2[ipt].Write()

        self.logger.info("\n\nBackground fitting finished, saved in %s\n\n", self.n_filemass_probscanfit)


    def fit_bkgparams_2_scan(self):
        """
        Function to fit+parametrise the bkg fit parameters versus ML probability

        FIXME: Make less chaotic (although it does the job)...
        """

        if self.p_period is "merged":
            self.logger.warning("Invalid option for ITSUpgrade analyser, skipping")
            return

        # Define limits first
        self.define_probscan_limits()

        fpars = TFile(self.n_filemass_probscanfit, "READ")
        fileout = TFile(self.n_file_params, "RECREATE")

        self.logger.info("Fitting background mass spectra for %s", self.p_period)
        TGaxis.SetMaxDigits(3)
        gStyle.SetOptStat(0000)

        fpar1 = []
        fpar1low = []
        fpar1high = []
        gpar1cent = []
        gpar1low = []
        gpar1high = []
        gpar0cent = []
        gpar0low = []
        gpar0high = []
        for ipt in range(self.p_nptfinbins):
            canv = TCanvas("c_%d" % ipt, "c_%d" % ipt, 1350, 400)
            canv.Divide(3)

            hbinwidth = TH1F("hbinwidth", "", self.p_num_bins,
                             self.p_mass_fit_lim[0], self.p_mass_fit_lim[1])
            hbinwidth.Rebin(self.rebins[ipt])
            self.binwidth[ipt] = hbinwidth.GetBinWidth(1)

            hpar0 = fpars.Get("hpar0_%d" % ipt)
            hpar1 = fpars.Get("hpar1_%d" % ipt)
            hentr = fpars.Get("hentr_%d" % ipt)

            canv.cd(3)
            hpar0.Draw("ep")
            gPad.SetTickx()
            gPad.SetTicky()
            hpar0.GetXaxis().SetRangeUser(self.probscan_min[ipt], self.probscan_max[ipt])

            canv.cd(2)
            hpar1.Draw("ep")
            fpar1.append(TF1("fpar1_%d" % ipt, self.fpar1func[ipt], self.probscan_min[ipt],
                             self.probscan_max[ipt]))
            hpar1.Fit("fpar1_%d" % ipt, "R,+")
            gPad.SetTickx()
            gPad.SetTicky()
            hpar1.GetXaxis().SetRangeUser(self.probscan_min[ipt], self.probscan_max[ipt])
            hpar1.GetYaxis().SetRangeUser(-2, 2)

            fpar1low.append(fpar1[ipt].Clone("%slow" % fpar1[ipt].GetName()))
            fpar1high.append(fpar1[ipt].Clone("%shigh" % fpar1[ipt].GetName()))
            fpar1low[ipt].SetLineStyle(2)
            fpar1high[ipt].SetLineStyle(2)
            for ipar in range(fpar1[ipt].GetNpar()):
                if fpar1[ipt].GetParameter(ipar) < 0:
                    fpar1low[ipt].SetParameter(ipar, fpar1[ipt].GetParameter(ipar) -
                                               fpar1[ipt].GetParError(ipar))
                    fpar1high[ipt].SetParameter(ipar, fpar1[ipt].GetParameter(ipar) +
                                                fpar1[ipt].GetParError(ipar))
                else:
                    fpar1low[ipt].SetParameter(ipar, fpar1[ipt].GetParameter(ipar) +
                                               fpar1[ipt].GetParError(ipar))
                    fpar1high[ipt].SetParameter(ipar, fpar1[ipt].GetParameter(ipar) -
                                                fpar1[ipt].GetParError(ipar))
            fpar1low[ipt].Draw("same")
            fpar1high[ipt].Draw("same")

            canv.cd(1)
            hentr.Draw("ep")
            gPad.SetTickx()
            gPad.SetTicky()
            hentr.GetXaxis().SetRangeUser(self.probscan_min[ipt], self.probscan_max[ipt])
            fentrall = TF1("fentr_%d" % ipt, "pol2", self.probscan_min[ipt],
                           self.probscan_max[ipt])
            fentrall.SetLineStyle(4)
            hentr.Fit("fentr_%d" % ipt, "R,+")

            fentr = []
            fentrfull = []
            fentrlow = []
            fentrhigh = []
            ranges_k = []
            minprob_forfit = self.probscan_max[ipt] - self.n_bkgentrfits * self.v_fitsplit[ipt]
            ranges_k.append(minprob_forfit)
            for k in range(self.n_bkgentrfits):
                fentr.append(TF1("fentr_%d_%d" % (ipt, k), "pol1",
                                 minprob_forfit + k * self.v_fitsplit[ipt],
                                 minprob_forfit + (k + 1) * self.v_fitsplit[ipt]))
                ranges_k.append(minprob_forfit + (k + 1) * self.v_fitsplit[ipt])
                fentr[k].SetLineColor(k+1)
                hentr.Fit("fentr_%d_%d" % (ipt, k), "R,+,0")

                #shift ranges a bit for when fit fails
                if fentr[k].GetParameter(0) == 0 and fentr[k].GetParameter(1) == 0:
                    fentr[k] = TF1("fentr_%d_%d" % (ipt, k), "pol1",
                                   minprob_forfit + (k - 0.1) * self.v_fitsplit[ipt],
                                   minprob_forfit + (k + 0.9) * self.v_fitsplit[ipt])
                    ranges_k[k+1] = minprob_forfit + (k + 0.9) * self.v_fitsplit[ipt]
                    hentr.Fit("fentr_%d_%d" % (ipt, k), "R,+,0")

                    if fentr[k].GetParameter(0) == 0 and fentr[k].GetParameter(1) == 0:
                        fentr[k] = TF1("fentr_%d_%d" % (ipt, k), "pol1",
                                       minprob_forfit + (k - 0.2) * self.v_fitsplit[ipt],
                                       minprob_forfit + (k + 0.8) * self.v_fitsplit[ipt])
                        ranges_k[k+1] = minprob_forfit + (k + 0.8) * self.v_fitsplit[ipt]
                        hentr.Fit("fentr_%d_%d" % (ipt, k), "R,+,0")
                fentr[k].Draw("same")

                fentrfull.append(TF1("%sfull" % fentr[k].GetName(), "pol1",
                                     self.probscan_min[ipt], self.probscan_max[ipt]))
                fentrlow.append(fentr[k].Clone("%slow" % fentr[k].GetName()))
                fentrhigh.append(fentr[k].Clone("%shigh" % fentr[k].GetName()))
                fentrfull[k].SetLineColor(k+1)
                fentrfull[k].SetLineStyle(2)
                fentrlow[k].SetLineStyle(2)
                fentrhigh[k].SetLineStyle(2)
                for ipar in range(fentr[k].GetNpar()):
                    fentrfull[k].SetParameter(ipar, fentr[k].GetParameter(ipar))
                    if fentr[k].GetParameter(ipar) < 0:
                        fentrlow[k].SetParameter(ipar, fentr[k].GetParameter(ipar) -
                                                 fentr[k].GetParError(ipar))
                        fentrhigh[k].SetParameter(ipar, fentr[k].GetParameter(ipar) +
                                                  fentr[k].GetParError(ipar))
                    else:
                        fentrlow[k].SetParameter(ipar, fentr[k].GetParameter(ipar) +
                                                 fentr[k].GetParError(ipar))
                        fentrhigh[k].SetParameter(ipar, fentr[k].GetParameter(ipar) -
                                                  fentr[k].GetParError(ipar))
                #fentrlow[k].Draw("same")
                #fentrhigh[k].Draw("same")
                fentrfull[k].Draw("same")

            gpar0cent.append(TGraph(0))
            gpar0low.append(TGraph(0))
            gpar0high.append(TGraph(0))
            gpar1cent.append(TGraph(0))
            gpar1low.append(TGraph(0))
            gpar1high.append(TGraph(0))
            minbc = 9999.
            minfit = 9999.
            maxfit = -1.
            for isc in range(self.n_probscan):
                mlcut = self.a_probscan[ipt][isc]
                k = -1
                for it_k in range(self.n_bkgentrfits):
                    if ranges_k[it_k] < mlcut:
                        k = k + 1

                fbc = hentr.GetBinContent(hentr.FindBin(mlcut))
                if fbc < minbc and fbc > 0: minbc = fbc
                if fbc <= 0: fbc = minbc

                v_fpar1 = fpar1[ipt].Eval(mlcut)
                v_fpar1low = fpar1low[ipt].Eval(mlcut)
                v_fpar1high = fpar1high[ipt].Eval(mlcut)
                par0new = self.calculate_par0(ipt, v_fpar1, fbc, self.bkg_fmap[self.p_bkgfunc[ipt]])
                if self.bkg_fmap[self.p_bkgfunc[ipt]] is "pol1":
                    checkneg = par0new + fpar1[ipt].Eval(mlcut) * self.p_mass_fit_lim[0]
                    if checkneg < 0:
                        v_fpar1 = 0
                        v_fpar1low = -0.1
                        v_fpar1high = 0.1
                        par0new = self.calculate_par0(ipt, v_fpar1, fbc, self.bkg_fmap[self.p_bkgfunc[ipt]])

                #Add here the code for defining lower and upper limits (a mess...)
                #if k >= 0:
                #    flow = hentr.GetBinContent(hentr.FindBin(mlcut))
                #    fhigh = hentr.GetBinContent(hentr.FindBin(mlcut))
                #    havezero = False
                #    havenallzero = True
                #    for z in range(k):
                #        if fentr[z].Eval(mlcut) < flow and fentr[z].Eval(mlcut) > 0:
                #            if fbc < self.minbc and (k - z) < 0.5 * self.n_bkgentrfits:
                #                flow = fentr[z].Eval(mlcut)
                #            elif fentr[z].Eval(mlcut) > 0.5 * fbc:
                #                flow = fentr[z].Eval(mlcut)
                #        if fentr[z].Eval(mlcut) > fhigh:
                #            if fbc < self.minbc and (k - z) < 0.5 * self.n_bkgentrfits:
                #                fhigh = fentr[z].Eval(mlcut)
                #            elif fentr[z].Eval(mlcut) > 1.5 * fbc:
                #                fhigh = fentr[z].Eval(mlcut)
                #        if fentr[z].Eval(mlcut) == 0:
                #            havezero = True
                #        if fentr[z].Eval(mlcut) > 0:
                #            havenallzero = False
                #    if flow < minfit and flow > 0: minfit = flow
                #    if flow <= 0 and minfit != 9999: flow = minfit
                #    if havezero and minfit < 9999: flow = minfit
                #    if flow <= 0 and minfit == 9999: flow = fbc
                #    if fhigh > 0: maxfit = fhigh
                #    if havenallzero and maxfit != -1: fhigh = maxfit
                #    par0low = self.calculate_par0(ipt, v_fpar1low, flow, self.bkg_fmap[self.p_bkgfunc[ipt]])
                #    par0high = self.calculate_par0(ipt, v_fpar1high, fhigh, self.bkg_fmap[self.p_bkgfunc[ipt]])
                #else:
                #    bcmin = hentr.GetBinContent(hentr.FindBin(mlcut)) - \
                #            hentr.GetBinError(hentr.FindBin(mlcut))
                #    bcmax = hentr.GetBinContent(hentr.FindBin(mlcut)) + \
                #            hentr.GetBinError(hentr.FindBin(mlcut))
                #    par0low = self.calculate_par0(ipt, v_fpar1low, bcmin, self.bkg_fmap[self.p_bkgfunc[ipt]])
                #    par0high = self.calculate_par0(ipt, v_fpar1high, bcmax, self.bkg_fmap[self.p_bkgfunc[ipt]])

                bccent = hentr.GetBinContent(hentr.FindBin(mlcut))
                bcmin = hentr.GetBinContent(hentr.FindBin(mlcut)) - \
                        hentr.GetBinError(hentr.FindBin(mlcut))
                bcmax = hentr.GetBinContent(hentr.FindBin(mlcut)) + \
                        hentr.GetBinError(hentr.FindBin(mlcut))
                if bccent < 50:
                    bcmin = bccent
                    bcmax = bccent
                fitcent = bccent
                fitbefore = bccent
                fitafter = bccent
                if k == 0:
                    if fentr[k].Eval(mlcut) > 0: fitcent = fentr[k].Eval(mlcut)
                    fitbefore = fitcent
                    fitafter = fitcent
                if k > 0 and k != self.n_bkgentrfits - 1:
                    if fentr[k].Eval(mlcut) > 0: fitcent = fentr[k].Eval(mlcut)
                    if fentr[k-1].Eval(mlcut) > 0: fitbefore = fentr[k-1].Eval(mlcut)
                    if fentr[k+1].Eval(mlcut) > 0: fitafter = fentr[k+1].Eval(mlcut)
                if k == self.n_bkgentrfits - 1:
                    if fentr[k].Eval(mlcut) > 0: fitcent = fentr[k].Eval(mlcut)
                    if fentr[k-1].Eval(mlcut) > 0: fitbefore = fentr[k-1].Eval(mlcut)
                    if bccent > 25 and fitbefore/bccent > 1.25: fitbefore = fitcent
                    fitafter = fitcent
                flow = min(bcmin, bccent, bcmax, fitbefore, fitcent, fitafter)
                fhigh = max(bcmin, bccent, bcmax, fitbefore, fitcent, fitafter)

                if flow == 0: flow = minfit
                if fhigh == 0: fhigh = maxfit
                minfit = flow
                maxfit = fhigh

                par0low = self.calculate_par0(ipt, v_fpar1, flow, self.bkg_fmap[self.p_bkgfunc[ipt]])
                par0high = self.calculate_par0(ipt, v_fpar1, fhigh, self.bkg_fmap[self.p_bkgfunc[ipt]])

                if par0low > par0new: par0low = par0new
                if par0high < par0new: par0high = par0new
                gpar0cent[ipt].SetPoint(isc, mlcut, par0new)
                gpar0low[ipt].SetPoint(isc, mlcut, par0low)
                gpar0high[ipt].SetPoint(isc, mlcut, par0high)
                gpar1cent[ipt].SetPoint(isc, mlcut, v_fpar1)
                gpar1low[ipt].SetPoint(isc, mlcut, v_fpar1low)
                gpar1high[ipt].SetPoint(isc, mlcut, v_fpar1high)

            canv.cd(3)
            gpar0cent[ipt].SetLineColor(2)
            gpar0low[ipt].SetLineColor(2)
            gpar0high[ipt].SetLineColor(2)
            gpar0cent[ipt].SetLineWidth(2)
            gpar0high[ipt].SetLineWidth(2)
            gpar0low[ipt].SetLineWidth(2)
            gpar0high[ipt].SetLineStyle(2)
            gpar0low[ipt].SetLineStyle(2)
            gpar0cent[ipt].Draw("same l")
            gpar0low[ipt].Draw("same l")
            gpar0high[ipt].Draw("same l")

            canv.cd(2)
            gpar1cent[ipt].SetLineColor(4)
            gpar1low[ipt].SetLineColor(4)
            gpar1high[ipt].SetLineColor(4)
            gpar1cent[ipt].SetLineWidth(2)
            gpar1high[ipt].SetLineWidth(2)
            gpar1low[ipt].SetLineWidth(2)
            gpar1high[ipt].SetLineStyle(2)
            gpar1low[ipt].SetLineStyle(2)
            gpar1cent[ipt].Draw("same l")
            gpar1low[ipt].Draw("same l")
            gpar1high[ipt].Draw("same l")

            fileout.cd()
            canv.Write("BkgFits_parametrised_%d" % ipt)

        for ipt in range(self.p_nptfinbins):
            fpar1[ipt].Write("fpar1cent_%d" % ipt)
            fpar1low[ipt].Write("fpar1low_%d" % ipt)
            fpar1high[ipt].Write("fpar1high_%d" % ipt)
            gpar0cent[ipt].Write("gpar0cent_%d" % ipt)
            gpar0low[ipt].Write("gpar0low_%d" % ipt)
            gpar0high[ipt].Write("gpar0high_%d" % ipt)
            gpar1cent[ipt].Write("gpar1cent_%d" % ipt)
            gpar1low[ipt].Write("gpar1low_%d" % ipt)
            gpar1high[ipt].Write("gpar1high_%d" % ipt)
        fileout.Close()

        self.logger.info("\n\nBackground parametrising finished, saved in %s\n\n", self.n_file_params)


    def calculate_par0(self, ipt, fpar1, fbc, function):
        if function == "expo":
            #Integrate[Exp[a + b x], {x, x1, x2}]
            #Solve[(Exp[a] (Exp[b x2] - Exp[b x1]))/b == c, a]
            if fpar1 == 0:
                return math.log(self.binwidth[ipt] * fbc)
            else:
                numerator = self.binwidth[ipt] * fpar1 * fbc
                denominator = math.exp(fpar1 * self.p_mass_fit_lim[1]) - math.exp(fpar1 * self.p_mass_fit_lim[0])
                return math.log(numerator / denominator)
        if function == "pol1":
            #Integrate[a + b x, {x, x1, x2}]
            #Solve[a (x2 - x1) + 0.5*b (x2^2 - x1^2) == c, a]
            numerator = self.binwidth[ipt] * fbc - 0.5 * fpar1 * (self.p_mass_fit_lim[1]*self.p_mass_fit_lim[1] - self.p_mass_fit_lim[0]*self.p_mass_fit_lim[0])
            denominator = self.p_mass_fit_lim[1] - self.p_mass_fit_lim[0]
            return numerator / denominator


    def parametrise_background_scan(self):
        """
        Function to run full background parametrisation in one go
        """

        if self.p_period is "merged":
            self.logger.warning("Invalid option for ITSUpgrade analyser, skipping")
            return
        self.fit_invmassbkg_scan()
        self.fit_bkgparams_2_scan()

    def expected_significance_print(self):
        """
        Calculate expected significance
        """

        if self.p_period is "merged":
            self.logger.warning("Invalid option for ITSUpgrade analyser, skipping")
            return

        # Define limits first
        self.define_probscan_limits()

        fmass = TFile(self.n_filemass_probscan, "READ")
        feff = TFile(self.n_fileeff_probscanfinal, "READ")
        fbkgcorr = TFile(self.n_filebkgcorrection, "READ")
        ffonll = TFile(self.n_filefonllinput, "READ")
        ftamu = TFile(self.n_filetamuinput, "READ")
        ffitpars = TFile(self.n_file_params, "READ")

        self.logger.info("Calculating expected significance for %s", self.p_period)
        TGaxis.SetMaxDigits(3)
        gStyle.SetOptStat(0000)

        hbkgcorr = fbkgcorr.Get("hCorrFacBs")
        hfonll = ffonll.Get("hFONLLcent")
        hfonllmin = ffonll.Get("hFONLLmin")
        hfonllmax = ffonll.Get("hFONLLmax")
        htamu = ftamu.Get("hTAMUcent")
        arrprint = []
        for ipt in range(self.p_nptfinbins):
            bkgcorr = hbkgcorr.GetBinContent(hbkgcorr.FindBin(self.lpt_finbinmax[ipt] - 0.1))
            fonll = hfonll.GetBinContent(hfonll.FindBin(self.lpt_finbinmax[ipt] - 0.1))
            fonllmin = hfonllmin.GetBinContent(hfonll.FindBin(self.lpt_finbinmax[ipt] - 0.1))
            fonllmax = hfonllmax.GetBinContent(hfonll.FindBin(self.lpt_finbinmax[ipt] - 0.1))
            tamu = htamu.GetBinContent(htamu.FindBin(self.lpt_finbinmax[ipt] - 0.1))
            for isc in range(self.n_probscan):
                if isc < self.n_probscan - self.v_ninterestingtrials:
                    continue
                heff = feff.Get("eff_%d" % isc)

                suffix = "%s%d_%d_%d" % (self.v_var_binning,
                                         self.lpt_finbinmin[ipt],
                                         self.lpt_finbinmax[ipt], isc)
                hsig = fmass.Get("hmass_sig" + suffix)
                hbkg = fmass.Get("hmass_bkg" + suffix)

                #Fit signal for 3sigma range
                fsig = TF1("fsig", "gaus", 5.32, 5.42)
                hsig.Fit("fsig", "R")
                mean = fsig.GetParameter(1)
                sigma = fsig.GetParameter(2)
                xmin = mean - 3 * sigma
                xmax = mean + 3 * sigma

                #Extract background
                hbkg.Rebin(self.rebins[ipt])
                fbkg = TF1("fbkg", self.bkg_fmap[self.p_bkgfunc[ipt]],
                           self.p_mass_fit_lim[0], self.p_mass_fit_lim[1])
                hbkg.Fit("fbkg", "R,E,+,0")

                grpar0cent = ffitpars.Get("gpar0cent_%d" % ipt)
                grpar0low = ffitpars.Get("gpar0low_%d" % ipt)
                grpar0up = ffitpars.Get("gpar0high_%d" % ipt)
                fpar1cent = ffitpars.Get("fpar1cent_%d" % ipt)
                fpar1low = ffitpars.Get("fpar1low_%d" % ipt)
                fpar1up = ffitpars.Get("fpar1high_%d" % ipt)

                fbkgcent = fbkg.Clone("fbkgcent")
                fbkgcent.SetParameter(0, grpar0cent.Eval(self.a_probscan[ipt][isc]))
                dpar1cent = fpar1cent.Eval(self.a_probscan[ipt][isc])
                if dpar1cent > 0: dpar1cent = 0
                fbkgcent.SetParameter(1, dpar1cent)

                fbkglow = fbkg.Clone("fbkglow")
                fbkglow.SetParameter(0, grpar0low.Eval(self.a_probscan[ipt][isc]))
                dpar1low = fpar1low.Eval(self.a_probscan[ipt][isc])
                if dpar1low > 0: dpar1low = 0
                fbkglow.SetParameter(1, dpar1low)

                fbkgup = fbkg.Clone("fbkgup")
                fbkgup.SetParameter(0, grpar0up.Eval(self.a_probscan[ipt][isc]))
                dpar1up = fpar1up.Eval(self.a_probscan[ipt][isc])
                if dpar1up > 0: dpar1up = 0
                fbkgup.SetParameter(1, dpar1up)

                bkgcentbc = hbkg.Integral(hbkg.FindBin(xmin), hbkg.FindBin(xmax))
                bkgcentorigfit = fbkg.Integral(xmin, xmax)/hbkg.GetBinWidth(1)
                bkgcentparafit = fbkgcent.Integral(xmin, xmax)/hbkg.GetBinWidth(1)

                #Extract efficiency
                effcent = heff.GetBinContent(heff.FindBin(self.lpt_finbinmax[ipt] - 0.1))
                effcenterr = heff.GetBinError(heff.FindBin(self.lpt_finbinmax[ipt] - 0.1))

                #Calculate expected significance
                expsigcent = 2 * (self.lpt_finbinmax[ipt] - self.lpt_finbinmin[ipt]) * 1 * \
                             self.brbs * self.brds * self.nevexpected * effcent * self.taa * \
                             fonll * tamu * self.gauss3sigmafac
                expbkgcentbc = bkgcorr * bkgcentbc * self.nevexpected / self.nevanalysed
                expbkgcentorigfit = bkgcorr * bkgcentorigfit * self.nevexpected / self.nevanalysed
                expbkgcentparafit = bkgcorr * bkgcentparafit * self.nevexpected / self.nevanalysed
                expsignfcentbc = expsigcent / math.sqrt(expsigcent + expbkgcentbc)
                expsignfcentorigfit = expsigcent / math.sqrt(expsigcent + expbkgcentorigfit)
                expsignfcentparafit = expsigcent / math.sqrt(expsigcent + expbkgcentparafit)

                arrprint.append(str("%.6f, %.4f, %.4f, %.4f, %.4f, %.4f" % (self.a_probscan[ipt][isc], expsigcent, expbkgcentbc, expsignfcentbc, expsignfcentorigfit, expsignfcentparafit)))
            arrprint.append(str("\n"))

        for ipr in range(len(arrprint)):
            print(arrprint[ipr])

    def bkgshapestudy_mass_histo(self):
        """
        Invariant mass histogram for injected Ds + HIJING pion
        """

        if self.p_period is "merged":
            self.logger.warning("Invalid option for ITSUpgrade analyser, skipping")
            return

        # Define limits first
        self.define_probscan_limits()
        for ipt in range(self.p_nptfinbins):
            self.a_probscan[ipt].insert(0, 0)

        myfile = TFile.Open(self.n_filemass_bkgshape, "recreate")

        self.logger.info("Doing bkg shape study for period %s", self.p_period)

        for ipt in range(self.p_nptfinbins):
            bin_id = self.bin_matching[ipt]
            df = pickle.load(openfile(self.lpt_recodecmerged_data[bin_id], "rb"))
            preselapplied = False

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

            df_dspr = df[df[self.v_dsprompt] == 1]
            df_dsfdbplus = df[df[self.v_dsfdbplus] == 1]
            df_dsfdbzero = df[df[self.v_dsfdbzero] == 1]
            df_dsfdlambdab = df[df[self.v_dsfdlambdab] == 1]
            df_dsfdbs = df[df[self.v_dsfdbs] == 1]
            unique_dspr = average([len(df_dspr["d_len_Ds"].unique()),
                                   len(df_dspr["imp_par_xy_Ds"].unique()),
                                   len(df_dspr["cos_PiDs_Ds"].unique())])
            unique_dsfdbplus = average([len(df_dsfdbplus["d_len_Ds"].unique()),
                                        len(df_dsfdbplus["imp_par_xy_Ds"].unique()),
                                        len(df_dsfdbplus["cos_PiDs_Ds"].unique())])
            unique_dsfdbzero = average([len(df_dsfdbzero["d_len_Ds"].unique()),
                                        len(df_dsfdbzero["imp_par_xy_Ds"].unique()),
                                        len(df_dsfdbzero["cos_PiDs_Ds"].unique())])
            unique_dsfdlambdab = average([len(df_dsfdlambdab["d_len_Ds"].unique()),
                                          len(df_dsfdlambdab["imp_par_xy_Ds"].unique()),
                                          len(df_dsfdlambdab["cos_PiDs_Ds"].unique())])
            unique_dsfdbs = average([len(df_dsfdbs["d_len_Ds"].unique()),
                                     len(df_dsfdbs["imp_par_xy_Ds"].unique()),
                                     len(df_dsfdbs["cos_PiDs_Ds"].unique())])

            myfile.cd()
            hNorm = TH1F("hNormUniqueDs" + ptstring, ";;Normalisation", 5, 0.5, 5.5)
            hNorm.GetXaxis().SetBinLabel(1, "prompt Ds")
            hNorm.GetXaxis().SetBinLabel(2, "feed-down Ds (B+)")
            hNorm.GetXaxis().SetBinLabel(3, "feed-down Ds (B0)")
            hNorm.GetXaxis().SetBinLabel(4, "feed-down Ds (Lb)")
            hNorm.GetXaxis().SetBinLabel(5, "feed-down Ds (Bs)")
            hNorm.SetBinContent(1, unique_dspr)
            hNorm.SetBinContent(2, unique_dsfdbplus)
            hNorm.SetBinContent(3, unique_dsfdbzero)
            hNorm.SetBinContent(4, unique_dsfdlambdab)
            hNorm.SetBinContent(5, unique_dsfdbs)
            hNorm.Write()

            for isc in range(len(self.a_probscan[ipt])):
                h_isc_match.Fill(isc, self.a_probscan[ipt][isc])
                selml_cv = "y_test_prob%s>%s" % (self.p_modelname, self.a_probscan[ipt][isc])
                df = df.query(selml_cv)

                #FIXME, for non-sequential case
                if isc > 0 and self.a_probscan[ipt][0] == 0 and preselapplied is False:
                    modelnamepref = self.p_modelname + "prefilter"
                    selml_prefilter = "y_test_prob%s>%s" % (modelnamepref, self.probscan_prefilter[ipt])
                    df = df.query(selml_prefilter)
                    preselapplied = True

                if self.lvar2_binmin is None:
                    suffix = "%s%d_%d_%d" % (self.v_var_binning,
                                             self.lpt_finbinmin[ipt],
                                             self.lpt_finbinmax[ipt], isc)
                    h_invmass_dspr = TH1F("hmass_DsPr" + suffix, "", self.p_num_bins,
                                          self.p_mass_fit_lim[0], self.p_mass_fit_lim[1])
                    h_invmass_dsfdbplus = TH1F("hmass_DsFDBplus" + suffix, "", self.p_num_bins,
                                               self.p_mass_fit_lim[0], self.p_mass_fit_lim[1])
                    h_invmass_dsfdbzero = TH1F("hmass_DsFDBzero" + suffix, "", self.p_num_bins,
                                               self.p_mass_fit_lim[0], self.p_mass_fit_lim[1])
                    h_invmass_dsfdlambdab = TH1F("hmass_DsFDLambdab" + suffix, "", self.p_num_bins,
                                                self.p_mass_fit_lim[0], self.p_mass_fit_lim[1])
                    h_invmass_dsfdbs = TH1F("hmass_DsFDBs" + suffix, "", self.p_num_bins,
                                            self.p_mass_fit_lim[0], self.p_mass_fit_lim[1])

                    df_dspr = df[df[self.v_dsprompt] == 1]
                    df_dsfdbplus = df[df[self.v_dsfdbplus] == 1]
                    df_dsfdbzero = df[df[self.v_dsfdbzero] == 1]
                    df_dsfdlambdab = df[df[self.v_dsfdlambdab] == 1]
                    df_dsfdbs = df[df[self.v_dsfdbs] == 1]

                    fill_hist(h_invmass_dspr, df_dspr.inv_mass)
                    fill_hist(h_invmass_dsfdbplus, df_dsfdbplus.inv_mass)
                    fill_hist(h_invmass_dsfdbzero, df_dsfdbzero.inv_mass)
                    fill_hist(h_invmass_dsfdlambdab, df_dsfdlambdab.inv_mass)
                    fill_hist(h_invmass_dsfdbs, df_dsfdbs.inv_mass)

                    h_norm_dspr = h_invmass_dspr.Clone("h_norm_dspr" + suffix)
                    h_norm_dspr.Rebin(self.rebins[ipt])
                    if h_norm_dspr.GetEntries() != 0:
                        h_norm_dspr.Scale(1. / h_norm_dspr.GetEntries())
                    h_norm_dspr.SetLineColor(kBlack)
                    fbkg_dspr = TF1("fbkg_dspr" + suffix, "pol2",
                                    self.p_mass_fit_lim[0], self.p_mass_fit_lim[1])
                    fbkg_dspr.SetLineColor(kBlack)
                    h_norm_dspr.Fit("fbkg_dspr" + suffix, "R,E,+,0")
                    fbkg_dspr1 = TF1("fbkg_dspr1" + suffix, "pol1",
                                     self.p_mass_fit_lim[0], self.p_mass_fit_lim[1])
                    fbkg_dspr1.SetLineColor(kBlack)
                    h_norm_dspr.Fit("fbkg_dspr1" + suffix, "R,E,+,0")

                    h_norm_dsfdbplus = h_invmass_dsfdbplus.Clone("h_norm_dsfdbplus" + suffix)
                    h_norm_dsfdbplus.Rebin(self.rebins[ipt])
                    if h_norm_dsfdbplus.GetEntries() != 0:
                        h_norm_dsfdbplus.Scale(1. / h_norm_dsfdbplus.GetEntries())
                    h_norm_dsfdbplus.SetLineColor(kRed+1)
                    fbkg_dsfdbplus = TF1("fbkg_dsfdbplus" + suffix, "pol2",
                                         self.p_mass_fit_lim[0], self.p_mass_fit_lim[1])
                    fbkg_dsfdbplus.SetLineColor(kRed+1)
                    h_norm_dsfdbplus.Fit("fbkg_dsfdbplus" + suffix, "R,E,+,0")
                    fbkg_dsfdbplus1 = TF1("fbkg_dsfdbplus1" + suffix, "pol1",
                                          self.p_mass_fit_lim[0], self.p_mass_fit_lim[1])
                    fbkg_dsfdbplus1.SetLineColor(kRed+1)
                    h_norm_dsfdbplus.Fit("fbkg_dsfdbplus1" + suffix, "R,E,+,0")

                    h_norm_dsfdbzero = h_invmass_dsfdbzero.Clone("h_norm_dsfdbzero" + suffix)
                    h_norm_dsfdbzero.Rebin(self.rebins[ipt])
                    if h_norm_dsfdbzero.GetEntries() != 0:
                        h_norm_dsfdbzero.Scale(1. / h_norm_dsfdbzero.GetEntries())
                    h_norm_dsfdbzero.SetLineColor(kBlue+1)
                    fbkg_dsfdbzero = TF1("fbkg_dsfdbzero" + suffix, "pol2",
                                         self.p_mass_fit_lim[0], self.p_mass_fit_lim[1])
                    fbkg_dsfdbzero.SetLineColor(kBlue+1)
                    h_norm_dsfdbzero.Fit("fbkg_dsfdbzero" + suffix, "R,E,+,0")
                    fbkg_dsfdbzero1 = TF1("fbkg_dsfdbzero1" + suffix, "pol1",
                                          self.p_mass_fit_lim[0], self.p_mass_fit_lim[1])
                    fbkg_dsfdbzero1.SetLineColor(kBlue+1)
                    h_norm_dsfdbzero.Fit("fbkg_dsfdbzero1" + suffix, "R,E,+,0")

                    h_norm_dsfdlambdab = h_invmass_dsfdlambdab.Clone("h_norm_dsfdlambdab" + suffix)
                    h_norm_dsfdlambdab.Rebin(self.rebins[ipt])
                    if h_norm_dsfdlambdab.GetEntries() != 0:
                        h_norm_dsfdlambdab.Scale(1. / h_norm_dsfdlambdab.GetEntries())
                    h_norm_dsfdlambdab.SetLineColor(kGreen+2)
                    fbkg_dsfdlambdab = TF1("fbkg_dsfdlambdab" + suffix, "pol2",
                                          self.p_mass_fit_lim[0], self.p_mass_fit_lim[1])
                    fbkg_dsfdlambdab.SetLineColor(kGreen+2)
                    h_norm_dsfdlambdab.Fit("fbkg_dsfdlambdab" + suffix, "R,E,+,0")
                    fbkg_dsfdlambdab1 = TF1("fbkg_dsfdlambdab1" + suffix, "pol1",
                                            self.p_mass_fit_lim[0], self.p_mass_fit_lim[1])
                    fbkg_dsfdlambdab1.SetLineColor(kGreen+2)
                    h_norm_dsfdlambdab.Fit("fbkg_dsfdlambdab1" + suffix, "R,E,+,0")

                    h_norm_dsfdbs = h_invmass_dsfdbs.Clone("h_norm_dsfdbs" + suffix)
                    h_norm_dsfdbs.Rebin(self.rebins[ipt])
                    if h_norm_dsfdbs.GetEntries() != 0:
                        h_norm_dsfdbs.Scale(1. / h_norm_dsfdbs.GetEntries())
                    h_norm_dsfdbs.SetLineColor(kOrange+2)
                    fbkg_dsfdbs = TF1("fbkg_dsfdbs" + suffix, "pol2",
                                      self.p_mass_fit_lim[0], self.p_mass_fit_lim[1])
                    fbkg_dsfdbs.SetLineColor(kOrange+2)
                    h_norm_dsfdbs.Fit("fbkg_dsfdbs" + suffix, "R,E,+,0")
                    fbkg_dsfdbs1 = TF1("fbkg_dsfdbs1" + suffix, "pol1",
                                       self.p_mass_fit_lim[0], self.p_mass_fit_lim[1])
                    fbkg_dsfdbs1.SetLineColor(kOrange+2)
                    h_norm_dsfdbs.Fit("fbkg_dsfdbs1" + suffix, "R,E,+,0")

                    self.binwidth[ipt] = h_norm_dspr.GetBinWidth(1)
                    h_norm_dspr.GetYaxis().SetTitle("Counts per %.2f MeV/#it{c}^{2} (norm.)" % self.binwidth[ipt])
                    h_norm_dspr.GetXaxis().SetTitle("#it{M}(KK#pi#pi) (GeV/#it{c}^{2})")

                    myfile.cd()
                    h_invmass_dspr.Write()
                    h_norm_dspr.Write()
                    fbkg_dspr.Write()
                    fbkg_dspr1.Write()
                    h_invmass_dsfdbplus.Write()
                    h_norm_dsfdbplus.Write()
                    fbkg_dsfdbplus.Write()
                    fbkg_dsfdbplus1.Write()
                    h_invmass_dsfdbzero.Write()
                    h_norm_dsfdbzero.Write()
                    fbkg_dsfdbzero.Write()
                    fbkg_dsfdbzero1.Write()
                    h_invmass_dsfdlambdab.Write()
                    h_norm_dsfdlambdab.Write()
                    fbkg_dsfdlambdab.Write()
                    fbkg_dsfdlambdab1.Write()
                    h_invmass_dsfdbs.Write()
                    h_norm_dsfdbs.Write()
                    fbkg_dsfdbs.Write()
                    fbkg_dsfdbs1.Write()
                else:
                    self.logger.fatal("2nd binning to be implemented (code exists)")
            h_isc_match.Write()
