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
main script for doing ml optimisation
"""
import os
import sys
import pickle
import time
from math import sqrt
from array import array
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from sklearn.model_selection import train_test_split
from sklearn.utils import shuffle
import xgboost as xgb
from hipe4ml import plot_utils
from hipe4ml.model_handler import ModelHandler
from ROOT import TFile, TCanvas, TH1F, TH2F, TF1, gROOT  # pylint: disable=import-error,no-name-in-module

from machine_learning_hep.logger import get_logger
from machine_learning_hep.utilities import createstringselection, openfile
from machine_learning_hep.utilities_selection import seldf_singlevar
from machine_learning_hep.utilities_selection import selectdfquery
from machine_learning_hep.utilities_selection import split_df_sigbkg
from machine_learning_hep.ml_functions import apply
from machine_learning_hep.ml_significance import calc_bkg, calc_signif, calc_eff, calc_sigeff_steps

class Optimiserhipe4ml:
    # Class Attribute
    species = "optimiser_hipe4ml"

    def __init__(self, data_param, case, typean, binmin, binmax, training_var, hyper_pars, raahp):

        self.logger = get_logger()

        self.do_mlprefilter = data_param.get("doml_asprefilter", None)
        self.v_bin = data_param["var_binning"]
        # directory
        dirmcml = data_param["multi"]["mc"]["pkl_skimmed_merge_for_ml_all"]
        dirdataml = data_param["multi"]["data"]["pkl_skimmed_merge_for_ml_all"]
        dirmcml_max = data_param["multi"]["mc"]["pkl_skimmed_merge_for_ml_all"] + "_max"
        dirdataml_max = data_param["multi"]["data"]["pkl_skimmed_merge_for_ml_all"] + "_max"
        self.v_max_ncand_merge = data_param["multi"]["max_ncand_merge"]
        if self.do_mlprefilter is False: #FIXME for multiple periods
            dirmcml = data_param["mlapplication"]["mc"]["pkl_skimmed_decmerged"][0] + "/prefilter"
            dirdataml = data_param["mlapplication"]["data"]["pkl_skimmed_decmerged"][0] + "/prefilter"
        dirdatatotsample = data_param["multi"]["data"]["pkl_evtcounter_all"]
        self.p_fracmerge = data_param["multi"]["data"]["fracmerge"][0] #FIXME for multiple periods
        # directory
        self.dirmlout = data_param["ml"]["mlout"]
        self.dirmlplot = data_param["ml"]["mlplot"]
        #if self.do_mlprefilter is True:
        #    self.dirmodel = self.dirmodel + "/prefilter"
        #    self.dirmlplot = self.dirmlplot + "/prefilter"
        #if self.do_mlprefilter is False:
        #    self.dirmodel = self.dirmodel + "/analysis"
        #    self.dirmlplot = self.dirmlplot + "/analysis"
        # ml file names
        self.n_reco = data_param["files_names"]["namefile_reco"]
        self.n_reco = self.n_reco.replace(".pkl", "_%s%d_%d.pkl" % (self.v_bin, binmin, binmax))
        self.n_evt = data_param["files_names"]["namefile_evt"]
        self.n_gen = data_param["files_names"]["namefile_gen"]
        self.n_gen = self.n_gen.replace(".pkl", "_%s%d_%d.pkl" % (self.v_bin, binmin, binmax))
        # ml files
        self.f_gen_mc = os.path.join(dirmcml, self.n_gen)
        self.f_reco_mc = os.path.join(dirmcml, self.n_reco)
        self.f_reco_data = os.path.join(dirdataml, self.n_reco)
        self.f_reco_mc_max = os.path.join(dirmcml_max, self.n_reco)
        self.f_reco_data_max = os.path.join(dirdataml_max, self.n_reco)
        self.f_evt_data = os.path.join(dirdataml, self.n_evt)
        self.f_evttotsample_data = os.path.join(dirdatatotsample, self.n_evt)
        # variables
        self.v_train = training_var
        if isinstance(self.v_train[0], list):
            self.logger.error("Training hipe4ml not yet suitable for pT dependent training variables")
        self.v_sig = data_param["variables"]["var_signal"]
        # parameters
        self.p_case = case
        self.p_typean = typean
        self.p_nbkg = data_param["ml"]["nbkg"]
        self.p_nbkgfd = data_param["ml"].get("nbkg", self.p_nbkg)
        self.p_nsig = data_param["ml"]["nsig"]
        self.v_fracbkgoversig = data_param["ml"].get("fracbkgoversig", None)
        self.p_tagsig = data_param["ml"]["sampletagforsignal"]
        self.p_tagbkg = data_param["ml"]["sampletagforbkg"]
        self.p_tagbkgfd = data_param["ml"].get("sampletagforbkgfd", 2)
        self.p_binmin = binmin
        self.p_binmax = binmax
        self.rnd_shuffle = data_param["ml"]["rnd_shuffle"]
        self.rnd_splt = data_param["ml"]["rnd_splt"]
        self.test_frac = data_param["ml"]["test_frac"]
        self.p_dofullevtmerge = data_param["dofullevtmerge"]

        self.p_evtsel = data_param["ml"]["evtsel"]
        self.p_triggersel_mc = data_param["ml"]["triggersel"]["mc"]
        self.p_triggersel_data = data_param["ml"]["triggersel"]["data"]

        # dataframes
        self.df_mc = None
        self.df_mcgen = None
        self.df_data = None
        self.df_data_signf = None
        self.df_sig = None
        self.df_bkg = None
        self.df_ml = None
        self.df_mltest = None
        self.df_mltrain = None
        self.df_sigtrain = None
        self.df_sigtest = None
        self.df_bkgtrain = None
        self.df_xtrain = None
        self.df_ytrain = None
        self.df_xtest = None
        self.df_ytest = None
        self.traintestdata = None
        self.ypredtrain_hipe4ml = None
        self.ypredtest_hipe4ml = None
        self.df_evt_data = None
        self.df_evttotsample_data = None
        # selections
        self.s_selbkgml = data_param["ml"]["sel_bkgml"]
        self.s_selbkgmlfd = data_param["ml"].get("sel_bkgmlfd", None)
        self.s_selsigml = data_param["ml"]["sel_sigml"]
        self.p_mltype = data_param["ml"]["mltype"]
        self.p_presel_gen_eff = data_param["ml"]["opt"]["presel_gen_eff"]
        # mask missing values
        data_maskmissingvalues = data_param.get("maskmissingvalues", None)
        if data_maskmissingvalues is not None:
            self.b_maskmissing = data_param["maskmissingvalues"].get("activate", False)
            self.v_varstomask = data_param["maskmissingvalues"].get("tomask", None)
        else:
            self.b_maskmissing = False
            self.v_varstomask = None

        self.s_suffix = None
        self.create_suffix()
        self.preparesample()

        self.p_hipe4ml_model = None
        self.p_hipe4ml_origmodel = None
        self.v_hipe4ml_pars = hyper_pars
        self.load_hipe4mlmodel()

        self.bayesoptconfig_hipe4ml = data_param["hipe4ml"]["hyper_par_opt"]["bayes_opt_config"]
        self.average_method_hipe4ml = data_param["hipe4ml"]["roc_auc_average"]
        self.nfold_hipe4ml = data_param["hipe4ml"]["hyper_par_opt"]["nfolds"]
        self.init_points = data_param["hipe4ml"]["hyper_par_opt"]["initpoints"]
        self.n_iter_hipe4ml = data_param["hipe4ml"]["hyper_par_opt"]["niter"]
        self.njobs_hipe4ml = data_param["hipe4ml"]["hyper_par_opt"]["njobs"]
        self.roc_method_hipe4ml = data_param["hipe4ml"]["roc_auc_approach"]
        self.raw_output_hipe4ml = data_param["hipe4ml"]["raw_output"]
        self.train_test_log_hipe4ml = data_param["hipe4ml"]["train_test_log"]

        #significance
        self.is_fonll_from_root = data_param["ml"]["opt"]["isFONLLfromROOT"]
        self.f_fonll = data_param["ml"]["opt"]["filename_fonll"]
        if self.is_fonll_from_root and "fonll_particle" not in data_param["ml"]["opt"]:
            self.logger.fatal("Attempt to read FONLL from ROOT file but field " \
                    "\"fonll_particle\" not provided in database")
        self.p_fonllparticle = data_param["ml"]["opt"].get("fonll_particle", "")
        self.p_fonllband = data_param["ml"]["opt"]["fonll_pred"]
        self.p_fragf = data_param["ml"]["opt"]["FF"]
        self.p_sigmamb = data_param["ml"]["opt"]["sigma_MB"]
        self.p_taa = data_param["ml"]["opt"]["Taa"]
        self.p_br = data_param["ml"]["opt"]["BR"]
        self.p_fprompt = data_param["ml"]["opt"]["f_prompt"]
        self.p_bkgfracopt = data_param["ml"]["opt"]["bkg_data_fraction"]
        self.p_nstepsign = data_param["ml"]["opt"]["num_steps"]
        self.p_bkg_func = data_param["ml"]["opt"]["bkg_function"]
        self.p_savefit = data_param["ml"]["opt"]["save_fit"]
        self.p_nevtml = None
        self.p_nevttot = None
        self.p_presel_gen_eff = data_param["ml"]["opt"]["presel_gen_eff"]
        self.p_mass_fit_lim = data_param["ml"]["opt"]['mass_fit_lim']
        self.p_bin_width = data_param["ml"]["opt"]['bin_width']
        self.p_num_bins = int(round((self.p_mass_fit_lim[1] - self.p_mass_fit_lim[0]) /
                                    self.p_bin_width))
        self.p_mass = data_param["mass"]
        self.p_raahp = raahp

        self.multiclass_labels = data_param["ml"].get("multiclass_labels", None)

        self.logger.info("Using the following training variables: %s", training_var)

    def create_suffix(self):
        string_selection = createstringselection(self.v_bin, self.p_binmin, self.p_binmax)
        self.s_suffix = f"{self.p_case}_{string_selection}"

    def preparesample(self):
        self.logger.info("Prepare Sample for hipe4ml")
        self.df_data = pickle.load(openfile(self.f_reco_data, "rb"))
        self.df_data_signf = pickle.load(openfile(self.f_reco_data, "rb"))
        if self.v_max_ncand_merge > 0:
            self.df_data = pickle.load(openfile(self.f_reco_data_max, "rb"))
        self.df_mc = pickle.load(openfile(self.f_reco_mc, "rb"))
        self.df_mcgen = pickle.load(openfile(self.f_gen_mc, "rb"))

        if self.b_maskmissing:
            self.df_data = self.df_data.replace(self.v_varstomask, value=np.nan)
            self.df_data_signf = self.df_data_signf.replace(self.v_varstomask, value=np.nan)
            self.df_mc = self.df_mc.replace(self.v_varstomask, value=np.nan)

        self.df_data = selectdfquery(self.df_data, self.p_evtsel)
        self.df_mc = selectdfquery(self.df_mc, self.p_evtsel)
        self.df_mcgen = selectdfquery(self.df_mcgen, self.p_evtsel)

        self.df_data = selectdfquery(self.df_data, self.p_triggersel_data)
        self.df_mc = selectdfquery(self.df_mc, self.p_triggersel_mc)
        self.df_mcgen = selectdfquery(self.df_mcgen, self.p_triggersel_mc)

        self.df_mcgen = self.df_mcgen.query(self.p_presel_gen_eff)
        arraydf = [self.df_data, self.df_mc, self.df_mc]
        self.df_mc = seldf_singlevar(self.df_mc, self.v_bin, self.p_binmin, self.p_binmax)
        self.df_mcgen = seldf_singlevar(self.df_mcgen, self.v_bin, self.p_binmin, self.p_binmax)
        self.df_data = seldf_singlevar(self.df_data, self.v_bin, self.p_binmin, self.p_binmax)

        self.df_sig, self.df_bkg, self.df_bkgfd = arraydf[self.p_tagsig], arraydf[self.p_tagbkg], arraydf[self.p_tagbkgfd]
        self.df_sig = seldf_singlevar(self.df_sig, self.v_bin, self.p_binmin, self.p_binmax)
        self.df_bkg = seldf_singlevar(self.df_bkg, self.v_bin, self.p_binmin, self.p_binmax)
        self.df_bkgfd = seldf_singlevar(self.df_bkgfd, self.v_bin, self.p_binmin, self.p_binmax)
        self.df_sig = self.df_sig.query(self.s_selsigml)
        self.df_bkg = self.df_bkg.query(self.s_selbkgml)
        if self.s_selbkgmlfd is not None:
            self.df_bkgfd = self.df_bkgfd.query(self.s_selbkgmlfd)
        else:
            self.df_bkgfd = pd.DataFrame()
        self.df_bkg["ismcsignal"] = 0
        self.df_bkg["ismcprompt"] = 0
        self.df_bkg["ismcfd"] = 0
        self.df_bkg["ismcbkg"] = 0

        self.p_nsig = min(len(self.df_sig), self.p_nsig)
        self.p_nbkg = min(len(self.df_bkg), self.p_nbkg)
        self.p_nbkgfd = min(len(self.df_bkgfd), self.p_nbkgfd)
        if self.p_nsig < self.p_nbkg and self.v_fracbkgoversig is not None:
            self.p_nbkg = self.v_fracbkgoversig * self.p_nsig
            if len(self.df_bkg) < self.p_nbkg:
                self.p_nbkg = len(self.df_bkg)
        if self.p_nsig < self.p_nbkgfd and self.v_fracbkgoversig is not None:
            self.p_nbkgfd = self.v_fracbkgoversig * self.p_nsig
            if len(self.df_bkgfd) < self.p_nbkgfd:
                self.p_nbkgfd = len(self.df_bkgfd)

        self.logger.info("Used number of signal events is %d", self.p_nsig)
        self.logger.info("Used number of background events is %d", self.p_nbkg)
        self.logger.info("Used number of background events (FD) is %d", self.p_nbkgfd)

        self.df_ml = pd.DataFrame()
        self.df_sig = shuffle(self.df_sig, random_state=self.rnd_shuffle)
        self.df_bkg = shuffle(self.df_bkg, random_state=self.rnd_shuffle)
        self.df_bkgfd = shuffle(self.df_bkgfd, random_state=self.rnd_shuffle)
        self.df_sig = self.df_sig[:self.p_nsig]
        self.df_bkg = self.df_bkg[:self.p_nbkg]
        self.df_bkgfd = self.df_bkgfd[:self.p_nbkgfd]
        self.df_sig[self.v_sig] = 1
        self.df_bkg[self.v_sig] = 0
        self.df_bkgfd[self.v_sig] = 2
        if len(self.df_bkgfd) == 0:
            self.n_classes = 2
            self.df_ml = pd.concat([self.df_sig, self.df_bkg])
        else:
            self.n_classes = 3
            self.df_ml = pd.concat([self.df_bkg, self.df_sig, self.df_bkgfd])
        self.df_mltrain, self.df_mltest = train_test_split(self.df_ml, test_size=self.test_frac,
                                                           random_state=self.rnd_splt)
        self.df_mltrain = self.df_mltrain.reset_index(drop=True)
        self.df_mltest = self.df_mltest.reset_index(drop=True)
        self.df_sigtrain, self.df_bkgtrain, self.df_bkgfdtrain = split_df_sigbkg(self.df_mltrain, self.v_sig)
        self.df_sigtest, self.df_bkgtest, self.df_bkgfdtest = split_df_sigbkg(self.df_mltest, self.v_sig)
        self.logger.info("Total number of candidates: train %d and test %d", len(self.df_mltrain),
                         len(self.df_mltest))
        self.logger.info("Number of signal candidates: train %d and test %d",
                         len(self.df_sigtrain), len(self.df_sigtest))
        self.logger.info("Number of bkg candidates: %d and test %d", len(self.df_bkgtrain),
                         len(self.df_bkgtest))
        if len(self.df_bkgfd) != 0:
            self.logger.info("Number of bkg FD candidates: %d and test %d", len(self.df_bkgfdtrain),
                             len(self.df_bkgfdtest))

        self.df_xtrain = self.df_mltrain[self.v_train]
        self.df_ytrain = self.df_mltrain[self.v_sig]
        self.df_xtest = self.df_mltest[self.v_train]
        self.df_ytest = self.df_mltest[self.v_sig]
        self.traintestdata = [self.df_xtrain, self.df_ytrain, self.df_xtest, self.df_ytest]
        #self.traintestdata = [self.df_mltrain, self.df_ytrain, self.df_mltest, self.df_ytest]

    def load_hipe4mlmodel(self):
        self.logger.info("Loading hipe4ml model")
        model_xgboost = xgb.XGBClassifier()
        self.p_hipe4ml_model = ModelHandler(model_xgboost, self.v_train, self.v_hipe4ml_pars)

    def set_hipe4ml_modelpar(self):
        self.logger.info("Setting hipe4ml hyperparameters")
        self.p_hipe4ml_model.set_model_params(self.v_hipe4ml_pars)

    def do_hipe4mlhyperparopti(self):
        self.logger.info("Optimising hipe4ml hyperparameters (Bayesian)")

        if self.n_classes > 2:
            if not (self.average_method_hipe4ml in ['macro', 'weighted'] and
                    self.roc_method_hipe4ml in ['ovo', 'ovr']):
                self.logger.fatal("Selected ROC configuration is not valid!")

            if self.average_method_hipe4ml == 'weighted':
                metric = f'roc_auc_{self.roc_method_hipe4ml}_{self.average_method_hipe4ml}'
            else:
                metric = f'roc_auc_{self.roc_method_hipe4ml}'
        else:
            metric = 'roc_auc'

        hypparsfile = f'{self.dirmlout}/HyperParOpt_pT_{self.p_binmin}_{self.p_binmax}.txt'
        outfilehyppars = open(hypparsfile, 'wt')
        sys.stdout = outfilehyppars
        self.p_hipe4ml_model.optimize_params_bayes(self.traintestdata, self.bayesoptconfig_hipe4ml,
                                                   metric, self.nfold_hipe4ml, self.init_points,
                                                   self.n_iter_hipe4ml, self.njobs_hipe4ml)
        outfilehyppars.close()
        sys.stdout = sys.__stdout__
        self.logger.info("Performing hyper-parameters optimisation: Done!")

    def do_hipe4mltrain(self):
        self.logger.info("Training + testing hipe4ml model")
        t0 = time.time()

        self.ypredtest_hipe4ml = self.p_hipe4ml_model.train_test_model(self.traintestdata, True,
                                                                       self.raw_output_hipe4ml,
                                                                       self.average_method_hipe4ml,
                                                                       self.roc_method_hipe4ml)
        self.ypredtrain_hipe4ml = self.p_hipe4ml_model.predict(self.traintestdata[0],
                                                               self.raw_output_hipe4ml)

        modelhandlerfile = f'{self.dirmlout}/ModelHandler_pT_{self.p_binmin}_{self.p_binmax}.pkl'
        self.p_hipe4ml_model.dump_model_handler(modelhandlerfile)
        modelfile = f'{self.dirmlout}/Model_pT_{self.p_binmin}_{self.p_binmax}.model'
        self.p_hipe4ml_model.dump_original_model(modelfile, True)

        self.p_hipe4ml_origmodel = self.p_hipe4ml_model.get_original_model()

        #filename_xtrain = self.dirmlout+"/xtrain_%s.pkl" % (self.s_suffix)
        #filename_ytrain = self.dirmlout+"/ytrain_%s.pkl" % (self.s_suffix)
        #pickle.dump(self.df_xtrain, openfile(filename_xtrain, "wb"), protocol=4)
        #pickle.dump(self.df_ytrain, openfile(filename_ytrain, "wb"), protocol=4)

        self.logger.info("Training + testing hipe4ml: Done!")
        self.logger.info("Time elapsed = %.3f", time.time() - t0)

    def get_hipe4mlmodel(self):
        self.logger.info("Getting already saved hipe4ml model")
        modelhandlerfile = f'{self.dirmlout}/ModelHandler_pT_{self.p_binmin}_{self.p_binmax}.pkl'
        self.p_hipe4ml_model = pickle.load(openfile(modelhandlerfile, 'rb'))
        self.p_hipe4ml_origmodel = self.p_hipe4ml_model.get_original_model()

    def do_hipe4mlplot(self, didtraining = True):
        self.logger.info("Plotting hipe4ml model")

        leglabels = ["Background", "Prompt signal"]
        outputlabels = ["Bkg", "SigPrompt"]
        listdf = [self.df_bkgtrain, self.df_sigtrain]
        if self.n_classes > 2:
            leglabels = ["Background", "Prompt", "Feed-down"]
            outputlabels = ["Bkg", "SigPrompt", "SigFeedDown"]
            listdf = [self.df_bkgtrain, self.df_sigtrain, self.df_bkgfdtrain]

        # _____________________________________________
        plot_utils.plot_distr(listdf, ["inv_mass", "pt_cand", "inv_mass_K0s", "pt_prong0", "pt_prong1", "pt_prong2"] + self.v_train, 100, leglabels, figsize=(12, 7),
                              alpha=0.3, log=True, grid=False, density=True)
        plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
        figname = f'{self.dirmlplot}/DistributionsAll_pT_{self.p_binmin}_{self.p_binmax}.pdf'
        plt.savefig(figname)
        plt.close('all')
        # _____________________________________________
        corrmatrixfig = plot_utils.plot_corr(listdf, ["inv_mass", "pt_cand", "inv_mass_K0s", "pt_prong0", "pt_prong1", "pt_prong2"] + self.v_train, leglabels)
        for figg, labb in zip(corrmatrixfig, outputlabels):
            plt.figure(figg.number)
            plt.subplots_adjust(left=0.2, bottom=0.25, right=0.95, top=0.9)
            figname = f'{self.dirmlplot}/CorrMatrix{labb}_pT_{self.p_binmin}_{self.p_binmax}.pdf'
            figg.savefig(figname)

        if didtraining:
            # _____________________________________________
            plt.rcParams["figure.figsize"] = (10, 7)
            mloutputfig = plot_utils.plot_output_train_test(self.p_hipe4ml_model, self.traintestdata,
                                                            80, self.raw_output_hipe4ml,
                                                            leglabels, self.train_test_log_hipe4ml,
                                                            density=True)
            if self.n_classes > 2:
                for figg, labb in zip(mloutputfig, outputlabels):
                    figname = f'{self.dirmlplot}/MLOutputDistr{labb}_pT_{self.p_binmin}_{self.p_binmax}.pdf'
                    figg.savefig(figname)
            else:
                figname = f'{self.dirmlplot}/MLOutputDistr_pT_{self.p_binmin}_{self.p_binmax}.pdf'
                mloutputfig.savefig(figname)
            # _____________________________________________
            plt.rcParams["figure.figsize"] = (10, 9)
            roccurvefig = plot_utils.plot_roc(self.traintestdata[3], self.ypredtest_hipe4ml,
                                              None, leglabels, self.average_method_hipe4ml,
                                              self.roc_method_hipe4ml)
            figname = f'{self.dirmlplot}/ROCCurveAll_pT_{self.p_binmin}_{self.p_binmax}.pdf'
            roccurvefig.savefig(figname)
            # _____________________________________________
            plt.rcParams["figure.figsize"] = (10, 9)
            roccurvettfig = plot_utils.plot_roc_train_test(self.traintestdata[3],
                                                           self.ypredtest_hipe4ml,
                                                           self.traintestdata[1],
                                                           self.ypredtrain_hipe4ml, None,
                                                           leglabels, self.average_method_hipe4ml,
                                                           self.roc_method_hipe4ml)
            figname = f'{self.dirmlplot}/ROCCurveTrainTest_pT_{self.p_binmin}_{self.p_binmax}.pdf'
            roccurvettfig.savefig(figname)
            # _____________________________________________
            precisionrecallfig = plot_utils.plot_precision_recall(self.traintestdata[3],
                                                                  self.ypredtest_hipe4ml,
                                                                  leglabels)
            figname = f'{self.dirmlplot}/PrecisionRecallAll_pT_{self.p_binmin}_{self.p_binmax}.pdf'
            precisionrecallfig.savefig(figname)
            # _____________________________________________
            plt.rcParams["figure.figsize"] = (12, 7)
            featuresimportancefig = plot_utils.plot_feature_imp(self.traintestdata[2][self.v_train],
                                                                self.traintestdata[3],
                                                                self.p_hipe4ml_model,
                                                                leglabels)
            n_plot = self.n_classes if self.n_classes > 2 else 1
            for i, figg in enumerate(featuresimportancefig):
                if i < n_plot:
                    figname = (f'{self.dirmlplot}/FeatureImportanceOpt{i}_'
                               f'pT_{self.p_binmin}_{self.p_binmax}.pdf')
                    featuresimportancefig[i].savefig(figname)
                else:
                    figname = (f'{self.dirmlplot}/FeatureImportanceAll_'
                               f'pT_{self.p_binmin}_{self.p_binmax}.pdf')
                    featuresimportancefig[i].savefig(figname)

    #pylint: disable=too-many-locals
    def do_significance(self, ptbin):
        self.logger.info("Doing significance optimisation")
        gROOT.SetBatch(True)
        gROOT.ProcessLine("gErrorIgnoreLevel = kWarning;")
        #first extract the number of data events in the ml sample
        self.df_evt_data = pickle.load(openfile(self.f_evt_data, 'rb'))
        self.p_nevtml = self.p_fracmerge[ptbin] * len(self.df_evt_data) / self.p_fracmerge[0]
        self.p_nevttot = self.p_nevtml
        if self.p_dofullevtmerge is True:
            self.df_evttotsample_data = pickle.load(openfile(self.f_evttotsample_data, 'rb'))
            self.p_nevttot = len(self.df_evttotsample_data)
        else:
            self.logger.warning("The total merged event dataframe was not merged for space limits")
        self.logger.debug("Number of data events used for ML: %d", self.p_nevtml)
        self.logger.debug("Total number of data events: %d", self.p_nevttot)
        #calculate acceptance correction. we use in this case all
        #the signal from the mc sample, without limiting to the n. signal
        #events used for training
        denacc = len(self.df_mcgen[self.df_mcgen["ismcprompt"] == 1])
        numacc = len(self.df_mc[self.df_mc["ismcprompt"] == 1])
        acc, acc_err = calc_eff(numacc, denacc)
        self.logger.debug("Acceptance: %.3e +/- %.3e", acc, acc_err)
        #calculation of the expected fonll signals
        ptmin = self.p_binmin
        ptmax = self.p_binmax
        delta_pt = ptmax - ptmin
        if self.is_fonll_from_root:
            df_fonll = TFile.Open(self.f_fonll)
            df_fonll_Lc = df_fonll.Get(self.p_fonllparticle+"pred_"+self.p_fonllband)
            prod_cross = df_fonll_Lc.Integral(ptmin*20, ptmax*20)* self.p_fragf * 1e-12 / delta_pt
            signal_yield = 2. * prod_cross * delta_pt * acc * self.p_taa \
                           / (self.p_sigmamb * self.p_fprompt)
            #now we plot the fonll expectation
            cFONLL = TCanvas('cFONLL', 'The FONLL expectation')
            df_fonll_Lc.GetXaxis().SetRangeUser(0, 16)
            df_fonll_Lc.Draw("")
            cFONLL.SaveAs("%s/FONLL_curve_%s.png" % (self.dirmlplot, self.s_suffix))
        else:
            df_fonll = pd.read_csv(self.f_fonll)
            df_fonll_in_pt = df_fonll.query('(pt >= @ptmin) and (pt < @ptmax)')[self.p_fonllband]
            prod_cross = df_fonll_in_pt.sum() * self.p_fragf * 1e-12 / delta_pt
            signal_yield = 2. * prod_cross * delta_pt * self.p_br * acc * self.p_taa \
                           / (self.p_sigmamb * self.p_fprompt)
            #now we plot the fonll expectation
            plt.figure(figsize=(20, 15))
            plt.subplot(111)
            plt.plot(df_fonll['pt'], df_fonll[self.p_fonllband] * self.p_fragf, linewidth=4.0)
            plt.xlabel('P_t [GeV/c]', fontsize=20)
            plt.ylabel('Cross Section [pb/GeV]', fontsize=20)
            plt.title("FONLL cross section " + self.p_case, fontsize=20)
            plt.semilogy()
            plt.savefig(f'{self.dirmlplot}/FONLL_curve_{self.s_suffix}.png')

        self.logger.debug("Expected signal yield: %.3e", signal_yield)
        signal_yield = self.p_raahp * signal_yield
        self.logger.debug("Expected signal yield x RAA hp: %.3e", signal_yield)

        df_data_sideband = self.df_data_signf.query(self.s_selbkgml)
        df_data_sideband = shuffle(df_data_sideband, random_state=self.rnd_shuffle)
        df_data_sideband = df_data_sideband.tail(round(len(df_data_sideband) * self.p_bkgfracopt))
        hmass = TH1F('hmass', '', self.p_num_bins, self.p_mass_fit_lim[0], self.p_mass_fit_lim[1])
        df_mc_signal = self.df_mc[self.df_mc["ismcsignal"] == 1]
        mass_array = df_mc_signal['inv_mass'].values
        for mass_value in np.nditer(mass_array):
            hmass.Fill(mass_value)

        gaus_fit = TF1("gaus_fit", "gaus", self.p_mass_fit_lim[0], self.p_mass_fit_lim[1])
        gaus_fit.SetParameters(0, hmass.Integral())
        gaus_fit.SetParameters(1, self.p_mass)
        gaus_fit.SetParameters(2, 0.02)
        self.logger.debug("To fit the signal a gaussian function is used")
        fitsucc = hmass.Fit("gaus_fit", "RQ")

        if int(fitsucc) != 0:
            self.logger.warning("Problem in signal peak fit")
            sigma = 0.

        sigma = gaus_fit.GetParameter(2)
        self.logger.debug("Mean of the gaussian: %.3e", gaus_fit.GetParameter(1))
        self.logger.debug("Sigma of the gaussian: %.3e", sigma)
        sig_region = [self.p_mass - 3 * sigma, self.p_mass + 3 * sigma]
        fig_signif_pevt = plt.figure(figsize=(20, 15))
        plt.xlabel('Threshold', fontsize=20)
        plt.ylabel(r'Significance Per Event ($3 \sigma$)', fontsize=20)
        plt.title("Significance Per Event vs Threshold", fontsize=20)
        plt.xticks(fontsize=18)
        plt.yticks(fontsize=18)
        fig_signif = plt.figure(figsize=(20, 15))
        plt.xlabel('Threshold', fontsize=20)
        plt.ylabel(r'Significance ($3 \sigma$)', fontsize=20)
        plt.title("Significance vs Threshold", fontsize=20)
        plt.xticks(fontsize=18)
        plt.yticks(fontsize=18)

        df_sig = self.df_mltest[self.df_mltest["ismcprompt"] == 1]

        df_sig = apply(self.p_mltype, ["xgboost"], [self.p_hipe4ml_origmodel],
                       df_sig, self.v_train, self.multiclass_labels)
        df_data_sideband = apply(self.p_mltype, ["xgboost"], [self.p_hipe4ml_origmodel],
                                 df_data_sideband, self.v_train, self.multiclass_labels)

        #hipe4ml only has xgboost, to keep same style as mlhep, use ["xgboost"]
        for name in ["xgboost"]:
            eff_array, eff_err_array, x_axis = calc_sigeff_steps(self.p_nstepsign, df_sig, name,
                                                                 self.multiclass_labels)
            bkg_array, bkg_err_array, _ = calc_bkg(df_data_sideband, name, self.p_nstepsign,
                                                   self.p_mass_fit_lim, self.p_bkg_func,
                                                   self.p_bin_width, sig_region, self.p_savefit,
                                                   self.dirmlplot, [self.p_binmin, self.p_binmax],
                                                   self.multiclass_labels)
            sig_array = [eff * signal_yield for eff in eff_array]
            sig_err_array = [eff_err * signal_yield for eff_err in eff_err_array]
            bkg_array = [bkg / (self.p_bkgfracopt * self.p_nevtml) for bkg in bkg_array]
            bkg_err_array = [bkg_err / (self.p_bkgfracopt * self.p_nevtml) \
                             for bkg_err in bkg_err_array]
            signif_array, signif_err_array = calc_signif(sig_array, sig_err_array, bkg_array,
                                                         bkg_err_array)

            if self.n_classes <= 2:
                plt.figure(fig_signif_pevt.number)
                plt.errorbar(x_axis, signif_array, yerr=signif_err_array, label=f'{name}',
                             elinewidth=2.5, linewidth=5.0)
                signif_array_ml = [sig * sqrt(self.p_nevtml) for sig in signif_array]
                signif_err_array_ml = [sig_err * sqrt(self.p_nevtml) for sig_err in signif_err_array]
                plt.figure(fig_signif.number)
                plt.errorbar(x_axis, signif_array_ml, yerr=signif_err_array_ml,
                             label=f'{name}_ML_dataset', elinewidth=2.5, linewidth=5.0)
                signif_array_tot = [sig * sqrt(self.p_nevttot) for sig in signif_array]
                signif_err_array_tot = [sig_err * sqrt(self.p_nevttot) for sig_err in signif_err_array]
                plt.figure(fig_signif.number)
                plt.errorbar(x_axis, signif_array_tot, yerr=signif_err_array_tot,
                             label=f'{name}_Tot', elinewidth=2.5, linewidth=5.0)
                plt.figure(fig_signif_pevt.number)
                plt.legend(loc="upper left", prop={'size': 30})
                plt.savefig(f'{self.dirmlplot}/Significance_PerEvent_{self.s_suffix}.png')
                plt.figure(fig_signif.number)
                plt.legend(loc="upper left", prop={'size': 30})
                plt.savefig(f'{self.dirmlplot}/Significance_{self.s_suffix}.png')
                with open(f'{self.dirmlplot}/Significance_{self.s_suffix}.pickle', 'wb') as out:
                    pickle.dump(fig_signif, out)
            else:
                hsignfvscut = TH2F(f'hsignfvscut_pT{self.p_binmin}_{self.p_binmax}',
                                   f';{self.multiclass_labels[0]};{self.multiclass_labels[1]};Expected Significance',
                                   len(x_axis) - 1, array("d", x_axis), len(x_axis) - 1, array("d", x_axis))
                for i0, thr0 in enumerate(x_axis):
                    for i1, thr1 in enumerate(x_axis):
                        binvar0 = hsignfvscut.GetXaxis().FindBin(thr0)
                        binvar1 = hsignfvscut.GetXaxis().FindBin(thr1)
                        hsignfvscut.SetBinContent(binvar0, binvar1, signif_array[i1 + i0 * len(x_axis)])

                csignf = TCanvas(f'csignf_pT{self.p_binmin}_{self.p_binmax}', '', 640, 540)
                hFrame = csignf.cd(1).DrawFrame(x_axis[0], x_axis[0], x_axis[-1], x_axis[-1],
                                                f';{self.multiclass_labels[0]};{self.multiclass_labels[1]};Expected Significance')
                hFrame.GetXaxis().SetNdivisions(505)
                hFrame.GetYaxis().SetNdivisions(505)
                hFrame.GetXaxis().SetDecimals()
                hFrame.GetYaxis().SetDecimals()
                hsignfvscut.DrawCopy('colzsame')
                csignf.SaveAs(f'{self.dirmlplot}/Significance2D_{self.s_suffix}.eps')
                fout = TFile(f'{self.dirmlplot}/Significance2D_{self.s_suffix}.root', "RECREATE")
                fout.cd()
                hsignfvscut.Write()
                fout.Close()
