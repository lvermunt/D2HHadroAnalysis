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
import sys
import time
import matplotlib.pyplot as plt

import xgboost as xgb
from hipe4ml import plot_utils
from hipe4ml.model_handler import ModelHandler
from hipe4ml.tree_handler import TreeHandler
from hipe4ml.analysis_utils import train_test_generator

from machine_learning_hep.logger import get_logger


class Optimiserhipe4mltree:
    # Class Attribute
    species = "optimiser_hipe4mltree"

    def __init__(self, data_param, binmin, binmax, training_var, bkg_sel, hyper_pars):

        self.logger = get_logger()

        # directory
        #self.do_mlprefilter = datap.get("doml_asprefilter", None)
        self.dirmlout = data_param["ml"]["mlout"]
        self.dirmlplot = data_param["ml"]["mlplot"]
        #if self.do_mlprefilter is True:
        #    self.dirmodel = self.dirmodel + "/prefilter"
        #    self.dirmlplot = self.dirmlplot + "/prefilter"
        #if self.do_mlprefilter is False:
        #    self.dirmodel = self.dirmodel + "/analysis"
        #    self.dirmlplot = self.dirmlplot + "/analysis"

        self.inputtreedata = "/Users/lvermunt/cernbox/Analyses/ML/input/hipe4mlTTree/data.root"
        self.inputtreemc = "/Users/lvermunt/cernbox/Analyses/ML/input/hipe4mlTTree/prompt.root"
        self.v_train = None
        self.p_binmin = binmin
        self.p_binmax = binmax

        self.s_selsigml = ""
        self.s_selbkgml = bkg_sel #"inv_mass < 1.82 or 1.92 < inv_mass < 2.00"
        self.v_bkgoversigfrac = 3
        self.v_sig = 1
        self.v_bkg = 0
        self.rnd_splt = data_param["ml"]["rnd_splt"]
        self.test_frac = data_param["ml"]["test_frac"]

        self.prompthandler = None
        self.datahandler = None
        self.bkghandler = None
        self.traintestdata = None
        self.ypredtrain_hipe4ml = None
        self.ypredtest_hipe4ml = None

        self.preparesample()

        self.p_hipe4ml_model = None
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

        self.multiclass_labels = data_param["ml"].get("multiclass_labels", None)

        self.logger.info("Using the following training variables: %s", self.v_train)

    def preparesample(self):
        self.logger.info("Prepare Sample for hipe4ml")

        self.signalhandler = TreeHandler(self.inputtreemc, 'treeMLDplus')
        nsigcand = self.signalhandler.get_n_cand()
        self.datahandler = TreeHandler(self.inputtreedata, 'treeMLDplus')
        self.bkghandler = self.datahandler.get_subset(self.s_selbkgml, size=nsigcand * self.v_bkgoversigfrac)
        self.traintestdata = train_test_generator([self.signalhandler, self.bkghandler],
                                                  [self.v_sig, self.v_bkg],
                                                  test_size=self.test_frac,
                                                  random_state=self.rnd_splt)

    def load_hipe4mlmodel(self):
        self.logger.info("Loading hipe4ml model")
        self.v_train = self.signalhandler.get_var_names()
        self.v_train.remove('inv_mass')
        self.v_train.remove('pt_cand')

        model_xgboost = xgb.XGBClassifier()
        self.p_hipe4ml_model = ModelHandler(model_xgboost, self.v_train)

    def set_hipe4ml_modelpar(self):
        self.logger.info("Setting hipe4ml hyperparameters")
        self.p_hipe4ml_model.set_model_params(self.v_hipe4ml_pars)

    def do_hipe4mlhyperparopti(self):
        self.logger.info("Optimising hipe4ml hyperparameters (Bayesian)")

        if not (self.average_method_hipe4ml in ['macro', 'weighted'] and
                self.roc_method_hipe4ml in ['ovo', 'ovr']):
            self.logger.fatal("Selected ROC configuration is not valid!")

        if self.average_method_hipe4ml == 'weighted':
            metric = f'roc_auc_{self.roc_method_hipe4ml}_{self.average_method_hipe4ml}'
        else:
            metric = f'roc_auc_{self.roc_method_hipe4ml}'

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

        self.p_hipe4ml_model.train_test_model(self.traintestdata, self.average_method_hipe4ml,
                                              self.roc_method_hipe4ml)
        self.ypredtrain_hipe4ml = self.p_hipe4ml_model.predict(self.traintestdata[0],
                                                               self.raw_output_hipe4ml)
        self.ypredtest_hipe4ml = self.p_hipe4ml_model.predict(self.traintestdata[2],
                                                              self.raw_output_hipe4ml)

        modelhandlerfile = f'{self.dirmlout}/ModelHandler_pT_{self.p_binmin}_{self.p_binmax}.pkl'
        self.p_hipe4ml_model.dump_model_handler(modelhandlerfile)
        modelfile = f'{self.dirmlout}/ModelHandler_pT_{self.p_binmin}_{self.p_binmax}.model'
        self.p_hipe4ml_model.dump_original_model(modelfile)

        self.logger.info("Training + testing hipe4ml: Done!")
        self.logger.info("Time elapsed = %.3f", time.time() - t0)

    def do_hipe4mlplot(self):
        self.logger.info("Plotting hipe4ml model")

        leglabels = ["Background", "Prompt signal"]
        outputlabels = ["Bkg", "SigPrompt"]

        # _____________________________________________
        plot_utils.plot_distr([self.bkghandler, self.signalhandler], self.v_train, 100, leglabels)
        plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
        figname = f'{self.dirmlplot}/DistributionsAll_pT_{self.p_binmin}_{self.p_binmax}.pdf'
        plt.savefig(figname)
        plt.close('all')
        # _____________________________________________
        corrmatrixfig = plot_utils.plot_corr([self.bkghandler, self.signalhandler],
                                             self.v_train, leglabels)
        for figg, labb in zip(corrmatrixfig, outputlabels):
            plt.figure(figg.number)
            plt.subplots_adjust(left=0.2, bottom=0.25, right=0.95, top=0.9)
            figname = f'{self.dirmlplot}/CorrMatrix{labb}_pT_{self.p_binmin}_{self.p_binmax}.pdf'
            figg.savefig(figname)
        # _____________________________________________
        plt.rcParams["figure.figsize"] = (10, 7)
        mloutputfig = plot_utils.plot_output_train_test(self.p_hipe4ml_model, self.traintestdata,
                                                        80, self.raw_output_hipe4ml,
                                                        leglabels, self.train_test_log_hipe4ml,
                                                        density=True)
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
        for i in range(0,len(featuresimportancefig)):
            figname = (f'{self.dirmlplot}/FeatureImportanceOpt{i}_'
                       f'pT_{self.p_binmin}_{self.p_binmax}.pdf')
            featuresimportancefig[i].savefig(figname)

