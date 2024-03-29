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
main script for doing data processing, machine learning and analysis
"""
import multiprocessing as mp
import os
import pickle
import random as rd
import sys

import numpy as np
import pandas as pd
import uproot

from machine_learning_hep.utilities import create_folder_struc, openfile
from machine_learning_hep.utilities import list_folders, createlist, appendmainfoldertolist
from machine_learning_hep.utilities import merge_method, merge_method_max2
from machine_learning_hep.utilities_selection import filter_bit_df, tag_bit_df, selectdfrunlist
from machine_learning_hep.utilities_selection import selectfidacc, selectdfquery, seldf_singlevar
from machine_learning_hep.ml_functions import apply, apply_mlprob

class Processer: # pylint: disable=too-many-instance-attributes
    # Class Attribute
    species = 'processer'

    # Initializer / Instance Attributes
    # pylint: disable=too-many-statements, too-many-arguments
    def __init__(self, case, datap, run_param, mcordata, p_maxfiles,
                 d_root, d_pkl, d_pklsk, d_pkl_ml, p_period, i_period,
                 p_chunksizeunp, p_chunksizeskim, p_maxprocess,
                 p_frac_merge, p_rd_merge, d_pkl_dec, d_pkl_decmerged,
                 checkiffileexist, doapply):

        self.datap = datap
        self.case = case
        self.first_check_if_file_exists = checkiffileexist
        self.doapply = doapply

        #directories
        self.d_root = d_root
        self.d_pkl = d_pkl
        self.d_pklsk = d_pklsk
        self.d_pkl_ml = d_pkl_ml

        #processing variables
        self.mcordata = mcordata
        self.p_frac_merge = p_frac_merge
        self.p_rd_merge = p_rd_merge
        self.period = p_period
        self.select_children = datap["multi"][mcordata].get("select_children", None)
        if self.select_children:
            # Make sure we have "<child>/" instead if <child> only. Cause in the latter case
            # "child_1" might select further children like "child_11"
            self.select_children = [f"{child}/" for child in self.select_children[i_period]]
        self.runlist = run_param.get(self.period, None)
        self.p_maxfiles = p_maxfiles
        self.p_chunksizeunp = p_chunksizeunp
        self.p_chunksizeskim = p_chunksizeskim
        self.p_maxprocess = p_maxprocess
        self.p_dofullevtmerge = datap["dofullevtmerge"]
        self.v_max_ncand_merge = datap["multi"]["max_ncand_merge"]
        self.p_max_frac_merge = datap["multi"][self.mcordata]["max_frac_merge"]

        #namefile root
        self.n_root = datap["files_names"]["namefile_unmerged_tree"]
        #troot trees names
        mcordata_int = 0
        if self.mcordata == "mc":
            mcordata_int = 1
        self.n_treereco = datap["files_names"]["treeoriginreco"][mcordata_int]
        self.n_treegen = datap["files_names"]["treeorigingen"][mcordata_int]
        self.n_treeevt = datap["files_names"]["treeoriginevt"][mcordata_int]
        #namefiles pkl
        self.n_reco = datap["files_names"]["namefile_reco"]
        self.n_evt = datap["files_names"]["namefile_evt"]
        self.n_evtorig = datap["files_names"]["namefile_evtorig"]
        self.n_gen = datap["files_names"]["namefile_gen"]

        #selections
        self.s_reco_unp = datap["sel_reco_unp"]
        self.s_good_evt_unp = datap["sel_good_evt_unp"]
        self.s_cen_unp = datap["sel_cen_unp"]
        if isinstance(self.s_cen_unp, list):
            if self.mcordata == "data":
                self.s_cen_unp = datap["sel_cen_unp"][0]
            else:
                self.s_cen_unp = datap["sel_cen_unp"][1]
        self.s_gen_unp = datap["sel_gen_unp"]
        self.s_reco_skim = datap["sel_reco_skim"]
        self.s_gen_skim = datap["sel_gen_skim"]
        # mask missing values
        data_maskmissingvalues = datap.get("maskmissingvalues", None)
        if data_maskmissingvalues is not None:
            self.b_maskmissing = datap["maskmissingvalues"].get("activate", False)
            self.v_varstomask = datap["maskmissingvalues"].get("tomask", None)
        else:
            self.b_maskmissing = False
            self.v_varstomask = None

        #bitmap
        self.b_trackcuts = datap["sel_reco_singletrac_unp"]
        self.b_std = datap["bitmap_sel"]["isstd"]
        self.b_mcsig = datap["bitmap_sel"]["ismcsignal"]
        self.b_mcsigprompt = datap["bitmap_sel"]["ismcprompt"]
        self.b_mcsigfd = datap["bitmap_sel"]["ismcfd"]
        self.b_mcbkg = datap["bitmap_sel"]["ismcbkg"]
        self.b_mcrefl = datap["bitmap_sel"]["ismcrefl"]
        self.b_dsprompt = datap["bitmap_sel"].get("isdsprompt", None)
        self.b_dsfdbplus = datap["bitmap_sel"].get("isdsfdbplus", None)
        self.b_dsfdbzero = datap["bitmap_sel"].get("isdsfdbzero", None)
        self.b_dsfdlambdab = datap["bitmap_sel"].get("isdsfdlambdab", None)
        self.b_dsfdbs = datap["bitmap_sel"].get("isdsfdbs", None)

        #variables name
        self.lpt_anbinmin = datap["sel_skim_binmin"]
        self.lpt_anbinmax = datap["sel_skim_binmax"]
        self.lpt_anbinmintr = datap["ml"].get("binmin", self.lpt_anbinmin)
        self.lpt_anbinmaxtr = datap["ml"].get("binmax", self.lpt_anbinmax)
        if len(self.lpt_anbinmintr) > len(self.lpt_anbinmin) and self.doapply:
        #if self.doapply:
            self.lpt_anbinmin = datap["ml"]["binminsk"]
            self.lpt_anbinmax = datap["ml"]["binmaxsk"]
        self.p_nptbins = len(self.lpt_anbinmin)
        self.v_all = datap["variables"]["var_all"]
        self.v_train = datap["variables"]["var_training"]
        if not isinstance(self.v_train[0], list):
            self.v_train = [self.v_train for _ in range(self.p_nptbins)]
        self.v_evt = datap["variables"]["var_evt"][self.mcordata]
        self.v_gen = datap["variables"]["var_gen"]
        self.v_evtmatch = datap["variables"]["var_evt_match"]
        self.v_bitvar = datap["bitmap_sel"]["var_name"]
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
        self.v_var_binning = datap["var_binning"]
        self.nprongs = datap["nprongs"]
        self.prongformultsub = datap.get("prongformultsub", [0] * self.nprongs)

        #list of files names
        self.l_path = None
        if os.path.isdir(self.d_root):
            self.l_path = list_folders(self.d_root, self.n_root, self.p_maxfiles,
                                       self.select_children)
        else:
            self.l_path = list_folders(self.d_pkl, self.n_reco, self.p_maxfiles,
                                       self.select_children)

        self.l_root = createlist(self.d_root, self.l_path, self.n_root)
        self.l_reco = createlist(self.d_pkl, self.l_path, self.n_reco)
        self.l_evt = createlist(self.d_pkl, self.l_path, self.n_evt)
        self.l_evtorig = createlist(self.d_pkl, self.l_path, self.n_evtorig)
        if self.mcordata == "mc":
            self.l_gen = createlist(self.d_pkl, self.l_path, self.n_gen)

        self.f_totevt = os.path.join(self.d_pkl, self.n_evt)
        self.f_evt_ml = os.path.join(self.d_pkl_ml, self.n_evt)
        self.f_totevtorig = os.path.join(self.d_pkl, self.n_evtorig)
        self.f_evtorig_ml = os.path.join(self.d_pkl_ml, self.n_evtorig)

        self.lpt_recosk = [self.n_reco.replace(".pkl", "_%s%d_%d.pkl" % \
                          (self.v_var_binning, self.lpt_anbinmin[i], self.lpt_anbinmax[i])) \
                          for i in range(self.p_nptbins)]
        self.lpt_gensk = [self.n_gen.replace(".pkl", "_%s%d_%d.pkl" % \
                          (self.v_var_binning, self.lpt_anbinmin[i], self.lpt_anbinmax[i])) \
                          for i in range(self.p_nptbins)]
        self.lpt_reco_ml = [os.path.join(self.d_pkl_ml, self.lpt_recosk[ipt]) \
                             for ipt in range(self.p_nptbins)]
        self.lpt_gen_ml = [os.path.join(self.d_pkl_ml, self.lpt_gensk[ipt]) \
                            for ipt in range(self.p_nptbins)]
        self.lpt_reco_mlmax = [os.path.join(self.d_pkl_ml + '_max', self.lpt_recosk[ipt]) \
                              for ipt in range(self.p_nptbins)]

        self.mptfiles_recosk = []
        self.mptfiles_gensk = []
        self.mptfiles_recosk = [createlist(self.d_pklsk, self.l_path,
                                           self.lpt_recosk[ipt]) for ipt in range(self.p_nptbins)]
        if self.mcordata == "mc":
            self.mptfiles_gensk = [createlist(self.d_pklsk, self.l_path,
                                              self.lpt_gensk[ipt]) for ipt in range(self.p_nptbins)]

        #Variables for ML applying
        self.do_mlprefilter = datap.get("doml_asprefilter", None)
        self.apply_w_pkl_layout = datap.get("apply_with_pkl_layout", None)
        self.overwrite_mlprob_mc = None
        if self.mcordata == "mc":
            self.overwrite_mlprob_mc = datap.get("overwrite_mlprob_mc", None)
        self.p_modelname = datap["mlapplication"]["modelname"]
        self.mltype = datap["ml"]["mltype"]
        self.lpt_model = datap["mlapplication"]["modelsperptbin"]
        self.lpt_modhandler_hipe4ml = datap["mlapplication"]["modelsperptbin_hipe4ml"]
        self.dirmodel = datap["ml"]["mlout"]
        #if self.do_mlprefilter is True:
        #    self.dirmodel = self.dirmodel + "/prefilter"
        #    self.p_modelname = self.p_modelname + "prefilter"
        #if self.do_mlprefilter is False:
        #    self.dirmodel = self.dirmodel + "/analysis"
        self.lpt_model = appendmainfoldertolist(self.dirmodel, self.lpt_model)
        self.lpt_modhandler_hipe4ml = appendmainfoldertolist(self.dirmodel, self.lpt_modhandler_hipe4ml)

        self.doml = datap["doml"]
        if self.do_mlprefilter is not None and self.doml is False:
            print("FATAL error: The ML prefilter feature cannot combine with rectangular cuts")
        if not self.doml:
            datap["mlapplication"]["probcutpresel"][self.mcordata] = [0 for _ in self.lpt_anbinmin]
            datap["mlapplication"]["probcutoptimal"] = [0 for _ in self.lpt_anbinmin]
            datap["mlapplication"]["ml_prefilter_probcut"] = [0 for _ in self.lpt_anbinmin]
        self.lpt_probcutfin = datap["mlapplication"]["probcutoptimal"]
        if self.do_mlprefilter is True or self.overwrite_mlprob_mc is True:
            self.lpt_probcutpre = datap["mlapplication"]["ml_prefilter_probcut"]
        else:
            self.lpt_probcutpre = datap["mlapplication"]["probcutpresel"][self.mcordata]
            if not isinstance(self.lpt_probcutpre[0], list):
                if self.lpt_probcutfin < self.lpt_probcutpre:
                    print("FATAL error: probability cut final must be tighter!")

        self.d_pkl_dec = d_pkl_dec
        self.d_pkl_decmerged = d_pkl_decmerged
        inputdir_forapply = None
        #if self.do_mlprefilter is True:
        #    self.d_pkl_dec = d_pkl_dec + "/prefilter"
        #    self.d_pkl_decmerged = d_pkl_decmerged + "/prefilter"
        #if self.do_mlprefilter is False:
        #    inputdir_forapply = d_pkl_dec.replace("/analysis", "/prefilter")
        #    self.d_pkl_dec = d_pkl_dec + "/analysis"
        #    self.d_pkl_decmerged = d_pkl_decmerged + "/analysis"

        self.mptfiles_recosk_forapply = []
        if self.do_mlprefilter is False:
            self.mptfiles_recosk_forapply = [createlist(inputdir_forapply, self.l_path,
                                                        self.lpt_recosk[ipt]) for ipt in range(self.p_nptbins)]
        else:
            self.mptfiles_recosk_forapply = [createlist(self.d_pklsk, self.l_path,
                                                        self.lpt_recosk[ipt]) for ipt in range(self.p_nptbins)]
        self.lpt_recodec = None
        if self.doml is True:
            if self.do_mlprefilter is True or self.overwrite_mlprob_mc is True or self.apply_w_pkl_layout is True:
                self.lpt_recodec = self.lpt_recosk
            else:
                if not isinstance(self.lpt_probcutpre[0], list):
                    self.lpt_recodec = [self.n_reco.replace(".pkl", "%d_%d_%.2f.pkl" % \
                                       (self.lpt_anbinmintr[i], self.lpt_anbinmaxtr[i], \
                                        self.lpt_probcutpre[i])) for i in range(self.p_nptbins)]
                else:
                    self.lpt_recodec = [self.n_reco.replace(".pkl", "%d_%d_%.2f%.2f.pkl" % \
                                       (self.lpt_anbinmintr[i], self.lpt_anbinmaxtr[i], \
                                        self.lpt_probcutpre[i][0], self.lpt_probcutpre[i][1])) \
                                        for i in range(self.p_nptbins)]
        else:
            self.lpt_recodec = [self.n_reco.replace(".pkl", "%d_%d_std.pkl" % \
                               (self.lpt_anbinmintr[i], self.lpt_anbinmaxtr[i])) \
                                                    for i in range(self.p_nptbins)]
        self.mptfiles_recoskmldec = [createlist(self.d_pkl_dec, self.l_path, \
                                     self.lpt_recodec[ipt]) for ipt in range(self.p_nptbins)]
        self.lpt_recodecmerged = [os.path.join(self.d_pkl_decmerged, self.lpt_recodec[ipt])
                                  for ipt in range(self.p_nptbins)]
        if self.mcordata == "mc":
            self.mptfiles_genskmldec = [createlist(self.d_pkl_dec, self.l_path, \
                                        self.lpt_gensk[ipt]) for ipt in range(self.p_nptbins)]
            self.lpt_gendecmerged = [os.path.join(self.d_pkl_decmerged, self.lpt_gensk[ipt])
                                     for ipt in range(self.p_nptbins)]

        self.multiclass_labels = datap["ml"].get("multiclass_labels", None)

    def unpack(self, file_index):

        print(self.l_reco[file_index])
        file_exists = True
        if self.first_check_if_file_exists is True:
            file_exists = self.check_if_file_exists(self.l_reco[file_index])
            if self.mcordata == "mc" and file_exists is True:
                file_exists = self.check_if_file_exists(self.l_reco[file_index])
            if file_exists is True:
                file_exists = self.check_if_file_exists(self.l_evt[file_index])
            if file_exists is True:
                file_exists = self.check_if_file_exists(self.l_evtorig[file_index])

            #Probably this is not needed and slows it down much. To be checked once
            if file_exists is True and self.mcordata == "mc":
                try:
                    pickle.load(openfile(self.l_gen[file_index], "rb"))
                except Exception as e:
                    print("Failed to load file", self.l_gen[file_index])
                    print("Possible error", str(e))
                    file_exists = False
            if file_exists is True and self.mcordata == "data":
                try:
                    pickle.load(openfile(self.l_reco[file_index], "rb"))
                except Exception as e:
                    print("Failed to load file", self.l_reco[file_index])
                    print("Possible error", str(e))
                    file_exists = False

            if file_exists is True:
                return

        treeevtorig = uproot.open(self.l_root[file_index])[self.n_treeevt]
        try:
            dfevtorig = treeevtorig.pandas.df(branches=self.v_evt)
        except Exception as e: # pylint: disable=broad-except
            print('Missing variable in the event root tree', str(e))
            print('Missing variable in the candidate root tree')
            print('I am sorry, I am dying ...\n \n \n')
            sys.exit()

        dfevtorig = selectdfrunlist(dfevtorig, self.runlist, "run_number")
        dfevtorig = selectdfquery(dfevtorig, self.s_cen_unp)
        dfevtorig = dfevtorig.reset_index(drop=True)
        pickle.dump(dfevtorig, openfile(self.l_evtorig[file_index], "wb"), protocol=4)
        dfevt = selectdfquery(dfevtorig, self.s_good_evt_unp)
        dfevt = dfevt.reset_index(drop=True)
        pickle.dump(dfevt, openfile(self.l_evt[file_index], "wb"), protocol=4)

        treereco = uproot.open(self.l_root[file_index])[self.n_treereco]
        try:
            dfreco = treereco.pandas.df(branches=self.v_all)
        except Exception as e: # pylint: disable=broad-except
            print('Missing variable in the candidate root tree', str(e))
            print('I am sorry, I am dying ...\n \n \n')
            #try:
            #    self.v_all.remove("armenteros_K0s")
            #    print("Try with removing armenteros_K0s, can have too large values")
            #    dfreco = treereco.pandas.df(branches=self.v_all)
            #except Exception as e: # pylint: disable=broad-except
            #    print('Didnt work... Exiting', str(e))
            #    sys.exit()
            #print("Worked, adding default value for 0.20")
            #arr_arm = [0.20] * len(dfreco)
            #dfreco["armenteros_K0s"] = arr_arm
            #pass
            try:
                self.v_all.remove("its_refit_prong2")
                print("Try with removing its_refit_prong2, can be truncated")
                dfreco = treereco.pandas.df(branches=self.v_all)
            except Exception as e: # pylint: disable=broad-except
                print('Didnt work... Exiting', str(e))
                sys.exit()
            print("Worked, adding default value for 0")
            arr_itsrefit2 = [0.] * len(dfreco)
            dfreco["its_refit_prong2"] = arr_itsrefit2
            pass


        dfreco = selectdfrunlist(dfreco, self.runlist, "run_number")
        dfreco = selectdfquery(dfreco, self.s_reco_unp)
        dfreco = pd.merge(dfreco, dfevt, on=self.v_evtmatch)

        isselacc = selectfidacc(dfreco.pt_cand.values, dfreco.y_cand.values)
        dfreco = dfreco[np.array(isselacc, dtype=bool)]

        if self.b_trackcuts is not None:
            dfreco = filter_bit_df(dfreco, self.v_bitvar, self.b_trackcuts)

        arraysub = [0 for ival in range(len(dfreco))]
        for iprong in range(self.nprongs):
            if self.prongformultsub[iprong] == 0:
                continue
            #print("considering prong %d for sub" % iprong)
            spdhits_thisprong = dfreco["spdhits_prong%s" % iprong].values
            ntrackletsthisprong = [1 if spdhits_thisprong[index] == 3 else 0 \
                                   for index in range(len(dfreco))]
            arraysub = np.add(ntrackletsthisprong, arraysub)
        if "n_tracklets" in self.v_evt:
            n_tracklets = dfreco["n_tracklets"].values
            n_tracklets_sub = None
            n_tracklets_sub = np.subtract(n_tracklets, arraysub)
            dfreco["n_tracklets_sub"] = n_tracklets_sub
        if "n_tracklets_corr" in self.v_evt:
            n_tracklets_corr = dfreco["n_tracklets_corr"].values
            n_tracklets_corr_sub = None
            n_tracklets_corr_sub = np.subtract(n_tracklets_corr, arraysub)
            dfreco["n_tracklets_corr_sub"] = n_tracklets_corr_sub
        if "n_tracklets_corr_shm" in self.v_evt:
            n_tracklets_corr_shm = dfreco["n_tracklets_corr_shm"].values
            n_tracklets_corr_shm_sub = None
            n_tracklets_corr_shm_sub = np.subtract(n_tracklets_corr_shm, arraysub)
            dfreco["n_tracklets_corr_shm_sub"] = n_tracklets_corr_shm_sub

        dfreco[self.v_isstd] = np.array(tag_bit_df(dfreco, self.v_bitvar,
                                                   self.b_std), dtype=int)
        if self.mcordata == "mc":
            dfreco[self.v_ismcsignal] = np.array(tag_bit_df(dfreco, self.v_bitvar,
                                                            self.b_mcsig), dtype=int)
            dfreco[self.v_ismcprompt] = np.array(tag_bit_df(dfreco, self.v_bitvar,
                                                            self.b_mcsigprompt), dtype=int)
            dfreco[self.v_ismcfd] = np.array(tag_bit_df(dfreco, self.v_bitvar,
                                                        self.b_mcsigfd), dtype=int)
            dfreco[self.v_ismcbkg] = np.array(tag_bit_df(dfreco, self.v_bitvar,
                                                         self.b_mcbkg), dtype=int)
        if self.mcordata == "data" and self.case.startswith('Bs'):
            dfreco[self.v_ismcbkg] = np.array(tag_bit_df(dfreco, self.v_bitvar,
                                                         self.b_mcbkg), dtype=int)

            # Assuming if you set one, you set all of them
            if self.v_dsprompt is not None:
                dfreco[self.v_dsprompt] = np.array(tag_bit_df(dfreco, self.v_bitvar,
                                                              self.b_dsprompt), dtype=int)
                dfreco[self.v_dsfdbplus] = np.array(tag_bit_df(dfreco, self.v_bitvar,
                                                               self.b_dsfdbplus), dtype=int)
                dfreco[self.v_dsfdbzero] = np.array(tag_bit_df(dfreco, self.v_bitvar,
                                                               self.b_dsfdbzero), dtype=int)
                dfreco[self.v_dsfdlambdab] = np.array(tag_bit_df(dfreco, self.v_bitvar,
                                                                 self.b_dsfdlambdab), dtype=int)
                dfreco[self.v_dsfdbs] = np.array(tag_bit_df(dfreco, self.v_bitvar,
                                                            self.b_dsfdbs), dtype=int)

        dfreco = dfreco.reset_index(drop=True)
        pickle.dump(dfreco, openfile(self.l_reco[file_index], "wb"), protocol=4)

        if self.mcordata == "mc":
            treegen = uproot.open(self.l_root[file_index])[self.n_treegen]
            dfgen = treegen.pandas.df(branches=self.v_gen)
            dfgen = selectdfrunlist(dfgen, self.runlist, "run_number")
            dfgen = pd.merge(dfgen, dfevtorig, on=self.v_evtmatch)
            dfgen = selectdfquery(dfgen, self.s_gen_unp)
            dfgen[self.v_isstd] = np.array(tag_bit_df(dfgen, self.v_bitvar,
                                                      self.b_std), dtype=int)
            dfgen[self.v_ismcsignal] = np.array(tag_bit_df(dfgen, self.v_bitvar,
                                                           self.b_mcsig), dtype=int)
            dfgen[self.v_ismcprompt] = np.array(tag_bit_df(dfgen, self.v_bitvar,
                                                           self.b_mcsigprompt), dtype=int)
            dfgen[self.v_ismcfd] = np.array(tag_bit_df(dfgen, self.v_bitvar,
                                                       self.b_mcsigfd), dtype=int)
            dfgen[self.v_ismcbkg] = np.array(tag_bit_df(dfgen, self.v_bitvar,
                                                        self.b_mcbkg), dtype=int)
            dfgen = dfgen.reset_index(drop=True)
            pickle.dump(dfgen, openfile(self.l_gen[file_index], "wb"), protocol=4)

    def skim(self, file_index):
        try:
            dfreco = pickle.load(openfile(self.l_reco[file_index], "rb"))
        except Exception as e: # pylint: disable=broad-except
            print('failed to open file', self.l_reco[file_index], str(e))
            sys.exit()
        #To (hopefully) fix double signd0 issue with database when unpacking
        dfreco = dfreco.loc[:,~dfreco.columns.duplicated()]

        for ipt in range(self.p_nptbins):
            dfrecosk = seldf_singlevar(dfreco, self.v_var_binning,
                                       self.lpt_anbinmin[ipt], self.lpt_anbinmax[ipt])
            dfrecosk = selectdfquery(dfrecosk, self.s_reco_skim[ipt])
            dfrecosk = dfrecosk.reset_index(drop=True)
            f = openfile(self.mptfiles_recosk[ipt][file_index], "wb")
            pickle.dump(dfrecosk, f, protocol=4)
            f.close()
            if self.mcordata == "mc":
                try:
                    dfgen = pickle.load(openfile(self.l_gen[file_index], "rb"))
                except Exception as e: # pylint: disable=broad-except
                    print('failed to open MC file', self.l_gen[file_index], str(e))
                dfgensk = seldf_singlevar(dfgen, self.v_var_binning,
                                          self.lpt_anbinmin[ipt], self.lpt_anbinmax[ipt])
                dfgensk = selectdfquery(dfgensk, self.s_gen_skim[ipt])
                dfgensk = dfgensk.reset_index(drop=True)
                pickle.dump(dfgensk, openfile(self.mptfiles_gensk[ipt][file_index], "wb"),
                            protocol=4)

    def applymodel(self, file_index):
        for ipt in range(self.p_nptbins):
            if os.path.exists(self.mptfiles_recoskmldec[ipt][file_index]):
                if os.stat(self.mptfiles_recoskmldec[ipt][file_index]).st_size != 0:
                    continue
            dfrecosk = pickle.load(openfile(self.mptfiles_recosk_forapply[ipt][file_index], "rb"))
            if self.doml is True:
                if os.path.isfile(self.lpt_model[ipt]) is False:
                    print("Model file not present in bin %d" % ipt)
                if self.b_maskmissing:
                    dfrecosk = dfrecosk.replace(self.v_varstomask, value=np.nan)
                mod = pickle.load(openfile(self.lpt_model[ipt], 'rb'))
                dfrecoskml = apply("BinaryClassification", [self.p_modelname], [mod],
                                   dfrecosk, self.v_train[ipt], None)
                probvar = "y_test_prob" + self.p_modelname
                dfrecoskml = dfrecoskml.loc[dfrecoskml[probvar] > self.lpt_probcutpre[ipt]]
            else:
                dfrecoskml = dfrecosk.query("isstd == 1")
            pickle.dump(dfrecoskml, openfile(self.mptfiles_recoskmldec[ipt][file_index], "wb"),
                        protocol=4)
            if self.do_mlprefilter is True and self.mcordata == "mc":
                dfgensk = pickle.load(openfile(self.mptfiles_gensk[ipt][file_index], "rb"))
                pickle.dump(dfgensk, openfile(self.mptfiles_genskmldec[ipt][file_index], "wb"),
                            protocol=4)

    def applymodel_hipe4ml(self, file_index):
        for ipt in range(self.p_nptbins):
            if os.path.exists(self.mptfiles_recoskmldec[ipt][file_index]):
                if os.stat(self.mptfiles_recoskmldec[ipt][file_index]).st_size != 0:
                    continue
            print(self.mptfiles_recosk[ipt][file_index])
            dfrecosk = pickle.load(openfile(self.mptfiles_recosk[ipt][file_index], "rb"))
            dfrecosk = seldf_singlevar(dfrecosk, self.v_var_binning,
                                       self.lpt_anbinmintr[ipt], self.lpt_anbinmaxtr[ipt])
            if self.doml is True:
                if os.path.isfile(self.lpt_modhandler_hipe4ml[ipt]) is False:
                    print("hipe4ml model file not present in bin %d" % ipt)
                if self.b_maskmissing:
                    dfrecosk = dfrecosk.replace(self.v_varstomask, value=np.nan)
                modhandler = pickle.load(openfile(self.lpt_modhandler_hipe4ml[ipt], 'rb'))
                mod = modhandler.get_original_model()

                #njobsgrid = {'n_jobs': -1}
                #mod.set_params(**njobsgrid)

                if self.mltype == "MultiClassification":
                    dfrecoskml = apply("MultiClassification", [self.p_modelname], [mod],
                                       dfrecosk, self.v_train[ipt], self.multiclass_labels)
                    probvar0 = 'y_test_prob' + self.p_modelname + self.multiclass_labels[0]
                    probvar1 = 'y_test_prob' + self.p_modelname + self.multiclass_labels[1]
                    dfrecoskml = dfrecoskml.loc[(dfrecoskml[probvar0] <= self.lpt_probcutpre[ipt][0]) &
                                                (dfrecoskml[probvar1] >= self.lpt_probcutpre[ipt][1])]
                else:
                    dfrecoskml = apply("BinaryClassification", [self.p_modelname], [mod],
                                       dfrecosk, self.v_train[ipt], None)
                    probvar = "y_test_prob" + self.p_modelname
                    dfrecoskml = dfrecoskml.loc[dfrecoskml[probvar] > self.lpt_probcutpre[ipt]]
            else:
                dfrecoskml = dfrecosk.query("isstd == 1")
            pickle.dump(dfrecoskml, openfile(self.mptfiles_recoskmldec[ipt][file_index], "wb"),
                        protocol=4)
            if (self.do_mlprefilter is True or self.apply_w_pkl_layout is True) and self.mcordata == "mc":
                dfgensk = pickle.load(openfile(self.mptfiles_gensk[ipt][file_index], "rb"))
                pickle.dump(dfgensk, openfile(self.mptfiles_genskmldec[ipt][file_index], "wb"),
                            protocol=4)

    def applymodel_hipe4ml_mlprob(self, file_index):
        for ipt in range(self.p_nptbins):
            if os.path.exists(self.mptfiles_recoskmldec[ipt][file_index]):
                if os.stat(self.mptfiles_recoskmldec[ipt][file_index]).st_size != 0:
                    continue
            print(self.mptfiles_recosk[ipt][file_index])
            dfrecosk = pickle.load(openfile(self.mptfiles_recosk[ipt][file_index], "rb"))
            dfrecosk = seldf_singlevar(dfrecosk, self.v_var_binning,
                                       self.lpt_anbinmintr[ipt], self.lpt_anbinmaxtr[ipt])
            if self.doml is True:
                if os.path.isfile(self.lpt_modhandler_hipe4ml[ipt]) is False:
                    print("hipe4ml model file not present in bin %d" % ipt)
                if self.b_maskmissing:
                    dfrecosk = dfrecosk.replace(self.v_varstomask, value=np.nan)
                modhandler = pickle.load(openfile(self.lpt_modhandler_hipe4ml[ipt], 'rb'))
                mod = modhandler.get_original_model()
                if self.mltype == "MultiClassification":
                    dfrecoskml = apply_mlprob("MultiClassification", [self.p_modelname], [mod],
                                              dfrecosk, self.v_train[ipt], self.multiclass_labels)
                    probvar0 = 'ml_prob_' + self.multiclass_labels[0]
                    probvar1 = 'ml_prob_' + self.multiclass_labels[1]
                    dfrecoskml = dfrecoskml.loc[(dfrecoskml[probvar0] <= self.lpt_probcutpre[ipt][0]) &
                                                (dfrecoskml[probvar1] >= self.lpt_probcutpre[ipt][1])]
                else:
                    dfrecoskml = apply_mlprob("BinaryClassification", [self.p_modelname], [mod],
                                              dfrecosk, self.v_train[ipt], None)
                    probvar = "ml_prob"
                    dfrecoskml = dfrecoskml.loc[dfrecoskml[probvar] > self.lpt_probcutpre[ipt]]
            else:
                dfrecoskml = dfrecosk.query("isstd == 1")
            pickle.dump(dfrecoskml, openfile(self.mptfiles_recoskmldec[ipt][file_index], "wb"),
                        protocol=4)
            if self.overwrite_mlprob_mc is True and self.mcordata == "mc":
                dfgensk = pickle.load(openfile(self.mptfiles_gensk[ipt][file_index], "rb"))
                pickle.dump(dfgensk, openfile(self.mptfiles_genskmldec[ipt][file_index], "wb"),
                            protocol=4)

    @staticmethod
    def callback(ex):
        print(ex)

    def parallelizer(self, function, argument_list, maxperchunk):
        chunks = [argument_list[x:x+maxperchunk] \
                  for x in range(0, len(argument_list), maxperchunk)]
        for chunk in chunks:
            print("Processing new chunck size=", maxperchunk)
            pool = mp.Pool(self.p_maxprocess)
            _ = [pool.apply_async(function, args=chunk[i],
                                  error_callback=self.callback) for i in range(len(chunk))]
            pool.close()
            pool.join()

    def process_unpack_par(self):
        print("doing unpacking", self.mcordata, self.period)
        create_folder_struc(self.d_pkl, self.l_path)
        arguments = [(i,) for i in range(len(self.l_root))]
        self.parallelizer(self.unpack, arguments, self.p_chunksizeunp)

    def process_skim_par(self):
        print("doing skimming", self.mcordata, self.period)
        create_folder_struc(self.d_pklsk, self.l_path)
        arguments = [(i,) for i in range(len(self.l_reco))]
        self.parallelizer(self.skim, arguments, self.p_chunksizeskim)
        if self.p_dofullevtmerge is True:
            merge_method(self.l_evt, self.f_totevt)
            merge_method(self.l_evtorig, self.f_totevtorig)

    def process_applymodel_par(self):
        print("doing apply model", self.mcordata, self.period)
        create_folder_struc(self.d_pkl_dec, self.l_path)
        arguments = [(i,) for i in range(len(self.mptfiles_recosk[0]))]
        self.parallelizer(self.applymodel, arguments, self.p_chunksizeskim)

    def process_applymodel_hipe4ml_par(self):
        print("doing apply model", self.mcordata, self.period)
        create_folder_struc(self.d_pkl_dec, self.l_path)
        arguments = [(i,) for i in range(len(self.mptfiles_recosk[0]))]
        self.parallelizer(self.applymodel_hipe4ml, arguments, self.p_chunksizeskim)

    def process_applymodel_hipe4ml_mlprob_par(self):
        print("doing apply model", self.mcordata, self.period)
        create_folder_struc(self.d_pkl_dec, self.l_path)
        arguments = [(i,) for i in range(len(self.mptfiles_recosk[0]))]
        self.parallelizer(self.applymodel_hipe4ml_mlprob, arguments, self.p_chunksizeskim)

    def process_mergeforml(self):
        print("doing merging", self.mcordata, self.period)
        nfiles = len(self.mptfiles_recosk[0])
        if nfiles == 0:
            print("increase the fraction of merged files or the total number")
            print(" of files you process")
        for ipt in range(self.p_nptbins):
            ntomerge = (int)(nfiles * self.p_frac_merge[ipt])
            rd.seed(self.p_rd_merge) #make sure we start with the same files for each pT bin
            filesel = rd.sample(range(0, nfiles), ntomerge)
            list_sel_recosk = [self.mptfiles_recosk[ipt][j] for j in filesel]
            merge_method(list_sel_recosk, self.lpt_reco_ml[ipt])
            if self.mcordata == "mc":
                list_sel_gensk = [self.mptfiles_gensk[ipt][j] for j in filesel]
                merge_method(list_sel_gensk, self.lpt_gen_ml[ipt])
            #only merge event for first pT bin (smallest)
            if ipt == 0:
                list_sel_evt = [self.l_evt[j] for j in filesel]
                list_sel_evtorig = [self.l_evtorig[j] for j in filesel]
                merge_method(list_sel_evt, self.f_evt_ml)
                merge_method(list_sel_evtorig, self.f_evtorig_ml)

    def process_mergeforml_max(self):
        print("doing merging_max", self.mcordata, self.period)
        nfiles = len(self.mptfiles_recosk[0])
        if nfiles == 0:
            print("increase the fraction of merged files or the total number")
            print(" of files you process")
        rd.seed(self.p_rd_merge)
        for ipt in range(self.p_nptbins):
            ntomerge = (int)(nfiles * self.p_max_frac_merge[ipt])
            filesel = rd.sample(range(0, nfiles), ntomerge)
            list_sel_recosk = [self.mptfiles_recosk[ipt][j] for j in filesel]
            merge_method_max2(list_sel_recosk, self.lpt_reco_mlmax[ipt], self.v_max_ncand_merge)

    def process_mergedec(self):
        for ipt in range(self.p_nptbins):
            merge_method(self.mptfiles_recoskmldec[ipt], self.lpt_recodecmerged[ipt])
            if self.mcordata == "mc":
                merge_method(self.mptfiles_gensk[ipt], self.lpt_gendecmerged[ipt])

    def process_mergedec_max(self):
        print("doing merging_max", self.mcordata, self.period)
        nfiles = len(self.mptfiles_recoskmldec[0])
        if nfiles == 0:
            print("increase the fraction of merged files or the total number")
            print(" of files you process")
        rd.seed(self.p_rd_merge)
        for ipt in range(self.p_nptbins):
            ntomerge = (int)(nfiles * self.p_max_frac_merge[ipt])
            filesel = rd.sample(range(0, nfiles), ntomerge)
            list_sel_recoskmldec = [self.mptfiles_recoskmldec[ipt][j] for j in filesel]
            merge_method_max2(list_sel_recoskmldec, self.lpt_recodecmerged[ipt], self.v_max_ncand_merge)
            if self.mcordata == "mc":
                list_sel_gensk = [self.mptfiles_gensk[ipt][j] for j in filesel]
                merge_method_max2(list_sel_gensk, self.lpt_gendecmerged[ipt], self.v_max_ncand_merge)

    def check_if_file_exists(self, filename):
        if os.path.exists(filename):
            if os.stat(filename).st_size != 0:
                return True
        return False
