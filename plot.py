#!/bin/env python

import plottery_wrapper as p

def plot(histname, xaxis_name):

    p.dump_plot(
        # reference frame (with lower value)
            fnames=["hadds/DYJETSbkg.root","hadds/ttdilep.root", "hadds/WWdilep.root"],
            #fnames=["hadds/WWdilep.root","hadds/WWinclusive.root"],
            # with larger value
            # data_fname="/home/users/joytzphysics/plots/VBSWWH_4p51.root",
            #sig_fnames=["hadds/VBSOSWWH_C2V_4.root", "hadds/VBSWZH_C2V_4.root", "hadds/VBSZZH_C2V_4.root"],
            sig_fnames=["hadds/VBSOSWWH_C2V_3.root", "hadds/VBSWZH_C2V_3.root", "hadds/VBSZZH_C2V_3.root"],
            legend_labels=["DYJETSbkg","ttdilep", "WWdilep"],
            #legend_labels=["WWdilep","WWinclusive"],
            #signal_labels=["VBSOSWWH C_{2V}=4", "VBSWZH C_{2V}=4", "VBSZZH C_{2V}=4"],
            signal_labels=["VBSOSWWH C_{2V}=3", "VBSWZH C_{2V}=3", "VBSZZH C_{2V}=3"],
            #signal_labels=["VBSZZH C_{2V}=3", "VBSZZH C_{2V}=4"],
            filter_pattern=histname,
            dogrep=False,
            dirname="allcuts",
            extraoptions={
                "nbins": 30,
                "lumi_value": 137,
                #"ratio_range": [0., 100.],
                #"legend_datalabel": "VBSWWH C_{2V}=4p5",
                "legend_scalex":2.0,
                #"ratio_ndivisions":503,
                "xaxis_ndivisions":503,
                #"ratio_name": "C_{2V}=4p5 / C_{2V}=1",
                #"ratio_xaxis_title":xaxis_name,
                "print_yield":True,
                "yield_prec":4,
                #"signal_scale":10,
                #"yaxis_log": True,
                },
            skip2d=True,
    )

if __name__ == "__main__":

    hist_list = ("softdropmass",
                 "massDiLep",
                 "ptDiLep",
                 "deltaEtaLep",
                 "dRLep",
                 "nJets",
                 "dRjets",
                 "deltaEtajets",
                 "dRVBF",
                 "h_cutflow",
                 "VBFjetMass",
                 "re_deltaEtaVBF",
                 "ptHbb",
                 "massHbb",
                 "wScore",
                 "hjetScore",
                 "zScore",
                 "ptMET",
                 "EtMET",
                 "ST",
                 "KT"
                 )

    for name in hist_list:
        plot(name, "GeV/c^2")

