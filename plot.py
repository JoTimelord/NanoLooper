#!/bin/env python

import plottery_wrapper as p

def plot(histname, xaxis_name):

    p.dump_plot(
        # reference frame (with lower value)
            fnames=["hadds/DYJETSbkg.root","hadds/ttdilep.root"],
            # fnames=["DYJETSbkg1.root"],
            # with larger value
            # data_fname="/home/users/joytzphysics/plots/VBSWWH_4p51.root",
            sig_fnames=["hadds/VBSWZH_C2V_3.root", "hadds/VBSOSWWH_C2V_3.root"],
            legend_labels=["DYJETS","ttdilep"],
            signal_labels=["VBSWZH C_{2V}=3", "VBSOSWWH C_{2V}=3"],
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
                 "nJets",
                 "ptDiLep",
                 "phiHbb",
                 "ptDijet",
                 "dRLep",
                 "dRjets",
                 "deltaEtajets",
                 "deltaEtaLep",
                 "dRVBF",
                 "dRVBF1FatJet",
                 "dRVBF2FatJet",
                 "hbbJetArea",
                 "mpRatioVBFjet1",
                 "mpRatioVBFjet2",
                 "tau21",
                 "tau31",
                 "n2b1",
                 "n3b1",
                 "h_cutflow",
                 "VBFjetMass",
                 "ptHbb",
                 "etaHbb",
                 "massHbb",
                 "wScore",
                 "hjetScore",
                 "zScore",
                 "ptMET",
                 "phiMET",
                 "EtMET",
                 "KT",
                 "ST",
                 "re_deltaEtaVBF",
                 "allN3b1",
                 )

    for name in hist_list:
        plot(name, "GeV/c^2")

