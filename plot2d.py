#!/bin/env python

import plottery_wrapper as p

def plot(histname):

    p.dump_plot(
        # reference frame (with lower value)
            fnames=["hadds/DYJETSbkg.root","hadds/WWdilep.root"],
            # fnames=["DYJETSbkg1.root"],
            # with larger value
            # data_fname="/home/users/joytzphysics/plots/VBSWWH_4p51.root",
            sig_fnames=["hadds/VBSOSWWH_C2V_3.root", "hadds/VBSWZH_C2V_3.root"],
            legend_labels=["DYJETS","WWdilep"],
            # signal_labels=["VBSOSWWH C_{2V}=3", "VBSWZH C_{2V}=3"],
            filter_pattern=histname,
            dogrep=False,
            dirname="allcuts",
            extraoptions={
                "bin_text_smart": True,
                "zaxis_noexponents": True,
                "xaxis_label": "hbb Score",
                "yaxis_label": "n3b1 Score",
                "zaxis_log": True,
                },
            skip2d=False,
            )

if __name__ == "__main__":
    hist_list = ("h2_hbb_n3b1",
                 "h2_z_n3b1",
                )

    for name in hist_list:
        plot(name)

