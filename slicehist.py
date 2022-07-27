#!/usr/local/bin/python
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mpl_toolkits
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import os
import numpy as np
import shutil

input_path = "/home/users/joytzphysics/NanoLooper/outputs/"
output_path = "/home/users/joytzphysics/NanoLooper/plots/"
process_type = ["DYJETSbkg", "VBSOSWWH_C2V_3","VBSWZH_C2V_3", "WWdilep", "WWinclusive"]
var1_type = "hbb Score"
var2_type = "N3B1"

class Process:
    typename = "notnamed"
    filelist = []
    scorelist = []
    var1name = "notnamed"
    var2name = "notnamed"

    def __init__(self, t, var1, var2):
        self.typename = t
        self.var1name = var1
        self.var2name = var2

    # input name is the directory of the .dat output
    def listfile(self, inputname):
        for file in os.listdir(inputname):
            if file.endswith("dat") and file.startswith(self.typename):
                filepath = inputname + file
                self.filelist.append(filepath)

    def getXY(self):
        score = []
        for file in self.filelist:
            data = pd.read_fwf(file)
            df = data[[self.var1name, self.var2name]]
            score.append(df)
        self.scorelist = pd.concat(score)

    def plot3DHist(self, outputname):
        fig = plt.figure(figsize = (10,8))
        ax = fig.add_subplot(111, projection='3d')
        xpos = (self.scorelist[[self.var1name]].to_numpy()).T
        ypos = (self.scorelist[[self.var2name]].to_numpy()).T
        H, xedges, yedges = np.histogram2d(xpos[0], ypos[0], bins=20, range=[[0,1],[-5,5]])
       # H, xedges, yedges = np.histogram2d(xpos[0], ypos[0], bins=20)
        colorlist = ['gold', 'turquoise', 'orange', 'aqua', 'violet', 'skyblue'] * (len(yedges)-1)
        for i in (0, len(yedges)-2):
            cs = colorlist[i]
            ax.bar(xedges[0:-1], H[i], zs=yedges[i], zdir='y', color=cs, alpha=0.6, width=1.0/20.0)
        ax.set_xlabel(self.var1name)
        ax.set_ylabel(self.var2name)
        ax.set_zlabel('raw # of events')

        ax.set_title(self.typename)

        plt.savefig(outputname + self.typename + ".png")

# --------------------------------------------------------------------------------------------
if __name__ == "__main__":
    dir = output_path
    if os.path.exists(dir):
        shutil.rmtree(dir)
    os.makedirs(dir)

    for processname in process_type:
        pro = Process(processname, var1_type, var2_type)
        pro.listfile(input_path)
        pro.getXY()
        pro.plot3DHist(output_path)


