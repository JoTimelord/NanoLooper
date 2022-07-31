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
process_type = ["DYJETSbkg", "VBSOSWWH_C2V_4", "VBSOSWWH_C2V_3","VBSWZH_C2V_4","VBSWZH_C2V_3", "WWdilep", "WWinclusive", "ttdilep"]
#process_type = ["VBSOSWWH_C2V_4", "VBSWZH_C2V_4"]
var1_type = "maxhbbscore" # plotted on the x-axis
var2_type = "N3B1" # plotted on the y-axis
xlim = [0,1]
ylim = [-5, 5]
zratio = 'k'
bins = 10

class Process:
    typename = "notnamed"
    filelist = []
    scorelist = []
    var1name = "notnamed"
    var2name = "notnamed"
    events = []
    xbins = []
    ybins = []

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
            column_headers = list(data.columns.values)
            df = data[[self.var1name, self.var2name]]
            # df = data.iloc[:,[1,2]]
            df = df.apply(pd.to_numeric, args=('coerce',))
            score.append(df)
        self.scorelist = pd.concat(score)

    def binning(self, x_range, y_range, bin_no):
        xpos = (self.scorelist[[self.var1name]].to_numpy()).T
        ypos = (self.scorelist[[self.var2name]].to_numpy()).T
        #H, xedges, yedges = np.histogram2d(xpos[0], ypos[0], bins=bin_no, range=[x_range,y_range])
        H, xedges, yedges = np.histogram2d(xpos[0], ypos[0], bins=bin_no)
        self.xbins = xedges
        self.ybins = yedges
        self.events = H.T

    def plotHeatMap(self, outputname, ratio):
        if ratio == 'k':
            zscales = self.events/1000.0
        elif ratio == 'M':
            zscales = self.events/1000000.0
        else:
            zscales = self.events
        fig = plt.figure(figsize = (10,8))
        ax = fig.add_subplot(111)
        ax.imshow(self.events, interpolation='none', aspect = 'auto')

        xticks = np.empty(len(self.xbins)-1)
        yticks = np.empty(len(self.ybins)-1)
        for i in range(0, len(self.xbins)-1):
            xticks[i] = np.mean([self.xbins[i], self.xbins[i+1]])
            yticks[i] = np.mean([self.ybins[i], self.ybins[i+1]])

        xticks = np.around(xticks, decimals=2)
        yticks = np.around(yticks, decimals=2)

        # Show all ticks and label them with the respective list entries
        ax.set_xticks(np.arange(len(xticks)))
        ax.set_yticks(np.arange(len(yticks)))
        ax.set_xticklabels(xticks)
        ax.set_yticklabels(yticks)

        # Rotate the tick labels and set their alignment.
#        plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
#            rotation_mode="anchor")

        # Loop over data dimensions and create text annotations.
        for i in range(len(self.xbins)-1):
            for j in range(len(self.ybins)-1):
                text = ax.text(j, i, zscales[i, j],
                    ha="center", va="center", color="w")

        ax.set_title(self.typename+"("+zratio+" events)")
        plt.xlabel(self.var1name)
        plt.ylabel(self.var2name)
        fig.tight_layout()
        plt.savefig(outputname + self.typename + "_heatmap.png")

    def plot3DHist(self, outputname):
        fig = plt.figure(figsize = (10,8))
        ax = fig.add_subplot(111, projection='3d')
        colorlist = ['gold', 'turquoise', 'orange', 'aqua', 'violet', 'skyblue'] * (len(self.ybins)-1)

        for i in range(len(self.ybins)-1):
            cs = colorlist[i]
            ax.bar(self.xbins[0:-1], self.events[i], zs=self.ybins[i], zdir='y', color=cs, alpha=0.6, width=self.xbins[-1]-self.xbins[-2])
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
        pro.binning(xlim, ylim, bins)
        pro.plotHeatMap(output_path, zratio)
        pro.plot3DHist(output_path)



