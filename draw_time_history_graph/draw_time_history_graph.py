"""
Todo
- auto-check data_stage from file name

"""
from pathlib import Path
import os
import numpy as np
import pandas as pd
import datetime
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import time
import gc
import concurrent.futures
from scipy import integrate

# magic code
mpl.rc('font',**{'family':'serif','serif':['Times New Roman']})
# mpl.rc('text', usetex=True)

# required parameters
data_stage = 1 # 0:precon 1:consolidation 2:liquefaction
LDT_installment = False
ClipGauge_installment = False

# variables
file_list_out = []
file_list_spe = []

# set data directory
dat_dir = r"F:\OneDrive\01_Research\01_Doctor\07_Triaxial_Test\07_Result\2020\20200106"
def main():
    # create result directory
    result_dir = Path(dat_dir) / "postprocessed result"
    Path(result_dir).mkdir(exist_ok=True)

    # listing all the .out file
    for all_file_list in os.listdir(dat_dir):
        base, ext = os.path.splitext(all_file_list)
        if ext == ".spe":
            file_list_spe.append(all_file_list)
        elif ext == ".out":
            file_list_out.append(all_file_list)
    
    # travarse file_list_out with iterator
    for file_name in file_list_out:
        # import file
        data_out = pd.read_csv(dat_dir + "\\" + file_name, delim_whitespace=True)

        ###########################################################################
        #                                 *CAUTION*                               #
        # Since some below lines are under revision, keep LDT_availability False  #
        # - don't consider ClipGauge_availability                                 #
        # - don't check which LDT parameter is used for calulation                #
        ###########################################################################
        
        # calculate parameters related to LDT and ClipGauge
        if LDT_installment == True:
            # calculate volume change based on clip gauge and LDT
            data_out["Specimen_Volume(mm3)"] = data_spe.iloc[3, data_stage + 2] * (1 - data_out["V-LDT2_(mm)"] / data_spe.iloc[8, data_stage + 2]) * (data_spe.iloc[0, data_stage + 2] + data_out["Clip-Gauge_(mm)"]) ** 2 * math.pi / 4

            # calculate volumetric strain 
            data_out["e_v-LDT-Clip(%)"] = - np.log(data_out["Specimen_Volume(mm3)"] / data_out["Specimen_Volume(mm3)"].iloc[0]) * 100

            # radial strain based on clip gauge
            data_out["e_r-Clip(%)"] = - np.log((data_out["Clip-Gauge_(mm)"] + data_spe.iloc[0, data_stage + 2]) / (data_out["Clip-Gauge_(mm)"].iloc[0] + data_spe.iloc[0, data_stage + 2])) * 100

            # baseline correlation on LDTs
            data_out["e_LDT2_base(%)"] = data_out["eLDT2(%)_"] - data_out["eLDT2(%)_"].iloc[0]

            # calcurate shear strain
            data_out["gamma(%)_based_LDT"] = (data_out["e_LDT2_base(%)"] - data_out["e_r-Clip(%)"]) / 2
        
        # calculate parameters under liquefaction
        if data_stage == 2:
            # baseline correlation on axial strain from EDT
            data_out["ea_base(%)"] = data_out["e(a)_(%)_"] - data_out["e(a)_(%)_"].iloc[0]

            # calcuration geometric mean effective stress
            data_out["Geometric_Mean_Eff._Str._(KPa)"] = data_out["s'(a)(kPa)"] ** (1/2) * data_out["s'(r)(kPa)"] ** (1/2) 

            # baseline correlation on radial strain from EDT and LCDPT
            data_out["gamma_based_EDT(%)"] = 3 * data_out["ea_base(%)"] / 4

            # baseline correlation on effective stress
            data_out["p'___(kPa)_base"] = data_out["p'___(kPa)"] - data_out["p'___(kPa)"].min()

            # calcurate normalized disspated energy
            data_out["Diss._Energy(kPa)"] = np.hstack(([0], integrate.cumtrapz(data_out["q____(kPa)"] / 2 / data_out["p'___(kPa)_base"] ,  data_out["gamma_based_EDT(%)"] / 100)))

        # draw time-dependent figures
        fig = plt.figure(figsize=(6, 10))
        fig.subplots_adjust(left=0.23, right=0.77, hspace=0)

        ax1 = fig.add_subplot(511)
        ax1.plot(data_out["Time(s)"], data_out["q____(kPa)"], color="k", linewidth = 0.75, label="Deviatoric Stress, " + r"$q$(kPa)")
        ax1.set_ylabel("Deviatoric Stress, " + r"$q$(kPa)",color="k",fontsize=10)
        plt.setp(ax1.get_xticklabels(), visible=False)
        ax1.tick_params(axis='x', which='both', length=0)
        ax1.grid(axis="x", linestyle=":", linewidth=0.5, color="k")
        plt.title(file_name.replace('_', '\_'))

        ax2 = fig.add_subplot(512)
        ax2.plot(data_out["Time(s)"], data_out["e(a)_(%)_"], color="b", linewidth = 0.75, label="Axial Strain, " + r"$\varepsilon_a$(%)")
        plt.setp(ax2.get_xticklabels(), visible=False)
        ax2.tick_params(axis='x', which='both', length=0)
        ax2.grid(axis="x", linestyle=":", linewidth=0.5, color="k")
        ax2.set_ylabel("Axial Strain, " + r"$\varepsilon_a$(\%)",color="b",fontsize=10)
        
        ax3 = fig.add_subplot(513)
        ax3.plot(data_out["Time(s)"], data_out["s(r)_(kPa)"], color="g", linewidth = 0.75, label="Total Stress, " + r"$\sigma_3$(kPa)")
        plt.setp(ax3.get_xticklabels(), visible=False)
        ax3.tick_params(axis='x', which='both', length=0)
        ax3.grid(axis="x", linestyle=":", linewidth=0.5, color="k")
        ax3.set_ylabel("Total Stress, " + r"$\sigma_3$(kPa)",color="g",fontsize=10)

        ax4 = fig.add_subplot(514)
        ax4.plot(data_out["Time(s)"], data_out["e(v)_(%)_"], color="r", linewidth = 0.75, label="Volumetric Strain, " + r"$\varepsilon_v$(%)")
        ax4.set_ylabel("Volumetric Strain, " + r"$\varepsilon_v$(\%)",color="r",fontsize=10)
        ax4.grid(axis="x", linestyle=":", linewidth=0.5, color="k")
        plt.setp(ax4.get_xticklabels(), visible=False)
        ax4.tick_params(axis='x', which='both', length=0)

        ax5 = fig.add_subplot(515)
        ax5.plot(data_out["Time(s)"], data_out["eLDT2(%)_"], color="c", linewidth = 0.75, label="LDT2 Strain, " + r"$\epsilon_{LDT2}$(%)")
        ax5.set_ylabel("LDT2 Strain, " + r"$\epsilon_{LDT2}$(%)",color="c",fontsize=10)
        ax5.grid(axis="x", linestyle=":", linewidth=0.5, color="k")
        ax5.set_xlabel("Ellapsed Time, " + r"$t$" + "(sec)")

        figure_path = result_dir / (Path(file_name).stem + ".png")
        plt.savefig(figure_path, dpi=300, bbox_inches='tight')

        # explicit garbage collection
        gc.collect

        


if __name__ == "__main__":
    main()