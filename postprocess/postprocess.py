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
from matplotlib import rc
from matplotlib import font_manager
import matplotlib.pyplot as plt
import time
import gc
import concurrent.futures
from scipy import integrate

# magic code
rc('font',**{'family':'serif','serif':['Times New Roman']})
# rc('text', usetex=True)

# required parameters
data_stage = 1 # 0:precon 1:consolidation 2:liquefaction
LDT_installment = False
ClipGauge_installment = False

# variables
file_list_out = []
file_list_spe = []

def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

# set data directory
dat_dir = r"F:\OneDrive\01_Research\01_Doctor\07_Triaxial_Test\07_Result\2020\20200108"
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
        fig = plt.figure(figsize=(8, 6))
        fig.subplots_adjust(left=0.23, right=0.77)
        
        ax1 = fig.add_subplot(111)
        ax1.plot(data_out["Time(s)"], data_out["q____(kPa)"], color="k", linewidth = 0.75, label="Deviatoric Stress, " + r"$q$(kPa)")
        ax1.set_ylabel("Deviatoric Stress, " + r"$q$(kPa)",color="k",fontsize=14)
        ax1.set_xlabel(r"Ellapsed Time, " + r"$t$(sec)",color="k",fontsize=14)
        ax1.yaxis.label.set_color("k")
        ax1.tick_params(axis='both', colors="k", direction='in')

        ax2=ax1.twinx()
        ax2.plot(data_out["Time(s)"], data_out["e(a)_(%)_"], color="b", linewidth = 0.75, label="Axial Strain, " + r"$\varepsilon_a$(%)")
        ax2.set_ylabel("Axial Strain, " + r"$\varepsilon_a$(\%)",color="b",fontsize=14)
        ax2.spines["right"].set_position(("axes", 1))
        ax2.yaxis.label.set_color("b")
        ax2.tick_params(axis='y', colors="b", direction='in')
        ax2.spines['right'].set_color("b")
        
        ax3=ax1.twinx()
        ax3.plot(data_out["Time(s)"], data_out["s(r)_(kPa)"], color="g", linewidth = 0.75, label="Total Stress, " + r"$\sigma_3$(kPa)")
        ax3.set_ylabel("Total Stress, " + r"$\sigma_3$(kPa)",color="g",fontsize=14)
        ax3.spines["right"].set_position(("axes", 1.2))
        make_patch_spines_invisible(ax3)
        ax3.spines["right"].set_visible(True)
        ax3.yaxis.label.set_color("g")
        ax3.tick_params(axis='y', colors="g", direction='in')
        ax3.spines['right'].set_color("g")

        ax4=ax1.twinx()
        ax4.plot(data_out["Time(s)"], data_out["e(v)_(%)_"], color="r", linewidth = 0.75, label="Volumetric Strain, " + r"$\varepsilon_v$(%)")
        ax4.set_ylabel("Volumetric Strain, " + r"$\varepsilon_v$(\%)",color="r",fontsize=14)
        ax4.spines["left"].set_position(("axes", -0.2))
        ax4.yaxis.set_label_position('left')
        ax4.yaxis.set_ticks_position('left')
        make_patch_spines_invisible(ax4)
        ax4.spines["left"].set_visible(True)
        ax4.yaxis.label.set_color("r")
        ax4.tick_params(axis='y', colors="r", direction='in')
        ax4.spines['left'].set_color("r")

        """
        ax5=ax1.twinx()
        ax5.plot(data_out["Time(s)"], data_out["eLDT2(%)_"], color="c", linewidth = 0.75, label="LDT2 Strain, " + r"$\epsilon_{LDT2}$(%)")
        ax5.set_ylabel("LDT2 Strain, " + r"$\epsilon_{LDT2}$(%)",color="c",fontsize=14)
        ax5.spines["left"].set_position(("axes", -0.3))
        ax5.yaxis.set_label_position('left')
        ax5.yaxis.set_ticks_position('left')
        make_patch_spines_invisible(ax5)
        ax5.spines["left"].set_visible(True)
        ax5.yaxis.label.set_color("c")
        ax5.tick_params(axis='y', colors="c", direction='in')
        ax5.spines['left'].set_color("c")
        """
        plt.title(file_name.replace('_', '\_'))
        
        figure_path = result_dir / (Path(file_name).stem + ".png")
        plt.savefig(figure_path, dpi=300)

        # explicit garbage collection
        gc.collect

        


if __name__ == "__main__":
    main()