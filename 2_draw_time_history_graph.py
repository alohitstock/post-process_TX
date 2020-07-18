"""
Todo
- auto-check parameter.data_stage from file name

"""
from pathlib import Path
import os
import numpy as np
import pandas as pd
import datetime
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager 
import time
import gc
import shutil
import concurrent.futures
from scipy import integrate
import parameter


# magic code
mpl.rcParams['font.family'] = 'serif'
# mpl.rcParams['text.latex.preamble']=[r'\usepackage{amsmath}']
# mpl.rcParams['mathtext.rm'] = 'serif'

# variables
file_list_out = []
file_list_spe = []

def main():
    # listing all the .out file in the directory and sub-directory
    dat_dir = Path(parameter.dat_dir)

    for temp_file in dat_dir.glob("**/*.out"):
        # import file
        data_out = pd.read_csv(temp_file, sep='\t')

        ###########################################################################
        #                                 *CAUTION*                               #
        # Since some below lines are under revision, keep LDT_availability False  #
        # - don't consider ClipGauge_availability                                 #
        # - don't check which LDT parameter is used for calulation                #
        ###########################################################################
        
        # calculate parameters related to LDT and ClipGauge
        if parameter.LDT_installment == True:
            # calculate volume change based on clip gauge and LDT
            data_out["Specimen_Volume(mm3)"] = data_spe.iloc[3, parameter.data_stage + 2] * (1 - data_out["V-LDT2_(mm)"] / data_spe.iloc[8, parameter.data_stage + 2]) * (data_spe.iloc[0, parameter.data_stage + 2] + data_out["Clip-Gauge_(mm)"]) ** 2 * math.pi / 4

            # calculate volumetric strain 
            data_out["e_v-LDT-Clip(%)"] = - np.log(data_out["Specimen_Volume(mm3)"] / data_out["Specimen_Volume(mm3)"].iloc[0]) * 100

            # radial strain based on clip gauge
            data_out["e_r-Clip(%)"] = - np.log((data_out["Clip-Gauge_(mm)"] + data_spe.iloc[0, parameter.data_stage + 2]) / (data_out["Clip-Gauge_(mm)"].iloc[0] + data_spe.iloc[0, parameter.data_stage + 2])) * 100

            # baseline correlation on LDTs
            data_out["e_LDT2_base(%)"] = data_out["eLDT2(%)_"] - data_out["eLDT2(%)_"].iloc[0]

            # calcurate shear strain
            data_out["gamma(%)_based_LDT"] = (data_out["e_LDT2_base(%)"] - data_out["e_r-Clip(%)"]) / 2
        
        # compute time difference
        # outファイルから更新日時を取得
        record_end_time_out = os.path.getmtime(temp_file)

        # 経過時間の最大値の取得
        ellapsed_time_out = data_out["Time(s)"].max()

        # PC内時間とオシロスコープ内の時間の差を計算
        delta_pc_oscillo = datetime.datetime.strptime(parameter.time_oscillo, "%H:%M:%S") - datetime.datetime.strptime(parameter.time_pc, "%H:%M:%S")

        # オシロスコープの時刻を表す列の追加
        data_out["Time_oscillo"] = record_end_time_out + delta_pc_oscillo.total_seconds() - ellapsed_time_out + data_out["Time(s)"]

        # calculate parameters under liquefaction
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

        # calcurate excess pore water pressure
        data_out["Delta u(kPa)"] = data_out["p'___(kPa)"].iloc[0] - data_out["p'___(kPa)"]

        # draw time-dependent figures
        fig1 = plt.figure(figsize=(6, 14))
        fig1.subplots_adjust(left=0.23, right=0.77, hspace=0)

        ax1 = fig1.add_subplot(711)
        ax1.plot(data_out["Time(s)"], data_out["q____(kPa)"], color="k", linewidth = 0.75, label="Deviatoric Stress, " + r"$q$(kPa)")
        ax1.set_ylabel("Deviatoric Stress, " + r"$q$(kPa)",color="k",fontsize=8,usetex=True)
        plt.setp(ax1.get_xticklabels(), visible=False)
        for tick in ax1.get_xticklabels():
            tick.set_fontname("Times")
        ax1.tick_params(axis='x', which='both', length=0)
        ax1.grid(axis="x", linestyle=":", linewidth=0.5, color="k")
        plt.title(str(temp_file.stem))

        ax2 = fig1.add_subplot(712)
        ax2.plot(data_out["Time(s)"], data_out["e(a)_(%)_"], color="b", linewidth = 0.75, label="Axial Strain, " + r"$\varepsilon_a$(%)")
        plt.setp(ax2.get_xticklabels(), visible=False)
        ax2.tick_params(axis='x', which='both', length=0)
        ax2.grid(axis="x", linestyle=":", linewidth=0.5, color="k")
        ax2.set_ylabel("Axial Strain, " + r"$\varepsilon_a$(\%)",color="b",fontsize=8,usetex=True)
        
        ax3 = fig1.add_subplot(713)
        ax3.plot(data_out["Time(s)"], data_out["s(r)_(kPa)"], color="g", linewidth = 0.75, label="Total Stress, " + r"$\sigma_3$(kPa)")
        plt.setp(ax3.get_xticklabels(), visible=False)
        ax3.tick_params(axis='x', which='both', length=0)
        ax3.grid(axis="x", linestyle=":", linewidth=0.5, color="k")
        ax3.set_ylabel("Total Stress, " + r"$\sigma_3$(kPa)",color="g",fontsize=8,usetex=True)

        ax4 = fig1.add_subplot(714)
        ax4.plot(data_out["Time(s)"], data_out["e(v)_(%)_"], color="r", linewidth = 0.75, label="Volumetric Strain, " + r"$\varepsilon_v$(%)")
        ax4.set_ylabel("Volumetric Strain, " + r"$\varepsilon_v$(\%)",color="r",fontsize=8,usetex=True)
        ax4.grid(axis="x", linestyle=":", linewidth=0.5, color="k")
        plt.setp(ax4.get_xticklabels(), visible=False)
        ax4.tick_params(axis='x', which='both', length=0)

        ax5 = fig1.add_subplot(715)
        ax5.plot(data_out["Time(s)"], data_out["eLDT2(%)_"], color="c", linewidth = 0.75, label="LDT2 Strain, " + r"$\epsilon_{LDT2}$(%)")
        ax5.set_ylabel("LDT2 Strain, " + r"$\varepsilon_{LDT2}$(\%)",color="c",fontsize=8,usetex=True)
        ax5.grid(axis="x", linestyle=":", linewidth=0.5, color="k")
        plt.setp(ax5.get_xticklabels(), visible=False)
        ax5.tick_params(axis='x', which='both', length=0)

        if "ClipGauge(%)" in data_out.columns:
            ax6 = fig1.add_subplot(716)
            ax6.plot(data_out["Time(s)"], data_out["ClipGauge(%)"], color="y", linewidth = 0.75, label="Clip Gauge Strain, " + "\n" + r"$\epsilon_{Clip}$(%)")
            ax6.set_ylabel("Clip Gauge Strain, " + "\n" + r"$\epsilon_{Clip}$(%)",color="y",fontsize=8,usetex=True)
            ax6.grid(axis="x", linestyle=":", linewidth=0.5, color="k")
            plt.setp(ax6.get_xticklabels(), visible=False)
            ax6.tick_params(axis='x', which='both', length=0)

        ax7 = fig1.add_subplot(717)
        ax7.plot(data_out["Time(s)"], data_out["p'___(kPa)"], color="m", linewidth = 0.75, label="Effective Stress, " + r"$p'$(kPa)")
        ax7.set_ylabel("Effective Stress, " + r"$p'$(kPa)", color="m", fontsize=8,usetex=True)
        ax7.grid(axis="x", linestyle=":", linewidth=0.5, color="k")
        ax7.set_xlabel("Ellapsed Time, " + r"$t$" + "(sec)",usetex=True)

        figure_path = str(temp_file.parent) + "\\" + temp_file.stem + ".png"
        if os.path.isfile(figure_path):
            os.remove(figure_path)   
        fig1.savefig(figure_path, dpi=300, bbox_inches='tight')
        plt.close(fig1)

        # explicit garbage collection
        gc.collect

        print(temp_file.stem + " | liq folder : " + str("csr" in temp_file.stem or "CSR" in temp_file.stem))
        if "csr" in temp_file.stem or "CSR" in temp_file.stem:

            fig2 = plt.figure(figsize=(5,10),dpi=200)
            plt.subplots_adjust(hspace=0.3)

            ax1 = fig2.add_subplot(211)
            ax1.plot(data_out["ea_base(%)"], data_out["q____(kPa)"], linewidth = 0.5, color="k")
            ax1.set_xlabel("Axial Strain, " + r"$\varepsilon_a$(\%)",usetex=True)
            ax1.set_ylabel("Deviatoric Stress, " + r"$q$(kPa)",usetex=True)
            # ax1.set_xlim([-10, 3])
            # ax1.set_ylim([-60, 60])
            ax1.grid(linestyle=":", linewidth=0.5, color="k")

            ax2 = fig2.add_subplot(212)
            ax2.plot(data_out["p'___(kPa)"], data_out["q____(kPa)"], linewidth = 0.5, color="k")
            ax2.set_xlabel("Mean Effective Stress, " + r"$p'$(kPa)",usetex=True)
            ax2.set_ylabel("Deviatoric Stress, " + r"$q$(kPa)",usetex=True)
            # ax2.set_xlim([-5, 120])
            # ax2.set_ylim([-60, 60])
            ax2.grid(linestyle=":", linewidth=0.5, color="k")

            """
            # insert the test specification
            textstr = '\n'.join(["TXCUC_Biho_01", "OCR1", r'$\rho_d =1.492g/cm^3$', r"$V_s$"])
            props = dict(boxstyle='square', facecolor='w', alpha=0.5)
            ax1.text(0.05, 0.95, textstr, transform=ax1.transAxes, fontsize=8,
                verticalalignment='top', bbox=props, usetex=True)
            """

            figure_path = str(temp_file.parent) + "\\" + temp_file.stem + "_liq-result.png"
            if os.path.isfile(figure_path):
                os.remove(figure_path)
            fig2.savefig(figure_path, dpi=300, bbox_inches='tight')

            plt.close(fig2)
        
        csv_file_path = str(temp_file.parent) + "\\" + temp_file.stem + "_postprocessed.csv"
        if os.path.isfile(csv_file_path):
            os.remove(csv_file_path)
        data_out.to_csv(csv_file_path, index=False)

if __name__ == "__main__":
    main()