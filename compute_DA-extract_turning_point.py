import pathlib
import pandas as pd
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

#required parameters
distance = 10
prominence = 1.5
n = 1 #smoothing

which_param = 1 # 0:ea, 1:q

def main():

    csv_path = pathlib.Path(r"F:\OneDrive\01_Research\01_Doctor\07_Triaxial_Test\07_Result\2020\20200122\20200122_15_liq_CSR015-OCR1-sc60\result\20200122_15_liq_CSR015-OCR1-sc60_postprocessed_simple_riv.csv")
    csv_dir = csv_path.parent
    if which_param == 0:
        export_csv_name = csv_path.stem + "_ea0.csv"
    elif which_param == 1:
        export_csv_name = csv_path.stem + "_q0.csv"
    export_csv_path = csv_dir / export_csv_name

    data_csv = pd.read_csv(csv_path)

    if which_param == 1:
        temp= np.abs(np.convolve(data_csv["q____(kPa)"], np.ones(n) / float(n), "same"))
    elif which_param == 0:
        temp= np.abs(np.convolve(data_csv["e(a)_(%)_"], np.ones(n) / float(n), "same"))

    max_id = signal.find_peaks(temp, distance=distance, prominence=prominence)
    min_id = signal.find_peaks(-temp, distance=distance, prominence=prominence)

    if which_param == 0:
        col_name = "e(a)_(%)_"
    elif which_param == 1:
        col_name = "q____(kPa)"
    
    plt.plot(data_csv["Time(s)"], data_csv[col_name])
    plt.plot(data_csv["Time(s)"].iloc[max_id[0]], data_csv[col_name].iloc[max_id[0]], "ro")
    plt.plot(data_csv["Time(s)"].iloc[min_id[0]], data_csv[col_name].iloc[min_id[0]], "bo")

    plt.show()
    if not "DA___(%)" in data_csv.columns:
        data_csv["DA___(%)"] = 0.0
    if not "Nc______" in data_csv.columns:
        data_csv["Nc______"] = 0.0

    ref_ea = data_csv["e(a)_(%)_"].iloc[0]
    ref_DA = 0.0
    temp_Nc = 0.0

    for index, row in data_csv.iterrows():
        temp_ea = row["e(a)_(%)_"]
        temp_DA = abs(temp_ea - ref_ea)
        if ref_DA < temp_DA:
            data_csv.at[index, "DA___(%)"] = temp_DA
            ref_DA = temp_DA
        else:
            data_csv.at[index, "DA___(%)"] = ref_DA
        temp_id = np.array(max_id[0]) - index
        if np.any(temp_id == 0):
            ref_ea = temp_ea
            temp_Nc += 0.5
        data_csv.at[index, "Nc______"] = temp_Nc
    
    plt.plot(data_csv["Nc______"], data_csv["e(a)_(%)_"], "o")
    plt.show()

    # search index of DA5%
    temp_DA = data_csv["DA___(%)"].values
    if np.where(temp_DA > 5)[0].size != 0:
        index_DA5 = np.where(temp_DA > 5)[0][0]
        data_csv["Nc_ratio"] = data_csv["Nc______"] / data_csv["Nc______"].iloc[index_DA5]
    else:
        data_csv["Nc_ratio"] = data_csv["Nc______"] / data_csv["Nc______"].max()
    
    
    # id = np.append(max_id[0], min_id[0])
    id = min_id[0]

    data_csv.to_csv(csv_path, index=False)
    data_csv_abs_q_peak = data_csv.iloc[id]
    data_csv_abs_q_peak.to_csv(export_csv_path, index=False)

if __name__ == "__main__":
    main()