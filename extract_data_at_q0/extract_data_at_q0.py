import pathlib
import pandas as pd
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

def main():
    csv_path = pathlib.Path(r"F:\OneDrive\01_Research\01_Doctor\07_Triaxial_Test\07_Result\2019\20190723_a\liq_result\20190723_06_sc30-csr0.25-ocr5_riv.csv")
    csv_dir = csv_path.parent
    export_csv_name = csv_path.stem + "_q0.csv"
    export_csv_path = csv_dir / export_csv_name

    data_csv = pd.read_csv(csv_path)

    n = 100

    temp= np.abs(np.convolve(data_csv["q____(kPa)"], np.ones(n) / float(n), "same"))

    max_id = signal.argrelmax(temp, order=500)
    min_id = signal.argrelmin(temp, order=500)

    plt.plot(data_csv["Time(s)"], data_csv["q____(kPa)"])
    plt.plot(data_csv["Time(s)"].iloc[max_id], data_csv["q____(kPa)"].iloc[max_id], "ro")
    plt.plot(data_csv["Time(s)"].iloc[min_id], data_csv["q____(kPa)"].iloc[min_id], "bo")

    plt.show()

    id = np.append(max_id, min_id)

    data_csv_abs_q_peak = data_csv.iloc[id]
    data_csv_abs_q_peak.to_csv(export_csv_path)

if __name__ == "__main__":
    main()