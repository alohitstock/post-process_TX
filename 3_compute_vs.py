"""
せん断波速度解析用プログラム

・必要なファイル
txt(せん断波速度のファイル)、spe(供試体データ)、out、dat

・ファイル構造
-[result]
--[vs_result]
---.png
--_riv.csv
-[parameter.dat_dir]
--[vs]
---.txt
--.spe
--.out
--.dat

・設定するべきパラメータ
parameter.accelerometer_distance : 圧密終了時における加速度計間距離
parameter.threshold_time : 電圧の基線計算に用いるデータの最大時刻
parameter.threshold_standard_deviation : 電圧の基線補正において偏差がこれ以内であれば解析続行
parameter.magnification_standard_deviation : せん断波が到達したとみなす閾値（基線の標準偏差の何倍かの電圧になった時点で波が到達したとする)
parameter.dat_dir : 解析データが入っているディレクトリ(フォルダ)

注意点
・変数のところは毎回正しい値を入力する！
・txt(せん断波速度のファイル)、spe(供試体データ)、out、datファイルを必ず入れること
・せん断波の読み取りはCloss-Closs法

追加したい機能
・次の計算方法の導入
    4. frequency domain
・液状化試験時と圧密試験時のspeファイルの読み取り値の変更(speファイルの整合性チェック、圧密後の修正密度計算)
・繰り返し回数、DAの計算
・入力波形で正規化
・累積損失エネルギー、1サイクルごとの損失エネルギーについての評価
"""

# Moduleのインポート
from pathlib import Path
import os
import numpy as np
import pandas as pd
import datetime
import math
import matplotlib.pyplot as plt
import time
import gc
import concurrent.futures
from scipy import integrate
import parameter

# function for txt(wave file) analysis
def vs_analysis(i, cur_dir, file_list_txt, data_csv, data_spe, result_png_dir):

    # Local Variables
    pass_flag = True
    measured_time_index = 0
    closs_shear_wave_velocity = 0.0
    rise_shear_wave_velocity = 0.0
    peak_shear_wave_velocity = 0.0
    corr_shear_wave_velocity = 0.0

    # 波形データの読み込み
    data_txt = pd.read_csv(file_list_txt[i], delimiter=",", header=None, names=[i for i in range(4)])

    # 波形データの抜き出し
    data_txt_wave = data_txt.iloc[9:, [0, 2, 3]].to_numpy().astype(np.float)

    # 時刻0未満のデータを抜き出し
    data_txt_wave_basis = data_txt_wave[np.where(data_txt_wave[:, 0] < parameter.threshold_time), 1:]

    # 時刻0未満のデータの電圧値の平均と偏差を計算する
    base_voltage = np.average(data_txt_wave_basis, axis=1)
    std_voltage = np.std(data_txt_wave_basis, axis=1)
    if std_voltage[0][0] > parameter.threshold_standard_deviation or std_voltage[0][1] > parameter.threshold_standard_deviation:
        # print(std_voltage)
        print(str(i + 1) + "/" + str(len(file_list_txt)) + " | File Name: " + file_list_txt[i].stem + "  Vs : N/A (too large std)")
        pass_flag = False
        return [pass_flag, measured_time_index, closs_shear_wave_velocity, rise_shear_wave_velocity, peak_shear_wave_velocity, corr_shear_wave_velocity]

    # 基線補正(個々のデータから平均値を引く)をかける
    data_txt_wave[:, 1:] -= base_voltage

    # 波形データの中で波の到達が検知された最初のインデックスを検索
    if parameter.detection_method_for_triggered_point == 0:
        # 基線補正後の波形データで最大値を検索
        maximum_value = np.max(np.abs(data_txt_wave), axis=0)
        trigger_time_index_1 = np.where(np.abs(data_txt_wave[:, 1]) > parameter.magnification_for_maximum_value * maximum_value[1])[0]
        trigger_time_index_2 = np.where(np.abs(data_txt_wave[:, 2]) > parameter.magnification_for_maximum_value * maximum_value[2])[0]
    elif parameter.detection_method_for_triggered_point == 1:
        trigger_time_index_1 = np.where(np.abs(data_txt_wave[:, 1]) > parameter.magnification_standard_deviation * std_voltage[0][0])[0]
        trigger_time_index_2 = np.where(np.abs(data_txt_wave[:, 2]) > parameter.magnification_standard_deviation * std_voltage[0][1])[0]

    if trigger_time_index_1.size == 0 or trigger_time_index_2.size == 0:
        print(str(i + 1) + "/" + str(len(file_list_txt)) + " | File Name: " + file_list_txt[i].stem + " | Vs : N/A (too small amp.)")
        pass_flag = False
        return [pass_flag, measured_time_index, closs_shear_wave_velocity, rise_shear_wave_velocity, peak_shear_wave_velocity, corr_shear_wave_velocity]

    trigger_time = np.array([trigger_time_index_1[0], trigger_time_index_2[0]])
    rise_transmission_time = abs(data_txt_wave[trigger_time[0], 0] - data_txt_wave[trigger_time[1], 0])

    # 上で検索したインデックスを基準として増やしていったときに、波形の正負が入れ替わるインデックスを検索
    clossing_time_index_1_negative = np.where(data_txt_wave[trigger_time[0]:, 1] < 0)[0]
    clossing_time_index_1_positive = np.where(data_txt_wave[trigger_time[0]:, 1] > 0)[0]
    clossing_time_index_2_negative = np.where(data_txt_wave[trigger_time[1]:, 2] < 0)[0]
    clossing_time_index_2_positive = np.where(data_txt_wave[trigger_time[1]:, 2] > 0)[0]

    if clossing_time_index_1_negative.size == 0 or clossing_time_index_1_positive.size == 0 or clossing_time_index_2_negative.size == 0 or clossing_time_index_2_positive.size == 0:
        print(str(i + 1) + "/" + str(len(file_list_txt)) + " | File Name: " + file_list_txt[i].stem + " | Vs : N/A (too large offset)")
        pass_flag = False
        return [pass_flag, measured_time_index, closs_shear_wave_velocity, rise_shear_wave_velocity, peak_shear_wave_velocity, corr_shear_wave_velocity]

    clossing_time = np.array([max(clossing_time_index_1_negative[0], clossing_time_index_1_positive[0]), max(clossing_time_index_2_negative[0], clossing_time_index_2_positive[0])])

    # clossing_timeのインデックスが最初から何番目からかを計算する
    clossing_time += trigger_time

    # clossing_timeのインデックスからせん断波の伝播時間を計算
    closs_transmission_time = abs(data_txt_wave[clossing_time[0], 0] - data_txt_wave[clossing_time[1], 0])

    # peak-peak法での計算
    # clossing_timeとtrigger_timeの中間で最も値の大きい地点のインデックスを検索する
    
    peak_time_index = np.array([np.argmax(np.abs(data_txt_wave[trigger_time[0]:clossing_time[0], 1])), np.argmax(np.abs(data_txt_wave[trigger_time[1]:clossing_time[1], 2]))])
    peak_time_index += trigger_time

    # peak-peak法でのせん断波の伝播時間を計算
    peak_transmission_time = abs(data_txt_wave[peak_time_index[0], 0] - data_txt_wave[peak_time_index[1], 0])

    ########################################
    #       Closs-corelation method        #
    ########################################
    # 相互相関関数の導出
    corr = np.correlate(data_txt_wave[:, 2], data_txt_wave[:, 1], "full")

    # ピーク値のインデックスの検索
    corr_xmin = -10000 * (data_txt_wave[1, 0] - data_txt_wave[0, 0]) 
    corr_xmax = -corr_xmin
    corr_x = np.linspace(corr_xmin, corr_xmax, 20001)
    dif_corr_x = np.abs(corr_x - closs_transmission_time)
    zero_index = math.floor((np.argmin(dif_corr_x) + dif_corr_x.shape[0] / 2) / 2)
    peak_corr_index = np.argmax(corr[zero_index:]) + zero_index

    # 遅延時間の計算
    corr_transmission_time = corr_x[peak_corr_index]

    # 波形データの取得時間を読み込み
    wave_record_time = data_txt.iloc[1, 1] + " " + data_txt.iloc[2, 1]
    wave_record_time = datetime.datetime.strptime(wave_record_time, "%m-%d-%Y %H:%M:%S.%f").timestamp()

    # 波形データに最も近い時刻のインデックスをoutファイルから検索
    time_dif = abs(data_csv["Time_oscillo"] - wave_record_time)
    measured_time_index = time_dif.idxmin()

    if parameter.axial_strain_parameter == 1:
        # 検索したインデックスにおけるeLDT2の変位量から、加速度計間距離を計算
        parameter.accelerometer_distance = parameter.accelerometer_distance * (1 - data_csv["V-LDT2_(mm)"].iloc[measured_time_index] / data_spe.iloc[8, parameter.data_stage + 2])
    else:
        parameter.accelerometer_distance = parameter.accelerometer_distance * math.exp(-1 * data_csv["e(a)_(%)_"].iloc[measured_time_index] / 100)

    # せん断波速度を計算
    closs_shear_wave_velocity = parameter.accelerometer_distance / closs_transmission_time / 1000
    rise_shear_wave_velocity = parameter.accelerometer_distance / rise_transmission_time / 1000
    peak_shear_wave_velocity = parameter.accelerometer_distance / peak_transmission_time / 1000
    corr_shear_wave_velocity = parameter.accelerometer_distance / corr_transmission_time / 1000 if corr_transmission_time > 0 else 0.0
    print(str(i + 1) + "/" + str(len(file_list_txt)) + " | File Name: " + file_list_txt[i].stem + " | Vs : " + "{:.1f}".format(closs_shear_wave_velocity))
    # フーリエ変換
    data_txt_wave_FFT = np.zeros((data_txt.iloc[9:, :].shape[0] * 10, data_txt.iloc[9:, :].shape[1]))
    data_txt_wave_FFT[:data_txt_wave.shape[0], :] = data_txt.iloc[9:, :].to_numpy().astype(np.float)
    fft = [np.abs(np.fft.fft(data_txt_wave_FFT[:, 1])) / (data_txt_wave_FFT.shape[0] / 2), np.abs(np.fft.fft(data_txt_wave_FFT[:, 2])) / (data_txt_wave_FFT.shape[0] / 2), np.abs(np.fft.fft(data_txt_wave_FFT[:, 3])) / (data_txt_wave_FFT.shape[0] / 2)]

    # 波形データをグラフで表示
    fig = plt.figure(figsize=(12, 18))
    plt.subplots_adjust(hspace=0.3, wspace=0.3)

    ax1 = fig.add_subplot(321)
    ax1.plot(data_txt_wave[:, 0], data_txt_wave[:, 1], color="b")
    ax1.plot(data_txt_wave[:, 0], data_txt_wave[:, 2], color="r")
    ax1.plot(data_txt_wave[trigger_time[0], 0], data_txt_wave[trigger_time[0], 1], marker="*", markersize=10, color="w", markerfacecolor="b")
    ax1.plot(data_txt_wave[trigger_time[1], 0], data_txt_wave[trigger_time[1], 2], marker="*", markersize=10, color="w", markerfacecolor="r")
    ax1.plot(data_txt_wave[clossing_time[0], 0], data_txt_wave[clossing_time[0], 1], marker=".", markersize=10, color="w", markerfacecolor="b")
    ax1.plot(data_txt_wave[clossing_time[1], 0], data_txt_wave[clossing_time[1], 2], marker=".", markersize=10, color="w", markerfacecolor="r")
    ax1.plot(data_txt_wave[peak_time_index[0], 0], data_txt_wave[peak_time_index[0], 1], marker="s", markersize=5, color="w", markerfacecolor="b")
    ax1.plot(data_txt_wave[peak_time_index[1], 0], data_txt_wave[peak_time_index[1], 2], marker="s", markersize=5, color="w", markerfacecolor="r")
    ax1.hlines(0, xmin=0, xmax=data_txt_wave[:, 0].max())

    ax1.set_xlabel("Time (sec)", fontname="Arial")

    ax1.set_ylim(-0.15, 0.15)
    ax1.set_ylabel("Volatage (V)", fontname="Arial")

    # 文字を追加
    analysis_result_label_closs_vs = "Closs Vs : " + '{:.1f}'.format(closs_shear_wave_velocity) + "m/s"
    analysis_result_label_rise_vs = "Rise Vs : " + '{:.1f}'.format(rise_shear_wave_velocity) + "m/s"
    analysis_result_label_peak_vs = "Peak Vs : " + '{:.1f}'.format(peak_shear_wave_velocity) + "m/s"
    analysis_result_label_corr_vs = "Corr Vs : " + '{:.1f}'.format(corr_shear_wave_velocity) + "m/s"
    analysis_result_label_acc_dis =  "Acc. Distance : " + '{:.3f}'.format(parameter.accelerometer_distance) + "mm"
    analysis_result_label = analysis_result_label_closs_vs + " " + analysis_result_label_rise_vs + "\n" + analysis_result_label_peak_vs + " " + analysis_result_label_corr_vs + "\n" + analysis_result_label_acc_dis
    ax1.text(0, -0.13, analysis_result_label, fontname="Arial", fontsize=8)

    # フーリエ変換のグラフを表示
    ax2 = fig.add_subplot(322)
    freq = np.linspace(0, 1 / (data_txt_wave[1, 0] - data_txt_wave[0, 0]), data_txt_wave_FFT.shape[0])
    ax2.plot(freq[:int(data_txt_wave_FFT.shape[0] / 2) + 1], fft[0][:int(data_txt_wave_FFT.shape[0] / 2) + 1], color="k")
    ax2.plot(freq[:int(data_txt_wave_FFT.shape[0] / 2) + 1], fft[1][:int(data_txt_wave_FFT.shape[0] / 2) + 1], color="b")
    ax2.plot(freq[:int(data_txt_wave_FFT.shape[0] / 2) + 1], fft[2][:int(data_txt_wave_FFT.shape[0] / 2) + 1], color="r")
    
    ax2.set_xlim(0, 10000)
    ax2.set_xlabel("Frequency (Hz)", fontname="Arial")

    ax2.set_ylim(0.000001, 1)
    ax2.set_ylabel("Amplitude (V)", fontname="Arial")
    ax2.set_yscale('log')

    # 相互相関関数のグラフを表示
    ax3 = fig.add_subplot(323)

    ax3.plot(corr_x, corr, color="k")
    ax3.plot(corr_x[peak_corr_index], corr[peak_corr_index], marker=".", markersize=15, markerfacecolor="r", markeredgecolor="w")

    ax3.set_xlabel("Delayed Time(s)", fontname="Arial")

    ax3.set_ylim(-15, 15)
    ax3.set_ylabel("C", fontname="Arial")

    # Show the current status in p-q plane
    ax4 = fig.add_subplot(324)
    ax4.plot(data_csv["p'___(kPa)"], data_csv["q____(kPa)"], color="k", linewidth=0.75)
    ax4.plot(data_csv["p'___(kPa)"].iloc[measured_time_index] , data_csv["q____(kPa)"].iloc[measured_time_index] ,marker=".", markersize=15, markerfacecolor="r", markeredgecolor="w")

    ax4.set_xlabel("Mean Principle Stress (kPa)", fontname="Arial")
    ax4.set_ylabel("Deviatoric Stress (kPa)", fontname="Arial")

    # Show the current status in ea-q plane
    ax5 = fig.add_subplot(325)
    ax5.plot(data_csv["e(a)_(%)_"], data_csv["q____(kPa)"], color="k", linewidth=0.75)
    ax5.plot(data_csv["e(a)_(%)_"].iloc[measured_time_index] , data_csv["q____(kPa)"].iloc[measured_time_index] ,marker=".", markersize=15, markerfacecolor="r", markeredgecolor="w")

    ax5.set_xlabel("Axial Strain (%)", fontname="Arial")
    ax5.set_ylabel("Deviatoric Stress (kPa)", fontname="Arial")

    # Show the current status in GMES-vs plane
    """
    ax6 = fig.add_subplot(326)
    ax6.plot(data_csv["Geometric Mean Eff. Str. (KPa)"], data_csv["Vs_closs(m/s)"], marker=".", markersize=15, markerfacecolor="k", markeredgecolor="w", alpha=0.3)
    ax6.plot(data_csv["Geometric Mean Eff. Str. (KPa)"].iloc[measured_time_index], data_csv["Vs_closs(m/s)"].iloc[measured_time_index], marker=".", markersize=15, markerfacecolor="r", markeredgecolor="w")

    ax6.set_xlabel("Geometric Mean Effective Stress (kPa)", fontname="Arial")
    ax6.set_ylabel("Vs based on closs-closs method (m/s)", fontname="Arial")
    """

    # 保存のファイルネームを設定
    png_name = file_list_txt[i].stem + "_{:.1f}".format(closs_shear_wave_velocity) + ".png"
    png_file_path = result_png_dir / png_name

    # 作成したVsのファイルの保存
    plt.savefig(png_file_path, dpi=200)
    plt.close()

    # explicit garbage collection
    gc.collect

    return [pass_flag, measured_time_index, closs_shear_wave_velocity, rise_shear_wave_velocity, peak_shear_wave_velocity, corr_shear_wave_velocity]

def main():

    dat_dir = Path(parameter.dat_dir)

    # find an directory including vs directory
    for cur_dir in dat_dir.glob("**/*/vs"):
        # local variables
        file_list_txt = [] # Vs測定ファイル
        file_list_spe = []
        file_list_csv = []
        file_list_dat = []

        # 解析データと同じ階層に解析結果のディレクトリを作成
        result_dir = Path(cur_dir).parent / "result"
        Path(result_dir).mkdir(exist_ok=True)
        result_png_dir = result_dir / "vs_result"
        Path(result_png_dir).mkdir(exist_ok=True)

        for file in result_png_dir.iterdir():
            if file.is_file():
                file.unlink()

        # find spe, csv and dat file in the directory
        for all_file_list in Path(cur_dir).parent.glob("*.*"):
            ext = all_file_list.suffix
            if ext == ".spe":
                file_list_spe.append(all_file_list)
            elif ext == ".csv":
                file_list_csv.append(all_file_list)
            elif ext == ".dat":
                file_list_dat.append(all_file_list)

        for all_file_list in Path(cur_dir).glob("*.*"):
            ext = all_file_list.suffix
            if ext == ".txt" or ext == ".TXT":
                file_list_txt.append(all_file_list)
    
        data_spe = pd.read_csv(file_list_spe[0], delim_whitespace=True, header=None)
        data_csv = pd.read_csv(file_list_csv[0])
        data_dat = pd.read_csv(file_list_dat[0], delim_whitespace=True)

        # make an column for Vs
        data_csv["Vs_closs(m/s)"] = 0.0
        data_csv["Vs_rise(m/s)"] = 0.0
        data_csv["Vs_peak(m/s)"] = 0.0
        data_csv["Vs_corr(m/s)"] = 0.0

        # check the directory for saving
        print("CSV : " + str(result_dir / (file_list_csv[0].stem + "_riv.csv")))

        # file_list_textの数の分だけforloop
        with concurrent.futures.ProcessPoolExecutor(max_workers=16) as executor:
            futures = [executor.submit(vs_analysis, i, cur_dir, file_list_txt, data_csv, data_spe, result_png_dir) for i in list(range(len(file_list_txt)))]
            for future in concurrent.futures.as_completed(futures):
                if future.result()[0] == True:
                    data_csv["Vs_closs(m/s)"].iloc[future.result()[1]] = future.result()[2]
                    data_csv["Vs_rise(m/s)"].iloc[future.result()[1]] = future.result()[3]
                    data_csv["Vs_peak(m/s)"].iloc[future.result()[1]] = future.result()[4]
                    data_csv["Vs_corr(m/s)"].iloc[future.result()[1]] = future.result()[5]

        # outファイルを出力
        csv_file_path = result_dir / (file_list_csv[0].stem + "_riv.csv")
        simple_csv_file_path = result_dir / (file_list_csv[0].stem + "_simple_riv.csv")
        data_csv.to_csv(csv_file_path)

        #簡略化されたoutファイルを出力(Vsが計算された行のみ)
        data_csv_simple = data_csv.iloc[np.where(data_csv["Vs_closs(m/s)"] != 0)[0], :]

        data_csv_simple.to_csv(simple_csv_file_path)

if __name__ == "__main__":
    main()