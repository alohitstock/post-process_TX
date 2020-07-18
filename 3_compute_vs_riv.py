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
from scipy.signal import find_peaks
import parameter

class ComputeVs():

    def __init__(self, path, skiprows=9, timecol=0, inputwavecol=1, outputwavecol1=2, outputwavecol2=3):
        
        self.file_path = Path(path)
        self.transmitted_time = {"closs": 0.,
                                 "rise" : 0.,
                                 "peak" : 0.,
                                 "closscorr" : 0.}
        self.isPassed = True

        # read txt file of input/output wave
        data_txt = pd.read_csv( self.file_path, 
                                delimiter=",", 
                                header=None, 
                                names=[i for i in range(4)])
        self.wave_cols = [timecol, inputwavecol, outputwavecol1, outputwavecol2]
        self.wave_cols_inoutput = [outputwavecol1, outputwavecol2]
        self.data_txt_wave = data_txt.iloc[skiprows:, self.wave_cols].to_numpy().astype(np.float)

        # baseline correction
        self.data_txt_wave = self.BaseCollection()

        # detect trigger index of waves
        if self.isPassed == True:
            self.trigger_time_index = self.TriggerTimeIndex()
        
        self.DrawFigure()

    def BaseCollection(self):
        # low-pass and high-pass filter (ongoing)


        # extract the rows when time is less than threshold_time
        data_txt_wave_basis = self.data_txt_wave[np.where(self.data_txt_wave[:, 0] < parameter.threshold_time)[0], :]

        # compute the average and the standard deviation of the extracted rows
        base_voltage = np.average(data_txt_wave_basis, axis=0)
        self.std_voltage = np.std(data_txt_wave_basis, axis=0)

        # if std_voltage is over threshold_standard_deviation, return error
        if (self.std_voltage > parameter.threshold_standard_deviation).sum() > 1:
            self.isPassed = False
            print("File Name: " + self.file_path.stem + " |  Vs : N/A (too large std of output)")
        
        self.data_txt_wave[:, self.wave_cols_inoutput] -= base_voltage[self.wave_cols_inoutput]
        return self.data_txt_wave


    def TriggerTimeIndex(self):

        maximum_value = np.max(np.abs(self.data_txt_wave), axis=0)
        
        # detect trigger index by using maximum value of waves
        if parameter.detection_method_for_triggered_point == 0:
            threshold_trigger = parameter.magnification_for_maximum_value * maximum_value
            trigger_time_index = self.TriggerTimeIndexValidation(threshold_trigger)
        
        # detect trigger index by using magnification_standard_deviation
        elif parameter.detection_method_for_triggered_point == 1:
            threshold_trigger = parameter.magnification_standard_deviation * self.std_voltage
            trigger_time_index = self.TriggerTimeIndexValidation(threshold_trigger)
        
        # detect trigger index by scipy peak search
        elif parameter.detection_method_for_triggered_point == 2:
            peak_index_1, _ = find_peaks(np.abs(self.data_txt_wave[:, self.wave_cols_inoutput[0]]),
                                         distance=200, height=0.1*maximum_value[self.wave_cols_inoutput[0]])
            peak_index_2, _ = find_peaks(np.abs(self.data_txt_wave[:, self.wave_cols_inoutput[1]]),
                                         distance=200, height=0.1*maximum_value[self.wave_cols_inoutput[1]])
            trigger_time_index_1 = 2*peak_index_1[0] - peak_index_1[1]
            trigger_time_index_2 = 2*peak_index_2[0] - peak_index_2[1]
            trigger_time_index = np.array([trigger_time_index_1, trigger_time_index_2])
        
        elif parameter.detection_method_for_triggered_point == 3:
            diff_data = np.diff(self.data_txt_wave[:, self.wave_cols_inoutput], axis=0)

            max_slope_index_1, _ = find_peaks(diff_data[:, 0], prominence=0.00004)
            max_slope_index_2, _ = find_peaks(diff_data[:, 1], prominence=0.00004)

            max_slope_index_1 = max_slope_index_1[0]
            max_slope_index_2 = max_slope_index_2[0]

            trigger_time_index_1 = int(max_slope_index_1 - self.data_txt_wave[max_slope_index_1 + 1, self.wave_cols_inoutput[0]] / diff_data[max_slope_index_1, 0])
            trigger_time_index_2 = int(max_slope_index_2 - self.data_txt_wave[max_slope_index_2 + 1, self.wave_cols_inoutput[1]] / diff_data[max_slope_index_2, 1])

            trigger_time_index = np.array([trigger_time_index_1, trigger_time_index_2])
        
        elif parameter.detection_method_for_triggered_point == 4:
            fft_data = np.fft.fftn(self.data_txt_wave[:, self.wave_cols_inoutput])
            timestep = self.data_txt_wave[1, self.wave_cols[0]] - self.data_txt_wave[0, self.wave_cols[0]]
            fft_data_freq = np.fft.fftfreq(self.data_txt_wave.shape[0], d=timestep)

            ifft_highpassed = np.fft.ifftn(fft_data[fft_data_freq < 5000])
            ifft_time = np.linspace(self.data_txt_wave[:, self.wave_cols[0]].min(),
                                    self.data_txt_wave[:, self.wave_cols[0]].max(),
                                    ifft_highpassed.shape[0])
            
            self.DrawTempFigure(ifft_time, ifft_highpassed)

            trigger_time_index_1 = 0
            trigger_time_index_2 = 0

            trigger_time_index = np.array([trigger_time_index_1, trigger_time_index_2])
        
        return trigger_time_index
    
    def TriggerTimeIndexValidation(self, threshold_trigger):
        trigger_time_index_1 = np.where(np.abs(self.data_txt_wave[:, self.wave_cols_inoutput[0]]) 
                                        > threshold_trigger[self.wave_cols_inoutput[0]])[0]
        trigger_time_index_2 = np.where(np.abs(self.data_txt_wave[:, self.wave_cols_inoutput[1]]) 
                                        > threshold_trigger[self.wave_cols_inoutput[1]])[0]
        
        if trigger_time_index_1.size == 0 or trigger_time_index_2.size == 0:
            self.isPassed = False
            print("File Name: " + self.file_path.stem + " | Vs : N/A (too small amp.)")
        
        trigger_time_index = np.array([trigger_time_index_1[0], trigger_time_index_2[0]])
        return trigger_time_index

    def DrawFigure(self):
        fig, axes = plt.subplots(1, 1, figsize=(6,9))
        fig.subplots_adjust(wspace=0.25)
        axes.plot(self.data_txt_wave[:, self.wave_cols[0]], self.data_txt_wave[:, self.wave_cols[1]])

        axes_sub = axes.twinx()
        axes_sub.plot(self.data_txt_wave[:, self.wave_cols[0]], self.data_txt_wave[:, self.wave_cols[2]])
        axes_sub.plot(self.data_txt_wave[self.trigger_time_index[0], self.wave_cols[0]], 
                                         self.data_txt_wave[self.trigger_time_index[0], self.wave_cols[2]],
                      marker="*")
        axes_sub.plot(self.data_txt_wave[:, self.wave_cols[0]], self.data_txt_wave[:, self.wave_cols[3]])
        axes_sub.plot(self.data_txt_wave[self.trigger_time_index[1], self.wave_cols[0]],
                                         self.data_txt_wave[self.trigger_time_index[1], self.wave_cols[3]],
                      marker="*")
    
        fig.savefig("20200421_300-200_test.png", dpi=600, bbox_inches='tight')
        plt.close(fig)
    
    def DrawTempFigure(self, arr_x, arr_y):
        fig, axes = plt.subplots(1, 1, figsize=(6,9))
        fig.subplots_adjust(wspace=0.25)
        axes.plot(arr_x, arr_y)
        # axes.set_ylim(0.00005, -0.00005)
        # axes.set_xscale("log")
        # axes.set_yscale("log")
        fig.savefig("temp.png", dpi=600, bbox_inches='tight')
        plt.close(fig)

def main():

    comVs = ComputeVs(r"D:\OneDrive\01_Research\01_Doctor\07_Triaxial_Test\07_Result\2020\20200421\20200421_06_B-con\vs\300-200.TXT")


if __name__ == "__main__":
    main()

