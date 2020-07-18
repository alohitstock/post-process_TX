import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
import numpy as np


class PrepareTrainingData():
    def __init__(self, path, skiprows=9):
        self.file_path = Path(path)
        self.skip_rows = skiprows
    
    def drawfigure(self):
        data_txt = pd.read_csv( self.file_path, 
                                delimiter=",", 
                                header=None, 
                                names=[i for i in range(4)])
        self.data_txt_wave = data_txt.iloc[self.skip_rows:, :].to_numpy().astype(np.float)
        fig, self.axes = plt.subplots(2, 1, figsize=(12,9))
        fig.subplots_adjust(wspace=0.25)

        for i in range(2):
            self.axes[i].plot(self.data_txt_wave[:, 0], self.data_txt_wave[:, i+2])
        
        self.ln_0, = self.axes[0].plot([], [], "o")
        self.ln_1, = self.axes[1].plot([], [], "o")

        fig.canvas.mpl_connect('button_press_event', self.onclick)
        fig.canvas.mpl_connect('motion_notify_event', self.movemouse)
        fig.canvas.mpl_connect('axes_enter_event', self.movemouse)
        fig.canvas.mpl_connect('axes_leave_event', self.movemouse)
        plt.show()
    
    def onclick(self, event):
        if event.dblclick:
            print(event.button)

    def movemouse(self, event):
        if event.xdata is None:
            return
        else:
            index_x = np.abs(self.data_txt_wave[:, 0] - event.xdata).argmin()
            if event.inaxes == self.axes[0]:
                print(self.ln_0)
                self.ln_0.set_data(self.data_txt_wave[index_x, 0], self.data_txt_wave[index_x, 2])
            else:
                self.ln_1.set_data(self.data_txt_wave[index_x, 0], self.data_txt_wave[index_x, 3])
            plt.draw()



def main():
    temp = PrepareTrainingData(r"D:\OneDrive\01_Research\01_Doctor\07_Triaxial_Test\07_Result\2020\20200421\20200421_06_B-con\vs\300-200.TXT")
    temp.drawfigure()


if __name__ == "__main__":
    main()
