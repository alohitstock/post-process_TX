"""
Todo
- ファイル、フォルダリストの選択+削除機能の追加
"""

import tkinter as tk
import tkinter.ttk as ttk
import os
from tkinter import filedialog
from pathlib import Path

class Frame(tk.Frame):

    def __init__(self, master=None):
        tk.Frame.__init__(self, master, height=540, width=960)
        self.master.title('Vs Analysis')

        self.iDirPath = __file__
        self.draw_widgets()

    def draw_widgets(self):
        # First Frame
        frame0 = tk.Frame(self, relief=tk.FLAT, bd=2)

        self.option_radio_var = tk.IntVar()
        self.option_radio_var.set(0)

        option_list = ["Align Files", "Draw Figures", "Compute Vs"]
        for i, option in enumerate(option_list):
            self.option_radio_buttons = tk.Radiobutton(frame0, text=option, variable=self.option_radio_var, value=i, command= lambda : self.list_directory())
            self.option_radio_buttons.grid(row=0, column=i)
        
        frame0.place(relx=0.025, rely=0)
        
        frame1 = tk.Frame(self, relief=tk.FLAT, bd=2)
        self.import_folder_button = tk.Button(frame1, text="Import...")
        self.import_folder_button["command"] = self.import_clicked
        self.import_folder_button.grid(row=0, column=0, columnspan=1, padx=5, pady=5)

        self.import_folder_entry = tk.Entry(frame1, text="None", width=50)
        self.import_folder_entry.grid(row=0, column=1, columnspan=1, padx=5, pady=5)

        frame1.place(relx=0.025, rely=0.05)

        self.import_directory_table = ttk.Treeview(frame1, columns=1, show="headings")
        self.import_directory_table.heading(1, text="Directory Name")
        self.import_directory_table.grid(row=1, column=0, padx=5, pady=5, columnspan=2, sticky=tk.W+tk.E)
        
        # Second Frame
        frame2 = tk.Frame(self, relief=tk.FLAT, bd=2)
        frame2.place(relx=0.6, rely=0.6)

    def import_clicked(self):
        iDir = os.path.abspath(os.path.dirname(__file__))
        self.iDirPath = filedialog.askdirectory(initialdir = iDir)
        self.import_folder_entry.delete(0, tk.END)
        self.import_folder_entry.insert(tk.END, self.iDirPath)
        self.list_directory()

    def list_directory(self):
        self.get_radio_values()
        print(self.option_num)
        self.iDirPath = Path(self.iDirPath)
        for i in self.import_directory_table.get_children():
            self.import_directory_table.delete(i)
        if self.option_num == 0:
            for out_file in self.iDirPath.glob("**/*.out"):
                self.import_directory_table.insert("", tk.END, values=(out_file.name))
        else:
            for out_file in self.iDirPath.glob("**/*.out"):
                self.import_directory_table.insert("", tk.END, values=(out_file.parent.name))


    def get_radio_values(self):
        self.option_num = self.option_radio_var.get()
        
if __name__ == '__main__':
    f = Frame()
    f.pack()
    f.mainloop()