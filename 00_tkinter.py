import tkinter as tk
from tkinter import filedialog
import tkinter.ttk as ttk
import os

class Application(tk.Frame):

    def __init__(self, master=None):
        super().__init__(master)
        self.pack()
        self.create_widgets()

    def create_widgets(self):
        self.import_folder = tk.Button(self)
        self.import_folder["text"] = "Import..."
        self.import_folder["command"] = self.import_clicked
        self.import_folder.grid(row=0, column=0, columnspan=1, padx=5, pady=5)

        self.import_folder_entry = tk.Entry(self)
        self.import_folder_entry["text"] = "None"
        self.import_folder_entry.grid(row=0, column=1, padx=5, pady=5)
        
        self.option_radio_frame = tk.Frame(self)
        self.option_radio_frame.grid(row=1, column=0, padx=5, pady=5)
        self.option_radio_frame.radio1 = tk.Radiobutton(self.option_radio_frame, 
                                                        text="Align Files", value=1)
        self.option_radio_frame.radio1 = tk.Radiobutton(self.option_radio_frame, 
                                                        text="Draw Figures", value=1)
        self.option_radio_frame.radio1 = tk.Radiobutton(self.option_radio_frame, 
                                                        text="Compute Vs", value=1)
        
        self.import_folder_table = ttk.Treeview(self)
        self.import_folder_table["columns"] = (1)
        self.import_folder_table["show"] = "headings"
        self.import_folder_table.heading(1, text="Folder Name")
        self.import_folder_table.grid(row=3, column=0, padx=5, pady=5)

    
    def import_clicked(self):
        iDir = os.path.abspath(os.path.dirname(__file__))
        iDirPath = filedialog.askdirectory(initialdir = iDir)
        self.import_folder_entry.insert(tk.END, iDirPath)


if __name__ == "__main__":
    root = tk.Tk()
    root.geometry("960x540")
    app = Application(master=root)
    app.mainloop()