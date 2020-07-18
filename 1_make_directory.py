import os
import shutil
from pathlib import Path
import parameter

# variables
file_list_out = []
file_list_spe = []

def main():
    # listing all the .out and .spe file
    for all_file_list in os.listdir(parameter.dat_dir):
        base, ext = os.path.splitext(all_file_list)
        if ext == ".spe":
            file_list_spe.append(all_file_list)
        elif ext == ".out":
            file_list_out.append(all_file_list)

    if len(file_list_out) != 0:
        # travarse file_list_out with iterator
        for file_name in file_list_out:

            # create result directory
            result_dir = parameter.dat_dir + "\\" + Path(file_name).stem
            Path(result_dir).mkdir(exist_ok=True)

            # cut&paste all the related files to new directory
            for file in Path(parameter.dat_dir).glob("*.spe"):
                shutil.copy(str(file.absolute()), str(result_dir))
            
            for file in Path(parameter.dat_dir).glob(Path(file_name).stem + "*"):
                shutil.move(str(file.absolute()), str(result_dir)) 
    
if __name__ == "__main__":
    main()