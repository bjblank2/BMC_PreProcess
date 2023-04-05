import numpy as np
import os
import shutil
from pathlib import Path
from distutils.dir_util import copy_tree
import compile_vasp_data as cvd
import matplotlib as mpl
import matplotlib.pyplot as plt
import postprocess_VASP as ppv
import support_fns as sfn
import pandas as pd
import mc_plot_fns as mcp


do_this = "COMP_VASP"
root_dir = "D:/NiMnIn_Vasp_Data/862" # Put the top directory with your vasp data here
vasp_data_file = './test' # The name of your output file
#vasp_data_file_pp = './NiMnIn_totpp'
#spin_style = ['threshold', 'threshold', 'threshold']  # options for spin_tol. Assuming [Ni Mn In]. choose 'threshold' or 'factor'
spin_tol = [.005, 2, 0]                         # insert spin parameters here, this assumes [Ni Mn In ]
use_spin_tol = True
species = ['Ni', 'Mn', 'In']  # this is the order that the post-processed data is reported, NEEDS TO BE Heusler format Ni2MnIn, Ni2FeGa.


if do_this == "COMP_VASP":
    cvd.import_vasp1(root_dir, vasp_data_file, species, use_spin_tol, spin_tol)
    print('\n')
    #M_structures = ppv.generate_m_structure(vasp_data_file, phase_tol, spin_style, spin_tol)
    #print('\n')
    #ppv.write_structures_processedvasp(M_structures, vasp_data_file_pp, vasp_data_A_pp, vasp_data_M_pp)


if do_this == "FIND_NAMES":
    for filename in Path('D:/All_VASP').rglob('OUTCAR'):
        #print(filename)
        file = filename.open('r')
        outcar = file.readlines()
        for line in outcar:
            if "TOTEN" in line:
                enrg = line.split()[4]
        print(str(enrg)+", "+str(filename))
