_author__ = 'brian'
# This file handles taking the input from vasp and calculating the fit.
# The functions that are important for the actual calculating the fit
# have comments

import numpy as np
import copy
import m_structure
import os
import matplotlib.pyplot as plt
#from scipy.optimize import least_squares
#from sklearn.linear_model import Ridge
#from sklearn import linear_model

def generate_m_structure(data_file, aust_tol, spin_style, spin_tol):
    m_struct_list = []
    data = open(data_file, 'r')
    lines = data.readlines()
    norms = []
    for i in range(len(lines)):
        if '#' in lines[i]:
            species = (lines[i].split())
            species.pop(0)
            print(lines[i+1])
            m_struct = m_structure.MStructureObj(lines[i+1], species, aust_tol)
            for j in range(m_struct.num_Atoms):
                atom_data = lines[i + j + 2]
                m_struct.set_atom_properties(j, atom_data, spin_style, spin_tol)
            #####################
            if check_zero_spin(m_struct.basis)=='False':
                if m_struct.phase_name != 'prem':
                    m_struct_list.append(m_struct)
            #################
            #norms.append(j_count)
    return m_struct_list

def scale_enrg(M_structures):
    e_comp0 = []
    e_comp50 = []
    for i in range(len(M_structures)):
        if M_structures[i].phase_name != "pmmmm":
            structure = M_structures[i]
            comp = structure.composition[2]/structure.composition[0]
            if comp == 0:
                e_comp0.append(structure.enrg)
            if comp == .5:
                e_comp50.append(structure.enrg)
    comp0_min = min(e_comp0)
    comp50_min = min(e_comp50)
    offset = (comp50_min-comp0_min)/.5
    #offset = 0
    enrg_list = []
    for i in range(len(M_structures)):
        if M_structures[i].phase_name != "pmmmm":
            comp = M_structures[i].composition[2]/M_structures[i].composition[0]
            old_enrg = M_structures[i].enrg
            M_structures[i].enrg = old_enrg-(offset*comp + comp0_min)
    emax = 0
    for i in range(len(M_structures)):
        if M_structures[i].enrg > emax:
            emax = M_structures[i].enrg
    #emax = emax*0
    for i in range(len(M_structures)):
        M_structures[i].enrg -= emax

def write_structures_processedvasp(structures,data_file_pp,data_file_A_pp,data_file_M_pp):
    file = open(data_file_pp, 'w')
    fileA = open(data_file_A_pp, 'w')
    fileM = open(data_file_M_pp, 'w')
    for i in range(len(structures)):
        file.write('# Ni Mn In \n')            # hard coding!
        mat = structures[i]
        out = [mat.composition, mat.name, mat.enrg, mat.phase_name, mat.LCs]
        if mat.phase_name == 'aust':
            file2 = fileA
        elif mat.phase_name == 'mart':
            file2 = fileM
        file2.write('# Ni Mn In \n')
        for j in range(len(out)):
            sums = str(out[j])
            file.write(sums.ljust(25))
            file2.write(sums.ljust(25))
        file.write("".ljust(10))
        file2.write("".ljust(10))
        #file.write(str(norms[i]).ljust(20))
        file.write("\n")
        file2.write("\n")
        for j in range(mat.num_Atoms):
            file.write("\t")
            file2.write("\t")
            out = [mat.basis[j].atom_index,mat.basis[j].species,mat.basis[j].spin,mat.basis[j].a_pos,mat.basis[j].b_pos,mat.basis[j].c_pos]
            for k in range(len(out)):
                sums = str(out[k])
                file.write(sums.ljust(17))
                file2.write(sums.ljust(17))
            file.write("\n")
            file2.write("\n")
    file.close()
    fileA.close()
    fileM.close()

def summarize_fitting_structures(structures,threshold):
   path = 'summary_fitting_structures'
   file = open(path, 'w')
   file.write("NAME".ljust(15) + "PHASE".ljust(7) + "LCONST".ljust(15) + "MAG".ljust(6) + "ENERG".ljust(17) + "SUMS->\n")

   Cluster_sum_length = len(structures[0].Cluster_sums)
   J_sum_length = len(structures[0].J_sums)
   C_sum_sum = np.array([0.0] * Cluster_sum_length)
   J_sum_sum = np.array([0.0] * J_sum_length)
   C_max = np.array([0.0] * Cluster_sum_length)
   J_max = np.array([0.0] * J_sum_length)
   for i in range(len(structures)):
       mat = structures[i]
       out = [mat.enrg, mat.Cluster_sums, mat.J_sums]
       np_C_sum = np.absolute(np.array(mat.Cluster_sums))
       np_J_sum = np.absolute(np.array(mat.J_sums))
       for j in range(len(np_C_sum)):
           if np_C_sum[j] > C_max[j]:
               C_max[j] = np_C_sum[j]
       for j in range(len(np_J_sum)):
           if np_J_sum[j] > J_max[j] :
               J_max[j] = np_J_sum[j]

       C_sum_sum = C_sum_sum + np_C_sum
       J_sum_sum = J_sum_sum + np_J_sum
       file.write(mat.name.ljust(15) + mat.phase_name.ljust(7) + str(round(mat.LCs[0],2)).ljust(5) + str(round(mat.LCs[1],2)).ljust(5) + str(round(mat.LCs[2],2)).ljust(5) + mat.mag_phase.ljust(7))
       for j in range(len(out)):
           sums = str(out[j])
           if j == 0:
               file.write(sums.ljust(17))
           else:
               file.write(sums.ljust(7))
       file.write("\n")
   file.write(str(C_sum_sum.tolist()) + str(J_sum_sum.tolist()))
   file.write("\n")

   C_max = C_max * len(structures)
   J_max = J_max * len(structures)
   # if C_max !=0:
   #     C_ratio = C_sum_sum/C_max
   #     J_ratio = J_sum_sum/J_max
   #     warning = False
   #     for i in range(len(C_ratio)):
   #         if np.isnan(C_ratio[i]):
   #             C_ratio[i] = 0.0
   #
   #         if C_ratio[i] < threshold:
   #             warning = True
   #     for i in range(len(J_ratio)):
   #         if np.isnan(J_ratio[i]):
   #             J_ratio[i] = 0.0
   #
   #         if J_ratio[i] < threshold:
   #             warning = True
   #     if warning:
   #         print("Warning! Check summary_fitting_structures")
   #     file.write(str(C_ratio.tolist()) + str(J_ratio.tolist()))
   file.write("\n")
   file.close()

def creat_aust(structure_list,set_point):
    for i in range(len(structure_list)):
        phase = structure_list[i].phase_name

def check_zero_spin(basis):
    zeros = 'True'
    for i in range(len(basis)):
        if (basis[i].spin != 0):
            zeros = 'False'
    return zeros

def set_KJ_ratio(m_structure_list,monomer_rule_list, cluster_rule_list, j_rule_list, KoJ):
    m_structure_KoJ = copy.deepcopy(m_structure_list[len(m_structure_list)-1])
    m_structure_KoJ.name = 'KoJ'
    m_structure_KoJ.enrg = 0.0
    m_structure_KoJ.original_enrg = 0.0
    for i in range(len(m_structure_KoJ.Monomer_sums)):
        m_structure_KoJ.Monomer_sums[i] = 0
    for i in range(len(m_structure_KoJ.Cluster_sums)):
        m_structure_KoJ.Cluster_sums[i]=0
    for i in range(len(m_structure_KoJ.J_sums)):
        m_structure_KoJ.J_sums[i] = 0
    for i in range(len(cluster_rule_list)):
        if cluster_rule_list[i].phase == 'mart':
            m_structure_KoJ.Cluster_sums[i] = -KoJ
        elif cluster_rule_list[i].phase == 'aust':
            m_structure_KoJ.Cluster_sums[i] = 1
    for i in range(len(j_rule_list)):
        if j_rule_list[i].phase == 'mart':
            m_structure_KoJ.J_sums[i] = -KoJ
        elif j_rule_list[i].phase == 'aust':
            m_structure_KoJ.J_sums[i] = 1
    m_structure_list.append(m_structure_KoJ)
