_author__ = 'brian'
# This file handles taking the input from vasp and calculating the fit.
# The functions that are important for the actual calculating the fit
# have comments

import numpy as np
from numpy.linalg import norm
import os


def import_vasp1(root_dir, output_dir,species):
    output = open(output_dir, 'w')
    enrg_list = []
    for subdir, dirs, files in os.walk(root_dir):
        contcar_lines = []
        outcar_lines = []
        flag = 0
        flag2 = 0
        for file in files:
            if 'CONTCAR' in files and 'OUTCAR' in files:
                name = subdir.replace(root_dir,"")
                name = name.replace("/","")
                flag = 1
                if file == "CONTCAR":
                    contcar = open(subdir + '/' + file, 'r')
                    contcar_lines = contcar.readlines()
                    contcar.close()
                if file == "OUTCAR":
                    outcar = open(subdir + '/' + file, 'r')
                    outcar_lines = outcar.readlines()
                    outcar_len = len(outcar_lines)
                    for i in range(outcar_len):
                        if 'Elapsed time (sec):' in outcar_lines[i]:
                            flag2 = 1
                    if flag2 == 0:
                        print(subdir)
                    outcar.close()
            else:
                flag = 0
        print(subdir)
        if len(contcar_lines) > 0 and len(outcar_lines) > 0 and flag == 1 and flag2 == 1:
            # species order is in contcar_line[0]
            order = contcar_lines[5].split()
            #print(order)
            # Need to change so it records "a,b,c" lattice vectors
            lc = float(contcar_lines[1])
            a_v = contcar_lines[2].split()
            a_v = [float(i) for i in a_v]
            a = np.sqrt(np.power(a_v[0] * lc, 2) + np.power(a_v[1] * lc, 2) + np.power(a_v[2] * lc, 2))
            b_v = contcar_lines[3].split()
            b_v = [float(i) for i in b_v]
            b = np.sqrt(np.power(b_v[0] * lc, 2) + np.power(b_v[1] * lc, 2) + np.power(b_v[2] * lc, 2))
            c_v = contcar_lines[4].split()
            c_v = [float(i) for i in c_v]
            c = np.sqrt(np.power(c_v[0] * lc, 2) + np.power(c_v[1] * lc, 2) + np.power(c_v[2] * lc, 2))
            cos_ = np.dot(b_v, c_v) / norm(b_v) / norm(c_v)
            alpha = np.arccos(np.clip(cos_, -1, 1))
            cos_ = np.dot(a_v, c_v) / norm(a_v) / norm(c_v)
            beta = np.arccos(np.clip(cos_, -1, 1))
            cos_ = np.dot(a_v, b_v) / norm(a_v) / norm(b_v)
            gamma = np.arccos(np.clip(cos_, -1, 1))
            comp = contcar_lines[6].split()
            composition = [[0,species[-1]]] *len(species)
            for i in range(len(comp)):
                composition[i]=[comp[i],order[i]]
            # for i in range(len(composition)):
            #     print(composition[i])
            #     print(files)
            composition_ordered = sorted(composition, key=lambda d: species.index(d[1]))

            for i in range(outcar_len):
                if "TOTEN" in outcar_lines[i]:
                    enrg = outcar_lines[i].split()
                    enrg = float(enrg[4])
            if enrg not in enrg_list:
                enrg_list.append(enrg)
                output.write("# ")
                for i in range(len(species)):
                    output.write(str(species[i])+ " ")
                output.write('\n')
                total_num = 0
                for i in range(len(composition_ordered)):
                    output.write(str(composition_ordered[i][0]) + "\t")
                    total_num += int(composition_ordered[i][0])
                mag_list = [[0,species[-1]]] * (total_num)
                for i in range(outcar_len):
                    if "magnetization (x)" in outcar_lines[i]:
                        for j in range(total_num):
                            mag = outcar_lines[i + j + 4].split()
                            if j < int(composition[0][0]):
                                kind = order[0]
                            elif j< int(composition[0][0]) + int(composition[1][0]):
                                kind = order[1]
                            else:
                                kind = order[2]
                            mag_list[j] = [mag[4],kind]
                #print(name)
                #print(mag_list)
                mag_list_ordered = sorted(mag_list, key=lambda d: species.index(d[1]))
                enrg = str(enrg)
                a = str(a)
                b = str(b)
                c = str(c)
                alpha = str(alpha)
                beta = str(beta)
                gamma = str(gamma)
                output_line = name + "\t" + enrg + "\t" + a + "\t" + b + "\t" + c + "\t" + alpha + "\t" + beta + "\t" + gamma +"\n"
                output.write(output_line)
                output_line = str(a_v[0]) + "\t" + str(a_v[1]) + "\t" + str(a_v[2]) + "\n"
                output.write(output_line)
                output_line = str(b_v[0]) + "\t" + str(b_v[1]) + "\t" + str(b_v[2]) + "\n"
                output.write(output_line)
                output_line = str(c_v[0]) + "\t" + str(c_v[1]) + "\t" + str(c_v[2]) + "\n"
                output.write(output_line)
                index = 0
                pos_list = []
                for i in range(8, 8 + total_num):
                    pos = contcar_lines[i].split()
                    if i-8 < int(composition[0][0]):
                        kind = order[0]
                    elif i-8 < int(composition[0][0]) + int(composition[1][0]):
                        kind = order[1]
                    else:
                        kind = order[2]
                    pos_list.append([pos,kind])

                # resort pos according to the order of species
                pos_list_ordered = sorted(pos_list, key=lambda d: species.index(d[1]))
                for i in range(total_num):
                    for j in range(len(species)):
                        if pos_list_ordered[i][1] == species[j]:
                            sp_index = j
                    output_line = "\t" + str(i ) + "\t" + str(sp_index) + "\t" + str(mag_list_ordered[index][0]) + "\t" + str(pos_list_ordered[i][0][0]) + "\t" + str(
                        pos_list_ordered[i][0][1]) + "\t" + str(pos_list_ordered[i][0][2]) + "\n"
                    output.write(output_line)
                    index += 1
    output.close()



def import_vasp(root_dir, output_dir,species):
    output = open(output_dir, 'w')
    for subdir, dirs, files in os.walk(root_dir):
        contcar_lines = []
        outcar_lines = []
        flag = 0
        for file in files:
            if 'CONTCAR' in files and 'OUTCAR' in files:
                name = subdir.replace(root_dir,"")
                name = name.replace("/","")
                flag = 1
                if file == "CONTCAR":
                    contcar = open(subdir + '/' + file, 'r')
                    contcar_lines = contcar.readlines()
                    contcar.close()
                if file == "OUTCAR":
                    outcar = open(subdir + '/' + file, 'r')
                    outcar_lines = outcar.readlines()
                    outcar_len = len(outcar_lines)
                    outcar.close()
            else:
                flag = 0
        if len(contcar_lines) > 0 and len(outcar_lines) > 0 and flag == 1:
            # species order is in contcar_line[0]
            order = contcar_lines[5].split()
            #print(order)
            lc = float(contcar_lines[1])                # ARE ALL OF THESE IN ORDER A, B, C???
            a = contcar_lines[2].split()
            a = float(a[0]) * lc
            b = contcar_lines[3].split()
            b = float(b[1]) * lc
            c = contcar_lines[4].split()
            c = float(c[2]) * lc
            composition = [0] * len(species)
            comp = contcar_lines[6].split()
            composition = [[0,species[-1]]] *len(species)
            for i in range(len(comp)):
                composition[i]=[comp[i],order[i]]

            composition_ordered = sorted(composition, key=lambda d: species.index(d[1]))


            output.write("# ")
            for i in range(len(species)):
               output.write(str(species[i])+ " ")
            output.write('\n')
            total_num = 0
            for i in range(len(composition_ordered)):
                output.write(str(composition_ordered[i][0]) + "\t")
                total_num += int(composition_ordered[i][0])
            mag_list = [[0,species[-1]]] * (total_num)
            for i in range(outcar_len):
                if "(sigma->0)" in outcar_lines[i]:
                    enrg = outcar_lines[i].split()
                    enrg = float(enrg[-1])
                if "TOTEN" in outcar_lines[i]:
                    enthalpy = outcar_lines[i].split()
                    enthalpy = float(enthalpy[4])
                if "pressure" in outcar_lines[i]:
                    pressure = outcar_lines[i].split()
                    ext_P = float(pressure[3])
                    pullay = float(pressure[8])
                    hess = outcar_lines[i-1]
                    hess = hess.split()
                    XX = hess[2]
                    YY = hess[3]
                    ZZ = hess[4]
                    XY = hess[5]
                    YZ = hess[6]
                    ZX = hess[7]
                if "length of vectors" in outcar_lines[i]:
                    lattice_const = outcar_lines[i+1]
                    latA = float(lattice_const[0])
                    latB = float(lattice_const[1])
                    latC = float(lattice_const[2])
                if "magnetization (x)" in outcar_lines[i]:
                    for j in range(total_num):
                        mag = outcar_lines[i + j + 4].split()
                        if j < int(composition[0][0]):
                            kind = order[0]
                        elif j < int(composition[0][0]) + int(composition[1][0]):
                            kind = order[1]
                        else:
                            kind = order[2]
                        mag_list[j] = [mag[4],kind]
            #print(name)
            #print(mag_list)
            mag_list_ordered = sorted(mag_list, key=lambda d: species.index(d[1]))
            enrg = str(enrg)
            a = str(a)
            b = str(b)
            c = str(c)
            output_line = name + "\t" + enrg + "\t" + a + "\t" + b + "\t" + c + "\t" + str(pullay) + "\t" + str(ext_P) + "\t" + XX + "\t" + YY + "\t" + ZZ + "\t" + XY + "\t" + YZ + "\t" + ZX + "\n"
            output.write(output_line)
            index = 0
            pos_list = []
            for i in range(8, 8 + total_num):
                pos = contcar_lines[i].split()
                if i-8 < int(composition[0][0]):
                    kind = order[0]
                elif i-8 < int(composition[0][0]) + int(composition[1][0]):
                    kind = order[1]
                else:
                    kind = order[2]
                pos_list.append([pos,kind])

            # resort pos according to the order of species
            pos_list_ordered = sorted(pos_list, key=lambda d: species.index(d[1]))
            for i in range(total_num):
                output_line = "\t" + str(i ) + "\t" + str(mag_list_ordered[index][0]) + "\t" + str(pos_list_ordered[i][0][0]) + "\t" + str(
                    pos_list_ordered[i][0][1]) + "\t" + str(pos_list_ordered[i][0][2]) + "\n"
                output.write(output_line)
                index += 1
    output.close()

