import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import pandas as pd
import matplotlib.lines as mlines

def get_file_data(rootdir, atom_numbs):
    sims = []
    raw_data = {}
    count = 0
    for subdir, dirs, files in os.walk(rootdir):
        for file in files:
            if "OUTPUT" in file:
                count += 1
                sim = pd.read_csv(os.path.join(subdir, file), header=8)
                sim['temp'] = sim['temp'].str.replace('#', '').astype(float)
                sims.append(sim)
                raw_data['temp_'     + str(count)] = np.flip(sim['temp'])
                raw_data['Cmag_'     + str(count)] = np.flip(sim[' Cmag'])/atom_numbs
                raw_data['enrg_'     + str(count)] = np.flip(sim[' enrg'])/atom_numbs
                raw_data['mag1_'     + str(count)] = np.flip(sim[' mag'])/atom_numbs
                raw_data['mag2_'     + str(count)] = np.flip(sim[' mag2'])/atom_numbs
                raw_data['mag3_'     + str(count)] = np.flip(sim[' mag3'])/atom_numbs
                raw_data['Xmag_'     + str(count)] = np.flip(sim[' Xmag'])/atom_numbs
                raw_data['var_e_'    + str(count)] = np.flip(sim[' var_e'])/atom_numbs
                raw_data['var_spin_' + str(count)] = np.flip(sim[' var_spin'])/atom_numbs
                sim_again = open(os.path.join(subdir, file))
                sim_lines = sim_again.readlines()
                raw_data["SRO_"      + str(count)] = float(sim_lines[6].split()[-1])
                raw_data["numbs_"    + str(count)] = [int(sim_lines[1].split()[1].strip(',')), int(sim_lines[1].split()[2].strip(',')), int(sim_lines[1].split()[3])]
                raw_data["Tenrg_"    + str(count)] = float(sim_lines[7].split(',')[0])
                raw_data["Menrg_"    + str(count)] = float(sim_lines[7].split(',')[1])
                raw_data["Cenrg_"    + str(count)] = raw_data["Tenrg_" + str(count)] - raw_data["Menrg_" + str(count)]
                sim_again.close()
    raw_data['count'] = count
    Cmag_avgs = []
    enrg_avgs = []
    mag1_avgs = []
    mag2_avgs = []
    mag3_avgs = []
    Xmag_avgs = []
    var_e_avgs = []
    var_spin_avgs = []
    enrg_Chem = np.mean([raw_data["Cenrg_" + str(i+1)] for i in range(count)])
    SRO = np.mean([raw_data["SRO_" + str(i+1)] for i in range(count)])

    for temp_ind in range(len(sims[0]['temp'])):
        Cmag_avg = 0
        enrg_avg = 0
        mag1_avg = 0
        mag2_avg = 0
        mag3_avg = 0
        Xmag_avg = 0
        var_e_avg = 0
        var_spin_avg = 0
        for i in range(len(sims)):
            Cmag_avg += sims[i][' Cmag'][temp_ind]/len(sims)/atom_numbs
            enrg_avg += sims[i][' enrg'][temp_ind]/len(sims)/atom_numbs
            mag1_avg += sims[i][' mag'][temp_ind]/len(sims)/atom_numbs
            mag2_avg += sims[i][' mag2'][temp_ind]/len(sims)/atom_numbs
            mag3_avg += np.power(sims[i][' mag3'][temp_ind]/len(sims)/atom_numbs,2)
            Xmag_avg += sims[i][' Xmag'][temp_ind]/len(sims)/atom_numbs
            var_e_avg += sims[i][' var_e'][temp_ind]/len(sims)/atom_numbs
            var_spin_avg += sims[i][' var_spin'][temp_ind]/len(sims)/atom_numbs
        Cmag_avgs.append(Cmag_avg)
        enrg_avgs.append(enrg_avg)
        mag1_avgs.append(mag1_avg)
        mag2_avgs.append(mag2_avg)
        mag3_avgs.append(np.sqrt(mag3_avg))
        Xmag_avgs.append(Xmag_avg)
        var_e_avgs.append(var_e_avg)
        var_spin_avgs.append(var_spin_avg)
    avg_data = {}
    avg_data['temp'] = np.flip([sims[0]['temp'][i] for i in range(len(sims[0]['temp']))])
    avg_data['Cmag'] = np.flip(Cmag_avgs)
    avg_data['enrg'] = np.flip(enrg_avgs)
    avg_data['mag1'] = np.flip(mag1_avgs)
    avg_data['mag2'] = np.flip(mag2_avgs)
    avg_data['mag3'] = np.flip(mag3_avgs)
    avg_data['Xmag'] = np.flip(Xmag_avgs)
    avg_data['var_e'] = np.flip(var_e_avgs)
    avg_data['var_spin'] = np.flip(var_spin_avgs)
    data = {'raw_data': raw_data, 'avg_data': avg_data, 'sro': SRO, 'enrg_chem': enrg_Chem}
    return data

def get_file_data_NEW(rootdir, atom_numbs):
    sims = []
    raw_data = {}
    count = 0
    for subdir, dirs, files in os.walk(rootdir):
        for file in files:
            if "OUTPUT" in file:
                count += 1
                sim = pd.read_csv(os.path.join(subdir, file), header=8)
                sim['temp'] = sim['temp'].str.replace('#', '').astype(float)
                sims.append(sim)
                raw_data['temp_'     + str(count)] = np.flip(sim['temp'])
                raw_data['Cmag_'     + str(count)] = np.flip(sim[' Cmag'])/atom_numbs
                raw_data['enrg_'     + str(count)] = np.flip(sim[' enrg'])/atom_numbs
                raw_data['mag1_'     + str(count)] = np.flip(sim[' mag'])/atom_numbs
                raw_data['mag2_'     + str(count)] = np.flip(sim[' magMn'])/atom_numbs
                #raw_data['mag3_'     + str(count)] = np.flip(sim[' mag3'])/atom_numbs
                raw_data['Xmag_'     + str(count)] = np.flip(sim[' Xmag'])/atom_numbs
                raw_data['var_e_'    + str(count)] = np.flip(sim[' var_e'])/atom_numbs
                raw_data['var_spin_' + str(count)] = np.flip(sim[' var_spin'])/atom_numbs
                sim_again = open(os.path.join(subdir, file))
                sim_lines = sim_again.readlines()
                raw_data["SRO_"      + str(count)] = float(sim_lines[6].split()[-1])
                raw_data["numbs_"    + str(count)] = [int(sim_lines[1].split()[1].strip(',')), int(sim_lines[1].split()[2].strip(',')), int(sim_lines[1].split()[3].strip(','))]
                raw_data["Tenrg_"    + str(count)] = float(sim_lines[7].split(',')[0])
                raw_data["Menrg_"    + str(count)] = float(sim_lines[7].split(',')[1])
                raw_data["Cenrg_"    + str(count)] = raw_data["Tenrg_" + str(count)] - raw_data["Menrg_" + str(count)]
                sim_again.close()
    raw_data['count'] = count
    Cmag_avgs = []
    enrg_avgs = []
    mag1_avgs = []
    mag2_avgs = []
    mag3_avgs = []
    Xmag_avgs = []
    var_e_avgs = []
    var_spin_avgs = []
    enrg_Chem = np.mean([raw_data["Cenrg_" + str(i+1)] for i in range(count)])
    SRO = np.mean([raw_data["SRO_" + str(i+1)] for i in range(count)])

    for temp_ind in range(len(sims[0]['temp'])):
        Cmag_avg = 0
        enrg_avg = 0
        mag1_avg = 0
        mag2_avg = 0
        mag3_avg = 0
        Xmag_avg = 0
        var_e_avg = 0
        var_spin_avg = 0
        for i in range(len(sims)):
            Cmag_avg += sims[i][' Cmag'][temp_ind]/len(sims)/atom_numbs
            enrg_avg += sims[i][' enrg'][temp_ind]/len(sims)/atom_numbs
            mag1_avg += sims[i][' mag'][temp_ind]/len(sims)/atom_numbs
            mag2_avg += sims[i][' magMn'][temp_ind]/len(sims)/atom_numbs
            Xmag_avg += sims[i][' Xmag'][temp_ind]/len(sims)/atom_numbs
            var_e_avg += sims[i][' var_e'][temp_ind]/len(sims)/atom_numbs
            var_spin_avg += sims[i][' var_spin'][temp_ind]/len(sims)/atom_numbs
        Cmag_avgs.append(Cmag_avg)
        enrg_avgs.append(enrg_avg)
        mag1_avgs.append(mag1_avg)
        mag2_avgs.append(mag2_avg)
        mag3_avgs.append(np.sqrt(mag3_avg))
        Xmag_avgs.append(Xmag_avg)
        var_e_avgs.append(var_e_avg)
        var_spin_avgs.append(var_spin_avg)
    avg_data = {}
    avg_data['temp'] = np.flip([sims[0]['temp'][i] for i in range(len(sims[0]['temp']))])
    avg_data['Cmag'] = np.flip(Cmag_avgs)
    avg_data['enrg'] = np.flip(enrg_avgs)
    avg_data['mag1'] = np.flip(mag1_avgs)
    avg_data['mag2'] = np.flip(mag2_avgs)
    avg_data['mag3'] = np.flip(mag3_avgs)
    avg_data['Xmag'] = np.flip(Xmag_avgs)
    avg_data['var_e'] = np.flip(var_e_avgs)
    avg_data['var_spin'] = np.flip(var_spin_avgs)
    data = {'raw_data': raw_data, 'avg_data': avg_data, 'sro': SRO, 'enrg_chem': enrg_Chem}
    return data

def get_file_data_2(rootdir, atom_numbs):
    sims = []
    raw_data = {}
    count = 0
    for subdir, dirs, files in os.walk(rootdir):
        for file in files:
            if "OUTPUT" in file:
                count += 1
                sim = pd.read_csv(os.path.join(subdir, file), header=8)
                sim['temp'] = sim['temp'].str.replace('#', '').astype(float)
                sims.append(sim)
                raw_data['temp_'     + str(count)] = np.flip(sim['temp'])
                raw_data['Cmag_'     + str(count)] = np.flip(sim[' Cmag'])
                raw_data['enrg_'     + str(count)] = np.flip(sim[' enrg'])
                raw_data['mag1_'     + str(count)] = np.flip(sim[' mag'])
                raw_data['mag2_'     + str(count)] = np.flip(sim[' magMn'])
                raw_data['Xmag_'     + str(count)] = np.flip(sim[' Xmag'])
                raw_data['var_e_'    + str(count)] = np.flip(sim[' var_e'])
                raw_data['var_spin_' + str(count)] = np.flip(sim[' var_spin'])
                sim_again = open(os.path.join(subdir, file))
                sim_lines = sim_again.readlines()
                raw_data["numbs_"    + str(count)] = [int(sim_lines[1].split()[1].strip(',')), int(sim_lines[1].split()[2].strip(',')), int(sim_lines[1].split()[3].strip(','))]
                raw_data["Tenrg_"    + str(count)] = float(sim_lines[7].split(',')[0])
                raw_data["Menrg_"    + str(count)] = float(sim_lines[7].split(',')[1])
                raw_data["Cenrg_"    + str(count)] = raw_data["Tenrg_" + str(count)] - raw_data["Menrg_" + str(count)]
                sim_again.close()
    raw_data['count'] = count
    Cmag_avgs = []
    enrg_avgs = []
    mag1_avgs = []
    mag2_avgs = []
    mag3_avgs = []
    Xmag_avgs = []
    var_e_avgs = []
    var_spin_avgs = []
    enrg_Chem = np.mean([raw_data["Cenrg_" + str(i+1)] for i in range(count)])

    for temp_ind in range(len(sims[0]['temp'])):
        Cmag_avg = 0
        enrg_avg = 0
        mag1_avg = 0
        mag2_avg = 0
        mag3_avg = 0
        Xmag_avg = 0
        var_e_avg = 0
        var_spin_avg = 0
        for i in range(len(sims)):
            Cmag_avg += sims[i][' Cmag'][temp_ind]/len(sims)
            enrg_avg += sims[i][' enrg'][temp_ind]/len(sims)
            mag1_avg += abs(sims[i][' mag'][temp_ind]/len(sims))*3456/2903.0
            mag2_avg += abs(sims[i][' magMn'][temp_ind]/len(sims))*3456/1175.0
            Xmag_avg += sims[i][' Xmag'][temp_ind]/len(sims)
            var_e_avg += sims[i][' var_e'][temp_ind]/len(sims)
            var_spin_avg += sims[i][' var_spin'][temp_ind]/len(sims)
        Cmag_avgs.append(Cmag_avg)
        enrg_avgs.append(enrg_avg)
        mag1_avgs.append(mag1_avg)
        mag2_avgs.append(mag2_avg)
        mag3_avgs.append(np.sqrt(mag3_avg))
        Xmag_avgs.append(Xmag_avg)
        var_e_avgs.append(var_e_avg)
        var_spin_avgs.append(var_spin_avg)
    avg_data = {}
    avg_data['temp'] = np.flip([sims[0]['temp'][i] for i in range(len(sims[0]['temp']))])
    avg_data['Cmag'] = np.flip(Cmag_avgs)
    avg_data['enrg'] = np.flip(enrg_avgs)
    avg_data['mag1'] = np.flip(mag1_avgs)
    avg_data['mag2'] = np.flip(mag2_avgs)
    avg_data['mag3'] = np.flip(mag3_avgs)
    avg_data['Xmag'] = np.flip(Xmag_avgs)
    avg_data['var_e'] = np.flip(var_e_avgs)
    avg_data['var_spin'] = np.flip(var_spin_avgs)
    data = {'raw_data': raw_data, 'avg_data': avg_data, 'enrg_chem': enrg_Chem}
    return data


def calc_entropy(data):
    s = 0
    S = []
    T = []
    for i in range(len(data['avg_data']['Cmag']) - 1):
        t2 = data['avg_data']['temp'][i + 1]
        t1 = data['avg_data']['temp'][i]
        C2 = data['avg_data']['Cmag'][i+1]
        C1 = data['avg_data']['Cmag'][i]
        s += (1 / 2 * (C2 / t2 - C1 / t1) * (t2 - t1) + (t2 - t1) * C1 / t1)
        S.append(s)
        T.append((t2+t1)/2)
    S_calc = {'temp': T, 'S': S}
    data['calc_vals'] = S_calc

def calc_free_enrg(data, enrg_offset):
    enrg = [(data['avg_data']['enrg'][i+1]+data['avg_data']['enrg'][i])/2 for i in range(len(data['avg_data']['enrg'])-1)]
    G = [enrg[i] + enrg_offset - data['calc_vals']['temp'][i]*data['calc_vals']['S'][i] for i in range(len(enrg))]
    data['calc_vals']['G'] = G

def calc_avg_e0(rootdir, atom_numbs):
    sros = []
    enrgs = []
    count = 0
    for subdir, dirs, files in os.walk(rootdir):
        for file in files:
            if "OUTPUT" in file:
                count += 1
                file = open(os.path.join(subdir, file), 'r')
                lines = file.readlines()
                for i in range(len(lines)):
                    if i > 10:
                        line = lines[i]
                        line = line.split()
                        sros.append(float(line[1].strip(',')))
                        enrgs.append(float(line[3]))
    data = [0, 0]
    data[0] = np.average(sros)
    data[1] = np.average(enrgs) / atom_numbs
    return data

