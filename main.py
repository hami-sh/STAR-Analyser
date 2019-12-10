import os
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np
import pandas as pd

pcawg_pe_sig = defaultdict(dict)
pcawg_pe_noise = defaultdict(dict)
pcawg_pe_overhang = defaultdict(dict)
pcawg_se_sig = defaultdict(dict)
pcawg_se_noise = defaultdict(dict)
pcawg_se_overhang = defaultdict(dict)

tcga_pe_sig = defaultdict(dict)
tcga_pe_noise = defaultdict(dict)
tcga_pe_overhang = defaultdict(dict)
tcga_se_sig = defaultdict(dict)
tcga_se_noise = defaultdict(dict)
tcga_se_overhang = defaultdict(dict)

arr_pcawg_pe_0 = []
arr_pcawg_se_0 = []
arr_tcga_pe_0 = []
arr_tcga_se_0 = []


len_pcawg_pe = 0
len_pcawg_se = 0
len_tcga_pe = 0
len_tcga_se = 0


def main():
    file_handler()
    grapher()

def grapher():
    # ----------PCAWG PE
    # sig
    graph(pcawg_pe_sig, "pcawg_pe_sig", "PCAWG PE Unique Maps", len_pcawg_pe)
    # noise
    graph(pcawg_pe_noise, "pcawg_pe_noise", "PCAWG PE MultiMaps", len_pcawg_pe)
    # overhang
    graph(pcawg_pe_overhang, "pcawg_pe_overhang", "PCAWG PE Maximum Overhang", len_pcawg_pe)
    # ----------PCAWG SE
    # sig
    graph(pcawg_se_sig, "pcawg_se_sig", "PCAWG SE Unique Maps", len_pcawg_se)
    # noise
    graph(pcawg_se_noise, "pcawg_se_noise", "PCAWG SE MultiMaps", len_pcawg_se)
    # overhang
    graph(pcawg_se_overhang, "pcawg_se_overhang", "PCAWG SE Maximum Overhang", len_pcawg_se)

    # TCGA PE
    # sig
    graph(tcga_se_sig, "tcga_se_sig", "TCGA SE Unique Maps", len_tcga_se)
    # noise
    graph(tcga_se_noise, "tcga_se_noise", "TCGA SE MultiMaps", len_tcga_se)
    # overhang
    graph(tcga_se_overhang, "tcga_se_overhang", "TCGA SE Maximum Overhang", len_tcga_se)

    # TCGA SE
    # sig
    graph(tcga_pe_sig, "tcga_pe_sig", "TCGA PE Unique Maps", len_tcga_pe)
    # noise
    graph(tcga_pe_noise, "tcga_pe_noise", "TCGA PE MultiMaps", len_tcga_pe)
    # overhang
    graph(tcga_pe_overhang, "tcga_pe_overhang", "TCGA PE Maximum Overhang", len_tcga_pe)


def graph(dictionary, name, title, nvalue):
    x = list(dict.fromkeys(dictionary.keys()))
    print("-------------------------------")

    one = list()
    two = list()
    other = list()
    keys = dictionary.keys()
    for item in keys:
        one.append(int(dictionary[item]['1']))
        two.append(int(dictionary[item]['2']))
        other.append(int(dictionary[item]['mid']))
    print(one)
    print(two)
    print(other)
    print("n = " + str(len(keys)) + " out of " + str(nvalue))
    df = pd.DataFrame(np.c_[one, two, other], index=x)
    ax = df.plot(kind='bar', title=title)
    ax.text(0.10, 1.1, 'n = ' + str(len(keys)) + " out of " + str(nvalue), horizontalalignment='center',
     verticalalignment='center',
     transform = ax.transAxes)
    plt.ylabel('Frequency')
    plt.legend(["1st Intron", "2nd Intron", "Other"])
    plt.savefig('plots/' + name + '.png', dpi=300)
    plt.show()

def file_handler():
    # handle those with transcript
    for filename in os.listdir('files'):
        # print("-----------------------")
        # print(filename)
        if filename != '0':
            fp = open('files/' + filename)
            for cnt, line in enumerate(fp):
                # print("Line {}: {}".format(cnt, line))
                # split line on space
                details = line.split()
                # print(details)
                details[0] = details[0].strip(".SJ.out.tab:chr8").split("/")[0]
                if 'pcawg' in filename:
                    # add to pcawg dict
                    if 'pe' in filename:
                        pcawg_pe_sig[details[0]][filename.split("_")[-1].split(".")[0]] = details[1]
                        pcawg_pe_noise[details[0]][filename.split("_")[-1].split(".")[0]] = details[2]
                        pcawg_pe_overhang[details[0]][filename.split("_")[-1].split(".")[0]] = details[3]

                    if 'se' in filename:
                        pcawg_se_sig[details[0]][filename.split("_")[-1].split(".")[0]] = details[1]
                        pcawg_se_noise[details[0]][filename.split("_")[-1].split(".")[0]] = details[2]
                        pcawg_se_overhang[details[0]][filename.split("_")[-1].split(".")[0]] = details[3]
                if 'tcga' in filename:
                    # add to tcga dict
                    if 'pe' in filename:
                        tcga_pe_sig[details[0]][filename.split("_")[-1].split(".")[0]] = details[1]
                        tcga_pe_noise[details[0]][filename.split("_")[-1].split(".")[0]] = details[2]
                        tcga_pe_overhang[details[0]][filename.split("_")[-1].split(".")[0]] = details[3]
                    if 'se' in filename:
                        tcga_se_sig[details[0]][filename.split("_")[-1].split(".")[0]] = details[1]
                        tcga_se_noise[details[0]][filename.split("_")[-1].split(".")[0]] = details[2]
                        tcga_se_overhang[details[0]][filename.split("_")[-1].split(".")[0]] = details[3]

    populate(pcawg_pe_sig)
    populate(pcawg_pe_noise)
    populate(pcawg_pe_overhang)
    populate(pcawg_se_sig)
    populate(pcawg_se_noise)
    populate(pcawg_se_overhang)
    populate(tcga_pe_sig)
    populate(tcga_pe_noise)
    populate(tcga_pe_overhang)
    populate(tcga_se_sig)
    populate(tcga_se_noise)
    populate(tcga_se_overhang)

    for filename in os.listdir('files/0'):
        fp = open('files/0/' + filename)
        for cnt, line in enumerate(fp):
            details = line.split()
            # print(details)
            details[0] = details[0].strip(".SJ.out.tab:chr8").split("/")[0]
            if 'pcawg' in filename:
                if 'pe' in filename:
                    arr_pcawg_pe_0.append(details[0])
                if 'se' in filename:
                    arr_pcawg_se_0.append(details[0])
            if 'tcga' in filename:
                if 'pe' in filename:
                    arr_tcga_pe_0.append(details[0])
                if 'se' in filename:
                    arr_tcga_se_0.append(details[0])

    total_pcawg = list(dict.fromkeys(arr_pcawg_pe_0)) + list(pcawg_pe_sig.keys())
    res = []
    for i in total_pcawg:
        if i not in res:
            res.append(i)
    global len_pcawg_pe
    len_pcawg_pe = len(res)
    print(len_pcawg_pe)

    print("+")
    se = list(pcawg_se_sig.keys()) + arr_pcawg_se_0
    res = []
    for i in se:
        if i not in res:
            res.append(i)
    global len_pcawg_se
    len_pcawg_se = len(res)
    print(len_pcawg_se)
    print("----------")
    tcga_pe = list(tcga_pe_sig.keys()) + list(dict.fromkeys(arr_tcga_pe_0))
    res = []
    for i in tcga_pe:
        if i not in res:
            res.append(i)
    global len_tcga_pe
    len_tcga_pe = len(res)
    print(len_tcga_pe)
    print("+")
    tcga_se = list(tcga_se_sig.keys()) + list(dict.fromkeys(arr_tcga_se_0))
    res = []
    for i in tcga_se:
        if i not in res:
            res.append(i)
    global len_tcga_se
    len_tcga_se = len(res)
    print(len_tcga_se)

def populate(dictionary):
    for key in dictionary.keys():
        if '1' not in dictionary[key]:
            dictionary[key]['1'] = '0'
        if '2' not in dictionary[key]:
            dictionary[key]['2'] = '0'
        if 'mid' not in dictionary[key]:
            dictionary[key]['mid'] = '0'


if __name__ == "__main__":
    main()