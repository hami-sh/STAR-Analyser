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

def main():
    file_handler()
    grapher()

def grapher():
    # ----------PCAWG PE
    # sig
    graph(pcawg_pe_sig, "pcawg_pe_sig", "PCAWG PE Unique Maps")
    # noise
    graph(pcawg_pe_noise, "pcawg_pe_noise", "PCAWG PE MultiMaps")
    # overhang
    graph(pcawg_pe_overhang, "pcawg_pe_overhang", "PCAWG PE Maximum Overhang")
    # ----------PCAWG SE
    # sig
    graph(pcawg_se_sig, "pcawg_se_sig", "PCAWG SE Unique Maps")
    # noise
    graph(pcawg_se_noise, "pcawg_se_noise", "PCAWG SE MultiMaps")
    # overhang
    graph(pcawg_se_overhang, "pcawg_se_overhang", "PCAWG SE Maximum Overhang")

    # TCGA PE
    # sig
    graph(tcga_se_sig, "tcga_se_sig", "TCGA SE Unique Maps")
    # noise
    graph(tcga_se_noise, "tcga_se_noise", "TCGA SE MultiMaps")
    # overhang
    graph(tcga_se_overhang, "tcga_se_overhang", "TCGA SE Maximum Overhang")

    # TCGA SE
    # sig
    graph(tcga_pe_sig, "tcga_pe_sig", "TCGA PE Unique Maps")
    # noise
    graph(tcga_pe_noise, "tcga_pe_noise", "TCGA PE MultiMaps")
    # overhang
    graph(tcga_pe_overhang, "tcga_pe_overhang", "TCGA PE Maximum Overhang")


def graph(dictionary, name, title):
    x = list(dict.fromkeys(dictionary.keys()))
    print("-------------------------------")
    print(len(x))
    print(x)

    one = list()
    two = list()
    other = list()
    keys = dictionary.keys()
    for item in keys:
        print(item)
        one.append(int(dictionary[item]['1']))
        two.append(int(dictionary[item]['2']))
        other.append(int(dictionary[item]['mid']))
    print(one)
    print(two)
    print(other)
    df = pd.DataFrame(np.c_[one, two, other], index=x)
    ax = df.plot(kind='bar', title=title)
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
                        print(pcawg_pe_sig)
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
    # print(pcawg_pe_sig)
    # pcawg_pe_noise
    # pcawg_pe_overhang
    # pcawg_se_sig
    # pcawg_se_noise
    # pcawg_se_overhang
    #
    # print(tcga_pe_sig)
    # tcga_pe_noise
    # tcga_pe_overhang
    # tcga_se_sig
    # tcga_se_noise
    # tcga_se_overhang


def populate(dictionary):
    for key in dictionary.keys():
        print(dictionary.get(key))
        if '1' not in dictionary[key]:
            dictionary[key]['1'] = '0'
        if '2' not in dictionary[key]:
            dictionary[key]['2'] = '0'
        if 'mid' not in dictionary[key]:
            dictionary[key]['mid'] = '0'


if __name__ == "__main__":
    main()