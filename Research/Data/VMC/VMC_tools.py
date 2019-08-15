##!/usr/bin/env python

import sys
import matplotlib.pyplot as plt
import pandas as pd
import os

def read_raw_data(where):
    files = os.listdir(where)
    for file in files:
        if file[-4:] == ".txt":
            energy = []
            noise = []
            gradient = []
            iteration = []
            time = []
            i = 0
            with open(where+'/'+file,'r') as f:
                lines = f.readlines()
                i += 1
                for line in lines:
                    line_split = line.split()
                    if len(line_split) > 6:
                        iteration.append(float(line_split[0]))
                        energy.append(float(line_split[1]))
                        n = line_split[2][1:-1]
                        noise.append(float(n))
                        gradient.append(float(line_split[3]))
                        time.append(float(line_split[-1]))
                pre_df = {"Iteration":iteration,"Energy":energy,"Noise":noise,"Gradient":gradient,"Time":time}
                df = pd.DataFrame(pre_df)
                df.to_csv("%s/%s.csv" % (where,file[:-4]))

def plot_data(where):
    files = os.listdir(where)
    plt.figure()
    step = 1
    max_iter = 50
    for file in files:
        if file[-4:] == ".csv":
            df = pd.read_csv("%s/%s" % (where,file))
            #plt.plot(df["Energy"])
            plt.ylabel("Energy")
            plt.xlabel("Iteration")
            plt.errorbar(df["Iteration"][:max_iter:step],df["Energy"][:max_iter:step],yerr = df["Noise"][:max_iter:step],capsize = 2,label = file[:-4])
    plt.legend()
    plt.savefig("energy.pdf")
    plt.show()

def plot_grad_norm(where):
    files = os.listdir(where)
    plt.figure()
    step = 1
    for file in files:
        if file[-4:] == ".csv":
            df = pd.read_csv("%s/%s" % (where,file))
            plt.plot(df["Gradient"],label = file[:-4])
            plt.ylabel("Grad_Norm")
            plt.xlabel("Iteration")
    plt.legend()
    plt.savefig('gradnorm.png')
    plt.show()

def conv_energy(where):
    files = os.listdir(where)
    for file in files:
        if file[-4:] == ".csv":
            df = pd.read_csv("%s/%s" % (where,file))
            e = df["Energy"].tolist()
            min_e = min(e)
            index = e.index(min_e)
            print( "%s min_e = %.8f +/- %.8f Iteration: %d | time: %f" % (file[:-4],e[index],df["Noise"][index],index,df["Time"][index]))
