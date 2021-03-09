#!/usr/bin/python3
import numpy as np
import glob
import matplotlib.pyplot as plt
import math
import scipy
from scipy.optimize import curve_fit

Omega=1.0

if not glob.glob('./data/*.npz'):

    resultList =glob.glob('/users/stud/ledwon/Documents/npzFiles/*.np[yz]')
    data=np.load(resultList[0])
    squaredQ=np.zeros(np.size(data['Q']))
    squaredP=np.zeros(np.size(data['P']))
    t=data['t']
    

    stdMat = np.zeros((len(t),len(resultList)))  

    i=0
    for file in resultList:
        results = np.load(file)
        stdMat[:,i] = np.power(results['Q'],2) - np.average(results[Q])
        varQ += np.power(results['Q'],2) - np.average(results['Q'])
        varP += np.power(results['P'],2) - np.average(results['P'])
        i+=1
        print(i) 

    std=np.std(stdMat, axis=1)
    norm=1.0/(float(len(resultList)))
    varQ *= norm
    varP *= norm
    std *= norm

    np.savez("./data/data", varQ=varQ,  t=t, std=std, varP=varP )


else:
    
    datafile=glob.glob('/users/stud/ledwon/Documents/heatbath/data/*.npz')
    data=np.load(datafile[0])
    varQ=data['varQ']
    varP=data['varP']
    std=data['std']
    t=data['t']


varQ = plt.figure(1)
plt.plot(t,varQ)
varQ.savefig("./img/varQ.pdf")
    



