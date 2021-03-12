#!/usr/bin/python3
import numpy as np
import glob
import matplotlib.pyplot as plt
import math
import scipy
from scipy.optimize import curve_fit

Omega=1.0

if not glob.glob('../data/*.npz'):

    resultList =glob.glob('/users/stud/ledwon/Documents/npzFiles/*.np[yz]')
    data=np.load(resultList[0])
    varQ=np.zeros(np.size(data['Q']))
    varP=np.zeros(np.size(data['P']))
    t=data['t']
    

#    stdMat = np.zeros((len(t),len(resultList)))  

    i=0
    for file in resultList:
        results = np.load(file)
#        stdMat[:,i] = np.power(results['Q'],2) - np.average(results[Q])
        varQ += np.power(results['Q'],2) - np.average(results['Q'])
        varP += np.power(results['P'],2) - np.average(results['P'])
#        results.close()
        i+=1
        print(i) 

#    std=np.std(stdMat, axis=1)
    norm=1.0/(float(len(resultList)))
    varQ *= norm
    varP *= norm
    
########
    std=0
    std *= norm


    np.savez("../data/data", varQ=varQ,  t=t, std=std, varP=varP )


else:
    
    datafile=glob.glob('/users/stud/ledwon/Documents/data/*.npz')
    data=np.load(datafile[0])
    varQ=data['varQ']
    varP=data['varP']
    std=data['std']
    t=data['t']


nSaved=5000
ds=(t[-1]-t[0])/nSaved
plotTimes=ds*range(1,len(varQ)+1)

print(len(plotTimes),len(varQ))

varQPlot = plt.figure(1)
plt.plot(varQ)
varQPlot.savefig("./img/varQ.pdf")
    



