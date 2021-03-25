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
    t=data['dts']*range(0,len(varQ)+1)
    gamma=data['gamma']
    

#    stdMat = np.zeros((len(t),len(resultList)))  

    i=0
    for file in resultList:
        results = np.load(file)
#        stdMat[:,i] = np.power(results['Q'],2) - np.average(results[Q])
        varQ += np.power(results['Q'],2) - np.average(results['Q'])
        varP += np.power(results['P'],2) - np.average(results['P'])
        i+=1
        print(i) 
        print(results['avgEnergyError'],results['maxEnergyError'])

#    std=np.std(stdMat, axis=1)
    norm=1.0/(float(len(resultList)))/4000.0
    varQ *= norm
    varP *= norm
    
    std=0
    std *= norm

    varQ=np.append(0,varQ)

    np.savez("../data/data", varQ=varQ,  t=t, std=std, varP=varP )


else:
    
    datafile=glob.glob('/users/stud/ledwon/Documents/data/*.npz')
    data=np.load(datafile[0])
    varQ=data['varQ']
    varP=data['varP']
    std=data['std']
    t=data['t']
    gamma=['gamma']

def theoDiff(x,a):
    return a*np.power(x,1.2)

startindex = int(math.floor(varQ.size*0.4))
fitTimes=t[startindex:-1]
fitVarQ=varQ[startindex:-1]
popt, pcov = curve_fit(theoDiff,fitTimes,fitVarQ)

varQPlot = plt.figure(1)
plt.loglog(t,varQ)
plt.loglog(fitTimes,fitVarQ)

varQPlot.savefig("./img/varQ.pdf")
    



