import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import re
import os

def readForces(fileName,skip):
    """Read forces or moment from a force file using the forces function object in openfoam"""

    f = open(fileName, "r")
    lines = f.readlines()
    f.close()

    x=[]

    for i in lines: #:Strip all (), , and [] from the file
        if '#' in i:
            continue
        else:
            i=(re.sub(' +', ',',re.sub('[\t]', ' ',re.sub('[()]', '',i.rstrip(os.linesep)))))
            x.append(i)

    del x[0:skip]   #: Delete the values below the skip defined in the gui

    data = np.empty((0,4))

    for i in x:
        d=np.fromstring(i, dtype=float, sep=',')
        #print(d)
        data = np.append(data,np.array([[d[0],d[1],d[2],d[3]]]),axis=0)

    return data

def getEvery(arr,rows,step):
    return arr.reshape(rows, step)

def averageLast10inEachRow(arr):
    lst = [] 
    for i in arr:
        me = np.mean(i[-11:])
        lst.append(me)
    return np.array(lst)

def getNumberFromFile(filename,searchString):
    f1 = open(filename, 'r')  # open the file for reading
    data = f1.readlines()  # read the entire file as a list of strings
    f1.close()
    for line in data:
        match = re.search(searchString+r' \s*(\d+(?:\.\d+)?)', " ".join(line.split()), re.IGNORECASE)
        if match:
            return float(match.group(1))
        pass


# input

rads = getNumberFromFile('./constant/MRFProperties','omega')
density = getNumberFromFile('./system/forces','rhoInf')
flowPoints = 10
ittPrFlow = 200

# Data for plotting
# Read the moment
m = readForces('./postProcessing/forces/0/moment.dat',0)
mSort = getEvery(m[:,3],flowPoints,ittPrFlow)
mFinal = averageLast10inEachRow(mSort)

# Read the flowRate
fr = np.loadtxt('./postProcessing/flowRatePatch/0/surfaceFieldValue.dat',comments='#')
frSort = getEvery(fr[:,1],flowPoints,ittPrFlow)
frFinal = averageLast10inEachRow(frSort)*3600 #m3/h
print(frFinal)

# Read the pressure difference
pd = np.loadtxt('./postProcessing/pressureDifferencePatch/0/fieldValueDelta.dat',comments='#')
pdSort = getEvery(pd[:,1],flowPoints,ittPrFlow)
pdFinal = averageLast10inEachRow(pdSort)*density*-1*0.00010197 # Convert to absolute pressure and mH2O
print(pdFinal)


# Shaft power
power = (np.abs(mFinal) * rads)/1000.0 #power in kW
print(power)

# Hydraulic power
hydPower = (frFinal * density * 9.81 * pdFinal) / 3.6e06

# Efficieny
eff = (hydPower/power)*100 # %
print(eff)

fig, (ax1, ax2, ax3) = plt.subplots(3,1,figsize=(11.69,8.27),sharex=True)
matplotlib.rc('xtick', labelsize=14)     
matplotlib.rc('ytick', labelsize=14)

ax1.plot(frFinal,pdFinal,label="Head (m)",linestyle='--', marker='o')
ax2.plot(frFinal,power*1000,label="Power (W)",linestyle='--', marker='o')
ax3.plot(frFinal,eff,label="Eff (%)",linestyle='--', marker='o')

ax1.set(ylabel='Head (m)')
ax2.set(ylabel='Power (W)')
ax3.set(xlabel='Flow (m3/h)', ylabel='Eff (%)')
ax1.minorticks_on()
ax2.minorticks_on()
ax3.minorticks_on()
ax1.grid(True, which='both')
ax2.grid(True, which='both')
ax3.grid(True, which='both')
ax1.legend(loc='upper right',ncol=1)
ax2.legend(loc='upper right',ncol=1)
ax3.legend(loc='upper right',ncol=1)
ax1.yaxis.label.set_size(16)
ax1.xaxis.label.set_size(16)
ax1.title.set_size(20)
ax2.yaxis.label.set_size(16)
ax2.xaxis.label.set_size(16)
ax2.title.set_size(20)
ax3.yaxis.label.set_size(16)
ax3.xaxis.label.set_size(16)
ax3.title.set_size(20)
#fig.savefig('pumpcurve.pdf', dpi = 300)
plt.show()