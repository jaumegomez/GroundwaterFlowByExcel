#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 10:31:53 2023

@author: jaime

This code implements the unconfined aquifer model described in the paper by Gómez-Hernández and Secci (2023)
"""

import numpy as np
import flopy

import matplotlib.pyplot as plt
import flopy.utils.binaryfile as bf

"""
The executable mf2005 must be accessible from within the script.
In my case, it is located in folder /usr/local/bin, which is not, by default included in the PATH environmental variable
The following lines check if mf2005 is available and if not it gives a warning and stops the script
To solve it, the folder path must be included in the PATH variable as is done here with /usr/local/bin 
This section can be removed if mf2005 is already in your PATH
"""

import os
import shutil

modflow_exe = 'mf2005'
path2add = '/usr/local/bin'
oldPath = os.environ.get('PATH')
newPath = f'{oldPath}:{path2add}'
os.environ['PATH']=newPath

mf2005_exe = shutil.which(modflow_exe)
if mf2005_exe is None:
    raise Exception("The executable mf2005 is not in your PATH. Make sure that it has been compiled and that its PATH is included as in the lines above")

"""
End of preliminary section to ensure that mf2005 is in the PATH
"""

modelname="UncPlanPaper"
model_ws = "temp"

mf = flopy.modflow.Modflow(modelname,exe_name=modflow_exe,model_ws=model_ws)


#Basic geometry
nrow = 19
ncol = 33
nlay = 1
Lx = 3300.0
Ly = 1900.0
delr = Lx / ncol
delc = Ly / nrow

#Load data (using the same variable names as the names of the spreadsheets)
#Active and bottom cells. In this case, both files are full so values are simply loaded as text
#   specifying the tab delimiter and reshaped into an np array of nlay, nrow and ncol
#   Values should have been written, layer by layer, then row by row
i = np.loadtxt("i.txt",delimiter='\t').reshape(nlay,nrow, ncol)
#Prescribed head values. This array contains only the values of heads at prescribed head boundaries
#   it is loaded using the genfromtxt function that will read a file with missing values
#   and that will replace the missing values by the filling_values; delimiter must be specified
hfix = np.genfromtxt("hfix.txt",delimiter='\t',filling_values=-99).reshape(nlay,nrow, ncol)
#Transmissivity
T = np.genfromtxt("T.txt",delimiter='\t',filling_values=0).reshape(nlay,nrow, ncol)
#Well extraction
W = np.genfromtxt("W.txt",delimiter='\t',filling_values=0).reshape(nlay,nrow, ncol)
#River stage, riverbed bottom elevation and river conductance
hR = np.genfromtxt("hR.txt",delimiter='\t',filling_values=0).reshape(nlay,nrow, ncol)
hB = np.genfromtxt("hB.txt",delimiter='\t',filling_values=0).reshape(nlay,nrow, ncol)
R = np.genfromtxt("R.txt",delimiter='\t',filling_values=0).reshape(nlay,nrow, ncol)
#Recharge (in m3/d)
QN = np.genfromtxt("QN.txt",delimiter='\t',filling_values=0).reshape(nlay,nrow, ncol)

#Heads depend on transmissivities and so do flows
#    only the water velocity will be affected by the thickness of the aquifer
bottom = 0  #array with nlay x nrow x ncol elements with the bottom elevation
ztop = 10 #arbitrarily set the top elevation 40 m above the bottom


#The model is at steady-state so the number of stress periods is 1, its lenght is 1, the number of steps is 1 and 
#   the steady-state flag is set to True
nper = 1
perlen = [1]
nstp = [1]
steady = [True]

#Specify space and time discretization
dis = flopy.modflow.ModflowDis(
    mf, nlay, nrow, ncol, delr=delr, delc=delc, top=ztop, botm=bottom,
    nper=nper, perlen=perlen,nstp=nstp,steady=steady
)

#Specify active and prescribed head boundary conditions in array ibound.
# set starting head for iterating at 100 and heads at noflowing cells at -99.00
# for ibound, first copy array A_i onto ibound, then modify ibound with the values in A_hfix that are not -99, which is the fillin_values used when loading it
ibound = np.copy(i)
ibound[hfix != -99] = -1
strt = 100
hnoflo = -99.00

bas = flopy.modflow.ModflowBas(mf,ibound=ibound,strt=strt,hnoflo=hnoflo)

#Read conductivities along rows and along columns. Vertical conductivity set to 1e-05 but never used since there is only one layer
hk = T/(ztop-bottom)

#set to 1, the values that are nan resulting from 0/0)
vka = 0.00001

#Note the need to specify that the layer is convertible laytyp=1 and wettable laywet=1 and the wetdry factor 
#   (see manual for the meaning of all the wetdry options)
#   notice that chani is set to a negative value meaning that the anisotropy is different for each cell and is specified
#   in the array hani, ipackcb is also defined to ensure that cell-by-cell fluxes are saved
lpf = flopy.modflow.ModflowLpf(mf,hk=hk,vka=vka,ipakcb=53,laytyp=0)

#WellData is not input as an array but as a dictionary, with key equal to the stress period
#   wel_data stacks a number of row vectors of equal length and transposes them into columns.
#   This can be done thanks to the versatily of np arrays that can be indexed by passing
#   vectors with the index values and returns a vector with the looked up values.
#   But first you must extract the layer, row and column indices where the wells are located

lays, rows, cols = np.where(W != 0)
#change sign to make extraction negative recharges as required by the wel package
wel_data = np.column_stack((lays, rows, cols, -W[lays, rows, cols]))
stress_period_data = {0:wel_data}
wel = flopy.modflow.ModflowWel(mf,stress_period_data=stress_period_data)

#RiverData
#   same approach for the river data

lays, rows, cols = np.where(R != 0)
riv_data = np.column_stack((lays, rows, cols, hR[lays, rows, cols],R[lays, rows, cols],hB[lays, rows, cols]))
stress_period_data= {0:riv_data}
riv = flopy.modflow.ModflowRiv(mf,stress_period_data=stress_period_data)

#RichargeData
#   Here you must specify infiltration rate, therefore the input flows in file D_QN have to 
#   be scaled down to infiltration rates
rech = QN/delr/delc
rch = flopy.modflow.ModflowRch(mf,nrchop=2,rech=rech)

#Output Control
#   Finally, the steps for which heads and flows must be saved is indicated
#   It also uses a dictionary as input, but now it must be specified the 
#   stress period and the stress step

stress_period_data = {(0,0): ["save head", "save drawdown","save budget",
                              "print head","print budget" ]}
oc = flopy.modflow.ModflowOc(
    mf, stress_period_data=stress_period_data, compact=True
)

#Solver
pcg = flopy.modflow.ModflowPcg(mf,mxiter=1000,relax=0.8)

#Write input files to the model_ws directory
mf.write_input()

#Run model
success, buff = mf.run_model()
if not success:
    raise Exception("MODFLOW did not terminate normally.")
    

# Extract the heads

hds = bf.HeadFile(f"{model_ws}/{modelname}.hds") #the f infront of the string allows 
                                       #replacing any variable in curly bracekts
                                       #by its value
head = hds.get_data(totim=1.0)

cbb = bf.CellBudgetFile(f"{model_ws}/{modelname}.cbc")
rec_list = cbb.list_records()
# get flow rates between adjacent cells
frf = cbb.get_data(text='FLOW RIGHT FACE')[0]
fff = cbb.get_data(text='FLOW FRONT FACE')[0]
coh = cbb.get_data(text='   CONSTANT HEAD')[0]
#flf = cbb.get_data(text='FLOW LOWER FACE')[0]

# Contour the heads

# extent defines horizontal min, horizontal max, vertical min, vertical max
# the way is written is to flip the image so that the first row of the head
# matrix is at the top

# The most common way to start a plot in matplotlib is with
# fig, ax = plt.subplots()


extent = (delr / 2.0, Lx - delr / 2.0, Ly - delc / 2.0, delc / 2.0)
fig = plt.figure(figsize=(10, 10))
#him=ax.imshow(head[0, :, :], vmin=87, vmax=100)   
#plt.colorbar(him)
ax = fig.add_subplot(1, 1, 1, aspect="equal")  #the subplot will take the 
                                               #"index " position (3rd arg) in 
                                               # grid with nrows (1st arg) and
                                               # ncols (2nd arg). Index starts
                                               # at 1
                                               # aspect "equal" means same 
                                               # scaling for x and y

                                               
h = head[0, :, :] 
h = np.ma.array(h, mask= h<0)

hmin = np.floor(np.min(h))
hmax = np.ceil(np.max(h))
nbcont = 15

ax.set_aspect(aspect="equal") #already specified when creating ax above
pmv = flopy.plot.PlotMapView(model=mf,ax=ax)
pmv.plot_ibound()
pmv.plot_grid()
cs=pmv.contour_array(head, levels=np.linspace(hmin, hmax,nbcont),masked_values=[hnoflo],linewidth=0.2)
plt.clabel(cs,fontsize=8)
#pmv.plot_bc(package=chd,color="red") 
pmv.plot_bc(package=riv,color="yellow")   

plt.show()

#Another contour map

fig = plt.figure(figsize=(10, 10*Ly/Lx))
ax = fig.add_subplot(1, 1, 1, aspect="equal") 
x = np.linspace(0, h.shape[1]-1, h.shape[1])

#y coordinates are flipped so that the first row is drawn on top
y = np.linspace(h.shape[0]-1, 0, h.shape[0])
X, Y = np.meshgrid(x, y)

# Create a banded colormap (you can customize this)
n_bands = 10
cmap = plt.get_cmap('viridis', n_bands)

# Create filled contour plot
contourf = ax.contourf(X, Y, h, levels=n_bands, cmap=cmap)

# Create contour lines
contour = ax.contour(X, Y, h, levels=n_bands, colors='black')

# Add labels to contour lines
ax.clabel(contour, inline=True, fontsize=10)

# Add colorbar
plt.colorbar(contourf)

# Add labels and title
ax.set_title('Contour Map')
ax.set_xlabel('X-axis')
ax.set_ylabel('Y-axis')

#plt.gca().set_aspect('auto')

plt.show()
    