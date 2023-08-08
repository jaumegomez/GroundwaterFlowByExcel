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

# 0. Load data (using the same variable names as the names of the spreadsheets)
#Active and bottom cells. In this case, both files are full so values are simply loaded as text
#   specifying the tab delimiter and reshaped into an np array of nlay, nrow and ncol
#   Values should have been written, layer by layer, then row by row


nlay, nrow, ncol = 51, 1, 113  # specify the number of layers, rows, and columns

A_i = np.genfromtxt("A_i.txt",delimiter='\t',filling_values=0).reshape(nlay,nrow, ncol)
A_Bot = np.loadtxt("A_Bot.txt",delimiter='\t').reshape(nlay,nrow, ncol)
#Prescribed head values. This array contains only the values of heads at prescribed head boundaries
#   it is loaded using the genfromtxt function that will read a file with missing values
#   and that will replace the missing values by the filling_values; delimiter must be specified
A_hfix = np.genfromtxt("A_hfix.txt",delimiter='\t',filling_values=-99).reshape(nlay,nrow, ncol)
#Conductivity along rows
A_Kx = np.genfromtxt("A_Kx.txt",delimiter='\t',filling_values=0).reshape(nlay,nrow, ncol)
#Conductivity along columns
A_Kz = np.genfromtxt("A_Kz.txt",delimiter='\t',filling_values=0).reshape(nlay,nrow, ncol)
#Well extraction
#A_W = np.genfromtxt("A_W.txt",delimiter='\t',filling_values=0).reshape(nlay,nrow, ncol)
#River stage, riverbed bottom elevation and river conductance
#A_hR = np.genfromtxt("A_hR.txt",delimiter='\t',filling_values=0).reshape(nlay,nrow, ncol)
#A_hB = np.genfromtxt("A_hB.txt",delimiter='\t',filling_values=0).reshape(nlay,nrow, ncol)
#A_R = np.genfromtxt("A_R.txt",delimiter='\t',filling_values=0).reshape(nlay,nrow, ncol)
#Recharge (in m3/d)
#D_QN = np.genfromtxt("D_QN.txt",delimiter='\t',filling_values=0).reshape(nlay,nrow, ncol)

# 1. Setup the MODFLOW modelmodelname="UncPlanPaper"
modelname = "UncXSectPaper"
model_ws = "temp"
mf = flopy.modflow.Modflow(modelname,exe_name=modflow_exe,model_ws=model_ws)

# 2. Define the Discretization (DIS) package
delr = 5.0
delc = 1000.0  # grid spacing along rows and columns
delz = 2.0
top = A_Bot[0,:,:] + delz  # top of the model
botm = A_Bot  # array of bottom elevations for each layer
Lx = ncol * delr
Ly = nlay * delz
#The model is at steady-state so the number of stress periods is 1, its lenght is 1, the number of steps is 1 and 
#   the steady-state flag is set to True
nper = 1
perlen = [1]
nstp = [1]
steady = [True]
dis = flopy.modflow.ModflowDis(mf, nlay, nrow, ncol, delr=delr, delc=delc, top=top, botm=botm)

# 3. Define the Basic (BAS) package
#Specify active and prescribed head boundary conditions in array ibound.
# set starting head for iterating at 100 and heads at noflowing cells at -99.00
# for ibound, first copy array A_i onto ibound, then modify ibound with the values in A_hfix that are not -99, which is the fillin_values used when loading it
ibound = np.copy(A_i)
ibound[A_hfix != -99] = -1 # boundary conditions, 1 is active, -1 is fixed head
strt = np.ones((nlay,nrow,ncol))*80.0 # used to set the prescribed at the boundaries and the initial head to start the iterations
strt[A_hfix != -99] = A_hfix[A_hfix != -99]
hnoflo = -99.00
bas = flopy.modflow.ModflowBas(mf, ibound=ibound, strt=strt, hnoflo=hnoflo)

# 4. Define the Layer-Property Flow (LPF) package
hk = A_Kx
vka = A_Kz  # hydraulic conductivity

ones = np.ones(nlay)  # all layers are convertible (1)  # all layers are rewetable
lpf = flopy.modflow.ModflowLpf(mf, hk=hk, vka=vka, laytyp=ones, laywet=ones, wetdry=ones, iwetit=1,wetfct=0.1,ihdwet=0,ipakcb=53)

#Output Control
#   Finally, the steps for which heads and flows must be saved is indicated
#   It also uses a dictionary as input, but now it must be specified the 
#   stress period and the stress step

stress_period_data = {(0,0): ["save head", "save drawdown","save budget",
                              "print head","print budget" ]}
oc = flopy.modflow.ModflowOc(
    mf, stress_period_data=stress_period_data, compact=True)

# Define the Preconditioned Conjugate-Gradient (PCG) package (solver)
# Set IWETIT to 1 to turn on wetting, and IHDWET to 1 to specify that wetting is done based on cell heads
pcg = flopy.modflow.ModflowPcg(mf,mxiter=50000,iter1=50,relax=0.8)

# 8. Write the MODFLOW model input files
mf.write_input()

# 9. Run the MODFLOW model
success, buff = mf.run_model()
# if not success:
#     raise Exception("MODFLOW did not terminate normally.")

success = True
if success:
    # Load outputs
    hds = bf.HeadFile(f'{model_ws}/{modelname}' + '.hds')
    head = hds.get_data(totim=1.0)

    # Create the plot
    fig = plt.figure(figsize=(10, 10*Ly/Lx))
    ax = fig.add_subplot(1, 1, 1, aspect='equal')

    # Using the PlotCrossSection module
    xsect = flopy.plot.PlotCrossSection(model=mf, line={'Row': 0})  # Cross-section along the first row
    patches = xsect.plot_ibound(color_noflow='red', color_ch='blue')
    linecollection = xsect.plot_grid(color='black')
    cs = xsect.plot_array(head, masked_values=[999.], head=head, alpha=0.5)
    
    # Add a color bar on the right of the plot
    plt.colorbar(cs, shrink=0.5)
    
    # Set plot title
    ax.set_title('Cross-section of head along Row 1')
    
    # Display the plot
    plt.show()
    
import flopy.utils.binaryfile as bf

# After running the model
if success:
    # Load outputs
    hds = bf.HeadFile(f'{model_ws}/{modelname}' + '.hds')
    head = hds.get_data(totim=1.0)

    # Load cell-by-cell budget file
    cbb = bf.CellBudgetFile(f'{model_ws}/{modelname}'+'.cbc')
    # get the list of records (you can print it to see what's available)
    rec_list = cbb.list_records()
    # get flow rates between adjacent cells
    frf = cbb.get_data(text='FLOW RIGHT FACE')[0]
 #   fff = cbb.get_data(text='FLOW FRONT FACE')[0]
    flf = cbb.get_data(text='FLOW LOWER FACE')[0]
    print("Flow Rates:")
    print("Right Face:\n", frf)
#    print("Front Face:\n", fff)
    print("Lower Face:\n", flf)
    
    # compute and print flow velocities
    dx = delr
    dy = delc
    dz = top - botm
    area = dy * dx
    v_right_face = frf / (dz[:, None, None] * dy)  # velocity at the right face of the cell
#    v_front_face = fff / (dz[:, None, None] * dx)  # velocity at the front face of the cell
    v_lower_face = flf / (dz[:, None, None] * area)  # velocity at the lower face of the cell
    print("Flow Velocities:")
    print("Right Face:\n", v_right_face)
#   print("Front Face:\n", v_front_face)
    print("Lower Face:\n", v_lower_face)
    
# Setup contour parameters
Lx = delr * ncol
Ly = delz * nlay
levels = np.linspace(67, 72, 9)
extent = (delr / 2.0, Lx - delr / 2.0, delc / 2.0, Ly - delc / 2.0)
print("Levels: ", levels)
print("Extent: ", extent)
    
fig = plt.figure(figsize=(10, 10*Ly/Lx))
mytimes = [1.0]
for iplot, time in enumerate(mytimes):
    print("*****Processing time: ", time)
    head = hds.get_data(totim=time)
    # Print statistics
    print("Head statistics")
    print("  min: ", head.min())
    print("  max: ", head.max())
    print("  std: ", head.std())

    # Extract flow right face and flow front face
    frf = cbb.get_data(text="FLOW RIGHT FACE", totim=time)[0]
    flf = cbb.get_data(text="FLOW LOWER FACE", totim=time)[0]

    # Create a map for this time

    ax = fig.add_subplot(1,1,1)
    ax.set_title(f"time {time}")

    # pmv = flopy.plot.PlotMapView(model=mf, layer=0, ax=ax)
    # qm = pmv.plot_ibound()
    # lc = pmv.plot_grid()
    # qm = pmv.plot_bc("GHB", alpha=0.5)
    # if head.min() != head.max():
    #     cs = pmv.contour_array(head, levels=levels, linewidths=3.0)
    #     plt.clabel(cs, inline=1, fontsize=10, fmt="%1.1f")
    #     quiver = pmv.plot_vector(frf, -fff)

    xsect = flopy.plot.PlotCrossSection(model=mf, line={'Row': 0})  # Cross-section along the first row
    qm = xsect.plot_ibound()
#    lc = xsect.plot_grid()
#    qm = xsect.plot_bc("GHB", alpha=0.5)
    if head.min() != head.max():
        cs = xsect.contour_array(head, masked_values=[-1.0e+30],levels=levels, linewidths=3.0)
        plt.clabel(cs, inline=1, fontsize=10, fmt="%1.1f")
        quiver = xsect.plot_vector(frf,0, -flf)
