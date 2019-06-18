# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 18:58:04 2017

@author: Razib
"""



############################################################
#Trial 



global E_field_new, ion_accel_length_, ion_drift_length
E_field_new =8770.62823
#ion_drift_length_new = 0.01494322
ion_accel_length = 27.0 *0.01
mass_amu = 1.66056e-27  
eV = 1.602176487e-19
global mass, charge # mass should be in kg and charge in Coulomb
#mass = m* mass_amu # in kg
#charge = q * eV # in Coulomb


testTOF = [6765, 4842, 7617]
testmq = [45, 19, 60]


##############################################################





import numpy as np
import matplotlib.pyplot as plt
from fast_histogram import histogram1d, histogram2d
from scipy.ndimage import center_of_mass
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.optimize as optimize
from scipy.interpolate import interp1d
from lmfit import minimize, Parameters, fit_report


def integrators(x, y, xlim):
    
    '''
    This definition return the index and the integrated y value setforth by
    xlim.
    Usage: index, ysum = integrators(xarray, yarray, xlim)
    
    * xarray and yarray must either be a list or a numpy array
    * xlim must be a list with LOWER RANGE and UPPER RANGE
    '''
    x = np.array(x)
    y = np.array(y)
    
    if x.shape != y.shape:
        print('x and y must be of same sized numpy array. No accumulation done!')
        return None, None
    
    if not isinstance(xlim, list):
        print('xlim must be a list' )
        return None, None

    else:
        idx = np.where((x > xlim[0]) & (x < xlim[1]))
        yIntegration = y[idx].sum()
        return idx, yIntegration
    
    
def pipico(xArray):
    xFirstSlice = []
    xSecondSlice = []
    start = 1
    for i in xArray:
        for j in xArray[start:]:
            xFirstSlice.append(i)
            xSecondSlice.append(j)
        start = start + 1 
    return np.array(xFirstSlice), np.array(xSecondSlice)


def tripico(array):
   list1 = []
   list2 = []
   list3 = []
   
      # c=1
   for i in range(0,len(array)-2):
       for j in range(i+1,len(array)-1):
           for k in range(i+2,len(array)):
               if k<=j:
                   continue
               list1.append(array[i])
               list2.append(array[j])
               list3.append(array[k])
     #  c=c+1
   x=np.asarray(list1)
   y=np.asarray(list2)
   z=np.asarray(list3)
   return x,y,z

def masked(arr, lessThan):
    arr = np.ma.masked_where(arr < lessThan, arr)
    return arr
    

def plotTushTush_old(arr2D):
    
    xAxis, yAxis, zAxis = arr2D.centerData
    zAxis += 1
    zAxis = np.transpose(np.log10(zAxis))
    zAxis_masked = masked(zAxis, np.log10(1.1))
    
    
    fig, ax = plt.subplots()
    
    cmaps = plt.cm.viridis
    cmaps.set_bad(color = 'white')
    ax.pcolormesh(xAxis, yAxis, zAxis_masked, cmap = cmaps)
    
    
def plotTushTush(xAxis, yAxis, zAxis):
        
    zAxis += 1
    zAxis = np.transpose(np.log10(zAxis))
    zAxis_masked = masked(zAxis, np.log10(1.1))
    
    
    fig, ax = plt.subplots()
    
    cmaps = plt.cm.jet
    cmaps.set_bad(color = 'white')
    img = ax.pcolormesh(xAxis, yAxis, zAxis_masked, cmap = cmaps)
    plt.colorbar(img)


def mqConversion(mq, testTOF, actualTOF):
    #for 170 eV
    #tof = [505, 2235, 6825, 20415]
    #testmq = [1, 18, 165, 1469]
    #
    tof = np.array(testTOF)
    sqrtMQ = np.sqrt(mq)
    k, t0 = np.polyfit(sqrtMQ, tof, 1)
    
    convMQ = ((np.array(actualTOF) - t0)/ k )**2
    return convMQ

def tofConversion(testmq, testTOF, mq):
    tof = np.array(testTOF)
    sqrtMQ = np.sqrt(testmq)
    k, t0 = np.polyfit(sqrtMQ, tof, 1)
    
    convTOF = k * np.sqrt(mq) + t0
    return convTOF


def plotproj(x, y, z, xlabel='xlabel', ylabel='ylabel', cmaps='jet', n=100):
    fig = plt.figure(n,figsize=(10, 10))

    grid = plt.GridSpec(8, 8, hspace=0., wspace=0.)
    main_ax = fig.add_subplot(grid[:-1, 1:])
    y_hist = fig.add_subplot(grid[:-1, 0], sharey=main_ax)
    x_hist = fig.add_subplot(grid[-1, 1:], sharex=main_ax)
    
    divider = make_axes_locatable(main_ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    
    # scatter points on the main axes
    im = main_ax.pcolormesh(x, y, z.T, norm=LogNorm(), cmap = cmaps)
    fig.colorbar(im, cax=cax, orientation='horizontal')
    #main_ax.set_xticklabels([])
    #main_ax.set_yticklabels([])
    
    sumy = z.sum(0)
    sumx = z.sum(1)
    
    x_hist.plot(x, sumx)
    x_hist.set_xlabel(xlabel)
     
    y_hist.plot(sumy, y)
    y_hist.set_xlim(sumy.max(), 0)
    y_hist.set_ylabel(ylabel)


class hist1d:
    def __init__(self, xmin, xmax, nbins = 10):
        
        self.nbins = nbins
        self.edges = np.linspace(xmin, xmax, nbins + 1)
        self.centers = (self.edges[:-1] + self.edges[1:])/2.0
        
        self.delta = 0.0
        self.range = (xmin, xmax + self.delta)
        self.hists = histogram1d([], nbins, self.range)
    
    def fill(self, arr):
        hists = histogram1d(arr, self.nbins, self.range)
        self.hists += hists
        
    @property
    def edgeData(self):
        return self.edges, self.hists
    
    @property
    def centerData(self):
        return self.centers, self.hists
    
    
class hist2d:
    def __init__(self, xmin, xmax, ymin, ymax, nxbins = 10, nybins = 10):
        self.nxbins = nxbins
        self.nybins = nybins
        
        self.xedges = np.linspace(xmin, xmax, nxbins + 1)
        self.xcenters = (self.xedges[:-1] + self.xedges[1:])/2.0
        
        self.yedges = np.linspace(ymin, ymax, nybins + 1)
        self.ycenters = (self.yedges[:-1] + self.yedges[1:])/2.0
        
        self.delta = 0.0
        self.xrange = (xmin, xmax + self.delta)
        self.yrange = (ymin, ymax + self.delta)
        self.hists = histogram2d([], [], [self.nxbins, self.nybins], [self.xrange, self.yrange])
        
        
    def fill(self, xarr, yarr):
        hists = histogram2d(xarr, yarr, [self.nxbins, self.nybins], [self.xrange, self.yrange] )
        self.hists += hists
    
    @property
    def edgeData(self):
        return self.xedges, self.yedges, self.hists
        
    @property
    def centerData(self):
        return self.xcenters, self.ycenters, self.hists
        

def init_hist1d(xmin, xmax, binSize):
    nBins = np.int64((xmax - xmin)/binSize)
    return xmin, xmax, nBins

def init_hist2d(xmin, xmax, xbinSize, ymin, ymax, ybinSize):
    xnBins = np.int64((xmax - xmin)/xbinSize)
    ynBins = np.int64((ymax - ymin)/ybinSize)
    return xmin, xmax, ymin, ymax, xnBins, ynBins
    
        
        
def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*np.exp(
                -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

def moments(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    total = data.sum()
    X, Y = np.indices(data.shape)
    x = (X*data).sum()/total
    y = (Y*data).sum()/total
    col = data[:, int(y)]
    width_x = np.sqrt(np.abs((np.arange(col.size)-y)**2*col).sum()/col.sum())
    row = data[int(x), :]
    width_y = np.sqrt(np.abs((np.arange(row.size)-x)**2*row).sum()/row.sum())
    height = data.max()
    return height, x, y, width_x, width_y

def fitgaussian(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    params = moments(data)
    errorfunction = lambda p: np.ravel(gaussian(*p)(*np.indices(data.shape)) -
                                 data)
    p, success = optimize.leastsq(errorfunction, params)
    return p

def center_by_moments(x_data, y_data, z_data, need_plot=True):
    """Returns xcenter and ycenter of the detector image"""
    data = z_data
    params = fitgaussian(data)
    fit = gaussian(*params)
    (height, x, y, width_x, width_y) = params
    
    if need_plot==True:
        plt.figure()
        fig=plt.figure(figsize=(15, 15), dpi= 200, facecolor='w', edgecolor='k')


        plt.pcolormesh(data.T, norm=LogNorm(), cmap=plt.cm.gist_earth_r)    
        plt.plot(x,y,'go',markersize=30)
        plt.contour(fit(*np.indices(data.shape)), cmap=plt.cm.copper)

        ax = plt.gca()
    
        
        plt.text(0.95, 0.05, """
        x : %.1f
        y : %.1f
        width_x : %.1f
        width_y : %.1f""" %(x, y, width_x, width_y),
                fontsize=16, horizontalalignment='right',
                verticalalignment='bottom', transform=ax.transAxes)
    
    indx, indy = data.shape
    x_index = np.linspace(0, indx-1, indx)
    y_index = np.linspace(0, indy-1, indy)

    fx = interp1d(x_index, x_data)
    fy = interp1d(y_index, y_data)
    
    x0 = fx(x)
    y0 = fy(y)
   
    return x0,y0



def cosineAnglefromMomenta(p1x, p1y, p1z, p2x, p2y, p2z):
    m1_p2D = np.vstack((p1x, p1y, p1z)).T
    m2_p2D = np.vstack((p2x, p2y, p2z)).T
    ang = ((m1_p2D * m2_p2D).sum(axis=1))/(np.linalg.norm(m1_p2D, axis=1) * np.linalg.norm(m2_p2D, axis=1))
    return ang


def singleListTwoConditions(ionTOF, cond1, cond2):

    chk1 = (ionTOF > cond1[0]) & (ionTOF < cond1[1])
    chk2 = (ionTOF > cond2[0]) & (ionTOF < cond2[1])
    
    if chk2.any() and chk1.any():
        cond = chk1 | chk2
    else: cond = chk1 & chk2
    return cond

def oneCondition(ionTOF, cond1):
    cond = (ionTOF > cond1[0]) & (ionTOF < cond1[1])
    return cond


def twoConditions(ionTOF, cond1, cond2):

    indices = ionTOF.size
    checkList = np.zeros(indices, dtype=bool)
    for i in range(0, indices):
        if cond1[0] < ionTOF[i] < cond1[1]:
            checkList[i] = True
            for j in range(i+1, indices):
                if cond2[0] < ionTOF[j] < cond2[1]:
                    checkList[j] = True
                    break
            break
        
    return checkList


def threeConditions(ionTOF, cond1, cond2, cond3):

    indices = ionTOF.size
    checkList = np.zeros(indices, dtype=bool)
    for i in range(0, indices):
        if cond1[0] < ionTOF[i] < cond1[1]:
            checkList[i] = True
            for j in range(i+1, indices):
                if cond2[0] < ionTOF[j] < cond2[1]:
                    checkList[j] = True
                    for k in range(j+1, indices):
                        if cond3[0] < ionTOF[k] < cond3[1]:
                            checkList[k] = True
                            break
                    break
            break
        
    return checkList




#####################################################################
#Jolly's additions: 
#The purpose of this is to return slope and time off_set 
def slope_t_offset(testTOF, mq):

    testTOF = np.array(testTOF)
    mqSqrt = np.sqrt(mq)
    
    k = p[0]
    t_offset = p[1]
    
   # allmq = ((allTOF-t_offset)/k)**2
   
    return k, t_offset

def mq2tof(testTOF, testmq, mq):
    k, t0 = slope_t_offset(testTOF, testmq)
    tof = k*np.sqrt(mq)+t0
    return tof
    
    



def GetZVelocity(particle_i_tof,m, charge,  E_field_new, ion_accel_length, ion_drift_length_new): # to calculate the initial velocity in z direction
    
    index = 0
    
    mass = m * mass_amu
    
    charge = charge * eV
    particle_i_tof = (particle_i_tof)*1e-9
    accel = abs(charge) * abs(E_field_new) / mass
 
    #need to figure out what this if statement is doing!  
    
    if abs(ion_drift_length_new < 1e-9):
        return ion_accel_length/particle_i_tof - accel*particle_i_tof*0.5
  
    # v here is final velocity, which starts with 0 and upon iterations acts as both initial and final velocities 
    v_step = 1e5
    
    #need to figure out what this if statement is doing!
    if mass > 1.5e-27:
        v_step = 1e3
    
    v = 0
    old_tof = CalculateTOFwithInitialV(v, mass, charge, E_field_new, ion_accel_length, ion_drift_length_new)
    #old_Ekin = 0. # what is this for? no other call references in dan.cpp too!
    
    for iteration in range(0, 11):
        v_new = v + v_step
        new_tof = CalculateTOFwithInitialV(v_new, mass, charge, E_field_new, ion_accel_length, ion_drift_length_new)
        step = (particle_i_tof - old_tof)/(new_tof - old_tof)
#        print ("step value before if statement: ", step)
        
        if step < 2 :
            step = step * v_step * 1.
            v_step = v_step *.5 
        else:
            step = step * v_step * 1.
            
        v = v + step
        old_tof = CalculateTOFwithInitialV(v, mass, charge, E_field_new, ion_accel_length, ion_drift_length_new)
        index = index + 1
#        print("Index is: ", index)
#        print("Step value after if statement: ", step)
        if abs(old_tof - particle_i_tof) < 0.005e-12:
            break
        
    return v

def CalculateTOFwithInitialV(v, mass, charge, E_field_new, ion_accel_length, ion_drift_length_new):
    
    accel = abs(charge)*abs(E_field_new)/mass
    #here b acts as final velocity : b^2 - v^2 = 2*a*ion_accel_length
    b = np.sqrt(v*v + 2. * ion_accel_length*accel)
    
    tof1 = (b - v) / accel
    tof2 = ion_drift_length_new / ((accel * tof1) + v)
    
    return tof1 + tof2


def fitEfield():
    
    params = Parameters()
    params.add('E_field', value = 8777.8, vary = True)
    params.add('ion_drift_length', value = 0.015, min=0.013, max=0.017)
    
    params.add('t0', value = 0.)
    
    args = (testTOF, testmq)

    result = minimize(obj_func_res, params, method='least_squares', args = args)   
    
    
    E_field_new = result.params['E_field']
    ion_drift_length_new = result.params['ion_drift_length']
    t0_new = result.params['t0']
    
    return E_field_new, ion_drift_length_new, t0_new



def obj_func_res(params, tof_reference, m_q_reference):# ion_accel_length = 0.27):
    
    paramsvals = params.valuesdict()  #unpack values 
    
    E_field = paramsvals['E_field']
    ion_drift_length = paramsvals['ion_drift_length']
    t0 = paramsvals['t0']
#    print(ion_drift_length)
    
  
    k_temp = 1e9 * np.sqrt(((4. * (ion_accel_length * (ion_accel_length + ion_drift_length) + ion_drift_length**2)/(2.* E_field * ion_accel_length)) * (mass_amu/eV)))
    
    tof_fit = k_temp * np.sqrt(m_q_reference) + t0
    
    residual_tof = (tof_fit - tof_reference)*(tof_fit - tof_reference)
    
    return residual_tof

def GetXVelocity(particle_i_tof,x_position, x_0):
    
    v_x = (x_position - x_0) / (particle_i_tof)
    return (v_x*1e6)


def GetYVelocity(particle_i_tof, y_position, y_0):
    v_y = (y_position - y_0) / (particle_i_tof)
    return (v_y*1e6) 


def auToSI_momentum(p_au):
   p_si = 2 * 1e-24 * p_au
   return (p_si)
    
def SIToau_momentum(p_si):
    p_au = 0.5 * 1e24 * p_si
    return(p_au)
    
def JtoeV(ker_J):
    ker_eV = 6.242*1e18 * ker_J
    return (ker_eV)


