import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import filters
import tqdm, random
from collections.abc import Iterable
from scipy.signal import argrelextrema, find_peaks, medfilt
from scipy.ndimage.filters import generic_filter
import matplotlib as mpl
from matplotlib import rc
from matplotlib.backends.backend_pdf import PdfPages
rc('font', size=10)
rc('font', family='Arial')
rc('pdf', fonttype=42)

##########################################################################################

def parabola(vi, di):

    x = np.linspace(vi,di,100)

    p1 = [di,0]
    p2 = [vi,0]
    p3 = [(vi+di)/2,np.abs(vi-di)/200]

    f1 = p1[1]*(x-p2[0])*(x-p3[0])/((p1[0]-p2[0])*(p1[0]-p3[0]))
    f2 = p2[1]*(x-p1[0])*(x-p3[0])/((p2[0]-p1[0])*(p2[0]-p3[0]))
    f3 = p3[1]*(x-p1[0])*(x-p2[0])/((p3[0]-p1[0])*(p3[0]-p2[0]))

    y = f1+f2+f3

    return x, y

def estimate_parameters(file, chromosome):
    
    # load data from file    
    data = []
    with open(file,'r') as f:
        data = f.readlines()
    data = [i.split('\t') for i in data]

    # filter the rows containing chromosome information
    data_other = []
    data_chr = []
    for i in data:
        if type(i[0]) is str:
            if chromosome != i[0].lower():
                data_other.append(float(i[3]))
            else:
                data_chr.append(float(i[3]))

    # estimate background
    data_other = np.array(data_other)
    perc = np.percentile(data_other, 99.9)
    background = np.mean(data_other[data_other>perc])

    # estimate prominence
    data_chr = np.array(data_chr)
    data_chr = np.clip(data_chr-background, 0, None)
    prominence = int(np.max(data_chr)/20)

    return background, prominence

def compute_4C_interactions(file, chromosome, vp_region, _lims,
                            close_contact_region=35000,
                            seed=1,
                            windowSize=31,
                            background=200,
                            peak_prominence=200,
                            peak_width=20,
                            pvalue=0.05,
                            verbose=False):
    '''
    file: text file in the form "chrNN   start   end   interaction_value" without headers
    vp_region: region of the viewpoint used in the 4C experiment
    ax: 
    '''

    # set random seed
    random.seed(seed)

    # load data from file    
    data = []
    with open(file,'r') as f:
        data = f.readlines()
    data = [i.split('\t') for i in data]

    # filter the rows containing chromosome information
    data_o = []
    for i in data:
        if type(i[0]) is str:
            if chromosome == i[0].lower():
                data_o.append([float(i[1]),float(i[2]),float(i[3])])
    data_o = np.array(data_o)

    # define the viewpoint spot center
    viewpoint = (vp_region[0]+vp_region[1])/2

    # print(data_o[:5])
    if len(data_o)==0:
        return False, False, False

    ### find the threshold with fitting a poisson
    from scipy.optimize import curve_fit
    from scipy.stats import norm

    # split data to left and right of the viewpoint
    data_right = data_o[data_o[:,0]>viewpoint]
    data_left = data_o[data_o[:,0]<viewpoint][::-1]

    if len(data_left[:,2]) < len(data_right[:,2]):
        data_right = data_right[:len(data_left[:,2]),:]
    else:
        data_left = data_left[:len(data_right[:,2]),:]
    _mean = (data_right[:,2]+data_left[:,2])/2
    data_left[:,2] = _mean
    data_right[:,2] = _mean

    # use running window
    data_right[:,2] = filters.gaussian_filter1d(data_right[:,2],sigma=1*windowSize)
    data_left[:,2] = filters.gaussian_filter1d(data_left[:,2],sigma=1*windowSize)

    # limit fitting data to region close to viewpoint
    data_right = data_right[(data_right[:,1]-viewpoint)<2*(_lims[1]-_lims[0])]
    data_left = data_left[(viewpoint-data_left[:,0])<2*(_lims[1]-_lims[0])]

    def exponential(k, N0, lamb):
        return N0*np.exp(-lamb*k)
    
    # compute right exponential and its standard deviation
    params_right, pcov = curve_fit(exponential,np.arange(len(data_right)), data_right[:,2])
    sigma_right = np.sqrt(np.diag(pcov))
    upper_right = np.array([
        exponential(np.arange(len(data_right)), params_right[0] + norm.ppf(1.-pvalue)*sigma_right[0], params_right[1] + norm.ppf(1.-pvalue)*sigma_right[1]),
        exponential(np.arange(len(data_right)), params_right[0] + norm.ppf(1.-pvalue)*sigma_right[0], params_right[1] - norm.ppf(1.-pvalue)*sigma_right[1]),
        exponential(np.arange(len(data_right)), params_right[0] - norm.ppf(1.-pvalue)*sigma_right[0], params_right[1] + norm.ppf(1.-pvalue)*sigma_right[1]),
        exponential(np.arange(len(data_right)), params_right[0] - norm.ppf(1.-pvalue)*sigma_right[0], params_right[1] - norm.ppf(1.-pvalue)*sigma_right[1]),
        ])
    stdr = np.std(upper_right, axis=0)
    yr = exponential(np.arange(len(data_right)), *params_right)

    # compute left exponential and its standard deviation
    params_left, pcov = curve_fit(exponential,np.arange(len(data_left)), data_left[:,2])
    sigma_left = np.sqrt(np.diag(pcov))
    upper_left = np.array([
        exponential(np.arange(len(data_left)), params_left[0] + norm.ppf(1.-pvalue)*sigma_left[0], params_left[1] + norm.ppf(1.-pvalue)*sigma_left[1]),
        exponential(np.arange(len(data_left)), params_left[0] + norm.ppf(1.-pvalue)*sigma_left[0], params_left[1] - norm.ppf(1.-pvalue)*sigma_left[1]),
        exponential(np.arange(len(data_left)), params_left[0] - norm.ppf(1.-pvalue)*sigma_left[0], params_left[1] + norm.ppf(1.-pvalue)*sigma_left[1]),
        exponential(np.arange(len(data_left)), params_left[0] - norm.ppf(1.-pvalue)*sigma_left[0], params_left[1] - norm.ppf(1.-pvalue)*sigma_left[1]),
        ])
    stdl = np.std(upper_left, axis=0)
    yl = exponential(np.arange(len(data_left)), *params_left)

    # combine the left and right in a single line
    exp = np.concatenate([yl[::-1],yr])
    xexp = np.concatenate([data_left[:,0][::-1],data_right[:,0]])
    expstd = np.concatenate([stdl[::-1],stdr])

    ### find maxima and minima
    interactions_gf = generic_filter(data_o[:,2], np.mean, size=windowSize)#filters.gaussian_filter1d(data_o[:,2],windowSize)
    data_gf = np.array([i for i in data_o])
    data_gf[:,2] = interactions_gf

    peaks, _ = find_peaks(np.clip(data_gf[:,2]-background,0,None), height=[None,None], prominence=[peak_prominence,None], width=peak_width)
    peaks = [p for p in peaks if np.abs(data_gf[p,0]-viewpoint)>close_contact_region]
    # print(peaks)

    data_peaks = []
    data_Zscore = []
    data_peaks_low = []
    for data in data_gf[peaks,:]:
        thr = (exp+norm.ppf(1.-pvalue)*expstd)[xexp==data[0]]
        if data[2] > thr:
            data_peaks.append(data)
            expidx = np.where(xexp==data[0])
            data_Zscore.append((data[2]-exp[expidx])/expstd[expidx])
        else:
            data_peaks_low.append(data)
    data_peaks = np.array(data_peaks)
    data_peaks_low = np.array(data_peaks_low)

    # limit data to the plot limits
    data_o = data_o[data_o[:,0]>_lims[0]]
    data_o = data_o[data_o[:,1]<_lims[1]]
    data_gf = data_gf[data_gf[:,0]>_lims[0]]
    data_gf = data_gf[data_gf[:,1]<_lims[1]]
    xexp_plot = xexp[(xexp>_lims[0])&(xexp<_lims[1])]
    exp_plot = exp[(xexp>_lims[0])&(xexp<_lims[1])]
    expstd_plot = expstd[(xexp>_lims[0])&(xexp<_lims[1])]

    if verbose:
        # plot gaussian filtered interactions
        fig = plt.figure()
        # plt.plot(data_right[:,0],data_right[:,2],'-g')
        # plt.plot(data_left[:,0],data_left[:,2],'-g')
        plt.plot(data_o[:,0],data_o[:,2], '-b', alpha=0.5)
        plt.plot(data_gf[:,0],data_gf[:,2], '-y', alpha=0.75)
        plt.plot(xexp_plot,exp_plot,'-r')
        plt.axvspan(viewpoint-close_contact_region, viewpoint+close_contact_region, alpha=0.2, color='black', lw=0)
        plt.fill_between(xexp_plot,exp_plot+norm.ppf(1.-pvalue)*expstd_plot,color='r',alpha=.2)
        if len(data_peaks)>0:
            plt.plot(data_peaks[:,0],data_peaks[:,2],'*r',ms=10)
        if len(data_peaks_low)>0:
            plt.plot(data_peaks_low[:,0],data_peaks_low[:,2],'*k',ms=10,alpha=.5)
        plt.xlim(_lims[0],_lims[1])
        plt.legend(['Interaction profile','Smoothened profile','Exponential fit'], loc=1)
        fig.show()

    return data_o, data_peaks, data_Zscore, viewpoint

def plot_4C_interactions(files, chromosome, viewpoints, _lims, 
                        close_contact_region=35000,
                        colors=None, 
                        colors_profile=None,
                        colorbar=False, 
                        seed=1,
                        windowSize=31,
                        backgrounds=None,
                        peaks_prominence=None,
                        peaks_width=None,
                        pvalue=0.00001,
                        verbose=False):

    '''
    '''
    N_ax = len(files)

    if backgrounds == None:
        backgrounds = [200 for i in N_ax]
    if peaks_prominence == None:
        peaks_prominence = [200 for i in N_ax]
    if peaks_width == None:
        peaks_width = [20 for i in N_ax]
    
    # set colors to default values
    if colors == None:
        cs = [
                mpl.cm.Blues,
                mpl.cm.Reds,
                mpl.cm.Greens,
                mpl.cm.YlOrBr,
                mpl.cm.RdPu,
            ]
        colors = []
        while len(colors)<=N_ax:
            for c in cs:
                colors.append(c)
    if colors_profile == None:
        cs = [
                'blue',
                'red',
                'green',
                'gold',
                'purple',
            ]
        colors_profile = []
        while len(colors_profile)<=N_ax:
            for c in cs:
                colors_profile.append(c)

    # setup figures
    fig, ax = plt.subplots(N_ax,1,figsize = (6,2*N_ax))
    if not isinstance(ax, Iterable):
        ax = [ax]
    plt.subplots_adjust(right=0.95,top=0.9,bottom=0.125,left=0.05,hspace=0)

    for i in tqdm.tqdm(range(N_ax)):
#        print(files[i])
        
        file = files[i]
        vp_region = viewpoints[i]
        color = colors[i]
        color_profile = colors_profile[i]
        
        # extract 4C profile and peak data
        data_o, data_peaks, data_Zscore, viewpoint = compute_4C_interactions(file, chromosome, vp_region, _lims,
                                    close_contact_region=close_contact_region,
                                    seed=seed,
                                    windowSize=windowSize,
                                    background=backgrounds[i],
                                    peak_prominence=peaks_prominence[i],
                                    peak_width=peaks_width[i],
                                    pvalue=pvalue,
                                    verbose=verbose)
        
        print(data_peaks[:,2])
        print(data_Zscore)
        if viewpoint==False:
            return False

        # prepare colormap for colorbar
        if len(data_peaks):
            norm = mpl.colors.Normalize(vmin=0, vmax=int(data_peaks[:,2].max()))
            cmap = mpl.cm.ScalarMappable(norm=norm, cmap=color)
            cmap.set_array([])

        # compute and plot parabola
        _max = np.max(data_o[:,2])
        if len(data_peaks):
            for d in data_peaks:
                x, y = parabola(viewpoint, (d[0]+d[1])/2)
                ax[i].plot(x,y,'-',lw=1.5,c=cmap.to_rgba(int(d[2])))
                _max = np.max([_max, np.max(y)])

        # plot 4C interaction profiles
        ax[i].fill_between(data_o[:,0],data_o[:,2],color=color_profile,alpha=0.5,lw=0)

        # set ax lims
        ax[i].set_xlim(_lims)
        ax[i].set_ylim(0, 1.5*_max)
        ax[i].ticklabel_format(axis="x", style="sci", scilimits=(0,2))
        ax[i].ticklabel_format(axis="y", style="sci", scilimits=(0,2))

        # add text on the current ax
        fname, _ = os.path.splitext(os.path.split(file)[-1])
        ax[i].text(_lims[0]+(_lims[1]-_lims[0])/50,1.5*_max*0.85,fname)

        # add colorbar in the current ax
        if colorbar:
            plt.colorbar(cmap, ax=ax[i], shrink=0.5)

        # mark all the viewpoints in the current ax
        for j in range(len(files)):
            v = (viewpoints[j][1]+viewpoints[j][0])/2
            ax[i].plot([v,v],[0,1.5*_max],'--k',lw=.5,dashes = (4,4))
            
    # remove all ax labels
    for i in range(len(files)-1):
        ax[i].get_xaxis().set_visible(False)
        ax[i].get_xaxis().set_ticks([])
        
    plt.show()

    return True

##########################################################################################

if __name__=='__main__':
    
    files = [ 
            os.path.join('..','..','..','..','..','hFOB',i) for i in [
            # os.path.join('..','..','..','..','..','hFOB',i) for i in [
            # os.path.join('..','..','..','..','..','hFOB',i) for i in [

                                        'CPED1_MP_hFOB.txt',
                                        'CPED1_R-11-12_hFOB.txt',
                                        'CPED1_R16_hFOB.txt',
                                        'WNT16_In2_hFOB.txt',
                                        'FAM3C_MP_hFOB.txt',
                                        ]
            ]

    chromosome = 'chr7'
    
    viewpoints = [
                    [120628556, 120630032],
                    [120771311, 120771670],
                    [120875513, 120877460],
                    [120970581, 120971655],
                    [121035619, 121036778],
                 ]

    peak_prominence = [
                        202,
                        280,
                        81,
                        95,
                        267
                        ]
    background = [
                    28,
                    228,
                    79,
                    15,
                    336
                    ]
    peak_width = [16 for i in files]
                
    _lims = (120400000,121750000)
    
    close_contact_region = 35000
    
    plot_4C_interactions(
                        files, chromosome, viewpoints, _lims, 
                        close_contact_region=close_contact_region,
                        colors=None, 
                        colors_profile=None,
                        colorbar=False, 
                        seed=1,
                        windowSize=30,
                        backgrounds=background,
                        peaks_prominence=peak_prominence,
                        peaks_width=peak_width,
                        pvalue=0.0005,
                        verbose=True
                        )

##########################################################################################
