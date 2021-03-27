import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import filters
import tqdm, random
from scipy.signal import argrelextrema, find_peaks
import viewer4C
import matplotlib as mpl
from matplotlib import rc
from matplotlib.backends.backend_pdf import PdfPages
rc('font', size=12)
rc('font', family='Arial')
# rc('font', serif='Times')
rc('pdf', fonttype=42)
# rc('text', usetex=True)


#############################################################################

def plot_4C_overview(files, chromosome, viewpoints, _lims, 
                    close_contact_region=35000,
                    strength_lower_limit=0.75,
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

    N_ax = len(files)

    if backgrounds == None:
        backgrounds = [200 for i in N_ax]
    if peaks_prominence == None:
        peaks_prominence = [200 for i in N_ax]
    if peaks_width == None:
        peaks_width = [20 for i in N_ax]

    # set colors to default values
    if colors is None:
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
    if colors_profile is None:
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

    if N_ax<2:
        return False

    fig, ax = plt.subplots(1,1,figsize = (6,1.5))
    plt.subplots_adjust(right=0.95,top=0.9,bottom=0.25,left=0.05,hspace=0)
    # axs = axs.flatten()

    all_lines = []
    all_gray_scales = []
    _ymax = 0
    for i in range(len(files)):

        file = files[i]
        vp_region = viewpoints[i]
        color = colors[i]
        color_profile = colors_profile[i]

        # extract 4C profile and peak data
        _, data_peaks, data_Zscore, viewpoint = viewer4C.compute_4C_interactions( file, chromosome, vp_region, _lims,
                                    close_contact_region=close_contact_region,
                                    seed=seed,
                                    windowSize=windowSize,
                                    background=backgrounds[i],
                                    peak_prominence=peaks_prominence[i],
                                    peak_width=peaks_width[i],
                                    pvalue=pvalue,
                                    verbose=verbose )


        d_min = int(data_peaks[:,2].max()*strength_lower_limit)
        d_max = int(data_peaks[:,2].max()+10)

        norm = mpl.colors.Normalize(vmin=0, vmax=d_max)
        cmap = mpl.cm.ScalarMappable(norm=norm, cmap=color)
        cmap.set_array([])

        lines = []
        gray_scales = []

        j=0
        for d in data_peaks:
            # print(d,d[2]/float(d_max))
            gray_scales.append(d[2]/float(d_max))
            x, y = viewer4C.parabola(viewpoint, d[0])
            lines.append([x,y])
            if d[2]>d_min:
                print('Int%d'%j,viewpoint,d[0],d[1],d[2],data_Zscore[j])
                _ymax = np.max([_ymax,np.max(y)])
                ax.plot(x,y,'-',lw=1.5,c=color_profile)

            j+=1
        all_lines.append(lines)
        all_gray_scales.append(gray_scales)
        

    norm = mpl.colors.Normalize(vmin=-1, vmax=1.)
    cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.Greys)
    cmap.set_array([])

    centers = [(v[1]+v[0])/2 for v in viewpoints]
    for i, lines in enumerate(all_lines):
        c1 = centers[i]
        gray_scales = all_gray_scales[i]
        for line1, gs1 in zip(lines,gray_scales):
            x1 = line1[0]
            y1 = line1[1]
            dist = [np.abs(x1[0]-c1),np.abs(x1[-1]-c1)]
            end1 = x1[-1]
            if dist[0]>dist[1]:
                end1 = x1[0]
            hit = False
            for j, c2 in enumerate(centers):
                if np.abs(end1-c2)<20000:
                    hit = True
                    jhit = j

            if hit:
                c2 = centers[jhit]
                gray_scales2 = all_gray_scales[jhit]
                hit_back = False
                for k,line2 in enumerate(all_lines[jhit]):
                    x2 = line2[0]
                    y2 = line2[1]
                    dist = [np.abs(x2[0]-c2),np.abs(x2[-1]-c2)]
                    # print(dist)
                    end2 = x2[-1]
                    if dist[0]>dist[1]:
                        end2 = x2[0]
                    if np.abs(end2-c1)<20000:
                        hit_back = True
                        xback = line2[0]
                        yback = line2[1]
                        gs2 = gray_scales2[k]
                
                if hit_back:
                    x,y = viewer4C.parabola(c1,c2)
                    _ymax = np.max([_ymax,np.max(y)])
                    print('Mutual',c1, c2, (gs1+gs2)/2)
                    # print((gs1+gs2)/2)
                    ax.plot(x,y,'-',lw=2,c=cmap.to_rgba((gs1+gs2)/2))

    for i in range(len(files)):
        file = files[i]
        viewpoint_l = viewpoints[i]
        color = colors[i]
        color_profile = colors_profile[i]
        viewpoint = (viewpoint_l[0]+viewpoint_l[1])/2

        ax.plot([viewpoint,viewpoint],[0,1.3*_ymax],'--k')

    ax.set_xlim(_lims)
    ax.set_ylim(0,1.3*_ymax)

    ax.get_yaxis().set_visible(False)
    ax.get_yaxis().set_ticks([])
        
    # add colorbar in the current ax
    if colorbar:
        plt.colorbar(cmap, ax=ax, shrink=0.5)
        
    plt.show()

    return True

######################################################################################################

if __name__=='__main__':
    import os
    
    celltypes = [
            [ os.path.join('..','..','..','..','..','hFOB',i) for i in [
                                        'CPED1_MP_hFOB.txt',
                                        'CPED1_R-11-12_hFOB.txt',
                                        'CPED1_R16_hFOB.txt',
                                        'WNT16_In2_hFOB.txt',
                                        'FAM3C_MP_hFOB.txt',
                                        ] ],
            [ os.path.join('..','..','..','..','..','saos2',i) for i in [
                                        'CPED1_MP_hFOB.txt',
                                        'CPED1_R-11-12_hFOB.txt',
                                        'CPED1_R16_hFOB.txt',
                                        'WNT16_In2_hFOB.txt',
                                        'FAM3C_MP_hFOB.txt',
                                        ] ],
            [ os.path.join('..','..','..','..','..','msc',i) for i in [
                                        # 'CPED1_MP_hFOB.txt',
                                        'CPED1_R-11-12_hFOB.txt',
                                        'CPED1_R16_hFOB.txt',
                                        'WNT16_In2_hFOB.txt',
                                        # 'FAM3C_MP_hFOB.txt',
                                        ] ],
        ]

    chromosome = 'chr7'

    viewpoints_all = [
        [
                [120628556, 120630032],
                [120771311, 120771670],
                [120875513, 120877460],
                [120970581, 120971655],
                [121035619, 121036778],
        ],
        [
                [120628556, 120630032],
                [120771311, 120771670],
                [120875513, 120877460],
                [120970581, 120971655],
                [121035619, 121036778],
        ],
        [
                # [120628556, 120630032],
                [120771311, 120771670],
                [120875513, 120877460],
                [120970581, 120971655],
                # [121035619, 121036778],
        ],
            ]

    colors_all = [
        [
                    mpl.cm.Blues,
                    mpl.cm.Reds,
                    mpl.cm.Greens,
                    mpl.cm.YlOrBr,
                    mpl.cm.RdPu,
        ],
        [
                    mpl.cm.Blues,
                    mpl.cm.Reds,
                    mpl.cm.Greens,
                    mpl.cm.YlOrBr,
                    mpl.cm.RdPu,
        ],
        [
                    # mpl.cm.Blues,
                    mpl.cm.Reds,
                    mpl.cm.Greens,
                    mpl.cm.YlOrBr,
                    # mpl.cm.RdPu,
        ],
            ]
    colors_profile_all = [
                ['blue','red','green','gold','purple'],
                ['blue','red','green','gold','purple'],
                [       'red','green','gold'         ]
                    ]

    peaks_prominence_all = [
                        [ 202, 280, 81, 95, 267 ],
                        [ 202, 280, 81, 95, 267 ],
                        [      280, 81, 95      ]
                    ]                            

    backgrounds_all = [
                [ 28, 228, 79, 15, 336 ],
                [ 28, 228, 79, 15, 336 ],
                [     228, 79, 15      ],
                ]

    peaks_width_all = [ [16]*5, [16]*5, [16]*3 ]
                
    _lims = (120400000,121750000)
    
    close_contact_region = 35000

    for m in range(len(celltypes)):
        # ax = axs[m]
        random.seed(1)
        
        files = celltypes[m]
        viewpoints = viewpoints_all[m]
        colors = colors_all[m]
        colors_profile = colors_profile_all[m]
        peaks_prominence = peaks_prominence_all[m]
        backgrounds = backgrounds_all[m]
        peaks_width = peaks_width_all[m]

        plot_4C_overview(files, chromosome, viewpoints, _lims, 
                        close_contact_region=close_contact_region,
                        strength_lower_limit=0.75,
                        colors=colors, 
                        colors_profile=colors_profile,
                        colorbar=False, 
                        seed=1,
                        windowSize=30,
                        backgrounds=backgrounds,
                        peaks_prominence=peaks_prominence,
                        peaks_width=peaks_width,
                        pvalue=0.00001,
                        verbose=False)
