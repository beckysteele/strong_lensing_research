
from pylab import *     # matplotlib
import math
import scipy
import numpy
from scipy import stats
from scipy.optimize import leastsq      # Levenberg-Marquandt algorithm
from scipy.stats.stats import spearmanr # Spearman correlation coefficient
from glob import glob

import ephem

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.cbook as cbook
import pyfits as py

#### readfile reads .txt files and returns an array of floats
def readfile (filename,splitchar):            # filename is the name (and location) of the file we wish to read
    seq = []                        # creates an empty string
    in_file = open(filename, "r")   # creates file object in_file, read-only
    lines = in_file.readlines()     # creates an array of the lines of the file

    for l in range(len(lines)):     # scans over all lines
        temp = lines[l].split()     # splits string into component words at white spaces
        if temp[0] != '#':
            print temp[0], len(temp)
            for i in range(len(temp)):  # scans over all elements in current line
                if temp[i] == 'nan':
                    temp[i] = -99.
                    # print 'found Naaaaan!'
                else:
                    try:                    
                        temp[i] = float(temp[i])    # try to make a float                    
                    except Exception:
                        #print temp[i],
                        #print "could not be converted to a float."
                        temp[i] = temp[i]       

            if len(temp) > 100:
                zrs = numpy.zeros(93)
                temp = numpy.concatenate(([float(temp[0]) ], zrs))
            seq.append(temp)
            print temp[0],
            print '\t',
            print len(temp)
    in_file.close()                 # closes the opened file
    output = seq
    return output                   # returns the output

def mkplots():


    # read in Ivan's merged catalog
    # /Users/bholwerd/Dropbox/BlendedSpectra

    can = py.open('/Users/bholwerd/Desktop/Science/Occulting-Galaxies/GAMA/BlendedSpectra/DoublezCandidates.fits')
    
    gama = can[1].data
    cols = can[1].columns
    colnames = can[1].columns.names
    print colnames

    gama_id        = gama['CATAID']
    ra             = gama['RA']
    dec            = gama['DEC']

    z1             = gama['Z']
    template1      = gama['TEMPLATE']
    cc_sigma1      = gama['CC_SIGMA']
    z2             = gama['Z2']
    template2      = gama['TEMPLATE2']
    cc_sigma2      = gama['CC_SIGMA2']
    specID         = gama['SPECID']
    is_dbest       = gama['IS_DBEST']
    vis_class      = gama['VISUALTYPE']
    spec_class     = gama['SPECTYPE']
    vis_comments   = gama['COMMENT']

    print is_dbest
    ok = np.where(is_dbest=='T')
    print ok
    gama_id        = gama_id[ok]
    ra             = ra[ok]
    dec            = dec[ok]
    z1             = z1[ok]
    template1      = template1[ok]
    cc_sigma1      = cc_sigma1[ok]
    z2             = z2[ok]
    template2      = template2[ok]
    cc_sigma2      = cc_sigma2[ok]
    specID         = specID[ok]
    vis_class      = vis_class[ok]
    spec_class     = spec_class[ok]
    vis_comments   = vis_comments[ok]

    print spec_class
       
    # numerelize the visual classifications
    # 1 - SE
    # 1.1 - B
    # 2 - SS
    # 2.1 - F
    # 2.2 - Phi
    # 2.3 - Q
    # 3 - ES
    # 4 - EE
    # 5 - Ef
    # 6 - M
    # 7 - E
    # 8 - S

    vc = np.zeros(len(vis_class),dtype='float')

    vc[np.where(vis_class=='SE')] = 1
    vc[np.where(vis_class=='SE/B')] = 1.1
    vc[np.where(vis_class=='SS')] = 2
    vc[np.where((vis_class=='SS/F')|(vis_class=='F'))] = 2.1
    vc[np.where((vis_class=='SS/Phi')|(vis_class=='Phi'))] = 2.2
    vc[np.where((vis_class=='SS/Q')|(vis_class=='Q'))] = 2.3

    vc[np.where(vis_class=='ES')] = 3
    vc[np.where(vis_class=='EE')] = 4
    vc[np.where((vis_class=='Ef')|(vis_class=='ES/Ef')|(vis_class=='EE/Ef'))] = 5
    vc[np.where((vis_class=='M')|(vis_class=='SE/M')|(vis_class=='SS/M')|(vis_class=='ES/M')|(vis_class=='S/M'))] = 6
    vc[np.where(vis_class=='E')] = 7
    vc[np.where(vis_class=='S')] = 8

    for i in range(len(vc)):
        print vc[i],
        print vis_class[i],
        print spec_class[i],
        print vis_comments[i]
    # exit()
    # collating ALL the different catalogs into a single one

    lens  = np.where(spec_class=='PG+ELG')
    ocg   = np.where(spec_class=='ELG+PG')
    uvocg   = np.where(spec_class=='ELG+ELG')
    EEocg   = np.where(spec_class=='PG+PG')


    figdir = '/Users/bholwerd/Desktop/Science/Occulting-Galaxies/GAMA/Catalog/Figs/'

    # read in the stellar mass catalog 
    can = py.open('/Users/bholwerd/Desktop/Science/GAMA_data/StellarMasses.fits')
    can.verify('fix')
    # stmass = can[0].data
    stmass = can[1].data
    cols = can[1].columns
    print cols
    colnames = can[1].columns.names
    print colnames

    # stellar mass

    stmass_id = stmass['CATAID']
    lgmass = stmass['logmstar']

    bs_lgmass = np.zeros(len(gama_id), dtype='float')
    # crosscorrelate based on GAMA id
    for i in range(len(gama_id)):
        for j in range(len(stmass_id)):
            if gama_id[i] == stmass_id[j]:
                bs_lgmass[i] = lgmass[j]

    # read in the SLACS information
    slacst3 = readfile('slacs9_t3.txt', ' ')
    t3 = np.asarray(slacst3)
    t3id = np.asarray(t3[:,0], dtype='string')
    t3z = np.asarray(t3[:,9], dtype='float')
    
    print 'Table 4'
    slacst4 = readfile('slacs9_t4.txt', ' ')
    t4 = np.asarray(slacst4)
    t4id = np.asarray(t4[:,0], dtype='string')
    t4m1 = np.asarray(t4[:,8], dtype='float') # Salpeter IMF
    t4m2 = np.asarray(t4[:,10], dtype='float') # Chabrier IMF

    slacs_z = np.zeros(len(t3id))
    slacs_m1 = np.zeros(len(t3id))
    slacs_m2 = np.zeros(len(t3id))

    for i in range(len(t3id)):
        for j in range(len(t4id)):
            if t3id[i] == t4id[j]:
                slacs_z[i] = t3z[i]
                slacs_m1[i] = t4m1[j]
                slacs_m2[i] = t4m2[j]

    slacst4 = readfile('S4TM_SMass.txt', ' ')
    t4 = np.asarray(slacst4)
    s4m_id = np.asarray(t4[:,0], dtype='string')
    s4m_z = np.asarray(t4[:,1], dtype='float') 
    s4m_m1 = np.asarray(t4[:,2], dtype='float') # Salpeter IMF
    s4m_m2 = np.asarray(t4[:,3], dtype='float') # Chabrier IMF



    # Redshift -stellar mass plot.

    fig0 = plt.figure()
    ax = fig0.add_subplot(111)
    ax.plot(slacs_m1, slacs_z, 'ko', alpha=0.5, markersize=15, label='SLACS')
    ax.plot(s4m_m1, s4m_z, 'g^', alpha=0.5, markersize=15, label='SLACS4MASSES')
    ax.plot(bs_lgmass, z1, 'r*', alpha=0.5, markersize=15, label='GAMA')

    ax.axis([8.5,13,0,0.6], prop={'size':25})
    ax.set_xlabel(r'Stellar Mass $log_{10}(M^*)$', fontsize=25)
    ax.set_ylabel(r'Redshift (z)', fontsize=25)
    ax.legend(loc='upper left', numpoints=1)
    savefig('./LENSES_mass_z.png')
    savefig('./LENSES_mass_z.pdf')
    savefig(figdir+'/LENSES_mass_z.png')
    savefig(figdir+'/LENSES_mass_z.pdf')


    plt.close("all")


    # Three-panel plot
    #     import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter
    nullfmt   = NullFormatter()         # no labels
    # definitions for the axes
    left, width = 0.1, 0.65
    lower, height = 0.1, 0.65
    lower_h = left_h = left+width+0.02
    
    rect_scatter = [left, lower, width, height]
    rect_histx = [left, lower_h, width, 0.2]
    rect_histy = [left_h, lower, 0.2, height]
    
    # start with a rectangular Figure
    plt.figure(1, figsize=(8,8))
    
    axScatter = plt.axes(rect_scatter)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)
    
    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)
    
    axHistx.yaxis.set_major_formatter(nullfmt)
    axHisty.xaxis.set_major_formatter(nullfmt)

    # axScatter.plot(mean_reff, H_36, 'ro', alpha=1, markersize=15, label='z=9-10 candidates')
    axScatter.plot(slacs_m1, slacs_z, 'ko', alpha=0.5, markersize=15, label='SLACS4MASSES')
    axScatter.plot(s4m_m1, s4m_z, 'g^', alpha=0.5, markersize=15, label='SLACS')
    axScatter.plot(bs_lgmass, z1, 'r*', alpha=1, markersize=15, label='GAMA')
    axScatter.legend(loc='upper left', numpoints=1)

    x = np.concatenate( (slacs_m1,s4m_m1,bs_lgmass), axis=0)
    y = np.concatenate( ( slacs_z,s4m_z,z1), axis=0)
    
    # now determine nice limits by hand:
    binwidth = 0.05
    xbinwidth = 0.3
    ybinwidth = 0.05
    xymax = np.max( [np.max(np.fabs(x)), np.max(np.fabs(y))] )
    lim = ( int(xymax/binwidth) + 1) * binwidth
    lim = 1.0
    
    axScatter.axis([8.5,12.5,0.,0.6], prop={'size':25})

    axScatter.set_xlabel(r'Stellar Mass $log_{10}(M^*)$', fontsize=20)
    axScatter.set_ylabel(r'Redshift (z)', fontsize=20)
    
    xbins = np.arange(8.5, 12.5, xbinwidth)
    ybins = np.arange(0., 0.6, ybinwidth)

    axHistx.hist(slacs_m1, bins=xbins, edgecolor='k', histtype='step', linewidth=2)
    axHistx.hist(s4m_m1, bins=xbins, edgecolor='g', histtype='step', linewidth=2, linestyle='dotted')
    axHistx.hist(bs_lgmass, bins=xbins, edgecolor='r', histtype='step', linewidth=4, linestyle='dashed')
    
    axHisty.hist(slacs_z, bins=ybins, orientation='horizontal',  edgecolor='k', histtype='step', linewidth=2)
    axHisty.hist(s4m_z, bins=ybins, orientation='horizontal',  edgecolor='g', histtype='step', linewidth=2, linestyle='dotted')
    axHisty.hist(z1, bins=ybins, orientation='horizontal',  edgecolor='r', histtype='step', linewidth=4, linestyle='dashed')

    
    axHistx.set_xlim( axScatter.get_xlim() )
    axHisty.set_ylim( axScatter.get_ylim() )
    axHisty.set_xlabel(r'#', fontsize=20)
    axHistx.set_ylabel(r'#', fontsize=20)
    # Referee asks for tick marks at unity along y-axes
    axHistx.set_yticks([0,10,20,30,40,50])
    axHisty.set_xticks([0,10,20,30,40,50])

    savefig('./LENSES_mass_z_threeplot.png')
    savefig('./LENSES_mass_z_threeplot.pdf')

    savefig('/Users/bholwerd/Desktop/Proposals/HST/C24/SNAP-GAMA-LENSES/Figures/LENSES_mass_z_threeplot.pdf')


    fig0 = plt.figure()
    ax = fig0.add_subplot(111)
    ax.hist(slacs_m1, 40, range=[8,12],normed=0, edgecolor='k', histtype='step', label='SLACS Lenses')
    ax.hist(bs_lgmass, 40, range=[8,12],normed=0, edgecolor='r', histtype='step', label='GAMA Lenses')

    ax.axis([8,12,0,20], prop={'size':25})
    ax.set_xlabel(r'Stellar Mass $log_{10}(M^*)$', fontsize=25)
    ax.set_ylabel(r'Number of Objects', fontsize=25)
    ax.legend(loc='upper left', numpoints=1)
    savefig('./LENSES_hist_mass.png')
    savefig('./LENSES_hist_mass.pdf')
    savefig(figdir+'/LENSES_hist_mass.png')
    savefig(figdir+'/LENSES_hist_mass.pdf')
    exit()


    # histogram of types (FG & BG in each pair).
    fig0 = plt.figure()
    ax = fig0.add_subplot(111)
    ax.hist(z1, 32, range=[0,0.8],normed=0, edgecolor='g', histtype='step', label='Foreground Galaxy')
    ax.hist(z2, 32, range=[0,0.8],normed=0, edgecolor='b', histtype='step', label='Background Galaxy')
    ax.set_xlabel(r'Redshift ($z$)', fontsize=25)
    ax.set_ylabel(r'Number of Galaxies', fontsize=25)
    ax.legend(loc='upper right')
    savefig(figdir+'/z-hist.png')
    savefig(figdir+'/z-hist.pdf')

    # B&W
    fig0 = plt.figure()
    ax = fig0.add_subplot(111)
    ax.hist(z1, 32, range=[0,0.8],normed=0, edgecolor='k', histtype='step', label='Foreground Galaxy')
    ax.hist(z2, 32, range=[0,0.8],normed=0, edgecolor='k', alpha=0.5, linewidth=3,histtype='step', label='Background Galaxy')
    ax.set_xlabel(r'Redshift ($z$)', fontsize=25)
    ax.set_ylabel(r'Number of Galaxies', fontsize=25)
    ax.legend(loc='upper left')
    savefig(figdir+'/z-hist_bw.png')
    savefig(figdir+'/z-hist_bw.pdf')



    # best template and visual classification
    fig0 = plt.figure()
    ax = fig0.add_subplot(111)
    ax.scatter(template1, template2, s=35, c=vc)

    ax.axis([0,10,0,50], prop={'size':25})

    ax.set_xlabel(r'Foreground Template', fontsize=25)
    ax.set_ylabel(r'Background Template', fontsize=25)
    savefig(figdir+'/template-vis.png')
    savefig(figdir+'/template-vis.pdf')


    # B&W



    # best template hist
    fig0 = plt.figure()
    ax = fig0.add_subplot(111)
    ax.hist(template1, 10, range=[40,50], normed=0, edgecolor='g',alpha=1, histtype='step', label='Foreground Galaxy')
    ax.hist(template2, 10, range=[40,50], normed=0, edgecolor='k', histtype='step', label='Background Galaxy')
    ax.axis([39,50,0,150], prop={'size':25})

    ax.set_xlabel(r'Template', fontsize=25)
    ax.set_ylabel(r'Number of Galaxies', fontsize=25)
    ax.legend(loc='upper right')
    savefig(figdir+'/templates-hist.png')
    savefig(figdir+'/templates-hist.pdf')


    fig0 = plt.figure()
    ax = fig0.add_subplot(111)
    ax.hist(template1, 10, range=[40,50], normed=0, edgecolor='k', alpha=1, histtype='step', label='Foreground Galaxy')
    ax.hist(template2, 10, range=[40,50], normed=0, edgecolor='k', alpha=0.5, linewidth=3, histtype='step', label='Background Galaxy')
    ax.axis([39,50,0,150], prop={'size':25})

    ax.set_xlabel(r'Template', fontsize=25)
    ax.set_ylabel(r'Number of Galaxies', fontsize=25)
    ax.legend(loc='upper right')
    savefig(figdir+'/templates-hist_bw.png')
    savefig(figdir+'/templates-hist_bw.pdf')



    print len(template1[np.where(template1==44.)])
    # exit()
    lensA  = np.where((template1<=44)&(template1>30)&(template2>44))
    ocgA   = np.where((template2<=44)&(template1>30)&(template1>44))
    uvocgA = np.where((template2<=44)&(template1>30)&(template1<=44)&(template1>30))
    EEocgA = np.where((template2>44)&(template1>44))
    


    lensB = np.where( ((template1==23)|(template1==24)|(template1==28))&(template2>24)&(template2<28))
    ocgB = np.where( ((template2==23)|(template2==24)|(template2==28))&(template1>24)&(template1<28))
    uvocgB = np.where(  (template1>24)&(template1<28) & (template2>24)&(template2<28))
    EEocgB = np.where( ((template1==23)|(template1==24)|(template1==28))&((template2==23)|(template2==24)|(template2==28)) )

    lens  = tuple(map(array, np.sort(np.concatenate((lensA,lensB), axis=1))))
    ocg   = tuple(map(array, np.sort(np.concatenate((ocgA,ocgB), axis=1))))
    uvocg = tuple(map(array, np.sort(np.concatenate((uvocgA,uvocgB), axis=1))))
    EEocg = tuple(map(array, np.sort(np.concatenate((EEocgA,EEocgB), axis=1))))




    fig0 = plt.figure()
    ax = fig0.add_subplot(111)
    ax.plot(z1[lens],  z2[lens],  'ko', alpha=0.2,markersize=10)
    ax.plot(z1[ocg],   z2[ocg],   'go', alpha=0.2,markersize=10)
    ax.plot(z1[uvocg], z2[uvocg], 'bo', alpha=0.2,markersize=10)
    ax.plot(z1[EEocg], z2[EEocg], 'ro', alpha=0.2,markersize=10)
    ax.axis([0,1,0,1])
    ax.set_xlabel(r'Foreground Galaxy Redshift ($z_{f}$)', fontsize=25)
    ax.set_ylabel(r'Background Galaxy Redshift ($z_{b}$)', fontsize=25)
    ax.legend(loc='upper right')
    savefig(figdir+'/z-z.png')
    savefig(figdir+'/z-z.pdf')


    print 'Number of Pair Types:'
    print 'ET+LT ',
    print len(z1[lens])
    print 'LT+ET ',
    print len(z1[ocg])
    print 'LT+LT ',
    print len(z1[uvocg])
    print 'ET+ET ',
    print len(z1[EEocg])
    print 'Sum: ',
    print len(z1[lens]) + len(z1[ocg]) + len(z1[uvocg]) + len(z1[EEocg])
    print 'Total: ',
    print len(z1)

    fig0 = plt.figure()
    ax = fig0.add_subplot(111)
    ax.hist(vc, 10, range=[0,10],normed=0, edgecolor='k',alpha=1, histtype='step', label='All')
    ax.hist(vc[lens], 10, range=[0,10],normed=0, edgecolor='m',alpha=1, histtype='step', label='E+S templates')
    ax.hist(vc[ocg], 10, range=[0,10],normed=0, edgecolor='g',alpha=1, histtype='step', label='S+E templates')
    ax.hist(vc[uvocg], 10, range=[0,10],normed=0, edgecolor='b',alpha=1, histtype='step', label='S+S templates')
    ax.hist(vc[EEocg], 10, range=[0,10],normed=0, edgecolor='r',alpha=1, histtype='step', label='E+E templates')
    ax.axis([0,9,0,100], prop={'size':25})
    ax.set_xticks([0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5])
    ax.set_xticklabels(['SE','SS','ES','EE','Ef','M','E','S'])
    ax.set_xlabel(r'Visual Classification', fontsize=25)
    ax.set_ylabel(r'Number of Pairs', fontsize=25)
    ax.legend(loc='upper right')
    savefig(figdir+'/vis-hist.png')
    savefig(figdir+'/vis-hist.pdf')

    fig0 = plt.figure()
    ax = fig0.add_subplot(111)
    h,b,p = ax.hist(vc, 10, range=[0,10],normed=0, edgecolor='k',alpha=1, histtype='step', label='All')
    h_lens,b_lens,p = ax.hist(vc[lens], 10, range=[0,10],normed=0, edgecolor='m',alpha=1, histtype='step', label='E+S templates')
    h_ocg,b_ocg,p = ax.hist(vc[ocg], 10, range=[0,10],normed=0, edgecolor='g',alpha=1, histtype='step', label='S+E templates')
    h_uvocg,b_uvocg,p = ax.hist(vc[uvocg], 10, range=[0,10],normed=0, edgecolor='b',alpha=1, histtype='step', label='S+S templates')
    h_ee,b_ee,p = ax.hist(vc[EEocg], 10, range=[0,10],normed=0, edgecolor='r',alpha=1, histtype='step', label='E+E templates')
    ax.axis([0,9,0,100], prop={'size':25})
    ax.set_xticks([0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5])
    ax.set_xticklabels(['SE','SS','ES','EE','Ef','M','E','S'])
    ax.set_xlabel(r'Visual Classification', fontsize=25)
    ax.set_ylabel(r'Number of Pairs', fontsize=25)
    ax.legend(loc='upper right')
    savefig(figdir+'/vis-hist.png')
    savefig(figdir+'/vis-hist.pdf')

    h = np.asarray(h,dtype='float')
    h_lens = np.asarray(h_lens,dtype='float')
    h_ocg = np.asarray(h_ocg,dtype='float')
    h_uvocg = np.asarray(h_uvocg,dtype='float')
    h_ee = np.asarray(h_ee,dtype='float')

    fig0 = plt.figure()
    ax = fig0.add_subplot(111)
    ax.bar(b[:-1], h_lens/h, 1, color='#ff3333')
    ax.bar(b[:-1], h_ocg/h, 1, color='#33ff33', bottom=h_lens/h)
    ax.bar(b[:-1], h_uvocg/h, 1, color='#3333ff', bottom=[(h_lens/h)[j] + (h_ocg/h)[j] for j in range(len(h_lens/h))])
    ax.bar(b[:-1], h_ee/h, 1, color='#33ffff', bottom=[(h_lens/h)[j] + (h_ocg/h)[j] + (h_uvocg/h)[j] for j in range(len(h_lens/h))])
    ax.axis([0,9,0,1], prop={'size':25})
    ax.set_xticks([0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5])
    ax.set_xticklabels(['SE','SS','ES','EE','Ef','M','E','S'])
    ax.set_xlabel(r'Visual Classification', fontsize=25)
    ax.set_ylabel(r'Fraction of Pairs', fontsize=25)
    ax.legend(loc='upper right')
    savefig(figdir+'/vis-bar.png')
    savefig(figdir+'/vis-bar.pdf')


    fig0 = plt.figure()
    ax = fig0.add_subplot(111)
    ax.bar(b[:-1], h_lens/h, 1, fill=False, hatch='/', rasterized=True)
    ax.bar(b[:-1], h_ocg/h, 1, fill=False, hatch='+', bottom=h_lens/h, rasterized=True)
    ax.bar(b[:-1], h_uvocg/h, 1, fill=False, hatch='.', bottom=[(h_lens/h)[j] + (h_ocg/h)[j] for j in range(len(h_lens/h))], rasterized=True)
    ax.bar(b[:-1], h_ee/h, 1, fill=False, bottom=[(h_lens/h)[j] + (h_ocg/h)[j] + (h_uvocg/h)[j] for j in range(len(h_lens/h))], rasterized=True)
    ax.axis([0,9,0,1], prop={'size':25})
    ax.set_xticks([0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5])
    ax.set_xticklabels(['SE','SS','ES','EE','Ef','M','E','S'])
    ax.set_xlabel(r'Visual Classification', fontsize=25)
    ax.set_ylabel(r'Fraction of Pairs', fontsize=25)
    ax.legend(loc='upper right')
    savefig(figdir+'/vis-bar-bw.png')
    savefig(figdir+'/vis-bar-bw.pdf', dpi=300)


    fig0 = plt.figure()
    ax = fig0.add_subplot(111)
    ax.bar(b[:-1], h_lens/h, 1, fill=True, color='k', alpha=0.6, rasterized=True)
    ax.bar(b[:-1], h_ocg/h, 1, fill=True, color='k', alpha=0.4, bottom=h_lens/h, rasterized=True)
    ax.bar(b[:-1], h_uvocg/h, 1, fill=True, color='k', alpha=0.2, bottom=[(h_lens/h)[j] + (h_ocg/h)[j] for j in range(len(h_lens/h))], rasterized=True)
    ax.bar(b[:-1], h_ee/h, 1, fill=True, color='k', alpha=0.0, bottom=[(h_lens/h)[j] + (h_ocg/h)[j] + (h_uvocg/h)[j] for j in range(len(h_lens/h))], rasterized=True)

    ax.annotate(r'P+EL', xy=(0.025, 0.9),  xycoords='axes fraction',
                horizontalalignment='left', verticalalignment='center')
    ax.annotate(r'EL+P', xy=(0.025, 0.7),  xycoords='axes fraction',
                horizontalalignment='left', verticalalignment='center')
    ax.annotate(r'EL+EL', xy=(0.025, 0.4),  xycoords='axes fraction',
                horizontalalignment='left', verticalalignment='center')
    ax.annotate(r'P+P', xy=(0.025, 0.1),  xycoords='axes fraction',
                horizontalalignment='left', verticalalignment='center')


    ax.axis([0,9,0,1], prop={'size':25})
    ax.set_xticks([0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5])
    ax.set_xticklabels(['SE','SS','ES','EE','Ef','M','E','S'])
    ax.set_xlabel(r'Visual Classification', fontsize=25)
    ax.set_ylabel(r'Fraction of Pairs', fontsize=25)
    ax.legend(loc='upper right')
    savefig(figdir+'/vis-bar-bw2.png')
    savefig(figdir+'/vis-bar-bw2.pdf', dpi=300)

    savefig('/Users/bholwerd/Dropbox/BlendedSpectra/PreparingPaper/Figures/vis-bar-bw2.pdf', dpi=300)

    fig0 = plt.figure()
    ax = fig0.add_subplot(111)
    ax.bar(b[:-1], h_lens, 1, fill=True, color='k', alpha=0.6, rasterized=True)
    ax.bar(b[:-1], h_ocg, 1, fill=True, color='k', alpha=0.4, bottom=h_lens, rasterized=True)
    ax.bar(b[:-1], h_uvocg, 1, fill=True, color='k', alpha=0.2, bottom=[(h_lens)[j] + (h_ocg)[j] for j in range(len(h_lens))], rasterized=True)
    ax.bar(b[:-1], h_ee, 1, fill=True, color='k', alpha=0.0, bottom=[(h_lens)[j] + (h_ocg)[j] + (h_uvocg)[j] for j in range(len(h_lens))], rasterized=True)
    ax.bar(b[:-1], h, 1, fill=False, edgecolor='k', rasterized=True)

    ax.annotate(r'P+EL', xy=(0.15, 53),  xycoords='data',
                horizontalalignment='left', verticalalignment='center')
    ax.annotate(r'EL+P', xy=(0.15, 40),  xycoords='data',
                horizontalalignment='left', verticalalignment='center')
    ax.annotate(r'EL+EL', xy=(0.15, 25),  xycoords='data',
                horizontalalignment='left', verticalalignment='center')
    ax.annotate(r'P+P', xy=(0.15, 7),  xycoords='data',
                horizontalalignment='left', verticalalignment='center')


    ax.axis([0,9,0,60], prop={'size':25})
    ax.set_xticks([0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5])
    ax.set_xticklabels(['SE','SS','ES','EE','Ef','M','E','S','Other'])
    ax.set_xlabel(r'Visual Classification', fontsize=25)
    ax.set_ylabel(r'Fraction of Pairs', fontsize=25)
    ax.legend(loc='upper right')
    savefig(figdir+'/vis-bar-bw3.png')
    savefig(figdir+'/vis-bar-bw3.pdf', dpi=300)

    # exit()
    
    fig0 = plt.figure()
    ax = fig0.add_subplot(111)
    ax.hist(template1, 30, range=[20,50],normed=0, edgecolor='g', histtype='step', label='Foreground Galaxy')
    ax.hist((template2+0.3),  30, range=[20,50], normed=0, edgecolor='b', histtype='step', label='Background Galaxy')
    ax.set_xlabel(r'Redshift ($z$)', fontsize=25)
    ax.set_ylabel(r'Number of Galaxies', fontsize=25)
    ax.legend(loc='upper left')
    savefig(figdir+'/template-hist-alt.png')
    savefig(figdir+'/template-hist-alt.pdf')



    fig0 = plt.figure()
    ax = fig0.add_subplot(111)
    # ax.hist(template1, 40, range=[20,50],normed=0, edgecolor='k', histtype='step', label='All Foreground Galaxy')
    ax.hist(template1[np.where(vc==1.)], 40, range=[20,50],normed=0, edgecolor='b', histtype='step', label='S-E')
    ax.hist(template1[np.where(vc==2.)], 40, range=[20,50],normed=0, edgecolor='g', histtype='step', label='S-S')
    ax.hist(template1[np.where(vc==4.)], 40, range=[20,50],normed=0, edgecolor='r', histtype='step', label='E-E')
    ax.hist(template1[np.where(vc==5.)], 40, range=[20,50],normed=0, edgecolor='m', histtype='step', label='Ef')
    ax.hist(template1[np.where(vc==3.)], 40, range=[20,50],normed=0, edgecolor='c', histtype='step', label='E-S')
    # 1 - SE
    # 1.1 - B
    # 2 - SS
    # 2.1 - F
    # 2.2 - Phi
    # 2.3 - Q
    # 3 - ES
    # 4 - EE
    # 5 - Ef
    # 6 - M
    # 7 - E
    # 8 - S

    ax.axis([20,50,0,50], prop={'size':25})

    ax.set_xlabel(r'Redshift ($z$)', fontsize=25)
    ax.set_ylabel(r'Number of Galaxies', fontsize=25)
    ax.legend(loc='upper left')
    savefig(figdir+'/template-fg-hist_vis.png')
    savefig(figdir+'/template_fg-hist_vis.pdf')

    #print template1[np.where(vc==2.)]
    #print template1[np.where(vc==3.)]
    #print template1[np.where(vc==4.)]
    # exit()

    # lenses
    # select passive fg & emission line background
    # plot histogram of redshifts
    fig0 = plt.figure()
    ax = fig0.add_subplot(111)
    ax.hist(z1[lens], 16, range=[0,0.8],normed=0, edgecolor='g', histtype='step', label='Foreground Galaxy')
    ax.hist(z2[lens], 16, range=[0,0.8],normed=0, edgecolor='b', histtype='step', label='Background Galaxy')
    ax.set_xlabel(r'Redshift ($z$)', fontsize=25)
    ax.set_ylabel(r'Number of Galaxies', fontsize=25)
    ax.legend(loc='upper right')
    savefig(figdir+'/Lens-z-hist.png')
    savefig(figdir+'/Lens-z-hist.pdf')

    fig0 = plt.figure()
    ax = fig0.add_subplot(111)
    ax.hist(z1[lens], 16, range=[0,0.8],normed=0, edgecolor='k', histtype='step', label='Foreground Galaxy')
    ax.hist(z2[lens], 16, range=[0,0.8],normed=0, edgecolor='k',alpha=0.5,linewidth=3, histtype='step', label='Background Galaxy')
    ax.set_xlabel(r'Redshift ($z$)', fontsize=25)
    ax.set_ylabel(r'Number of Galaxies', fontsize=25)
    ax.legend(loc='upper right')
    savefig(figdir+'/Lens-z-hist_bw.png')
    savefig(figdir+'/Lens-z-hist_bw.pdf')

    # perfect ocg
    # select emission line FG and passive background
    fig0 = plt.figure()
    ax = fig0.add_subplot(111)
    ax.hist(z1[ocg], 16, range=[0,0.8],normed=0, edgecolor='g', histtype='step', label='Foreground Galaxy')
    ax.hist(z2[ocg], 16, range=[0,0.8],normed=0, edgecolor='b', histtype='step', label='Background Galaxy')
    ax.set_xlabel(r'Redshift ($z$)', fontsize=25)
    ax.set_ylabel(r'Number of Galaxies', fontsize=25)
    ax.legend(loc='upper right')
    savefig(figdir+'/ocg-z-hist.png')
    savefig(figdir+'/ocg-z-hist.pdf')

    fig0 = plt.figure()
    ax = fig0.add_subplot(111)
    ax.hist(z1[ocg], 16, range=[0,0.8],normed=0, edgecolor='k', histtype='step', label='Foreground Galaxy')
    ax.hist(z2[ocg], 16, range=[0,0.8],normed=0, edgecolor='k', alpha=0.5, linewidth=3, histtype='step', label='Background Galaxy')
    ax.set_xlabel(r'Redshift ($z$)', fontsize=25)
    ax.set_ylabel(r'Number of Galaxies', fontsize=25)
    ax.legend(loc='upper right')
    savefig(figdir+'/ocg-z-hist_bw.png')
    savefig(figdir+'/ocg-z-hist_bw.pdf')


    fig0 = plt.figure()
    ax = fig0.add_subplot(111)
    ax.hist(z1[lens], 20, normed=0, edgecolor='g', linestyle='dashed', histtype='step', label='Lensing  Galaxy')
    ax.hist(z2[lens], 20, normed=0, edgecolor='b', linestyle='dashed', histtype='step', label='Lensed Galaxy')
    ax.hist(z1[ocg], 20, normed=0, edgecolor='g', histtype='step', label='Occulting Galaxy')
    ax.hist(z2[ocg], 20, normed=0, edgecolor='b', histtype='step', label='Occulted Galaxy')
    ax.set_xlabel(r'Redshift ($z$)', fontsize=25)
    ax.set_ylabel(r'Number of Galaxies', fontsize=25)
    ax.legend(loc='upper right')
    savefig(figdir+'/gama-z-hist.png')
    savefig(figdir+'/gama-z-hist.pdf')


    # mergers?
    # Delta-z < 600 km/s

    # histogram of delta-z
    fig0 = plt.figure()
    ax = fig0.add_subplot(111)
    ax.hist(z2-z1, 25, range=[0,0.5],normed=0, edgecolor='k', histtype='step')
    # fg_ba = 1.-fg_gal_e
    
    y = np.arange(0,40,0.01)
    x = 3./1000. +0.*y
    # ax.plot(x,y,'k--')
    ax.set_xlabel(r'Redshift ($\Delta z$)', fontsize=25)
    ax.set_ylabel(r'Number of Galaxies', fontsize=25)
    ax.axis([0,0.5,0,30], prop={'size':25})
    savefig(figdir+'/Delta-z-hist.png')
    savefig(figdir+'/Delta-z-hist.pdf')

    print 'number of mergers:'
    # print len(z1(np.where( (z2-z1) <= 3./1000.)))

    # make the catalog file
    pairtype = np.empty(len(z1), dtype=object)
    pairtype[lens] = 'ET+LT'
    pairtype[ocg] = 'LT+ET'
    pairtype[uvocg] = 'LT+LT'
    pairtype[EEocg] = 'ET+ET'

    # print pairtype

    field = np.empty(len(z1), dtype=object)

    g09 = np.where((ra > 129.0)&(ra<141.0))
    g12 = np.where((ra > 174.0)&(ra<186))
    g15 = np.where((ra > 211.5)&(ra<223.5))

    field[g09] = 'G09'
    field[g12] = 'G12'
    field[g15] = 'G15'

    

    cat = open('./cat.tex','w')

    for i in range(len(z1)):
        cat.write('%d \t & %f & %f & %4.2f & %d & %d \t & %4.2f & %d & %d & %s & %s & %s \\\\ \n' % (gama_id[i],        
                                                                                             ra[i],             
                                                                                             dec[i],            
                                                                                             z1[i],             
                                                                                             template1[i],      
                                                                                             cc_sigma1[i],     
                                                                                             z2[i],             
                                                                                             template2[i],      
                                                                                             cc_sigma2[i],      
                                                                                             field[i],
                                                                                             pairtype[i],
                                                                                                vis_class[i]))
    cat.close()


    cat = open('./HST_targets.txt','w')

    # potential HST follow-up targets:
    for i in range(len(z1)):
        if z1[i] < 0.05:
            cat.write('GAMA-%d \t %f %f %4.2f %d %d \t %4.2f %d %d %s %s \n' % (gama_id[i],        
                                                                            ra[i],             
                                                                            dec[i],            
                                                                            z1[i],             
                                                                            template1[i],      
                                                                            cc_sigma1[i],     
                                                                            z2[i],             
                                                                            template2[i],      
                                                                            cc_sigma2[i],      
                                                                            field[i],
                                                                            pairtype[i]))
    cat.close()


    #for idd in gama_id[ocg]:
    #    print '\includegraphics[width=0.24\textwidth]{./cutouts/%d_i}' % (idd)


    for idd in gama_id[uvocg]:
        print '\includegraphics[width=0.24\\textwidth]{./cutouts/%d_i}' % (idd)

mkplots()
