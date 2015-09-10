#!/usr/bin/env python

print('\n\nPLOTFITS: \n\tUsage: specify .fits (binary fits file) or .dat file as argv[1].\n'+\
        '\tThis program plots the Lyman series at any given redshift (z_sel)\n'+\
        '\tAND also: SEARCHES for metal lines\n'+\
        '\tand plots all lines for given metal, sorted with the strongest lines first (oscillator strengths)\n\n')

from astropy.io import fits
from astropy.io.votable import parse
import os.path
from numpy import *
import matplotlib.pylab as plt
from sys import argv,exit

from bokeh.plotting import figure, output_notebook, show, gridplot

cmap = plt.get_cmap('Set1')

class Plotfits:
    def __init__(self, file):
        self.datafileloc = str(file)
        self.datafile, self.extension = os.path.splitext(self.datafileloc)


    # Use running mean
    def runningMeanFast(self, x, N):
        # Input:    x   array to find running mean of
        #           N   window in points
        # NB! Gives fringing at the edges!!!!
        return convolve(x, ones((N,))/N)[(N-1):]

    def smooth(self, x,beta,window_len):
        """ kaiser window smoothing """
        # ref: http://glowingpython.blogspot.no/2012/02/convolution-with-numpy.html
        #window_len=11
        if (window_len%2 == 0):
            window_len += 1
            spacing_arr = int(window_len)/2
        else:
            spacing_arr = int(window_len)/2
        # extending the data at beginning and at the end
        # to apply the window at the borders
        s = r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
        w = kaiser(window_len,beta)
        y = convolve(w/w.sum(),s,mode='valid')
    #   return y[5:len(y)-5]
        print spacing_arr, window_len
        return y[spacing_arr:len(y)-spacing_arr]

    def load(self):
        if (self.extension == '.fits'):
            print 'Loading .fits file...'
            hdudata = fits.open(self.datafileloc)
            print hdudata.info()
            #raw_input('PRESS ENTER TO CONTINUE')

            # Obtain info necc to produce wavelength array
            start_pixel = float(hdudata[0].header['CRPIX1'])
            cntr_wl     = 10**float(hdudata[0].header['CRVAL1'])
            disp        = float(hdudata[0].header['CDELT1'])
            self.N           = int(hdudata[0].header['NAXIS1'])

            self.wls = zeros(self.N)
            for i in xrange(self.N):
                self.wls[i] = cntr_wl * 10 ** ((i+1-start_pixel) * disp)
            self.data = hdudata[0].data[0]
            # Thats it
        else:
            print 'Loading text file on form: [wavelength, data]'
            data_file   = loadtxt(self.datafileloc)
            self.N           = len(data_file)
            self.wls         = data_file[:,0]
            self.data        = data_file[:,1]

    def find_averages(self):
        #dataMEANlong = runningMeanFast(data, len(data)/10)
        #dataMEANshort= runningMeanFast(data, len(data)/3000)
        self.dataMEANlong = self.smooth(self.data, 4, len(self.data)/200)
        self.dataMEANshort= self.smooth(self.data, 32, len(self.data)/3000)
        print 'Done.'

    def define_redshifts(self, z_sel):
        # Eval z_sel: either list (vector) or float
        self.z_sel = array(map(float, z_sel.split(',')))

    def find_redshifted_lines(self):
        # next: read atomic.xml and find relevant lab wavelengths
        mypath = os.path.dirname(os.path.realpath(__file__))
        atoms = parse(os.path.expanduser(mypath+'/atomic.xml'))
        atomstable = atoms.get_first_table()
        self.atomstable = atomstable
        
        promels = ['Si', 'C', 'N', 'O', 'Al', 'Fe', 'Ca', 'B', 'Be']
        ionallow= ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'IX', 'X']
        isoallow= [28, 12, 14, 16, 27, 56, 40, 10, 9]
        wlshift = 5.0   # WL displacement in A
        wls     = self.wls
        z_sel   = self.z_sel

        lymaninf = where((atomstable.array['Element'] == 'H') \
                & (atomstable.array['Ion'] == 'I') \
                & (atomstable.array['MassNumber'] == 1) )
        deutinf  = where((atomstable.array['Element'] == 'H') \
                & (atomstable.array['Ion'] == 'I') \
                & (atomstable.array['MassNumber'] == 2) )
        promlinesBOOL = in1d(atomstable.array['Element'], promels)
        promlinesIonBOOL = in1d(atomstable.array['Ion'], ionallow)
        promlinesIsoBOOL = in1d(atomstable.array['MassNumber'], isoallow)
        promlines = where((promlinesBOOL) \
                & (promlinesIonBOOL) \
                & (promlinesIsoBOOL) \
                & (atomstable.array['fval'] >= 0.05) \
                & (atomstable.array['RestWave'] <= max(wls)/(1+z_sel[0]) )\
                & (atomstable.array['RestWave'] >= min(wls)/(1+z_sel[0]) ) )
        implines = where((atomstable.array['fval'] >= 0.05) \
                & (atomstable.array['RestWave'] <= max(wls)/(1+z_sel[0]) )\
                & (atomstable.array['RestWave'] >= min(wls)/(1+z_sel[0]) ) )

        # Sort prominent lines, most important first
        promlines_id = atomstable.array['fval'][promlines].argsort()
        self.promlines1 = promlines[0][promlines_id]  # sort by strenght
        self.promlines1 = self.promlines1[::-1]

        self.promlines_wl = atomstable.array['RestWave'][promlines].argsort()
        self.promlines2  = promlines[0][promlines_id]

        self.lymanseries = outer((1 + z_sel),  atomstable.array['RestWave'][lymaninf])
        self.lymanserieslab = atomstable.array['RestWave'][lymaninf]
        self.deutseries  = outer((1 + z_sel),  atomstable.array['RestWave'][deutinf])
        self.promseries1  = outer((1 + z_sel), atomstable.array['RestWave'][self.promlines1])
        self.promseries2  = outer((1 + z_sel), atomstable.array['RestWave'][self.promlines2])

        self.promels  =  promels 
        self.ionallow =  ionallow
        self.isoallow =  isoallow
        self.wlshift  =  wlshift 

    def find_important_lines(self):
        # Choose
        # lines where the running mean is smaller than the larger running mean
        # AND
        # lines that are in the IMPLINES array
        wlsINT      = self.wls.astype(int)
        linesINT    = self.promseries2.astype(int)
        #linesINT    = ((1+z_sel) * atomstable.array['RestWave']).astype(int)
        linesCandidates= in1d(wlsINT,linesINT)
        # Choose those who have absorption features
        self.find_averages()
        absorption  = where((linesCandidates) \
                #& (dataMEANshort < dataMEANlong))
                & (self.dataMEANshort/self.dataMEANlong < 0.9))
        linesPresent = in1d(linesINT[0],wlsINT[absorption]) # BOOL on linesINT
        self.seriesPresent = zeros((len(self.z_sel),sum(linesPresent)))  # must: sum as only valid elements are incl in fin arr
        for i in range(len(self.z_sel)):
            self.seriesPresent[i] = self.promseries2[i][linesPresent]    # wavelengths
        self.linesPresentWL=self.promlines2[linesPresent]     # index

        self.prom_els1= min(24, len(self.promseries1[0]))
        self.prom_els2= min(24, len(self.seriesPresent[0]))

    def plot_hydrogen(self):
        # Link
        prom_els1 = self.prom_els1 
        prom_els2 = self.prom_els2
        wls         = self.wls
        z_sel       = self.z_sel
        data        = self.data

        datafileloc = self.datafileloc

        lymanseries = self.lymanseries
        deutseries  = self.deutseries
        lymanserieslab  = self.lymanserieslab
        wlshift     = self.wlshift
        seriesPresent   = self.seriesPresent
        atomstable      = self.atomstable

        
        

        fig2, axes      = plt.subplots(12,1, figsize=(10,20))
        fig2b           = figure()  # Bokeh
        fig3, axes3     = plt.subplots(12,1, figsize=(10,20))
        fig3b, axes3b   = plt.subplots(8,1, figsize=(10,20))
        #fig3, axes3= plt.subplots(prom_els1/2,2)
        fig4, axes4     = plt.subplots((prom_els2+1)/2,2, figsize=(10,20))

        # Make axes iterable, going from start to end (no higher dim)
        axes = axes.ravel()
        axes3= axes3.ravel()
        axes4= axes4.ravel()

        zcolor  = [(cmap(i*255/len(z_sel))) for i in range(len(z_sel))]

        output_notebook();
        # Bokeh
        b_fig = figure(x_range=[lymanseries[0][3] - wlshift,\
                                      lymanseries[0][3] + wlshift], \
                             y_range=[-0.1, 1.])
        b_fig.line(wls, data)
        for k in range(len(z_sel)):
            b_fig.line([lymanseries[k][3],lymanseries[k][3]], [-1,1],\
                    color=zcolor[k])
            b_fig.line([deutseries[k][0], deutseries[k][0]], [-1,1], color='black')
        #axes_int[i].set_xlim(lymanseries[0][i+3] - wlshift, lymanseries[0][i+3] + wlshift)
        #axes_int[i].set_ylim(-0,1.1)
        #axes_int[i].tick_params(which='both', labelsize=8)
        #axes_int[i].locator_params(nbins=3, axis='y')
        #axes_int[i].text(1.07, 0.50, (r'Ly-%i: %.2f $\AA$' % (i+1, lymanserieslab[i+3],)), \
        #        verticalalignment='center', horizontalalignment='center', \
        #        size='small',transform=axes[i].transAxes)
        show(b_fig)

        for i in range(12):
            # Standard 
            axes[i].plot(wls, data)
            for k in range(len(z_sel)):
                zlines, = axes[i].plot([lymanseries[k][i+3],lymanseries[k][i+3]], [-1,1],\
                        color=zcolor[k])
            axes[i].plot([deutseries[0][i], deutseries[0][i]], [-1,1], 'k--')
            axes[i].set_xlim(lymanseries[0][i+3] - wlshift, lymanseries[0][i+3] + wlshift)
            axes[i].set_ylim(-0,1.1)
            axes[i].tick_params(which='both', labelsize=8)
            axes[i].locator_params(nbins=3, axis='y')
            axes[i].text(1.07, 0.50, (r'Ly-%i: %.2f $\AA$' % (i+1, lymanserieslab[i+3],)), \
                    verticalalignment='center', horizontalalignment='center', \
                    size='small',transform=axes[i].transAxes)

            # Ly 13-24
            axes3[i].plot(wls, data)
            for k in range(len(z_sel)):
                axes3[i].plot([lymanseries[k][i+3+12],lymanseries[k][i+3+12]], [-1,1], \
                        color=zcolor[k])
            if (i+12 < 18):
                axes3[i].plot([deutseries[0][i+12], deutseries[0][i+12]], [-1,1], 'k--')
            axes3[i].set_xlim(lymanseries[0][i+3+12] - wlshift, lymanseries[0][i+3+12] + wlshift)
            axes3[i].set_ylim(-0,1.1)
            axes3[i].tick_params(which='both', labelsize=8)
            axes3[i].locator_params(nbins=3, axis='y')
            axes3[i].text(1.07, 0.50, (r'Ly-%i: %.2f $\AA$' % (i+1+12, lymanserieslab[i+3+12],)), \
                    verticalalignment='center', horizontalalignment='center', \
                    size='small',transform=axes3[i].transAxes)

            # Ly 25-33
            if (i+12+12 < 32):
                axes3b[i].plot(wls, data)
                for k in range(len(z_sel)):
                    axes3b[i].plot([lymanseries[k][i+3+12+12],lymanseries[k][i+3+12+12]], [-1,1], \
                            color=zcolor[k])
                axes3b[i].set_xlim(lymanseries[0][i+3+12+12] - wlshift, lymanseries[0][i+3+12+12] + wlshift)
                axes3b[i].set_ylim(-0,1.1)
                axes3b[i].tick_params(which='both', labelsize=8)
                axes3b[i].locator_params(nbins=3, axis='y')
                axes3b[i].text(1.07, 0.50, (r'Ly-%i: %.2f $\AA$' % (i+1+12+12,lymanserieslab[i+3+12+12],)), \
                        verticalalignment='center', horizontalalignment='center', \
                        size='small',transform=axes3b[i].transAxes)
        for m in range(len(z_sel)):
            axes[0].text(-0.05,-m*0.8, 'z= '+str(z_sel[m]), color='white', \
                    verticalalignment='bottom', horizontalalignment='right', \
                    bbox={'facecolor':zcolor[m], 'alpha':1.0, 'pad':10}, \
                    size='small', transform=axes[0].transAxes)


        fig2.subplots_adjust(bottom=0.03)
        fig2.subplots_adjust(top=0.93)
        fig2.subplots_adjust(hspace=0.55)
        axes[0].set_title("H (red) and D (black) Lyman series 1-12, z="+str(z_sel[0]))
        fig3.subplots_adjust(bottom=0.03)
        fig3.subplots_adjust(top=0.93)
        fig3.subplots_adjust(hspace=0.55)
        axes3[0].set_title("H (red) and D (black) Lyman series 13-24, z="+str(z_sel[0]))
        fig3b.subplots_adjust(bottom=0.03)
        fig3b.subplots_adjust(top=0.93)
        fig3b.subplots_adjust(hspace=0.55)
        axes3b[0].set_title("H (red)Lyman series 25-32, z="+str(z_sel[0]))
        
        for i in range(prom_els2):
            print i
            # Other relevant lines?
            axes4[i].plot(wls, data)
            for k in range(len(z_sel)):
                axes4[i].plot([seriesPresent[k][i], seriesPresent[k][i]], [-1,1], \
                        color=zcolor[k])
            axes4[i].set_xlim(seriesPresent[0][i] - wlshift, seriesPresent[0][i] + wlshift)
            axes4[i].set_ylim(0,1.1)
            elementText = \
                    str(atomstable.array['MassNumber'][self.linesPresentWL[i]])+\
                    atomstable.array['Element'][self.linesPresentWL[i]]+' '+\
                    atomstable.array['Ion'][self.linesPresentWL[i]] + ', wl='+\
                    str(atomstable.array['RestWave'][self.linesPresentWL[i]])
            axes4[i].text(0.95, 0.01, elementText, \
                    verticalalignment='bottom', horizontalalignment='right', \
                    bbox={'facecolor':'white', 'alpha':0.7, 'pad':10}, \
                    transform=axes4[i].transAxes)
            axes4[i].tick_params(which='both', labelsize=8)
            plt.subplots_adjust(bottom=0.03)
            plt.subplots_adjust(top=0.93)
        axes4[0].set_title("Prominent metal lines, z="+str(z_sel))
        axes4[1].set_title(datafileloc, size='small')
        for m in range(len(z_sel)):
            axes4[0].text(-0.05,-m*0.8, 'z= '+str(z_sel[m]), color='white', \
                    verticalalignment='bottom', horizontalalignment='right', \
                    bbox={'facecolor':zcolor[m], 'alpha':1.0, 'pad':10}, \
                    size='small', transform=axes4[0].transAxes)

        plt.show()


    def find_metal_lines(self, sel_el):
        """ Plots metal lines, 
            only transitions in sel_el:
            0: 'Si',\n1: 'C',\n2: 'N',\3: 'O',\n4: 'Al',\n5: 'Fe',\n6: 'Ca',\7: 'B',\n8: 'Be' """

        atomstable = self.atomstable
        ionallow    = self.ionallow
        isoallow    = self.isoallow
        promels     = self.promels

        wls         = self.wls
        z_sel       = self.z_sel


        # Find the elements in atomic.xml that fulfills reqs:
        sel_el_ion = in1d(atomstable.array['Ion'], ionallow)
        sel_el_iso = in1d(atomstable.array['MassNumber'], isoallow)
        sel_el_inf = where(sel_el_ion \
                & (sel_el_iso) \
                & (atomstable.array['Element'] == promels[sel_el])  \
                & (atomstable.array['RestWave'] <= max(wls)/(1+z_sel[0]) )\
                & (atomstable.array['RestWave'] >= min(wls)/(1+z_sel[0]) ) )

        # Sort these elements based on oscillator strength
        sel_el_osc = atomstable.array['fval'][sel_el_inf].argsort()
        sel_el_osc = sel_el_osc[::-1]

        # Have now found the location of the elements, next: find relev. wavelengths
        sel_el_series_ = outer((1 + z_sel), atomstable.array['RestWave'][sel_el_inf[0]] )
        sel_el_series  = zeros((len(z_sel), len(sel_el_osc)))
        for i in range(len(z_sel)):
            sel_el_series[i] = sel_el_series_[i][sel_el_osc]
        self.sel_el_inf    = sel_el_inf[0][sel_el_osc]

        self.sel_el_series = sel_el_series

    def plot_metal_lines(self):
        # Link
        #prom_els1 = self.prom_els1 
        #prom_els2 = self.prom_els2
        wls         = self.wls
        z_sel       = self.z_sel
        data        = self.data

        datafileloc = self.datafileloc

        #lymanseries = self.lymanseries
        #deutseries  = self.deutseries
        #lymanserieslab  = self.lymanserieslab
        wlshift     = self.wlshift
        seriesPresent   = self.seriesPresent
        atomstable      = self.atomstable

        sel_el_series = self.sel_el_series
        sel_el_inf    = self.sel_el_inf

        # Plot em
        plotwindows_N = len(sel_el_series[0])/12
        plotwindows   = []
        plotwindows  += [12] * plotwindows_N
        if (len(sel_el_series[0])%12 != 0):
            plotwindows.append(len(sel_el_series[0])%12)
        if (len(sel_el_series[0]) == 0):
            print 'No transitions to plot.'

        for j in range(len(plotwindows)):
            curplot, curaxes = plt.subplots(plotwindows[j],1, figsize=(10,20))
            curplot, curaxes = [curplot], [curaxes]
            #curaxes = curaxes.ravel()
            zcolor  = [(cmap(i*255/len(z_sel))) for i in range(len(z_sel))]
            for i in range(plotwindows[j]):
                if plotwindows[j] > 1:
                    curaxes[0][i].plot(wls, data)
                    for k in range(len(z_sel)):
                        zline, = curaxes[0][i].plot([sel_el_series[k][i+12*j]]*2, [-1,1], color=zcolor[k])
                    if (i == 0):
                        # Add z-legend at left-hand side:
                        for m in range(len(z_sel)):
                            curaxes[0][0].text(-0.05,-m*0.8, 'z= '+str(z_sel[m]), \
                                    verticalalignment='bottom', horizontalalignment='right', \
                                    bbox={'facecolor':zcolor[m], 'alpha':1.0, 'pad':10}, \
                                    size='small', transform=curaxes[0][0].transAxes)
                    curaxes[0][i].set_xlim(sel_el_series[0][i+12*j] - wlshift, \
                            sel_el_series[0][i+12*j] + wlshift)
                    curaxes[0][i].set_ylim(0, 1.1)
                    curaxes[0][i].tick_params(which='both', labelsize=8)
                    curaxes[0][i].locator_params(nbins=3, axis='y')
                    elementText = ('%d%s %s, wl=%.3f shifted=%.1f' % \
                            (atomstable.array['MassNumber'][sel_el_inf[i+12*j]], \
                            atomstable.array['Element'][sel_el_inf[i+12*j]], \
                            atomstable.array['Ion'][sel_el_inf[i+12*j]], \
                            atomstable.array['RestWave'][sel_el_inf[i+12*j]], \
                            atomstable.array['RestWave'][sel_el_inf[i+12*j]] * (1 + z_sel[0]) ) )
                    curaxes[0][i].text(0.95, 0.01, elementText, \
                            verticalalignment='bottom', horizontalalignment='right', \
                            bbox={'facecolor':'white', 'alpha':0.7, 'pad':10}, \
                            transform=curaxes[0][i].transAxes)
                    curTempColor = {'I':'cyan', 'II':'green', 'III':'yellow', 'IV':'orange','V':'red','VI':'purple','VII':'blue'}
                    curIon = atomstable.array['Ion'][sel_el_inf[i+12*j]]
                    curaxes[0][i].text(1.05, 0.50, curIon, \
                            verticalalignment='center', horizontalalignment='center', \
                            bbox={'facecolor':curTempColor[curIon], 'alpha':0.7, 'pad':10}, \
                            transform=curaxes[0][i].transAxes)
                    curaxes[0][i].tick_params(which='both', labelsize=8)
                else:
                    curaxes[0].plot(wls, data)
                    for k in range(len(z_sel)):
                        curaxes[0].plot([sel_el_series[k][i+12*j]]*2, [-1,1])
                    curaxes[0].set_xlim(sel_el_series[0][i+12*j] - wlshift, \
                            sel_el_series[0][i+12*j] + wlshift)
                    curaxes[0].set_ylim(0, 1.1)
                    curaxes[0].tick_params(which='both', labelsize=8)
                    curaxes[0].locator_params(nbins=3, axis='y')
                    elementText = \
                        str(atomstable.array['MassNumber'][sel_el_inf[i+12*j]])+\
                        atomstable.array['Element'][sel_el_inf[i+12*j]]+' '+\
                        atomstable.array['Ion'][sel_el_inf[i+12*j]] + ', wl='+\
                        str(atomstable.array['RestWave'][sel_el_inf[i+12*j]])
                    curaxes[0].text(0.95, 0.01, elementText, \
                        verticalalignment='bottom', horizontalalignment='right', \
                        bbox={'facecolor':'white', 'alpha':0.7, 'pad':10}, \
                        transform=curaxes[0].transAxes)
                    curaxes[0].tick_params(which='both', labelsize=8)
                # Remove some spacing
                curplot[0].subplots_adjust(bottom=0.03)
                curplot[0].subplots_adjust(top=0.91)
                curplot[0].subplots_adjust(hspace=0.55)
                # Add title
                curplot[0].suptitle(datafileloc+" z="+str(z_sel)+'\n'+ \
                        'Sorted: highest oscillator strenghts first')
        plt.show()

    def plot_spectrum(self):
        ##fig1, ax1 = plt.subplots(1, figsize=(20,10))
        #ax1.plot(self.wls, self.data, self.wls, self.dataMEANlong, self.wls, self.dataMEANshort)
        #plt.show()

        output_notebook()
        p = figure()
        p.line(self.wls, self.data)
        p.line(self.wls, self.dataMEANlong, legend='Continuum Mean')
        p.line(self.wls, self.dataMEANshort, legend='Line mean')
        show(p)




    ## PROJECT: FIND SYSTEMS
    ## Current approach: use number of hits on transitions
##    z_pos = linspace(0,4,4001)
##    sys_confirmed = []
##    for i in xrange(len(z_pos)):
##        tmp_promlines = where((promlinesBOOL) \
##            & (promlinesIonBOOL) \
##            & (promlinesIsoBOOL) \
##            & (atomstable.array['fval'] >= 0.05) \
##            & (atomstable.array['RestWave'] <= max(wls)/(1+z_pos[i]) )\
##            & (atomstable.array['RestWave'] >= min(wls)/(1+z_pos[i]) ) )
##        promseries3 = (1 + z_pos[i]) * atomstable.array['RestWave'][tmp_promlines]
##        linesINT3 = promseries3.astype(int)
##        linesCandidates = in1d(wlsINT,linesINT3)
##        absorption = where((linesCandidates) \
##                & (dataMEANshort / dataMEANlong < 0.9 ))
##        linesPresent = in1d(linesINT3,wlsINT[absorption])
##        seriesPresent = zeros(1,sum(linesPresent))
##
##        seriesPresent = promseries3[linesPresent]
##        print len(seriesPresent)
##        if len(seriesPresent) > 13:
##            sys_confirmed.append([z_pos[i],len(seriesPresent)])
