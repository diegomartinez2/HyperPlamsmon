#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Excel_read.py
#
#  Copyright 2024 Diego Martinez Gutierrez <diego.martinez@ehu.eus>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#
# ---------------------------
# Importación de los módulos
# ---------------------------

# -------
# Clases
# -------
class Plasmon_analysis(object):
    """Code for the analysis of Silkin's plasmon data."""

    def __init__(self, arg, namefile):
        super(Plasmon_analysis, self).__init__()
        self.arg = arg
        self.namefile = namefile
        self.pars = np.zeros((51,3))
        self.pars2 = np.zeros((51,3))
        self.perr_lorentz = np.zeros((51,3))
        self.perr_lorentz2 = np.zeros((51,3))

    def load_data(self):
        with open(self.namefile) as file:
             lines = [line.rsplit() for line in file]
        Spec = np.zeros(255051)
        for i in range(255051):
             Spec[i]=lines[i][6]
        #data = np.resize(Spec,(51,int(len(Spec)/51)))
        data = np.resize(Spec,(51,5001))
        Frequency = np.zeros(5001)
        for i in range(5001):
             Frequency[i]=lines[i][2]
        print ("q_x=",lines[0][0])
        return data, Frequency, lines[0][0]

    def load_big_file(self, index, filename="A7_EPS.dat", diagonal=False):
        q_x, q_y, w, epsilon = np.loadtxt(filename,usecols=(0,1,2,6), unpack=True)
        print (len(np. unique(q_x))) #np. unique(my_array, return_counts=True)
        print (len(np. unique(q_y)))
        print (len(np. unique(w)))
        self.data = np.resize(epsilon,(len(np. unique(q_x)),
                                        len(np. unique(q_y)),
                                        len(np. unique(w))))
        """
        direccion q_y -> data[i]
        direccion q_x=q_y=i -> data = data[i][i] for i in {0..51}
        direccion
        """
        data = np.zeros((len(np. unique(q_x)),len(np. unique(w))))
        print ("+++++")
        if (diagonal == True):
            for i in range(len(np. unique(q_x))):
                print ("[" + "=" * i +  " " * ((len(np. unique(q_x)) - i)) + "]" +  str(i) + "%",end='\r', flush=True)
                for j in range(len(np. unique(w))):
                    data[i,j] = self.data[i,i,j]

        print ("*************")
        print (data.shape)
        print ("--------------")

        print (np.shape(self.data[int(index)]),"=(51,5001)?")
        if (diagonal == True):
            return data, w[:5001]
        else:
            return self.data[int(index)], w[:5001]

    def plot_contour(self,data):
        plt.style.use('_mpl-gallery-nogrid')
        fig, ax1 = plt.subplots(1,1)
        fig.set_size_inches(10, 5)
        fig.set_dpi(100)
        cax = ax1.imshow(data.T,
        #	vmin = 0.0 , vmax = 0.004,
        	vmin = 0.0 , vmax = 0.2,
        	cmap=plt.colormaps['jet'], origin='lower',
        	interpolation='gaussian', aspect='auto')
        ax1.set_ylabel(r'Frequency (cm$^{-1}$)', fontsize=12)

        cbar = fig.colorbar(cax)
        plt.tight_layout()
        plt.show()
        pass

    def plot_contour2(self,data,data2):
        plt.style.use('_mpl-gallery-nogrid')
        fig, ax1 = plt.subplots(2,1)
        fig.set_size_inches(10, 5)
        fig.set_dpi(100)
        cax = ax1[0].imshow(data.T,
        #	vmin = 0.0 , vmax = 0.004,
        #	vmin = 0.0 , vmax = 0.175,
        	vmin = 0.0 , vmax = 0.12,
        	cmap=plt.colormaps['jet'], origin='lower',
        	interpolation='gaussian', aspect='auto')
        cax2 = ax1[1].imshow(data2.T,
        #	vmin = 0.0 , vmax = 0.004,
        #	vmin = 0.0 , vmax = 0.3,
        	cmap=plt.colormaps['jet'], origin='lower',
        	interpolation='gaussian', aspect='auto')
        #ax1[1].set_ylabel(r'Frequency (cm$^{-1}$)', fontsize=12)
        print (self.pars[:,1])
        cax3 = ax1[1].scatter(x = range(len(self.pars[:,1])), y = (self.pars[:,1])*5001,c='k',marker='x',s=10)
        cax3 = ax1[1].scatter(x = range(len(self.pars[:,1])), y = (self.pars[:,1]+self.pars[:,2])*5001,c='k',marker='1',s=10)
        cax3 = ax1[1].scatter(x = range(len(self.pars[:,1])), y = (self.pars[:,1]-self.pars[:,2])*5001,c='k',marker='2',s=10)
        cax3 = ax1[1].axvline(x = len(self.pars[:,1]), linestyle = "--", color = "red")

        cbar = fig.colorbar(cax)
        cbar2 = fig.colorbar(cax2)
        if (len( sys.argv) == 3):
            plt.savefig("Ajuste_{}".format("original_vs_fitted"))
        else:
            plt.savefig("Ajuste_{}_{}".format(self.namefile,"original_vs_fitted"))
        plt.tight_layout()
        plt.show()
        pass

    def plot_contour3(self,data,data2):
        plt.style.use('_mpl-gallery-nogrid')
        fig, ax1 = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=True)
        fig.set_size_inches(10, 5)
        fig.set_dpi(100)
        cax = ax1[0].imshow(data.T,
        #	vmin = 0.0 , vmax = 0.004,
        #	vmin = 0.0 , vmax = 0.175,
        	vmin = 0.0 , vmax = 0.12,
        	cmap=plt.colormaps['jet'], origin='lower',
        	interpolation='gaussian', aspect='auto')
        #plt.tick_params(direction='in', length=6, width=2, colors='k', right=True, labelright='on')
        cax2 = ax1[1].imshow(data2.T,
        #	vmin = 0.0 , vmax = 0.004,
        #	vmin = 0.0 , vmax = 0.3,
            vmin = 0.0 , vmax = 0.12,
        	cmap=plt.colormaps['jet'], origin='lower',
        	interpolation='gaussian', aspect='auto')
        x_max_range = len(self.pars[:,1])
        for xi in range(x_max_range):
            cax3 = ax1[1].scatter(x = x_max_range-xi, y = self.pars[xi,1]*5001,c='k',marker='x',s=10)
            y_min =(self.pars[xi,1]-self.pars[xi,2])*5001
            y_max =(self.pars[xi,1]+self.pars[xi,2])*5001
            cax3 = ax1[1].scatter(x = x_max_range-xi, y = y_min,c='k',marker='1',s=10)
            cax3 = ax1[1].scatter(x = x_max_range-xi, y = y_max,c='k',marker='2',s=10)

            cax3 = ax1[1].scatter(x = x_max_range+xi, y = self.pars2[xi,1]*5001,c='k',marker='x',s=10)
            y_min =(self.pars2[xi,1]-self.pars2[xi,2])*5001
            y_max =(self.pars2[xi,1]+self.pars2[xi,2])*5001
            cax3 = ax1[1].scatter(x = x_max_range+xi, y = y_min,c='k',marker='1',s=10)
            cax3 = ax1[1].scatter(x = x_max_range+xi, y = y_max,c='k',marker='2',s=10)
        cax3 = ax1[1].axvline(x = len(self.pars[:,1]), linestyle = "--", color = "red")

        cbar = fig.colorbar(cax)
        cbar2 = fig.colorbar(cax2)
        ax1[0].set_title('Original data')
        ax1[0].set_ylabel('$\omega$ (eV)')
        ax1[1].set_ylabel('$\omega$ (eV)')
        ax1[0].set_xlabel('$q_x$ (a.u.)')
        ax1[1].set_xlabel('$q_x$ (a.u.)')
        ax1[1].set_title('Lorentz fitting')
        ax1[0].set_xticks([0,51,102])
        ax1[0].set_yticks([0,5001])
        ax1[0].set_xticklabels(["($\pi$,0)","(0,0)","($\pi$,$\pi$)"])
        ax1[0].set_yticklabels(["0","1"])

        plt.tight_layout()
        if (len( sys.argv) == 3):
            plt.savefig("Ajuste_d_{}".format("original_vs_fitted"))
        else:
            plt.savefig("Ajuste_d_{}_{}".format(self.namefile,"original_vs_fitted"))
        plt.tight_layout()
        plt.show()
        pass

    def locate_1Lorenztian(self, Frequency, Spec, index, big=False):
        def _1Lorentzian(x, amp1, cen1, wid1):
            return amp1*wid1**2/((x-cen1)**2+wid1**2)
        peaks, _ = find_peaks(Spec,width=5,rel_height=0.3)
        amp = 0.009
        cen = 0.009
        wid = 0.001
        popt_lorentz, pcov_lorentz = scipy.optimize.curve_fit(_1Lorentzian, Frequency, Spec, p0=[amp, cen, wid])

        perr_lorentz = np.sqrt(np.diag(pcov_lorentz))
        pars = popt_lorentz[:]
    #---------------------BIG---------------------------------------------------
        if big:
            if (Frequency[0] != 0):
                Frequency_big = np.append(Frequency,Frequency[:len(Frequency)//3]+Frequency[-1])
            else:
                Frequency_big = np.append(Frequency,Frequency[1:len(Frequency)//3]+Frequency[-1])
            lorentz_peak = _1Lorentzian(Frequency_big, *pars)
        else:
            lorentz_peak = _1Lorentzian(Frequency, *pars)
    #-------------------END_BIG-------------------------------------------------
        print ("-------------Peak ",index,"-------------")
        print ("amplitude = %0.2f (+/-) %0.2f" % (pars[0], perr_lorentz[0]))
        print ("center = %0.2f (+/-) %0.2f" % (pars[1], perr_lorentz[1]))
        print ("width = %0.2f (+/-) %0.2f" % (pars[2], perr_lorentz[2]))
        print ("area = %0.2f" % np.trapz(lorentz_peak))
        print ("--------------------------------")

        self.pars[index] = pars
        self.perr_lorentz[index] = perr_lorentz
        return lorentz_peak

    def fitting_Lorentz(self,frequencies,data,size,big=False):
        #self.Fitted_data = np.zeros((size,5001))
        self.Fitted_data = np.zeros((size,len(frequencies)))
        for i in range(size):
            self.Fitted_data[i] = self.locate_1Lorenztian(frequencies,data[i],i)
        #---------------------------new-----------------------------------------
        if big:
            if (frequencies[0] != 0):
                frequencies_big = np.append(frequencies,frequencies[:len(frequencies)//3]+frequencies[-1])
            else:
                frequencies_big = np.append(frequencies,frequencies[1:len(frequencies)//3]+frequencies[-1])
            self.Fitted_data = np.zeros((size,len(frequencies_big)))
            for i in range(size):
                self.Fitted_data[i] = self.locate_1Lorenztian(frequencies,data[i],i,big=True)
        #-----------------------------------------------------------------------
        pass

    def fitting_Lorentz2(self,frequencies,data,size,big=False):
        #self.Fitted_data2 = np.zeros((size,5001))
        self.Fitted_data2 = np.zeros((size,len(frequencies)))
        for i in range(size):
            self.Fitted_data2[i] = self.locate_1Lorenztian2(frequencies,data[i],i)
        #---------------------------new-----------------------------------------
        if big:
            if (frequencies[0] != 0):
                frequencies_big = np.append(frequencies,frequencies[:len(frequencies)//3]+frequencies[-1])
            else:
                frequencies_big = np.append(frequencies,frequencies[1:len(frequencies)//3]+frequencies[-1])
            self.Fitted_data2 = np.zeros((size,len(frequencies_big)))
            for i in range(size):
                self.Fitted_data2[i] = self.locate_1Lorenztian2(frequencies,data[i],i,big=True)
        #-----------------------------------------------------------------------
        pass

    def locate_1Lorenztian2(self, Frequency, Spec, index, big=False):
        def _1Lorentzian(x, amp1, cen1, wid1):
            return amp1*wid1**2/((x-cen1)**2+wid1**2)
        peaks, _ = find_peaks(Spec,width=5,rel_height=0.3)
        amp = 0.009
        cen = 0.009
        wid = 0.001
        popt_lorentz, pcov_lorentz = scipy.optimize.curve_fit(_1Lorentzian, Frequency, Spec, p0=[amp, cen, wid])

        perr_lorentz = np.sqrt(np.diag(pcov_lorentz))
        pars = popt_lorentz[:]
    #---------------------BIG---------------------------------------------------
        if big:
            if (Frequency[0] != 0):
                Frequency_big = np.append(Frequency,Frequency[:len(Frequency)//3]+Frequency[-1])
            else:
                Frequency_big = np.append(Frequency,Frequency[1:len(Frequency)//3]+Frequency[-1])
            lorentz_peak = _1Lorentzian(Frequency_big, *pars)
        else:
            lorentz_peak = _1Lorentzian(Frequency, *pars)
    #-------------------END_BIG-------------------------------------------------
        print ("-------------Peak ",index,"-------------")
        print ("amplitude = %0.2f (+/-) %0.2f" % (pars[0], perr_lorentz[0]))
        print ("center = %0.2f (+/-) %0.2f" % (pars[1], perr_lorentz[1]))
        print ("width = %0.2f (+/-) %0.2f" % (pars[2], perr_lorentz[2]))
        print ("area = %0.2f" % np.trapz(lorentz_peak))
        print ("--------------------------------")

        self.pars2[index] = pars
        self.perr_lorentz2[index] = perr_lorentz

        return lorentz_peak



# ----------
# Funciones
# ----------


def main(args):
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
