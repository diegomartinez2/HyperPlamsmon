#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Plasmon.py
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
import Eliashberg
import numpy as np
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


def main(arg):
    #namefile = "xbv"
    #print (len(sys.argv))
    if (len( sys.argv ) > 1):
        print ("nombre de fichero:",arg[1])
        namefile = arg[1]
        print ("Creando el objeto plasmon")
        plasmon = Plasmon_analysis(arg,namefile)
        print ("Leyendo datos")

        if (len( sys.argv) == 3):
            index = arg[2]
            #data, frequencies = plasmon.load_big_file(index, arg[1])
            data, frequencies = plasmon.load_big_file(index, arg[1],diagonal=False)
            data_d, frequencies_d = plasmon.load_big_file(index, arg[1], diagonal=True)
            #data.append(data_d)
        #    data = np.vstack((np.flip(data, axis=1), data_d))
            #frequencies.append(frequencies_d)
        #    frequencies = np.vstack((np.flip(frequencies), frequencies_d))
        else:
            data, frequencies, qx= plasmon.load_data()
#--------------------------diagonal-----------------------------
#        if (len(sys.argv[3]) == 4):
#            for i in range(51):
#                data[i,:] = plasmon.data[i,i,:]
#---------------------------------end---diagonal----------------
#-------------------------fix to remove extra polaron data--------------
        if (len( sys.argv) == 3):
            for i in range(5):
                for j in range(2000,5001,1):
                    if (data[i][j]!=0.0):
                        print("data[{}][{}]=".format(i,j),data[i][j])
            # if i < 5:
            #     if j > 2000:
                        data[i][j] = 0
                    if (data_d[i][j]!=0.0):
                        print("data_d[{}][{}]=".format(i,j),data_d[i][j])
            # if i < 5:
            #     if j > 2000:
                        data_d[i][j] = 0
        else:
            for i in range(5):
                for j in range(2000,5001,1):
                    if (data[i][j]!=0.0):
                        print("data[{}][{}]=".format(i,j),data[i][j])
            # if i < 5:
            #     if j > 2000:
                        data[i][j] = 0
#------------------------end-removing data------------------------------
        print ("data[49][2110]=",data[49][2110]) # aqui hay datos
        print ("data[50][2110]=",data[50][2110]) # aqui no hay datos???
        print("***********************")
        print (np.shape(data),"=(51,5001)?")
        print ("Frequencies length=",len(frequencies))
        print ("dibuja")
        plasmon.plot_contour(data)
        print("calcula ajuste lorentzian")
        #none = plasmon.locate_1Lorenztian(frequencies,data[30])
        #f=open('file_data_fittings.txt','a')
        #f.write("-amplitude--(+/-)error--center--(+/-)error--width--(+/-)error-\n")
        plasmon.fitting_Lorentz(frequencies,data, 51, big=True) #size = 51*2 = 102
        #f.close()
        if (len( sys.argv) == 3):
            np.savetxt("data_fitting_amplitudes_{}.txt".format(namefile),np.c_[plasmon.pars[:,0], plasmon.perr_lorentz[:,0]],header='#-----amplitude--(+/-)error---', footer='-------------')
            np.savetxt("data_fitting_center_{}.txt".format(namefile),np.c_[plasmon.pars[:,1], plasmon.perr_lorentz[:,1]],header='#-----center(e.V.)--(+/-)error---', footer='-------------')
            np.savetxt("data_fitting_width_{}.txt".format(namefile),np.c_[plasmon.pars[:,2], plasmon.perr_lorentz[:,2]],header='#-----width(e.V.)--(+/-)error---', footer='-------------')
        else:
            np.savetxt("data_fitting_amplitudes_{}.txt".format(namefile),np.c_[plasmon.pars[:,0], plasmon.perr_lorentz[:,0]],header='#-----amplitude--(+/-)error---for--qx={}'.format(qx), footer='-------------')
            np.savetxt("data_fitting_center_{}.txt".format(namefile),np.c_[plasmon.pars[:,1], plasmon.perr_lorentz[:,1]],header='#-----center(e.V.)--(+/-)error---for--qx={}'.format(qx), footer='-------------')
            np.savetxt("data_fitting_width_{}.txt".format(namefile),np.c_[plasmon.pars[:,2], plasmon.perr_lorentz[:,2]],header='#-----width(e.V.)--(+/-)error---for--qx={}'.format(qx), footer='-------------')


        #f.close()
        if (len( sys.argv) == 3):
            plasmon.fitting_Lorentz2(frequencies_d,data_d, 51, big=True) #size = 51*2 = 102
            np.savetxt("data_fitting_amplitudes_d_{}.txt".format(namefile),np.c_[plasmon.pars2[:,0], plasmon.perr_lorentz2[:,0]],header='#-----amplitude--(+/-)error---', footer='-------------')
            np.savetxt("data_fitting_center_d_{}.txt".format(namefile),np.c_[plasmon.pars2[:,1], plasmon.perr_lorentz2[:,1]],header='#-----center(e.V.)--(+/-)error---', footer='-------------')
            np.savetxt("data_fitting_width_d_{}.txt".format(namefile),np.c_[plasmon.pars2[:,2], plasmon.perr_lorentz2[:,2]],header='#-----width(e.V.)--(+/-)error---', footer='-------------')
        # else:
        #     np.savetxt("data_fitting_amplitudes_d_{}.txt".format(namefile),np.c_[plasmon.pars2[:,0], plasmon.perr_lorentz2[:,0]],header='#-----amplitude--(+/-)error---for--qx={}'.format(qx), footer='-------------')
        #     np.savetxt("data_fitting_center_d_{}.txt".format(namefile),np.c_[plasmon.pars2[:,1], plasmon.perr_lorentz2[:,1]],header='#-----center(e.V.)--(+/-)error---for--qx={}'.format(qx), footer='-------------')
        #     np.savetxt("data_fitting_width_d_{}.txt".format(namefile),np.c_[plasmon.pars2[:,2], plasmon.perr_lorentz2[:,2]],header='#-----width(e.V.)--(+/-)error---for--qx={}'.format(qx), footer='-------------')

        print ("dibuja")
        plasmon.plot_contour(plasmon.Fitted_data)
        plasmon.plot_contour2(data,plasmon.Fitted_data)
        if (len( sys.argv) == 3):
            print ("dibuja_doble")
            Fitted_data = np.vstack((np.flip(plasmon.Fitted_data, axis=0), plasmon.Fitted_data2))
            #plasmon.plot_contour(Fitted_data)
            #plasmon.plot_contour2(data,Fitted_data)
            plasmon.plot_contour3(np.vstack((np.flip(data, axis=0), data_d)),Fitted_data)

            # np.savetxt('Lambda.txt', (lambda_1))
            # np.savetxt('Lambda_from_a2F.txt', np.array((frequencies[1:],superconductor.lambda_2)).T, header='frequencies,Lambda')

        file1 = './pars_data.txt'
        file2 = './frequencies_data.txt'
        if (len( sys.argv) == 3):
            if not (os.path.isfile(file1) and os.path.isfile(file2)):
                data, frequencies = plasmon.load_big_file(0, arg[1],diagonal=False)
                plasmon.fitting_Lorentz(frequencies,data, 51, big=False) #TEST
                pars = plasmon.pars
                for index in range(0,51): #51 eller 50 -> 2601 eller 2550 : antes index in range(1,51)
                    data_t, frequencies_t = plasmon.load_big_file(index, arg[1],diagonal=False)
    #-------------------------fix to remove extra polaron data--------------
                    if (index<6):
                        for i in range(5):
                            for j in range(2000,5001,1):
                                if (data[i][j]!=0.0):
                                    print("data[{}][{}]=".format(i,j),data[i][j])
                                    data[i][j] = 0
    #------------------------end-removing data------------------------------
                    #data = np.vstack((data, data_t))
                    #frequencies = np.vstack((frequencies, frequencies_t))
                    plasmon.fitting_Lorentz(frequencies,data_t, 51, big=False) #TEST
                    pars = np.vstack((pars,plasmon.pars))
                #print ("frequencies=",frequencies)
                #print ("shape(data)=",np.shape(data),"=(51,5001)?")
                #print ("shape(frequencies)=",np.shape(frequencies),"=(51,5001)?")
                #plasmon.fitting_Lorentz(frequencies,data, 51)
                #plasmon.fitting_Lorentz(np.tile(frequencies,51),data, 2601)
            #    print ("shape(pars)=",np.shape(pars))
                #print (plasmon.pars2[:,0])
            #    print('Min(pars[0])=',np.amin(pars[:,0]))
                #print ('Min(pars2[0])=',np.amin(plasmon.pars2[:,0]))

                np.savetxt(file1, pars) #no negative values
                np.savetxt(file2, frequencies) #no negative values
            else:
                pars = np.loadtxt(file1) #no negative values,
                frequencies = np.loadtxt(file2) #no negative values
            superconductor = Eliashberg(pars)
            superconductor.read_Ne()
        #    print("Omega range:",np.min(superconductor.pars[:,1]-np.abs(np.max(superconductor.pars[:,2]))),'::',np.max(superconductor.pars[:,1])+np.abs(np.max(superconductor.pars[:,2])))
            lambda_1 = superconductor.Lambda_new(frequencies)
            print('Lambda_1=',lambda_1,'[]?') # Lambda calculated from Lambda_q
            print('Lambda_2=',superconductor.lambda_2,'[]?') #Lambda calculated fron Eliashberg function
            #print("shape(lambda_q)=",np.shape(superconductor.lambda_q_lista))
        #    print("lambda_q=",superconductor.lambda_q_lista[:50])

            fig_lambda_q = plt.figure(figsize=(10,6))
            ax = fig_lambda_q.add_subplot(1, 1, 1)
            ax.plot(superconductor.lambda_w_lista,superconductor.w[1:])
            ax.set_title('$\lambda$ vs. $\omega$')
            ax.set_xlabel('$\lambda(\omega)$')
            ax.set_ylabel('$\omega$')
            #ax.set_xticks([0,len(superconductor.w)])
            #ax.set_yticks([0,5001])
            #ax.set_xticklabels(["0","$\pi$"])
            #ax.set_xticklabels([superconductor.w[0],superconductor.w[-1]])
            #ax.set_yticklabels(["0","1"])
            plt.tight_layout()
            plt.show()
            fig_lambda_q.savefig("Ajuste_d_{}".format("lambda_w"))
            a2F_lista = []
            #print("frequencies=",frequencies)
            #for w in frequencies:
            #    a2F_lista.append(superconductor.a2F_new(w))
            a2F_lista = superconductor.a2F_new(frequencies)

            #print ("a2F(w)=",a2F_lista)
            fig_a2F = plt.figure(figsize=(10,6))
            ax = fig_a2F.add_subplot(1, 1, 1)
            ax.plot(a2F_lista,frequencies)
            ax.set_title('a2F vs. $\omega$')
            ax.set_ylabel('$\omega$ (eV)')
            ax.set_xlabel('a2F')

            plt.plot ()
            plt.tight_layout()
            plt.show()
            fig_a2F.savefig("Ajuste_d_{}".format("a2F"))

            plot_all(np.vstack((np.flip(data, axis=0), data_d)),
                np.vstack((np.flip(plasmon.Fitted_data, axis=0), plasmon.Fitted_data2)),
                plasmon.pars, plasmon.pars2,
                superconductor.lambda_w_lista,superconductor.w,
                a2F_lista,frequencies)
            T_c = superconductor.T_c(mu_par=0.1) #mu*=0.1 y mu*=0.15. Son los valores típicos.
            print("T_c=",T_c,"eV:: T_c=",T_c*11604,"K")
            print("T_c=",T_c,"K")
            np.savetxt("lambda_and_T_C.txt",(superconductor.lambda_2,T_c), header='Lambda, T_c (K)')
            #np.savetxt("lambda_and_T_C.txt",(superconductor.lambda_2,T_c,T_c*11604), header='Lambda, T_C (eV), T_c (K)')
        else:
            np.savetxt('./pars_data_{}.txt'.format(arg[1]), plasmon.pars)
    else:
        print ("Arguments are namefile and the index of q_x as second argument if you want the BIG FILE")
    pass

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
