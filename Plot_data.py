#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Plot_data.py
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
import numpy as np
import matplotlib.pyplot as plt
import Excel_read
import Plasmon

def plot_excel_data(qx,qy,Omega,Gamma,file_HP):
    for index_q in range(0,901,50):
        mask = qx == index_q
        print("mask=",mask)
        plt.errorbar(qy[mask], Omega[mask], yerr=Gamma[mask], fmt='o', capsize=5, ecolor='red', markerfacecolor='blue')
        # Etiquetas
        plt.xlabel('qx')
        plt.ylabel('Omega')
        plt.title('Fitting')
        # Establecer el rango de los ejes
        plt.xlim(0, 1400) # Rango para el eje x
        plt.ylim(0, 0.5) # Rango para el eje y
        plt.savefig("Figura_Slava_fitting_{0}_{1}".format(index_q,file_HP))
        plt.show()

    exit()

def plot_plasmon(qx,qy,Omega,Gamma,data,mask_value=50, diagonal = True, index = 0, plot_error = False, save_fig = False):
    print ("Leyendo datos")

    for index in range(0,51):
        plasmon = Plasmon.Plasmon_analysis(index,"A7_EPS.dat")
        if diagonal:
            mask = qx == qy
            data, frequencies = plasmon.load_big_file(index, diagonal=True)
        else:
            mask = qx == mask_value
            data, frequencies = plasmon.load_big_file(index, diagonal=False)

        plt.style.use('_mpl-gallery-nogrid')
        fig, ax1 = plt.subplots(1,1)
        fig.set_size_inches(10, 5)
        fig.set_dpi(100)
        cax = ax1.imshow(data.T[:2500,:],
        #	vmin = 0.0 , vmax = 0.004,
        #	vmin = 0.0 , vmax = 0.2,
        	vmin = 0.0 , vmax = 0.1,
        	cmap=plt.colormaps['jet'], origin='lower',
        	interpolation='gaussian', aspect='auto')
        ax1.set_ylabel(r'Frequency (eV)', fontsize=12)
        ax1.tick_params(axis='both',direction='inout', length=8)
        ax1.set_xticks([0,10,20,25,30,40,50])
        ax1.set_yticks([0,500,1000,1500,2000,2500])
#        ax1.set_xticklabels(["0",5*np.pi/50,10*np.pi/50,15*np.pi/50,20*np.pi/50,np.pi/2,"$\pi$"])
        ax1.set_xticklabels(["0","$\pi/5$","$2\pi/5$","$\pi/2$","$3\pi/5$","$4\pi/5$","$\pi$"])
        ax1.set_yticklabels(["0","0.1","0.2","0.3","0.4","0.5"])
        if diagonal:
               ax1.title("$q_x=q_y")
        print (mask)

        if plot_error:
            cax3 = ax1.errorbar(x = qy[mask]/50, y = Omega[mask]*5001, yerr = Gamma[mask]*5001/2, fmt = 'o') #¿y*5001? para la escala
        cbar = fig.colorbar(cax)
        print("qx=",self.qx[1:10])
        plt.tight_layout()
        if save_fig:
            plt.savefig("Figura_e_{}".format(index))
        else:
            plt.show()
        pass

def plot_all(mask_value=50, diagonal = True, index = 1, plot_error = True, save_fig = False):

    print ("Leyendo datos")

    plasmon = Plasmon.Plasmon_analysis(index,"A7_EPS.dat")
    if diagonal:
        data, frequencies = plasmon.load_big_file(index, diagonal=True)
    else:
        data, frequencies = plasmon.load_big_file(index, diagonal=False)

    plt.style.use('_mpl-gallery-nogrid')
    fig, ax1 = plt.subplots(1,1)
    fig.set_size_inches(10, 5)
    fig.set_dpi(100)
    cax = ax1.imshow(data.T[:2500,:],
    #	vmin = 0.0 , vmax = 0.004,
    #	vmin = 0.0 , vmax = 0.2,
    	vmin = 0.0 , vmax = 0.1,
    	cmap=plt.colormaps['gist_stern'], origin='lower', #plt.colormaps['jet']; plt.colormaps['gnuplot']
    	interpolation='gaussian', aspect='auto')
    ax1.set_ylabel(r'Frequency (eV)', fontsize=12)
    ax1.set_xlabel(r'$q_x$', fontsize=12)
    #ax1.set_xticks([0,25,50])
    ax1.tick_params(axis='both',direction='inout', length=8)
    ax1.set_yticks([0,500,1000,1500,2000,2500])
    #ax1.set_xticklabels(["0","\frac{\pi}{2}","$\pi$"])
    ax1.set_xticks([0,10,20,25,30,40,50])
    ax1.set_xticklabels(["0","$\pi/5$","$2\pi/5$","$\pi/2$","$3\pi/5$","$4\pi/5$","$\pi$"])
    ax1.set_yticklabels(["0","0.1","0.2","0.3","0.4","0.5"])

    if diagonal:
        if plot_error:
            qx,qy,Omega,Gamma,Ratio = Excel_read.Excel_data(filename="1DP_c3")
            mask = qx == qy
            cax3 = ax1.errorbar(x = qy[mask]/50 -1, y = Omega[mask]*5001, yerr = Gamma[mask]*5001/2, fmt = 'o') #¿y*5001? para la escala
            qx,qy,Omega,Gamma,Ratio = Excel_read.Excel_data(filename="HPI_c2")
            mask = qx == qy
            cax2 = ax1.errorbar(x = qy[mask]/50 -1, y = Omega[mask]*5001, yerr = Gamma[mask]*5001/2, fmt = 's') #¿y*5001? para la escala
            qx,qy,Omega,Gamma,Ratio = Excel_read.Excel_data(filename="HPII_c2")
            mask = qx == qy
            cax1 = ax1.errorbar(x = qy[mask]/50 -1, y = Omega[mask]*5001, yerr = Gamma[mask]*5001/2, fmt = 'v') #¿y*5001? para la escala
    else:
        if plot_error:
            qx,qy,Omega,Gamma,Ratio = Excel_read.Excel_data(filename="1DP_c3")
            mask = qx == mask_value
            cax3 = ax1.errorbar(x = qy[mask]/50 -1, y = Omega[mask]*5001, yerr = Gamma[mask]*5001/2, fmt = 'o') #¿y*5001? para la escala; -1 because the array in python begin in '0' and the original data do not have data in qx=0
            qx,qy,Omega,Gamma,Ratio = Excel_read.Excel_data(filename="HPI_c2")
            mask = qx == mask_value
            cax2 = ax1.errorbar(x = qy[mask]/50 -1, y = Omega[mask]*5001, yerr = Gamma[mask]*5001/2, fmt = 's') #¿y*5001? para la escala
            qx,qy,Omega,Gamma,Ratio = Excel_read.Excel_data(filename="HPII_c2")
            mask = qx == mask_value
            cax1 = ax1.errorbar(x = qy[mask]/50 -1, y = Omega[mask]*5001, yerr = Gamma[mask]*5001/2, fmt = 'v') #¿y*5001? para la escala


    cbar = fig.colorbar(cax)
    #print("qx=",self.qx[1:10])
    #plt.tight_layout()
    if diagonal:
            plt.title("$q_x=q_y$")
    else:
            if index==0:
                 plt.title("$q_y=0$")
            else:
                 plt.title("$q_y=$"+str(index/50)+"*$\pi$")
    plt.tight_layout()
    if save_fig:
        if diagonal:
            plt.savefig("Figura_diagonal_{}".format(index))
        else:
            plt.savefig("Figura_q_{}".format(index))
    else:
        plt.show()
    pass

plot_all(mask_value=50, diagonal = True, index = 1, plot_error = True, save_fig = False)
pass
