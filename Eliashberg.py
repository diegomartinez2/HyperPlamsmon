#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Eliashberg.py
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
from scipy import integrate
import scipy as scipy
from ase.units import create_units
import os.path
import os

import Excel_read
# -------
# Clases
# -------
class Eliashberg(object):
    """docstring for Eliashberg.
    This object calculates the Lambda from two methods:
    From the Eliashberg function amd fron the lambda at the q points.
    For more info: DOI: 10.1103/PhysRevB.103.134305
    """

    def __init__(self, qx,qy,Omega,Gamma,Ratio):    ###, arg):
        super(Eliashberg, self).__init__()
        #self.arg = arg
        """
        Set parameters into object and transformation units.
        """
        units = create_units('2014')   #new way of defining units
        self.from_cm1_to_Hartree = units.invcm/units.Hartree #0.0000045563352812122295
        self.from_GHz_to_Hartree = self.from_cm1_to_Hartree /29.9792458 # from GHz to Hartree
        self.from_Ry_to_Hartree = units.Ry/units.Hartree #0.5 # from Ry to Hartree
        """
        cm^-1 to Hartree: Hartree = cm^-1 / 219474.6
        Hartree to cm^-1: cm^-1 = Hartree * 219474.6

        GHz to Hartree: Hartree = GHz / 6.5796839 × 10^9
        Hartree to GHz: GHz = Hartree * 6.5796839 × 10^9

        eV to Hartree: Hartree = eV*0,0367493
        """
        self.from_cm1_to_eV = units.invcm/units.eV #0.00012398425731484318
        self.from_GHz_to_eV = 0.000004135669661004052
        self.qx = qx
        self.qy = qy
        self.Omega = Omega
        self.Gamma =  Gamma     #los datos proporcionados son del full widths
        self.Ratio = Ratio

        self.indice_zeros = 0
        self.N_qs = ((np.max(qx)*np.max(qy))+np.max([np.max(qx),np.max(qy)]))/(2*50)
        self.N = 50*50 #test for 3 and 5 (maybe *4 for the 4 cuadrants and *2 for the HPII AND HPI :: and  the 1DP has 3 plasmons)
        print("factor N_qs=",self.N_qs,"::",self.N,"::",len(Omega))

    def a2F(self,x,method = 0):
        """
        Calculates the Eliashberg spectral functions
        a2F(w)=1/(2*pi*N(E_F)*Nq)Sum{(gamma_q/omega_q)*delta(omega-omega_q)}
        ---input---
        x: coordiate to calculate the Eliashberg function.
        method: Method to use in the calculation.
        ---output---
        a2F = factor1*summa: The Eliashberg function at "x"
        """

        center = self.Omega[:] #*put units correctly...
        width = self.Gamma[:] #*put units the same as center
        gauss_width = 5*self.from_cm1_to_eV#(units.invcm/units.Hartree) should be aprox. 5 cm-1 (1, 5 or 10)
        if (method == 1):
            #---------method1-------vvvv---
            summa = 0
            for i in range(len(self.qx)):
                simmetry_factor =  2 #8
                if self.qx[i] ==  0:
                    simmetry_factor =  1 #4
                    if self.qy[i] ==  0:
                        simmetry_factor =  0 #1
                elif self.qy[i] ==  0:
                    simmetry_factor =  1 #4
                elif self.qx[i] == self.qy[i]:
                    simmetry_factor =  1 #4
                summa += self.Lambda_q(i)*self.Omega[i] * self.gaussian(x,self.Omega[i],gauss_width) * simmetry_factor
            #return summa/(2*(len(self.Omega)))
            return summa/(2*self.N)
        else:
            #---------method2------vvvv-
            summa = 0

            for i in range(len(self.qx)):
                simmetry_factor =  2 #8
                if self.qx[i] ==  0:
                    simmetry_factor =  1 #4
                    if self.qy[i] ==  0:
                        simmetry_factor = 0 #1
                elif self.qy[i] ==  0:
                    simmetry_factor =  1 #4
                elif self.qx[i] == self.qy[i]:
                    simmetry_factor =  1 #4

                summa += (self.Ratio[i]) * self.gaussian(x,self.Omega[i],gauss_width) * simmetry_factor
            #return summa/(2*np.pi*self.N_ef*len(self.Omega))
            #return summa/(2*np.pi*self.N_ef*self.N_qs)
            return summa/(2*np.pi*self.N_ef*self.N)

    def Lambda_q(self,i):
        """
        Calculates the Lambda(q) functions
        Lambda_q=1/(pi*N(E_F))*gamma_q/omega_q²
        ---input---
        i: index of the frequencies and widths of the lorentzian fitting of the plasmon
        ---output---
        Lamb_q: Lambda(q)
        """
        if (self.Omega[i]!=0):
            Lamb_q=(1/(np.pi*self.N_ef)) * (self.Ratio[i])/(self.Omega[i]) #fix from omega to omega²
        else:
            print("error, Omega<0")
            exit()
        return Lamb_q

    def read_Ne(self,filename="out_DOS.dat"):
        """
        Read the number of particles under the Fermi level (???) ->'N(E_F)
        Set the numbre of elements.
        """
        self.energy, self.Ne = np.loadtxt(filename,usecols=(0,1), unpack=True)
        numero_de_elementos = np.trapz(self.Ne[:(np.where(self.energy==0.0)[0][0])],dx=np.absolute(self.energy[0]-self.energy[-1])/len(self.energy)) # maybe this is better??
        numero_de_elementos2 = integrate.simpson(self.Ne[:(np.where(self.energy==0.0)[0][0])],dx = np.absolute(self.energy[0]-self.energy[-1])/len(self.energy))
        total_states = np.trapz(self.Ne,dx=np.absolute(self.energy[0]-self.energy[-1])/len(self.energy))
        print("Número de elementos:",numero_de_elementos,"|total=",total_states)
        print("Número de elementos2:",numero_de_elementos2,"|total=",total_states)
        print("N(eF)=",self.Ne[(np.where(self.energy==0.0)[0][0])]) #states/spin/eV/unit cell

        #------^--------
        print ("Density of states at Fermi level per cell=",self.Ne[np.where(self.energy==0.0)])
        print ("Number of elements (for the factor)=",numero_de_elementos)
        self.N_ef = self.Ne[(np.where(self.energy==0.0)[0][0])] #1.8855775
        pass

    def Lambda(self,Frequencies):
        """
        Calculates the Lambda by two methods, notice that it must calculate the integral
        in a range that takes the Lorenztian obtained by the Plasmon_analysis object.
        Lambda_1=1/N * Sum{Lambda_q}
        Lambda_2=2*Int_0^Inf{a²F(w)/w dw}  <---DOI: 10.1103/PhysRevB.103.134305
        """
        center = self.Omega
        width = self.Gamma
        #method 1 ----------------------------

        summa1 = 0
        self.lambda_q_lista = np.array([])
        for i in range(len(center)):  #summa in q (6x6x6 or q_x*q_y)
            symmetry_factor =  2 #8
            if self.qx[i] ==  0:
                symmetry_factor =  1 #4
                if self.qy[i] ==  0:
                    symmetry_factor =  0 #1 #test taking out the qx=qy=0
            elif self.qy[i] ==  0:
                symmetry_factor =  1 #4
            elif self.qx[i] == self.qy[i]:
                symmetry_factor =  1 #4

            summa1 += self.Lambda_q(i)*symmetry_factor
            self.lambda_q_lista = np.append(self.lambda_q_lista, self.Lambda_q(i)*symmetry_factor)

        #Lambda_1=summa1/(len(center)-self.indice_zeros)
        Lambda_1=summa1/self.N
        #print("Test Lambda with lambda_q=",summa1/(len(center)),":test lambda_q[N]=",summa1/self.N,":lambda_q[N_qs]=",summa1/self.N_qs)
        #method 2 -------------------------------
        self.lambda_2=[]
        mask = Frequencies >  0
        self.w = Frequencies[mask]
        # w = Frequencies[Frequencies != 0]
        # mask = w >  0
        # self.w = w[mask]
        res = np.absolute(self.a2F(self.w))
        a2F_x = np.absolute(np.divide(res, self.w))
        #res = np.absolute(self.a2F(w))
        #a2F_x = np.absolute(np.divide(res, w))
        #self.lambda_2 = 2*integrate.simpson(a2F_x,w)
        self.lambda_2 = 2*np.trapz(a2F_x, self.w)
        print("plot a2F(w)")
        self.plot_lambda(res,self.w)
        print("plot Lambda(w)")
        #self.lambda_w_lista = []
        #-------removin negatives w's---
        #mask = w >=  0
        #mask = self.w >  0
        self.w_0 = self.w #w[mask]
        self.lambda_w_lista = np.zeros(len(self.w_0))
        diff_w_0 = np.diff(self.w_0)
        print("w_0 always growing?=",np.all(diff_w_0 > 0)) # this test is OK
        for i in range(1,len(self.w_0)): #Frequencies[Frequencies != 0]:
            w_1 = self.w_0[:i]
            #func_w = np.divide(self.a2F(w_1), w_1)
            func_w = np.absolute(self.a2F(w_1)/w_1)
            #partial_value=2*integrate.simpson(func_w,w_1)
            partial_value = 2* np.trapz(func_w,w_1)
            self.lambda_w_lista[i] = partial_value #test to see the evolution of Lambda
        diff_lambda_w = np.diff(self.lambda_w_lista)
        print("lambda(w) always growing?=",np.all(diff_lambda_w > 0)) # this test fails!!!
        self.plot_lambda(self.lambda_w_lista,self.w_0)
        return Lambda_1

    def T_c(self,mu_par,lambda_t):
        """
        Allen-Dynes formula (Advanced McMillan's equation and its application for the analysis of highly-compressed superconductors)
        Tc=W_log/1.2*exp(-1.04*(1+self.lambda_2)/(self.lambda_2-mu_*(1+0.62*self.lambda_2)))
        Note this is an approximation formula and is not accurate for some superconductors
        formula from DOI 10.1007/s10948-017-4295-y
        """
        AllenDynes = True
        out2=-1.04*(1+lambda_t)/(lambda_t-mu_par*(1+0.62*lambda_t))
        if (AllenDynes):
            print ("f1,f2=",self.f_1(mu_par,lambda_t),self.f_2(mu_par,lambda_t))
            out=self.f_1(mu_par,lambda_t)*self.f_2(mu_par,lambda_t)*self.w_log(lambda_t)/1.2 * np.exp(out2)
        else:
            out=self.w_log(lambda_t)/1.2 * np.exp(out2)

        return out/8.617333262e-5 #Boltzman constant in eV/K

    def w_log(self,lambda_t):
        """
        w_log calculated from Eliashberg
        """
        #eV_to_K = 11604.5250061657

        #mask = self.w >=  0
        #mask = self.w >  0
        w = self.w#[mask]#*eV_to_K
        print("max(w)=",max(w))
        w_log = np.exp(
            (2.0/lambda_t)*
            #integrate.simpson(
            #(np.divide(self.a2F(self.w), self.w)*np.log(self.w)),self.w))
            np.trapz((np.divide(self.a2F(w), w)*np.log(w)),w)
            )
        return w_log

    def w_2(self,lambda_t):
        #eV_to_K=11604
        w = self.w#*eV_to_K
        w_2 = np.sqrt(
            (2.0/lambda_t)*
            #integrate.simpson(
            np.trapz((self.a2F(w)*w),w)
            )
        return w_2

    def f_1(self,mu_par,lambda_t):
        LAMBDA_temp = 2.46*(1+3.8*mu_par)
        return np.power(1+np.power(lambda_t/LAMBDA_temp,3/2),1/3)

    def f_2(self,mu_par,lambda_t):
        LAMBDA_temp =1.82*(1+6.3*mu_par)*(self.w_2(lambda_t)/self.w_log(lambda_t))
        print("mu_par=",mu_par)#OK
        print ("LAMBDA_temp=",LAMBDA_temp) #ERROR
        print ("w_2(lambda)=",self.w_2(lambda_t)) #ERROR
        print ("w_log(lambda)=",self.w_log(lambda_t))
        print ("return f2=",1 + (( self.w_2(lambda_t)/self.w_log(lambda_t) - 1) * lambda_t**2)/(
        (lambda_t**2) + (LAMBDA_temp**2)))
        return 1 + (( self.w_2(lambda_t)/self.w_log(lambda_t) - 1) * lambda_t**2)/(
        (lambda_t**2) + (LAMBDA_temp**2))


    def gaussian(self,x, center,gauss_width):
        """
        Gaussian distribution
        ---input---
        x: x position of the gaussian graph, in this case the frequency
        center: center of the gaussian distribution, in this case the frequency of the plasmon Lorenztian approximation
        gauss_width: width of the gaussian distribution, in this case the witdh of the plasmon lorentzian aproximation
        ---output---
        g: gaussian value in 'x'
        """
        mu = center #self.pars[:1] #center
        sigma = gauss_width #self.pars[:,2] #width
        #d = 1/(2*np.sqrt(2*np.pi)*sigma)
        d = 1/(sigma*np.sqrt(2*np.pi))
        g = d*np.exp(-(x-mu)**2 / ( 2.0 * sigma**2 )  )
        #exp(-(w_aux-w(j,l))**2.0d0/(2.0d0*broad**2.0d0))/ (broad*sqrt(twopi)) #copy from qe-5.1.0_elph
        return g

    def test_gaussian(self,w,w_q,N_q):
        """
        Test for the gaussian, the result must be the DOS.
        ---input---
        w: frequencies to test
        w_q: frequencies of the plasmon grid
        N_q: number of terms in w_q
        ---output---
        test_dos:
        """
        width = 5#*self.from_cm1_to_eV #self.from_cm1_to_Hartree*2 #1,5,10 cm-1
        test_dos = 1/N_q # is this equal to len(W_q)?
        for i in range(len(w_q)): #q=6x6x6 (???) for the test data. q=50x50=250 for the plasmon
            test_dos += self.gaussian(w,w_q[i],width) #This must be the DOS, integrate this to check the gaussian.
        print("Test of the gaussian: this must be the integral of the DOS=",test_dos)
        return test_dos

    def plot_lambda(self,x,w):
        plt.figure(figsize=(10,6))
        #plt.figure().set_figwidth(15)
        #plt.figure().set_figheight(2)
        plt.plot(x,w)
        plt.tight_layout()
        plt.show()
        pass
# ----------
# Funciones
# ----------

"""
Alternate functions
"""
def T_c(w,mu_par,lambda_t,a2F):
        """
        Allen-Dynes formula (Advanced McMillan's equation and its application for the analysis of highly-compressed superconductors)
        Tc=W_log/1.2*exp(-1.04*(1+self.lambda_2)/(self.lambda_2-mu_*(1+0.62*self.lambda_2)))
        Note this is an approximation formula and is not accurate for some superconductors
        formula from DOI 10.1007/s10948-017-4295-y
        """
        AllenDynes = True
        out2=-1.04*(1+lambda_t)/(lambda_t-mu_par*(1+0.62*lambda_t))
        if (AllenDynes):
            print ("f1,f2=",f_1(mu_par,lambda_t),f_2(mu_par,w,lambda_t,a2F))
            out=f_1(mu_par,lambda_t)*f_2(mu_par,w,lambda_t,a2F)*w_log(w,lambda_t,a2F)/1.2 * np.exp(out2)
        else:
            out=w_log(w,lambda_t,a2F)/1.2 * np.exp(out2)

        return out/8.617333262e-5 #Boltzman constant in eV/K

def w_log(w,lambda_t,a2F):
        """
        w_log calculated from Eliashberg
        """
        #eV_to_K=11604
        #eV_to_K = 11604.5250061657

        #mask = self.w >=  0
        #mask = self.w >  0
        #w = self.w#[mask]#*eV_to_K
        #print("max(w)=",max(w))
        for i in range(len(w)):
          if w[i]<0:
            print("error en w<0",i)
            exit()
        a2F_w = a2F/w #np.divide(a2F, w)
        log_w = np.log(w)
        w_log = np.exp((2.0/lambda_t)*np.trapz(a2F_w*log_w,w))
            #integrate.simpson((np.divide(self.a2F(self.w), self.w)*np.log(self.w)),self.w))
            #np.trapz(a2F_w*np.log(w)),w)
        #print(w_log,"(eV)")
        return w_log

def w_2(w,lambda_t,a2F):
        #eV_to_K=11604
        #w = w#*eV_to_K
        a2F_w = a2F*w
        w_2 = np.sqrt((2.0/lambda_t)*
            #integrate.simpson(
            np.trapz(a2F_w,w))
        return w_2

def f_1(mu_par,lambda_t):
        LAMBDA_temp = 2.46*(1+3.8*mu_par)
        return np.power(1+np.power(lambda_t/LAMBDA_temp,3/2),1/3)

def f_2(mu_par,w,lambda_t,a2F):
        LAMBDA_temp =1.82*(1+6.3*mu_par)*(w_2(w,lambda_t,a2F)/w_log(w,lambda_t,a2F))
        return 1 + (( w_2(w,lambda_t,a2F)/w_log(w,lambda_t,a2F) - 1) * lambda_t**2)/(
        (lambda_t**2) + (LAMBDA_temp**2))



def main(args):
    """
    -----
    Flags
    -----
    fix_frequencies: Set frequencies to a set range to compare the same outputs. If false it set the range to a optimized one.
    plot_excel_data: Plots the data.
    """
    fix_frequencies = False
    plot_excel_data = False
    if len(args)<2:
        print("Arguments are 1DP, HPI or HPII:")
        exit()
    if args[1]=="1DP":
        file_HP = "1DP"
        qx,qy,Omega,Gamma,Ratio = Excel_read.Excel_data(filename="{}_c3".format(file_HP),flag_1DP=False)
    elif args[1]=="HPI":
        file_HP = "HPI"
        qx,qy,Omega,Gamma,Ratio = Excel_read.Excel_data(filename="{}_c2".format(file_HP))
    elif args[1]=="HPII":
        file_HP = "HPII"
        qx,qy,Omega,Gamma,Ratio = Excel_read.Excel_data(filename="{}_c2".format(file_HP))
    else:
        print("Arguments are 1DP, HPI or HPII:",args[1])
        exit()
    if len(args)==3:
        if args[2]=="fix":
            fix_frequencies = True

    if plot_excel_data:
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
    superconductor = Eliashberg(qx,qy,Omega,Gamma,Ratio)
    superconductor.read_Ne()

    print("Omega range:",np.min(superconductor.Omega-np.abs(np.max(superconductor.Gamma))),'::',np.max(superconductor.Omega+np.abs(np.max(superconductor.Gamma))))
    if fix_frequencies:
        frequencies = np.linspace(0,0.5,10000) #test, range of frequencies
    else:
        frequencies = np.linspace(0,np.max(superconductor.Omega+np.abs(np.max(superconductor.Gamma))),10000)

    lambda_1 = superconductor.Lambda(frequencies)
    print('Lambda_1[{}]='.format(file_HP),lambda_1,'[Lambda calculated from Lambda_q for {}]?'.format(file_HP)) # Lambda calculated from Lambda_q
    print('Lambda_2[{}]='.format(file_HP),superconductor.lambda_2,'[Lambda calculated fron Eliashberg function for {}]?'.format(file_HP)) #Lambda calculated fron Eliashberg function
    np.savetxt('salida_{}.txt'.format(file_HP),np.c_[superconductor.qx,superconductor.qy,np.full(len(superconductor.Omega), superconductor.N_ef),superconductor.Omega,superconductor.Gamma,superconductor.lambda_q_lista],header='#---q_x---q_y---N_ef[States/Spin/eV/Unit Cell]---Omega[eV]---Gamma[eV]---Lambda_q')
#-----plot-start
    fig_lambda_q = plt.figure(figsize=(10,6))
    ax = fig_lambda_q.add_subplot(1, 1, 1)
    ax.plot(superconductor.lambda_w_lista,superconductor.w_0)
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
    fig_lambda_q.savefig("{0}_Ajuste_d_{1}".format(args[1],"lambda_w"))
    a2F_lista = []
    a2F_lista = superconductor.a2F(frequencies)

    fig_a2F = plt.figure(figsize=(10,6))
    ax = fig_a2F.add_subplot(1, 1, 1)
    ax.plot(a2F_lista,frequencies)
    ax.set_title('a2F vs. $\omega$')
    ax.set_ylabel('$\omega$ (eV)')
    ax.set_xlabel('a2F')

    plt.plot ()
    plt.tight_layout()
    plt.show()
    fig_a2F.savefig("{0}_Ajuste_d_{1}".format(args[1],"a2F"))
#---plot-end
    np.savetxt("Lambda_lista{}".format(file_HP),np.vstack((superconductor.w_0,superconductor.lambda_w_lista)).T)
    np.savetxt("a2F_{}".format(file_HP),np.vstack((frequencies, a2F_lista)).T)
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
