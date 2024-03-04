#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  untitled.py
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
class Eliashberg(object):
    """docstring for NombredeClase."""

    def __init__(self, arg):
        super(NombredeClase, self).__init__()
        self.arg = arg
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
        self.Gamma =  Gamma
        self.Ratio = Ratio

        self.indice_zeros = 0
        self.N_qs = ((np.max(qx)*np.max(qy))+np.max([np.max(qx),np.max(qy)]))/(2*50)
        self.N = 50*50 #test for 3 and 5 (maybe *4 for the 4 cuadrants and *2 for the HPII AND HPI :: and  the 1DP has 3 plasmons)
        print("factor N_qs=",self.N_qs,"::",self.N,"::",len(Omega))

    def a2F(self,x,method = 0):
        """
        Calculates the Eliashberg spectral functions
        ---input---
        x: coordiate to calculate the Eliashberg function.
        method: Method to use in the calculation.
        ---output---
        a2F = factor1*summa: The Eliashberg function at "x"
        """

        center = self.Omega[:] #*put units correctly...
        width = self.Gamma[:] #*put units the same as center
        gauss_width = 5*self.from_cm1_to_eV#(units.invcm/units.Hartree) #0.00002 # test the units of this... should be aprox. 5 cm-1 (1, 5 or 10)
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
                summa += self.Lambda_q_new(i)*self.Omega[i] * self.gaussian(x,self.Omega[i],gauss_width) * simmetry_factor
            #return summa/(2*(len(self.Omega)-self.indice_zeros))
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
            #return summa/(2*np.pi*self.N_ef*(len(self.Omega)-self.indice_zeros))
            #return summa/(2*np.pi*self.N_ef*len(self.Omega))
            #return summa/(2*np.pi*self.N_ef*self.N_qs)
            return summa/(2*np.pi*self.N_ef*self.N)

    def Lambda_q(self,i):
        """
        Calculates the Lambda(q) functions
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
        #     Lamb_q=0
        #     self.indice_zeros += 1
        #     print("Indice_zeros=",self.indice_zeros)
        return Lamb_q

    def read_Ne(self,filename="out_DOS.dat"):
        """
        Read the number of particles under the Fermi level (???)
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

            summa1 += self.Lambda_q_new(i)*symmetry_factor
            self.lambda_q_lista = np.append(self.lambda_q_lista, self.Lambda_q_new(i)*symmetry_factor)

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
        res = np.absolute(self.a2F_new(self.w))
        a2F_x = np.absolute(np.divide(res, self.w))
        #res = np.absolute(self.a2F_new(w))
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
            #func_w = np.divide(self.a2F_new(w_1), w_1)
            func_w = np.absolute(self.a2F_new(w_1)/w_1)
            #partial_value=2*integrate.simpson(func_w,w_1)
            partial_value = 2* np.trapz(func_w,w_1)
            self.lambda_w_lista[i] = partial_value #test to see the evolution of Lambda
        diff_lambda_w = np.diff(self.lambda_w_lista)
        print("lambda(w) always growing?=",np.all(diff_lambda_w > 0)) # this test fails!!!
        self.plot_lambda(self.lambda_w_lista,self.w_0)
        return Lambda_1
# ----------
# Funciones
# ----------
def NombredeFuncion(arg):
    pass



def main(args):
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
