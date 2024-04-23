#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Calculate_T_c.py
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
import Excel_read
import Eliashberg
import numpy as np
import matplotlib.pyplot as plt

# ----------
# Funciones
# ----------

def plotter(w, Lambda_1DP, Lambda_HPI, Lambda_HPII, Lambda_Total, a2F_1DP, a2F_HPI, a2F_HPII, a2F_Total):
    plt.plot(w, (Lambda_1DP+Lambda_HPI+Lambda_HPII), label="Summa lambdas")
    plt.plot(w, Lambda_Total, label="Lambda Total")
    plt.legend()
    plt.show()

def plotter2(w, Lambda_1DP, Lambda_HPI, Lambda_HPII, Lambda_Total, a2F_1DP, a2F_HPI, a2F_HPII, a2F_Total):
    plt.plot(w, (a2F_1DP+a2F_HPI+a2F_HPII), label="Summa a2F")
    plt.plot(w, a2F_Total, label="a2F Total")
    plt.legend()
    plt.show()

def main(args):
    """
    Calculates the critical superconductivity temperature using the Eliashberg function and the Lambda.
    """
    w_a, a2F_1DP, a2F_HPI, a2F_HPII = np.loadtxt("a2F_combinados.txt", usecols=(0, 1, 2, 3), unpack=True)
    w, Lambda_1DP, Lambda_HPI, Lambda_HPII = np.loadtxt("Lambda_combinados.txt", usecols=(0, 1, 2, 3), unpack=True)
    a2F_Total = a2F_1DP+a2F_HPI+a2F_HPII
    Lambda_Total = Lambda_1DP+Lambda_HPI+Lambda_HPII
    plotter(w, Lambda_1DP, Lambda_HPI, Lambda_HPII, Lambda_Total, a2F_1DP, a2F_HPI, a2F_HPII, a2F_Total)
    # Lambda = Lambda_1DP+Lambda_HPI+Lambda_HPII
    Lambda = Lambda_Total
    a2F = a2F_Total[1:]
    mask = Lambda > 0
    Tc = Eliashberg.T_c(w[mask],0.1,Lambda[-1],a2F[mask])
    print("T_c=",Tc,"(K) sumando 1DP+HPI+HPII")


    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
