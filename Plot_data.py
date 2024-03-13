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
