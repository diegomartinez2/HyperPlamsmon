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
import numpy as np

# Leer los datos de los tres archivos
datos1 = np.loadtxt('a2F_1DP', delimiter=' ')
datos2 = np.loadtxt('a2F_HPI', delimiter=' ')
datos3 = np.loadtxt('a2F_HPII', delimiter=' ')

# Asegurarse de que todos los archivos tienen la misma longitud
assert len(datos1) == len(datos2) == len(datos3), "Los archivos deben tener la misma longitud"

# Combinar los datos en un solo array
datos_combinados = np.column_stack((datos1[:, 0], datos1[:, 1], datos2[:, 1], datos3[:, 1]))

# Guardar el array combinado en un nuevo archivo
np.savetxt('a2F_combinados.txt', datos_combinados, delimiter=' ', header='#w(eV)-------------------a2F(w)[1DP]--------------a2F(w)[HPI]---------------a2F(w)[HPII]')

#----repetimos-----------
# Leer los datos de los tres archivos
datos1 = np.loadtxt('Lambda_lista1DP', delimiter=' ')
datos2 = np.loadtxt('Lambda_listaHPI', delimiter=' ')
datos3 = np.loadtxt('Lambda_listaHPII', delimiter=' ')

# Asegurarse de que todos los archivos tienen la misma longitud
assert len(datos1) == len(datos2) == len(datos3), "Los archivos deben tener la misma longitud"

# Combinar los datos en un solo array
datos_combinados = np.column_stack((datos1[:, 0], datos1[:, 1], datos2[:, 1], datos3[:, 1]))

# Guardar el array combinado en un nuevo archivo
np.savetxt('Lambda_combinados.txt', datos_combinados, delimiter=' ', header='#w(eV)-------------------Lambda(w)[1DP]-----------Lambda(w)[HPI]-----------Lambda(w)[HPII]')
