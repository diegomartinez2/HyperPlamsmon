
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
np.savetxt('a2F_combinados.txt', datos_combinados, delimiter=' ')

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
np.savetxt('Lambda_combinados.txt', datos_combinados, delimiter=' ')
