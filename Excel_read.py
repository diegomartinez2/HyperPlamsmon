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
import os
import glob
import xlrd

# -------
# Clases
# -------


# ----------
# Funciones
# ----------
def excel_read(file_directory=".",header="1DP"):
    """
    Function that reads the data from the excel files in a directory.
    Header options are "1DP", "HPI", "HPII"
    """
    if (file_directory=="."):
        current_directory = os.getcwd()
    else:
        current_directory = file_directory
    # Find all .xlsx files in the current directory
    excel_files = glob.glob(os.path.join(current_directory, header+"*.xlsx"))
    data=[]
    # Iterate over each Excel file
    for file_name in excel_files:
        # Open the workbook
        workbook = xlrd.open_workbook(file_name)

        # Open the worksheet (assuming the first sheet is the one to be processed)
        worksheet = workbook.sheet_by_index(0)

        out=np.zeros((worksheet.nrows,worksheet.ncols))
        # Iterate the rows and columns
        for i in range(worksheet.nrows):
            row_values = []
            for j in range(worksheet.ncols):
                # Check if the cell is empty and fill with zero if true
                cell_value = worksheet.cell_value(i, j) if worksheet.cell_type(i, j) != xlrd.XL_CELL_EMPTY else  0
                # Add the cell value to the output array
                out[i][j] = cell_value

        data.append(out)
    #return qx,qy,out
    return data
def read_1_excel_file(file_directory=".",filename="1DP"):#filenames=('1DP','HPI','HPII')
    """
    This function only reads one excel file in qx,qy,w,e,r... format and outputs the data inside
    out=[[qx,qy,w,e,r,w2,e2,r2,w3,e3,r3],...]
    """
    if (file_directory=="."):
        current_directory = os.getcwd()
    else:
        current_directory = file_directory
    #excel_file = glob.glob(os.path.join(current_directory, filename+".xlsx"))
    # Open the workbook
    #workbook = xlrd.open_workbook(excel_file)
    workbook = xlrd.open_workbook(os.path.join(current_directory, filename+".xlsx"))
    # Open the worksheet (assuming the first sheet is the one to be processed)
    worksheet = workbook.sheet_by_index(0)

    out=np.zeros((worksheet.nrows,worksheet.ncols))
    # Iterate the rows and columns
    for i in range(worksheet.nrows):
        row_values = []
        for j in range(worksheet.ncols):
            # Check if the cell is empty and fill with zero if true
            cell_value = worksheet.cell_value(i, j) if worksheet.cell_type(i, j) != xlrd.XL_CELL_EMPTY else  0
            # Add the cell value to the output array
            out[i][j] = cell_value
    return out
    pass
def Excel_data_parser(data):
    qx=np.zeros(len(data))
    qy=np.zeros(len(data))
    Omega=np.zeros(len(data))
    Gamma=np.zeros(len(data))
    Ratio=np.zeros(len(data))
    for i in range(len(data)):
            qx[i]=data[i][0]
            qy[i]=data[i][1]
            Omega[i]=data[i][2]
            Gamma[i]=data[i][3]
            Ratio[i]=data[i][4]
    return qx,qy,Omega,Gamma,Ratio
    pass

def Excel_data(filename,flag_1DP = False):
    """
    Read the data for the plasmons in cuprate from the excel files.
    If the file is for the 1DP...
    """
    out = Eliashberg.read_1_excel_file(filename=filename) #filenames=('1DP','HPI','HPII');filenames=('1DP_c','HPI_c','HPII'_c)
    #qx,qy,Omega,Gamma,Ratio = Eliashberg.Excel_data_parser(out)
    qx=out[:,0]
    qy=out[:,1]
    Omega=out[:,2]
    Gamma=out[:,3]
    Ratio=out[:,4]
    if flag_1DP: #make true for "1DP" or "1DP_c"
        qx=np.append(qx,out[:,0])
        qy=np.append(qy,out[:,1])
        Omega=np.append(Omega,out[:,5])
        Gamma=np.append(Gamma,out[:,6])
        Ratio=np.append(Ratio,out[:,7])
        qx=np.append(qx,out[:,0])
        qy=np.append(qy,out[:,1])
        Omaga=np.append(Omega,out[:,8])
        Gamma=np.append(Gamma,out[:,9])
        Ratio=np.append(Ratio,out[:,10])
    return qx,qy,Omega/1000,Gamma/1000,Ratio

def main(args):
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
