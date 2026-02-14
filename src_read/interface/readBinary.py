# ======================================================================
# -*- name  : readBinary.py
# -*- coding: utf-8 -*-
# ======================================================================

# import
import sys
import numpy as np

# ======================================================================
# getBinary(f, max_len):
# ----------------------------------
# Only the original data is extracted in order from the beginning.
# Returns the extracted data and a list of headers and footers
# ======================================================================
def getBinary(f, max_len):
    index=0
    headerList=[]
    footerList=[]
    data=b''

    # Repeat until all bytes are read
    while index<max_len:
        
        # read header
        f.seek(index)
        [header]=np.frombuffer(f.read(4), '>i')
        index+=4
        headerList.append(header)
                          
        # Extract the contents
        f.seek(index)
        data+=f.read(header)
        index+=header
        
        # read foote
        f.seek(index)
        [header]=np.frombuffer(f.read(4), '>i')
        index+=4
        footerList.append(header)
    
    return data, headerList, footerList

# ======================================================================
# readBinary(binary, index, count, endian):
# ----------------------------------
# Reads the number of bytes (count) given from the read start position (index),
#    and returns the data decoded by the specified endian and the next read start
#    position (new_index = index + count).
# ======================================================================
def readBinary(binary, index, count, endian):
    
    # read
    data=np.frombuffer(binary[index:index+count], endian)
    new_index=index+count
    
    return data, new_index

# ======================================================================
# resize(b, s):
# ----------------------------------
# Convert a one-dimensional array (b) to a multidimensional array of size (s) (for Fortran)
# ======================================================================
def resize(b, s):
    a=np.reshape(b, s, order='F')

    return a

# ======================================================================
# running program
# ======================================================================
if __name__ == '__main__':
    
    # program description
    print('readBinary.py')
    print('--------------------')
    print(' > function [getBinary]  get binaryData , from datafile(fortran : header + Data + footer ... )')
    print(' > function [readBinary] read binaryData')
    print(' > function [resize]     resize data , by fortran rule')
    print('======================================================================')

    
