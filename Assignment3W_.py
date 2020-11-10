#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 09:26:14 2020

@author: maritutheim
"""

import matplotlib.pyplot as plt
import h5py as hp 
import numpy as np
import scipy.signal as sp

H1 = hp.File("H-H1_LOSC_4_V2-1126259446-32.hdf5", "r")
H1data = H1["strain/Strain"]

L1 = hp.File("L-L1_LOSC_4_V2-1126259446-32.hdf5", "r")
L1data = L1["strain/Strain"]

samples = len(H1data) # same for both H1 and L1

samplerate = 4096 # in Hz - same value for L1 and H1

time = samples/samplerate # in seconds - same for both H1 and L1

time_vec = np.linspace(0,time,samples) # time vector


#Hanford

H1_FT = np.fft.rfft(H1data) # FFT of Hanford signal
N = len(H1_FT) # number of elements in Hanford dataset

f_Hz = np.linspace(0,samplerate/2,N) # frequency array for xaxis

w = sp.hann(samples) # window function
WFT_H1 = np.fft.rfft(w*H1data) # window function of FFT of Hanford 


# Livingston

L1_FT = np.fft.rfft(L1data) #FFT of Livingston signal
N_L1 = len(L1_FT) # number of elements in Livingston dataset

w_L1 = sp.hann(samples) # window function
WFT_L1 = np.fft.rfft(w*L1data) # window function of FFT of Livingston 

#8

WFH = 1/ (np.abs(WFT_H1)) 
WFL = 1/ (np.abs(WFT_L1))

WH = WFT_H1*WFH
WL = WFT_L1*WFL

#8.4


ir_WH = np.fft.irfft(WH)

plt.plot(time_vec, ir_WH)
plt.xlabel("Time (sec)")
plt.ylabel("Whitened strain y[n]")
plt.title("Hanford")
plt.show()

ir_WL = np.fft.irfft(WL)
plt.plot(time_vec, ir_WL)
plt.xlabel("Time (sec)")
plt.ylabel("Whitened strain y[n]")
plt.title("Livingston")
plt.show()

plt.plot(time_vec, ir_WH)
plt.xlim(right=16.5)
plt.xlim(left=16.2)
plt.show()




#def H(L):
#    return 1/L * ((1-np.exp(-1j*om*L))/(1-np.exp(-1j*om)))

#9.1
om300 = 2*np.pi*300 / 4096
om1 = np.linspace(0,np.pi,len(ir_WH))
def H_omega(om1,L):
    return 1/L * ((1-np.exp(-1j*om1*L))/(1-np.exp(-1j*om1)))


L = np.arange(1,13)

h = np.zeros(len(L))
for i in range(0, len(h)-1):
    h[i] = np.abs(H_omega(om300,L[i]))

plt.plot(L, H_omega(om300,L))
plt.plot(L,h)
plt.hlines(0.5,0,12)
plt.show()

#9.6
om = np.linspace(0,300,len(ir_WH))
def H_test(om,L):
    return 1/L * ((1-np.exp(-1j*om1*L))/(1-np.exp(-1j*om1)))

L=8

FIR = H_omega(om,L)*ir_WH

#plt.plot(FIR)

plt.show()

#9.2

omtest = np.linspace(0,np.pi,1000)
test = 10*np.log10(abs(H_omega(omtest,L))**2)

N2 = len(test)
f_Hz2 = np.linspace(0,samplerate/2,N2)


plt.plot(f_Hz2,test)
plt.show()
#plt.plot(L, h)
#plt.plot(H(L))

# test 9.2
#H1_FT = np.fft.rfft(H1data) 

#y= np.zeros(len(H1data))

#for k in range(0,8):
#    y += 1/8 * (H1data[i]-[k])
    
#plt.plot(y)
h = np.ones(8) /8
s = sp.convolve(ir_WH,h, mode = "same")

plt.plot(time_vec,s)
plt.show()


