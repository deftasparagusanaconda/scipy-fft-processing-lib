import scipy.fft as fft, scipy.io, math, matplotlib.pyplot as plt, numpy as np

infreq = []
for i in range(31):
    infreq.append(1000*2**(i/3-17/3))

inlog2freq = [] #based on freq. not independently calculated.
for i in infreq:
    inlog2freq.append(math.log2(i))

rate = 44100 #in samples per second
length = 8192 #in samples

list = [0,0,0,0,0]

fft.irfft.()

