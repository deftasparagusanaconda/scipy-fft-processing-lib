import scipy.fft, scipy.io, math, matplotlib.pyplot as plt, numpy as np

# input = scipy.io.wavfile.read("/mnt/files_2/Coding Stuff/Loudness Analysis/middle sine -18dB.wav") #assign the wav file to a variable
# input = scipy.io.wavfile.read("/home/chonga2520/EQAPO config files/impulse.wav")
input = scipy.io.wavfile.read("/mnt/files_2/Coding Stuff/Loudness Analysis/NUU$HI - Sakura.wav")

rate  = input[0]      #the sampling rate of the wav file
sound = input[1]      #the actual numbers from the wav file

left  = sound[:,0]    #left channel of the wav file
right = sound[:,1]    #right channel of the wav file
mono = left/2+right/2 #mono channel of the wav file
side = left/2-right/2 #side channel of the wav file

freq = []             #a list for the freqs that the fft will represent
for i in range(len(sound)//2+1): 
    freq.append(rate/len(sound)*(i))
log2freq = []         #freqs list but in log2 values, excluding the 0Hz component
for i in freq[1:]: 
    log2freq.append(math.log2(i))

leftrfft  = scipy.fft.rfft(left,norm="forward")
rightrfft = scipy.fft.rfft(right,norm="forward")
monorfft  = scipy.fft.rfft(mono,norm="forward")
siderfft  = scipy.fft.rfft(side,norm="forward")

#now we have the fourier transform of the L,R,M channels. now comes the data part.
#this section converts cartesian complex coordinates to polar magnitude/phase coordinates

leftmag  = abs(leftrfft)
rightmag = abs(rightrfft)
monomag  = np.array(abs(monorfft))
sidemag  = np.array(abs(siderfft))

leftphase  = []
rightphase = []
monophase  = []
sidephase  = []
for i in leftrfft:
    leftphase.append(math.atan2(i.imag,i.real))
for i in rightrfft:
    rightphase.append(math.atan2(i.imag,i.real))
for i in monorfft:
    monophase.append(math.atan2(i.imag,i.real))
for i in siderfft:
    sidephase.append(math.atan2(i.imag,i.real))

#this section converts magnitudes to decibels and radians to degrees
leftmagdB  = []
rightmagdB = []
monomagdB  = []
sidemagdB  = []
# for i in leftmag:
#     leftmagdB.append(20*math.log10(i))
# for i in rightmag:
#     rightmagdB.append(20*math.log10(i))
# for i in monomag:
#     monomagdB.append(20*math.log10(i))
# for i in sidemag:
#     sidemagdB.append(20*math.log10(i))

leftphasedeg  = []
rightphasedeg = []
monophasedeg  = []
sidephasedeg  = []
for i in leftphase:
    leftphasedeg.append(i*180/math.pi)
for i in rightphase:
    rightphasedeg.append(i*180/math.pi)
for i in monophase:
    monophasedeg.append(i*180/math.pi)
for i in sidephase:
    sidephasedeg.append(i*180/math.pi)

#plt.semilogx(freq,leftmagdB)
#plt.semilogx(freq,rightmagdB)
#plt.semilogx(freq,monomagdB)
#plt.semilogx(freq,leftphasedeg)
#plt.semilogx(freq,rightphasedeg)
#plt.semilogx(freq,monophasedeg)
#plt.loglog(freq,leftmag)
#plt.loglog(freq,rightmag)
#plt.plot(sound)
#plt.show()

