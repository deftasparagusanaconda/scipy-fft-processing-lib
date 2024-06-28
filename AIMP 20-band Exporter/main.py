# for now this file only supports stereo wav files, and i dont really think i need to support anything else so, the "for now" was probably "forever-ly" lol
# and yes i did comment on everything i need to know how everything works easily
# the EQ derived for AIMP's EQ traces between points using lerp in decibel magnitude and log2 frequency space

import scipy.fft, scipy.io, matplotlib.pyplot as plt, numpy as np, csv, math, random as rand, statistics as stat #import stuffs
import time #speed testing package

input = scipy.io.wavfile.read("/home/chonga2520/EQAPO config files/MC2 EQ 18.wav") #assign the wav file to a variable

rate  = input[0]   #the sampling rate of the wav file
sound = input[1]   #the actual numbers from the wav file

left  = sound[:,0] #left channel of the wav file
right = sound[:,1] #right channel of the wav file
mono = left/2+right/2 # mono channel of the wav file

freq = []         #a list for the freqs that the fft will represent
for i in range(len(sound)//2+1): 
    freq.append(rate/len(sound)*(i))
log2freq = []     #freqs list but in log2 values, excluding the 0Hz component
for i in freq[1:]: 
    log2freq.append(math.log2(i))

leftrfft  = scipy.fft.rfft(left,norm="backward")
rightrfft = scipy.fft.rfft(right,norm="backward")
monorfft  = scipy.fft.rfft(mono,norm="backward")

#now we have the fourier transform of the L,R,M channels. now comes the data part.
#this section converts cartesian complex coordinates to polar magnitude/phase coordinates

leftmag  = abs(leftrfft)
rightmag = abs(rightrfft)
monomag  = np.array(abs(monorfft))

leftphase  = []
rightphase = []
monophase  = []
for i in leftrfft:
    leftphase.append(math.atan2(i.imag,i.real))
for i in rightrfft:
    rightphase.append(math.atan2(i.imag,i.real))
for i in monorfft:
    monophase.append(math.atan2(i.imag,i.real))

#this section converts magnitudes to decibels
leftmagdB  = []
rightmagdB = []
monomagdB  = []
for i in leftmag:
    leftmagdB.append(20*math.log10(i))
for i in rightmag:
    rightmagdB.append(20*math.log10(i))
for i in monomag:
    monomagdB.append(20*math.log10(i))

#setting things up for AIMP export

AIMPfreq = []
for i in range(20):
    AIMPfreq.append(1000*2**(i/2-5))
AIMPlog2freq = []
for i in AIMPfreq:
    AIMPlog2freq.append(math.log2(i))
AIMPmagdB = [] #arbitrary "all 0dB" definition
for i in range(20):
    AIMPmagdB.append(0)

# AIMPmagdB = [1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0]
AIMPmagdB = [3.6, 2.5, 0.6, -1.6, -3, -3.8, -4.2, -4.1, -3.6, -2.8, -1.7, -0.1, 2.4, -0.1, -0.5, 4.4, 3.2, 5.1, -10.9, -15]
# AIMPmagdB = [2.2205527445602335, 1.716746279100155, -1.2380418943245608, -3.5543645013361402, -5.315312157945168, -6.3691357935633075, -6.807803447277956, -6.705433834541214, -6.100335830557888, -5.04722723311267, -3.8081726584004785, -1.88920948239678, 1.0379262910371823, -1.8244366133365735, -1.0824767937795177, 3.0440029270621825, 2.629489656230579, 5.462670545047224, -29.91524223798182, 0]

def twonear(x,y): #returns the index of the two nearest numbers to x in y
    #x is number
    #y is array or list
    if y[0] > x:
        return(0,1)
    elif y[-1] < x:
        return(-2,-1)
    else:
        l = 0
        h = -1
        for i in range(len(y)):
           if y[i] <= x:
               l = i
           else:
               break
        for i in range(len(y)):
            if y[::-1][i] >= x:
                h = len(y)-i-1
            else:
                break
        return(l,h)

def calcEQmagdB(arr): #returns the EQmagdB of the final EQ from AIMPmagdB and DC component from monomagdB
    output = [monomagdB[0]]
    for i in range(len(log2freq)): #append AC components into EQmagdB. derived from AIMPmagdB using lerp of points, in decibel magnitude and log2 frequency space.
        (l,h) = twonear(log2freq[i],AIMPlog2freq)
        factor = (log2freq[i]-AIMPlog2freq[l])/(AIMPlog2freq[h]-AIMPlog2freq[l])
        output.append(arr[h]*factor+arr[l]*(1-factor))
    return(output)

def calcEQmag(arr): #eats EQmagdB and poops EQmag
    output = []
    for i in arr:
        output.append(10**(i/20))
    return(output)

EQmagdB = calcEQmagdB(AIMPmagdB)
EQmag = calcEQmag(EQmagdB)

plt.semilogx(freq,EQmagdB)
plt.show()

audbound = [29,15000] #desired audible frequency range in Hertz
audboundlog2 = [math.log2(audbound[0]),math.log2(audbound[1])] #audfreqbound in log2 frequency
audboundlog2nearest = [twonear(audboundlog2[0],log2freq)[0],twonear(audboundlog2[1],log2freq)[1]] #log2freq index of frequency bounds. the numbers are the closest numbers just outside of the range but still in log2freq.
l = audboundlog2nearest[0] #redefinition to make naming easier
h = audboundlog2nearest[1] #redefinition to make naming easier

log2scaler = [] #just a hadamard product scaler the auderr function uses to scale frequencies according to log2 scale. higher freqs weigh less as they are denser. it is 1/frequency.
for i in freq[1:]:
    log2scaler.append(1/i)

def auderr(x,y): #return the audible error between two arrays, as a sum across all discrete frequencies. higher frequencies are weighted less according to log2 scale
    audx = x[l:h]*np.array(log2scaler)[l:h]
    audy = y[l:h]*np.array(log2scaler)[l:h]
    return(sum(abs(audx-audy)))

def audMAD(arr): #returns the audible mean absolute deviation of an array. frequencies are scaled according to their density (1/freq)
    audarr = arr[l:h]
    mean = sum(audarr*log2scaler[l:h])/sum(log2scaler[l:h])
    output = sum(abs((np.full(len(audarr),mean)-audarr))*log2scaler[l:h])/sum(log2scaler[l:h])
    return output

    

error = audMAD(EQmag-monomag)

print("starting EQ   : ", AIMPmagdB)
print("starting error:", error)

# print(freq)
# print(2**log2freq[l],2**log2freq[h]) #prints the freq bounds the program uses, in Hertz
plt.plot(log2freq[l:h],EQmag[l+1:h+1])
plt.plot(log2freq[l:h],monomag[l+1:h+1])
plt.plot(log2freq[l:h],(EQmag-monomag)[l+1:h+1])

plt.show()

##iterative error minimization algorithm. takes steps of scaled random size in one dependent variable at a time.
iter = 500 #number of iterations 
step_factor = 10 #how the step_size scales with error
step_size = error*step_factor #how large of a step it takes

list = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
EQchangesdB = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
hit_count = 0

for i in range(iter):
    AIMPmagdBnew = AIMPmagdB.copy()
    rndnum = (rand.random())*2-1 #random number in [-1, 1)
    rndint = round(rand.random()*19) #random integer in [0, 19]
    AIMPmagdBnew[rndint] += rndnum*step_size #changes random element in AIMPmagdBnew and adds random number scaled by step size
    EQmagdBnew = calcEQmagdB(AIMPmagdBnew)
    EQmagnew = calcEQmag(EQmagdBnew)
    errornew = audMAD(EQmagnew-monomag)
    if errornew < error:
        error = errornew
        hit_count += 1
        AIMPmagdB[rndint] += rndnum*step_size
        EQchangesdB[rndint] += rndnum*step_size
        step_size = errornew*step_factor

print("EQ changes: " + str(EQchangesdB))
print("final error: " + str(error))
print("final EQ: " + str(AIMPmagdB))
print("hit_count:", hit_count)

EQmagdB = calcEQmagdB(AIMPmagdB)
EQmag = calcEQmag(EQmagdB)
EQaudmagmean = sum(np.array(EQmag[l:h])*log2scaler[l:h])/sum(log2scaler[l:h])
diffaudmagmean = sum((EQmag-monomag)[l:h]*log2scaler[l:h])/sum(log2scaler[l:h])

print(EQaudmagmean,diffaudmagmean)

for i in range(len(AIMPmagdB)): #audible magnitude correction
    AIMPmagdB[i] -= 20*math.log10(EQaudmagmean)

EQmagdB = calcEQmagdB(AIMPmagdB)
EQmag = calcEQmag(EQmagdB)
EQaudmagmean = sum(np.array(EQmag[l:h])*log2scaler[l:h])/sum(log2scaler[l:h])
diffaudmagmean = sum((EQmag-monomag)[l:h]*log2scaler[l:h])/sum(log2scaler[l:h])

print(EQaudmagmean,diffaudmagmean)

plt.plot(log2freq[l:h],EQmag[l+1:h+1])
plt.plot(log2freq[l:h],monomag[l+1:h+1])
plt.plot(log2freq[l:h],(EQmag-monomag)[l+1:h+1]-diffaudmagmean)
plt.show()
# plt.plot(freq[l:h],EQmag[l:h])
# plt.plot(freq[l:h],monomag[l:h])
# plt.plot(freq[l:h],(EQmag-monomag)[l:h])




##check whether fft shows a 0dB sine as 0dB
##           try to work purely in log2 freqs and magnitudes

# outputfile = open("/mnt/files_2/Coding Stuff/fft stufs/output.csv","w+")
# outputobj = csv.writer(outputfile, delimiter=',')
# outputobj.writerow(dft2)