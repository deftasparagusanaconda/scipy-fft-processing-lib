import scipy.fft, scipy.io, matplotlib.pyplot as plt, numpy, csv, math #import stuffs

input = scipy.io.wavfile.read("/home/chonga2520/EQAPO config files/m72.wav") #assign the wav file to a variable

rate = input[0]  #the sampling rate of the wav file
sound = input[1] #the actual numbers from the wav file
freqs = []       #a list for the freqs that the fft will represent
for i in range(len(sound)//2+1): #appends values for freqs
    freqs.append(rate/len(sound)*(i))
log2freqs = []      #freqs list but in log2 values, excluding the 0Hz component
for i in freqs[1:]: #appends values for log2freqs
    log2freqs.append(math.log2(i))

left = sound[:,0]  #left channel of the wav file
right = sound[:,1] #right channel of the wav file

rfft2 = scipy.fft.rfft2(sound,norm="forward")
# rightrfft = scipy.fft.rfft(right,norm="forward")

mag = abs(rfft2)

phase = []
for i in rfft2: #appends values for phase
    phase.append(math.atan2(i.imag,i.real)*180/math.pi)

decibelmag = []
for i in mag: #appends values for decibelmag
    decibelmag.append(20*math.log10(i))


plt.plot(log2freqs,decibelmag)
plt.show()

# outputfile = open("/mnt/files_2/Coding Stuff/fft stufs/output.csv","w+")
# outputobj = csv.writer(outputfile, delimiter=',')
# outputobj.writerow(dft2)