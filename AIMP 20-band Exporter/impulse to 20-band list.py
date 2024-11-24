# specialized python script to convert impulse wav to AIMP 20-band EQ dB values
# supports only 2-channel wav files

import sys, scipy.io, scipy.fft, numpy, math, random
import matplotlib.pyplot as plt, csv

if len(sys.argv) == 1:
	filepath = input("filepath: ")
else:
	filepath = sys.argv[1]

print("importing", filepath)
input = scipy.io.wavfile.read(filepath)
print("done! now setting up variables...")

freq_low_limit = 29
freq_high_limit = 15000
sampling_rate  = input[0]
sound = input[1]
mono = sound[:,0]/2+sound[:,1]/2

print("\nfreq_low_limit =", freq_low_limit)
print("freq_high_limit =", freq_high_limit)
print("sampling_rate =", sampling_rate)
print("len(sound) =", len(sound))

print("done! now performing real-domain fast fourier transform...")
mono_rfft = scipy.fft.rfft(mono,norm="backward")
print("done! now generating arrays...")
mono_mag = numpy.array(abs(mono_rfft))
mono_mag_dB = []
for i in mono_mag:
	mono_mag_dB.append(20*math.log10(i))

freq = []
for i in range(len(sound)//2+1): 
	freq.append(sampling_rate/len(sound)*(i))
freq_log2 = []
for i in freq[1:]: 
	freq_log2.append(math.log2(i))

# variables for AIMP EQ curve
AIMP_mag_dB = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

AIMP_freq = []
for i in range(20):
	AIMP_freq.append(1000*2**(i/2-5))
AIMP_freq_log2 = []
for i in AIMP_freq:
	AIMP_freq_log2.append(math.log2(i))

## predefinitions before approximation of EQ curve -----------------------------

print("done! now getting ready for approximation...")

# index of the two nearest successive numbers in sorted list
def nearest_two(target,list):
	if list[0] > target:
		return(0,1)
	elif list[-1] < target:
		return(-2,-1)
	else:
		low = 0
		high = -1
		for i in range(len(list)):
			if list[i] <= target:
				low = i
			else:
				break
		for i in range(len(list)):
			if list[::-1][i] >= target:
				high = len(list)-i-1
			else:
				break
		return(low,high)

# calculate l & h, the index of the nearest freqs just outside of freq limits
print("done! now calculating indices of frequency limits...")
freq_low_limit_log2 = math.log2(freq_low_limit)
freq_high_limit_log2 = math.log2(freq_high_limit)
l = nearest_two(freq_low_limit_log2,freq_log2)[0]+1
h = nearest_two(freq_high_limit_log2,freq_log2)[1]+1

print("\nl =", l, "\t\t", freq[l])
print("h =", h, "\t", freq[h])
print("\nnow setting up conversion functions...")
# convert a 20-value list to an EQ curve
# uses linear interpolation
def list_to_curve(arr):
	output = [mono_mag_dB[0]]	# DC component
	for i in range(len(freq_log2)):
		(l,h) = nearest_two(freq_log2[i], AIMP_freq_log2)
		factor = (freq_log2[i]-AIMP_freq_log2[l])/(AIMP_freq_log2[h]-AIMP_freq_log2[l])
		output.append( arr[l]*(1-factor) + arr[h]*factor )
	return(output)

# convert a list of dB values to linear values
def dB_to_lin(arr):
	output = []
	for i in arr:
		output.append(10**(i/20))
	return output

# the weighting for frequencies, as data points are denser at higher frequencies
print("done! now initializing frequency weightage...")
freq_weight = []
for i in freq[1:]:
	freq_weight.append(1/math.log2(i+1))

# audible absolute sum of an array
# sum of absolute of all values in [l, h] weighted by freq_weight
print("done! now setting up error calculator...")
def audible_error(arr):
	return sum(abs(arr[l:h]*freq_weight[l:h]))

## first approximation ---------------------------------------------------------

# rough approximation using linear interpolation between two nearest points
# basically a curve_to_list function

print("done! approximation may commence")
print("\ngenerating initial approximation...")

for i in range(19):
	a, b = nearest_two(AIMP_freq_log2[i], freq_log2)
	factor = (AIMP_freq_log2[i]-freq_log2[a])/(freq_log2[b]-freq_log2[a])
	change = mono_mag[a]*(1-factor) + mono_mag[b]*factor
	AIMP_mag_dB[i] += 20*math.log10(change)

# fix 22KHz outlier
AIMP_mag_dB[19] = 2*AIMP_mag_dB[18] - AIMP_mag_dB[17]

print("done!")
print("\nfreq\t dB")
for i in range(20):
	print(round(AIMP_freq[i]),"\t",AIMP_mag_dB[i])

EQ_mag_dB = list_to_curve(AIMP_mag_dB)
EQ_mag = dB_to_lin(EQ_mag_dB)
error = audible_error(numpy.subtract(EQ_mag,mono_mag))

print("\ninitial error:", error)

# plot initial approximation
plt.semilogx(AIMP_freq, AIMP_mag_dB)

# plot initial error curve
plt.semilogx(freq[l:h+1],numpy.subtract(list_to_curve(AIMP_mag_dB),mono_mag_dB)[l:h+1])

## iterative error minimization algorithm --------------------------------------

# takes steps of pre-determined size in one value at a time

print("\nsetting up iterative error minimization algorithm...")

step_size = 0.05

print("stepping in sizes of", step_size)

EQ_changes_dB = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
hit_count = 0
try_count = 0
sweep_count = 0
try_limit = 3000

print("done! now sweeping...\n")

while try_count < try_limit:
	changes = 0
	for j in range(19):
		try_count += 1
		AIMP_mag_dB_new = numpy.add(AIMP_mag_dB,EQ_changes_dB)
		AIMP_mag_dB_new[j] += step_size
		EQ_mag_dB_new = list_to_curve(AIMP_mag_dB_new)
		EQ_mag_new = dB_to_lin(EQ_mag_dB_new)
		error_new = audible_error(numpy.subtract(EQ_mag_new,mono_mag))
		if error_new < error:
			error = error_new
			hit_count += 1
			changes += 1
			EQ_changes_dB[j] += step_size
	for j in range(19):
		try_count += 1
		AIMP_mag_dB_new = numpy.add(AIMP_mag_dB,EQ_changes_dB)
		AIMP_mag_dB_new[j] -= step_size
		EQ_mag_dB_new = list_to_curve(AIMP_mag_dB_new)
		EQ_mag_new = dB_to_lin(EQ_mag_dB_new)
		error_new = audible_error(numpy.subtract(EQ_mag_new,mono_mag))
		if error_new < error:
			error = error_new
			hit_count += 1
			changes += 1
			EQ_changes_dB[j] -= step_size
	sweep_count += 1
	if changes == 0:
		print("no new hits :)")
		print("stopping...")
		break
	else:
		print("sweep      :", sweep_count)
		print("error      :", error_new)
		print("new hits   :", changes,"/ 38")
		print("total hits :", hit_count,"/", try_count, "\n")

# fix 22KHz outlier
AIMP_mag_dB[19] = 2*AIMP_mag_dB[18] - AIMP_mag_dB[17]

EQ_mag_dB = list_to_curve(AIMP_mag_dB)
EQ_mag = dB_to_lin(EQ_mag_dB)
error = audible_error(numpy.subtract(EQ_mag,mono_mag))
print("\ninitial error = ", error,)
EQ_mag_dB = list_to_curve(numpy.add(AIMP_mag_dB,EQ_changes_dB))
EQ_mag = dB_to_lin(EQ_mag_dB)
error = audible_error(numpy.subtract(EQ_mag,mono_mag))
print("  final error = ", error)
print("    hit count = ", hit_count, "/", try_count)

## output ----------------------------------------------------------------------

print("\nfreq\t dB\t first\tchanges")
for i in range(20):
        print(round(AIMP_freq[i]), "\t", round(numpy.add(AIMP_mag_dB,EQ_changes_dB)[i],1), "\t", round(AIMP_mag_dB[i],1), "\t", round(EQ_changes_dB[i],1))

# plot a straight line
plt.semilogx([freq[l], freq[h]], [0,0])

# plot impulse curve
plt.semilogx(freq[l:h+1], mono_mag_dB[l:h+1])

# plot EQ curve
plt.semilogx(AIMP_freq, numpy.add(AIMP_mag_dB,EQ_changes_dB))

# plot EQ curve approximation
plt.semilogx(freq[l:h+1], EQ_mag_dB[l:h+1])

# plot error curve
plt.semilogx(freq[l:h+1],numpy.subtract(list_to_curve(numpy.add(AIMP_mag_dB,EQ_changes_dB)),mono_mag_dB)[l:h+1])

plt.show()
