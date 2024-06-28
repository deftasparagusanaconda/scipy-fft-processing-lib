import math, numpy as np, matplotlib.pyplot as plt

samples = 8192 #how many frequencies to take as points
startfreq = 44100/8192 #the first frequency to take. also interpreted as resolution

freq = []
a = startfreq #running var for loop
for i in range(samples):
    freq.append(a)
    a += startfreq

log2freq = []
for i in freq:
    log2freq.append(math.log2(i))

def intervalcount(arr,l,u):
    count = 0
    for i in arr:
        if (i>=l) & (i<=u):
            count += 1
    return(count)

lower = math.log2(startfreq) #running var for loop
delta = 0.004

list = []
for i in range(4096):
    list.append(intervalcount(log2freq,lower,lower+delta))
    lower += delta

list2 = []
for i in list:
    if i == 0:
        list2.append(1)
    else:
        list2.append(math.log2(i))

plt.plot(list2)
plt.show()