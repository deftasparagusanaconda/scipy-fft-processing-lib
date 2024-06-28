import math, numpy as np, matplotlib.pyplot as plt

samples = 8192 #how many frequencies to take as points
startfreq = 5 #the first frequency to take. also interpreted as resolution

a = startfreq #running var for loop
log2freq = []
for i in range(samples):
    log2freq.append(a)
    a += startfreq

def intervalcount(arr,l,u):
    count = 0
    for i in arr:
        if (i>=l) & (i<u):
            count += 1
    return(count)

lower = startfreq #running var for loop
delta = 40

list = []
for i in range(4096):
    list.append(intervalcount(log2freq,lower,lower+delta))
    lower += delta

list2 = []
for i in list:
    if i == 0:
        list2.append(0)
    else:
        list2.append(math.log2(i))

plt.plot(list2)
plt.show()