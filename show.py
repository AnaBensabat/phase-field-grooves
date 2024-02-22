import matplotlib.pylab as plt
import numpy as np

info = open("volume_vs_time.txt", 'r')
t = []
v = []
for line in info:
    x = line.split('\t')
    t.append(float(x[0]))
    v.append(float(x[1]))

print(min(v))
print(max(v))
plt.plot(v[10:])
plt.show()


info = open("area_vs_time.txt", 'r')
t = []
a = []
for line in info:
    x = line.split('\t')
    t.append(float(x[0]))
    a.append(float(x[1]))

print(a[:20])
plt.plot(t,a)
plt.show()
