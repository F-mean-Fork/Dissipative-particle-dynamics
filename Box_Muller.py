import matplotlib.pyplot as plt
import numpy as np
import random
import time 
import ctypes as ct

fortlib = ct.CDLL(r'C:\Users\ivanl\Desktop\ITMO\Course_3\Polymer_Brushes\dpd\myflib_1.so')

# fortlib.GaussRand.argtypes = [ct.POINTER(ct.c_float)]
# fortlib.GaussRand.restypes = ct.c_float
# x = 10
# fortlib.GaussRand(x)
# print(x)

# n = 1000000
# y = [random.random() for i in range(n)]

# def box_muller(r, phi):
#     z0 = np.cos(2*np.pi*phi)*np.sqrt(-2*np.log(r))
#     z1 = np.sin(2*np.pi*phi)*np.sqrt(-2*np.log(r))
#     return z0, z1

# def rand_norm(label = False):
#     while True:
#         label = not label
#         if label:
#             r, phi  = np.random.rand(), np.random.rand()
#             z0 = np.cos(2*np.pi*phi)*np.sqrt(-2*np.log(r))
#             z1 = np.sin(2*np.pi*phi)*np.sqrt(-2*np.log(r))
#             yield z0
#         else:
#             yield z1
 
# t1 = time.time()     
# r = rand_norm(False)    
# y_n = [next(r) for i in range(n)]

# for i in range(0, len(y), 2):
#     z0, z1 = box_muller(y[i], y[i+1]) 
#     y_n.append(z0)  
#     y_n.append(z1)  

# t2 = time.time()
# print(t2-t1)

# hist, bin_edges = np.histogram(y_n, bins=20, density=True)
# bin_centers = [0.5*(bin_edges[i-1]+bin_edges[i]) for i in range(1,len(bin_edges))]
# plt.xlabel('x')
# plt.ylabel('weight')
# plt.ylim(0.0, 1.1)
# plt.plot(bin_centers,hist)
# plt.show()
