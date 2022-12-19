import numpy as np
import matplotlib.pyplot as plt
import random

def print_picture(X, Y, Z, types, bonds, namefile = 'test'):
    file = open(namefile + '.ent', 'w')
    n = 0
    for x, y, z, t in zip(X, Y, Z, types): 
        n += 1        
        file.write(f'HETATM%5d  {t}%12d' % (n, n) +'%12.3f%8.3f%8.3f \n' % (x, y, z))
    
    
    for i in bonds:
        dx = X[i[0]-1] - X[i[1]-1] 
        dy = Y[i[0]-1] - Y[i[1]-1] 
        dz = Z[i[0]-1] - Z[i[1]-1] 
        r = dx**2 + dy**2 + dz**2
        if r < 2.0:
            file.write(f'CONECT%5d%5d' % (i[0], i[1]) +'\n')
        
    file.close()

def periodic_test(x, box):
    if np.abs(x) > 0.5 * box:
        return x - x/np.abs(x) * box
    else:
        return x
        
def random_vector():
    vector = np.random.sample(3) - 0.5
    vector[2] = np.abs(vector[2]) * 2.0
    s = np.sum(vector[:]**2)
    return vector / np.sqrt(s)

class MatchFlat:
    half_ball = 0.288675
    def __init__(self, num_points, num_primary, box, r = 0.288675):
        # Major parameters
        self.num_points = num_points
        self.num_primary = num_primary
        self.box = box
        # Coordinates and accelerations
        self.x = (np.random.sample(num_points) - 0.5) * self.box[0]
        self.y = (np.random.sample(num_points) - 0.5) * self.box[1]
        self.z = np.array([-self.box[2] * 0.5 + MatchFlat.half_ball] * self.num_points)
        self.ax = np.zeros(num_points, dtype=float)
        self.ay = np.zeros(num_points, dtype=float)
        # Vector

        
        self.r_cut = self.box[0] * 0.5
        self.types = ['N'] * (self.num_points - self.num_primary) + ['O'] * self.num_primary
        self.dt = 0.01
        self.dbl = self.double_list()
        self.r_min = self.box[0]
        self.repulsive = 1.0
        
    def double_list(self):
        ub1 = []
        ub2 = []
        for i in range(self.num_points-1):
            for j in range(i+1, self.num_points):
                ub1.append(i)
                ub2.append(j)
        return list(zip(ub1, ub2))
    
    def soft_force(self, x, a = 1.0):
        if x < self.r_cut:
            return (self.r_cut- x) / self.r_cut * a
        else:
            return 0.0
    
    def move(self):
        while True:
            self.new_acceleration()
            self.new_coordinates()
            yield True
        
    def new_coordinates(self):
            for i in range(self.num_points):
                self.x[i]  += self.ax[i] * self.dt
                self.y[i] += self.ay[i] * self.dt
                self.x[i] = periodic_test(self.x[i], self.box[0])
                self.y[i] = periodic_test(self.y[i], self.box[1])
            
    def new_acceleration(self):
        self.ax[:] = 0.0
        self.ay[:] = 0.0   
        self.r_min = self.box[0]
        for i, j in self.dbl:
            dx = self.x[i] - self.x[j] 
            dy = self.y[i] - self.y[j]             
            dx = periodic_test(dx, self.box[0])
            dy = periodic_test(dy, self.box[1])     
            r = np.sqrt(dx**2 + dy**2)
            if r < self.r_min:
                self.r_min = r  
            if self.types[i] == self.types[j]:
               f = self.soft_force(r, a = self.repulsive) 
            else:
                f = self.soft_force(r)     
            self.ax[i] += f * dx / r
            self.ay[i] += f * dy / r
            self.ax[j] -= f * dx / r
            self.ay[j] -= f * dy / r 
    
def random_surface(M, N, box):
    MF = MatchFlat(N+M, M, [box, box, box])
    # print_picture(MF.x, MF.y, MF.z, MF.types)
    label = True
    one_step = MF.move()
    result = [0.0]
    EPS = 1e-5
    MF.repulsive = 5.0
    
    while True:
        label = next(one_step)
        result.append(MF.r_min)
        if np.abs(result[-1] - result[-2]) < EPS:
            break
        
    result[-1] = 0.0  
    MF.repulsive = 1.0
    while True:
        label = next(one_step)
        result.append(MF.r_min)
        if np.abs(result[-1] - result[-2]) < EPS:
            break
    
    return MF.x, MF.y, MF.z, MF.types

def chess_surface(M, N, box):
    assert N > 0 and M > 0
    assert N == M 
    assert int((np.sqrt(N+M))**2) == N+M
    
    k = int(np.sqrt(N+M))
    coordinates = []
    for i in range(k):
        for j in range(k):
            coordinates.append((i,j))

    x = []
    y = []
    z = [MatchFlat.half_ball - box * 0.5] * (N+M)
    types = []
    delta = box / k 
    for cord in coordinates:
        x.append(delta*(0.5 + cord[0]) - box*0.5)
        y.append(delta*(0.5 + cord[1]) - box*0.5)
        if (cord[0]+cord[1]) % 2 == 0:
            types.append('N')
        else:
            types.append('O')
    return x, y, z, types

def brush_create(m1, m2, n1, n2, box, segment_lenght = (1/3)**(1/3)):
    assert box[0] == box[1]
    
    try:
        x0, y0, z0, t = chess_surface(m1, m2, box[0])
    except:
        x0, y0, z0, t = random_surface(m1, m2, box[0])
    
    x, y, z, types = [], [], [], []
    
    counter = 0
    bonds = []
    for X0, Y0, Z0, T in zip(x0, y0, z0, t):
        print(X0, Y0, Z0, T)
        if T == 'N':
            chain_lenght = n1
        else:
            chain_lenght = n2
        
        x.append(X0)
        y.append(Y0)
        z.append(Z0)
        types.append(T)      
        counter += 1
        # print(T, counter)

        for _ in range(chain_lenght-1):
            label = True
            while label:
                
                r = random_vector()*segment_lenght
                xt = x[-1] + r[0]
                yt = y[-1] + r[1]
                zt = z[-1] + r[2]
                if -box[2] * 0.5 < zt < box[2]*0.5:
                    xt = periodic_test(xt, box[0])
                    yt = periodic_test(yt, box[1])
                    x.append(xt)
                    y.append(yt)
                    z.append(zt)
                    types.append(T)
                    counter += 1
                    bonds.append([counter-1, counter])
                    label = False
                    
    return x, y, z, types, bonds
    
if __name__ == '__main__':
    M = 8
    N = 7
    box = 8.0

    x, y, z, t, bonds = brush_create(m1 = 20, m2 = 20, n1 = 15, n2 = 10, box = [8.0, 8.0, 10.0])
    print_picture(x, y, z, t, bonds)
