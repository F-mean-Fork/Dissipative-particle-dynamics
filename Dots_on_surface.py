import numpy as np
import matplotlib.pyplot as plt
import random

def print_picture(X, Y, Z, types, namefile = 'test'):
    file = open(namefile + '.ent', 'w')
    n = 0
    for x, y, z, t in zip(X, Y, Z, types): 
        n += 1        
        file.write(f'HETATM%5d  {t}%12d' % (n, n) +'%12.3f%8.3f%8.3f \n' % (x, y, z))
    file.close()

def periodic_test(x, box):
    if np.abs(x) > 0.5 * box:
        return x - x/np.abs(x) * box
    else:
        return x
        
def random_vector():
    vector = np.random.sample(3)- 0.5
    s = np.sum(vector[:]**2)
    return vector / np.sqrt(s)

class MatchFlat:
    def __init__(self, num_points, num_primary, N, box, r = 0.288675):
        # Major parameters
        self.num_points = num_points
        self.num_primary = num_primary
        self.box = box
        self.N = N                     # Chains` amount
        self.l = 1.0                   # Segment` lenth
        self.H = self.box[2]           # Box height
        # Coordinates and accelerations
        self.x = (np.random.sample(num_points) - 0.5) * self.box[0]
        self.y = (np.random.sample(num_points) - 0.5) * self.box[1]
        self.z = np.array([-self.box[2] * 0.5 + 0.288675] * self.num_points)
        self.ax = np.zeros(num_points, dtype=float)
        self.ay = np.zeros(num_points, dtype=float)
        # Vector
        # self.vec_x = random.uniform(-1, 1)
        # self.vec_y = random.uniform(-1, 1)
        # self.vec_z = random.uniform(-1, 1)
        
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
    
    # def chain_construct(self):
    #     Z = [self.z]
    #     for i in range(1, self.N):
    #         l_x = random.uniform(-1, 1)
    #         l_y = random.uniform(-1, 1)
    #         l_z = random.uniform(-1, 1)
    #         s = np.sqrt(l_x**2 + l_y**2 + l_z**2)
    #         z_new = Z[-1] + s
    #         if z_new < self.H:
    #             Z.append(z_new)
    #         else:
                
        # # Initial chains` lenths
        # vector_lenth_x = self.x
        # vector_lenth_y = self.y
        # vector_lenth_z = self.z
        # # Chain` lenth check
        # #for i in range(self.N):
        # if vector_lenth_x < self.box[0]:
        #     vector_lenth_x += self.l * self.r_x / s
        # else:
        #     vecor_lenth_x -= self.l * self.r_x / s
            
        # if vector_lenth_y < self.box[1]:
        #     vector_lenth_y += self.l * self.r_y / s
        # else:
        #     vector_lenth_y -= self.l * self.r_y / s
            
        # if vector_lenth_z < self.H:
        #     vector_lenth_z += self.l * self.r_z / s
        # else:
        #     vector_lenth_z -= self.l * self.r_z / s           
            
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
    
    return MF.x, MF.y, MF.types

# plt.plot(result, color = 'black')
# plt.show()
    
# print_picture(MF.x, MF.y, MF.z, MF.types, namefile = 'end')
print(random_vector())