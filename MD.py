import numpy as np

class MD:
    def __init__(self, N, box, dt):
        self.N = N
        self.box = box
        self.x = np.zeros(N, dtype=float)
        self.y = np.zeros(N, dtype=float)
        self.z = np.zeros(N, dtype=float)
        self.vx = np.zeros(N, dtype=float)
        self.vy = np.zeros(N, dtype=float)
        self.vz = np.zeros(N, dtype=float)
        self.ax = np.zeros(N, dtype=float)
        self.ay = np.zeros(N, dtype=float)
        self.az = np.zeros(N, dtype=float)
        self.step = 0
        self.type = ['N'] * self.N
        self.r_cut = 1.0
        self.dbl = self.double_list()
        self.dt = dt
    
    def initial_coordinates (self):
        self.x = (np.random.sample(self.N) - 0.5) * self.box
        self.y = (np.random.sample(self.N) - 0.5) * self.box
        self.z = (np.random.sample(self.N) - 0.5) * self.box
        
    def print_picture(self, namefile = 'test'):
        file = open(namefile + '.ent', 'w')
        for i in range (self.N):         
            file.write(f'HETATM%5d  {self.type[i]}%12d' % (i+1, i+1) +'%12.3f%8.3f%8.3f \n' % (self.x[i], self.y[i], self.z[i]))
        file.close()
        
    def soft_force(self, a, x):
        if x < self.r_cut:
            return (1.0 - x)*a
        else:
            return 0.0
    
    def double_list(self):
        ub1 = []
        ub2 = []
        for i in range(self.N-1):
            for j in range(i+1, self.N):
                ub1.append(i)
                ub2.append(j)
        return zip(ub1, ub2)
    
    def new_acceleration(self):
        self.ax[:] = 0.0
        self.ay[:] = 0.0
        self.az[:] = 0.0
        r_min = self.box
        
        for i, j in self.dbl:
            dx = self.x[i] - self.x[j] 
            dy = self.y[i] - self.y[j]
            dz = self.z[i] - self.z[j]    
            
            dx = self.periodic_test(dx)
            dy = self.periodic_test(dy)
            dz = self.periodic_test(dz)
            
            r = np.sqrt(dx**2 + dy**2 + dz**2)
            if r < r_min:
                r_min = r
            # print(r)
            f = self.soft_force(1.0, r)            
            self.ax[i] += f * dx / r
            self.ay[i] += f * dy / r
            self.az[i] += f * dz / r
            self.ax[j] += - f * dx / r
            self.ay[j] += - f * dy / r
            self.az[j] += - f * dz / r
        print(r_min, 'R minimal`noe')
        
    def new_coordinates(self):
        self.x += self.ax**2 * self.dt * 0.5
        self.y += self.ay**2 * self.dt * 0.5
        self.z += self.az**2 * self.dt * 0.5
        for i in range(self.N):              
            self.x[i] = self.periodic_test(self.x[i])
            self.y[i] = self.periodic_test(self.y[i])
            self.z[i] = self.periodic_test(self.z[i])
        
    def periodic_test(self, x):
        if np.abs(x) > 0.5 * self.box:
            return x - x/np.abs(x) * self.box
        else:
            return x
        
    
System1 = MD(10, 6, 0.001) 
print(System1.x)
System1.initial_coordinates()
print(System1.x)
# System1.print_picture()

for i in range(10):
    System1.new_acceleration()
    System1.new_coordinates()
