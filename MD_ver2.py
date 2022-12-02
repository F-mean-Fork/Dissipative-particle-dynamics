import numpy as np
import matplotlib.pyplot as plt
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
        self.r_cut = 2.4
        self.dbl = self.double_list()
        self.dt = dt
        self.r_min = self.box
        self.E_k = 0.0
        self.T = 0.0
        self.U = 0.0
    
    def initial_coordinates (self):
        self.x = (np.random.sample(self.N) - 0.5) * self.box
        self.y = (np.random.sample(self.N) - 0.5) * self.box
        self.z = (np.random.sample(self.N) - 0.5) * self.box
        
    def initial_velocity(self):
        self.vx = 2*np.random.sample(self.N) - 1
        self.vy = 2*np.random.sample(self.N) - 1
        self.vz = 2*np.random.sample(self.N) - 1
        Sum_vx = np.sum(self.vx) / self.N
        Sum_vy = np.sum(self.vy) / self.N
        Sum_vz = np.sum(self.vz) / self.N
        self.vx -= Sum_vx
        self.vy -= Sum_vy
        self.vz -= Sum_vz
    
    def print_picture(self, namefile = 'test'):
        file = open(namefile + '.ent', 'w')
        for i in range (self.N):         
            file.write(f'HETATM%5d  {self.type[i]}%12d' % (i+1, i+1) +'%12.3f%8.3f%8.3f \n' % (self.x[i], self.y[i], self.z[i]))
        file.close()
        
    def soft_force(self, x, a=1.0):
        if x < self.r_cut:
            return (1.0 - x)*a
        else:
            return 0.0
    
    def soft_energy(self, a, r):
        if r < self.r_cut:
            return 0.5*a*(1-r)**2
        else:
            return 0.0
    
    def double_list(self):
        ub1 = []
        ub2 = []
        for i in range(self.N-1):
            for j in range(i+1, self.N):
                ub1.append(i)
                ub2.append(j)
        return list(zip(ub1, ub2))
    
    def new_acceleration(self, force):
        self.U = 0.0
        self.ax[:] = 0.0
        self.ay[:] = 0.0
        self.az[:] = 0.0
        self.r_min = self.box      
        for i, j in self.dbl:
            dx = self.x[i] - self.x[j] 
            dy = self.y[i] - self.y[j]
            dz = self.z[i] - self.z[j]               
            dx = self.periodic_test(dx)
            dy = self.periodic_test(dy)
            dz = self.periodic_test(dz)            
            r = np.sqrt(dx**2 + dy**2 + dz**2)
            if r < self.r_min:
                self.r_min = r
            f = force(r) 
            self.U += self.LJ_energy(r)         
            self.ax[i] += f * dx / r
            self.ay[i] += f * dy / r
            self.az[i] += f * dz / r
            self.ax[j] -= f * dx / r
            self.ay[j] -= f * dy / r
            self.az[j] -= f * dz / r
        
    def new_coordinates(self):
        for i in range(self.N):
            self.x[i] += 0.5*self.ax[i] * self.dt**2 + self.vx[i] * self.dt
            self.y[i] += 0.5*self.ay[i] * self.dt**2 + self.vy[i] * self.dt
            self.z[i] += 0.5*self.az[i] * self.dt**2 + self.vz[i] * self.dt
            self.x[i] = self.periodic_test(self.x[i])
            self.y[i] = self.periodic_test(self.y[i])
            self.z[i] = self.periodic_test(self.z[i])
        
    def periodic_test(self, x):
        if np.abs(x) > 0.5 * self.box:
            return x - x/np.abs(x) * self.box
        else:
            return x
    
    def new_velocity(self):
        self.vx += self.ax * self.dt*0.5
        self.vy += self.ay * self.dt*0.5
        self.vz += self.az * self.dt*0.5
        
    def temperature(self):
        self.E_k = 0.5 * np.sum(self.vx**2 + self.vy**2 + self.vz**2) / self.N
        self.T = 2/3 * self.E_k
        
    def termostat(self):
        betta = np.sqrt(np.sum(self.vx**2 + self.vy**2 + self.vz**2) / (3*self.N))
        self.vx /= betta
        self.vy /= betta
        self.vz /= betta 
    
    def LJ_force(self, r, sigma=1.0, eps=1.0):
        R = (sigma / r)**6
        if r < self.r_cut:
            return 4*eps*R/r*(R-1)
        else:
            return 0.0
        
    def LJ_energy(self, r, sigma=1.0, eps=1.0):
        R = (sigma / r)**6
        if r < self.r_cut:
            return 4*eps*(R**2-R)
        else:
            return 0.0
        
N_step = 400
System1 = MD(108, 6, 0.1) 
System1.initial_coordinates()
System1.print_picture('begin')
System1.r_cut = 1 
System1.new_acceleration(System1.soft_force)
for i in range(N_step):    
    System1.new_coordinates()
    System1.new_acceleration(System1.soft_force)
    System1.temperature()

N_step = 4000
System1.initial_velocity()
System1.termostat()
System1.r_cut = 2.4
System1.dt = 0.01
U_list = []
E_k_list = []
System1.new_acceleration(System1.LJ_force)
for i in range(N_step):    
    System1.new_coordinates()
    System1.new_velocity()
    System1.new_acceleration(System1.LJ_force)
    System1.new_velocity()
    System1.temperature()
    # System1.termostat()
    U_list.append(System1.U/System1.N)
    E_k_list.append(System1.E_k)
    print(i, System1.r_min, System1.U/System1.N, System1.E_k)
System1.print_picture('end')     

plt.plot(E_k_list, color = "red", linestyle = "--", label = "E_k")
plt.plot(U_list, color = 'black', linestyle = "-", label = "U" )
plt.legend()
plt.draw()
plt.show()

def print_distribution(left, right, N, arr):
  step = (right - left) / N
  x = np.zeros(N)
  count = 0
  for i in arr:
    if left < i < right:
      n = int((i - left) / step)
      x[n] +=1
      count +=1
  return np.arange(left+step*0.5, right, step), x / count
  
v = np.sqrt(System1.vx**2 + System1.vx**2 + System1.vx**2)
x, y = print_distribution(-3, 3, 10, v)
plt.plot(x, y)
plt.show()
