import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np 
import glob
import time

#todo:
#http://matplotlib.org/examples/mplot3d/scatter3d_demo.html

speed = 12
mass_multiplicator = 2 
def read_file(filename):
    
    x = []
    y = []
    masses = []
    
    for line in open(filename):
        values = line.split(',')
        x.append(float(values[1]))
        y.append(float(values[2]))
        masses.append(float(values[0]))
    return x,y

def get_masses(filename):
    
    masses = []
    
    for line in open(filename):
        values = line.split(',')
        masses.append(float(values[0]))
    return masses

file_pattern = "spec/output/nbody_sim*"
files = sorted(glob.glob(file_pattern))

x_min = 0
x_max = 0
y_min = 0
y_max = 0

for f in files:
    x,y = read_file(f)
    x_min = min(min(x),x_min)
    y_min = min(min(y),y_min)

    x_max = max(max(x),x_max)
    y_max = max(max(y),y_max)

width = x_max - x_min
height = y_max - y_min
print(str([x_min,y_min,width,height ]))
print("found: " + str(files))
    
masses = get_masses(files[0])
min_mass = min(masses)
for i in range(0,len(masses)):
    masses[i] = mass_multiplicator*masses[i]/min_mass
               
x,y = read_file(files.pop(0))

fig = plt.figure()

#ax = fig.add_axes([x_min,y_min,width,height ])
#ax = fig.add_axes([0,-100,20,200])
scat = plt.scatter(x, y, c='blue',label='start', edgecolors='none') #,s=masses)
plt.xlim(x_min - abs(0.1*x_min),x_max + abs(0.1*x_max))
plt.ylim(y_min - abs(0.1*y_min),y_max + abs(0.1*y_max))
plt.draw()

#plt.show()

# axes = plt.gca()
# axes.set_xlim([x_min,x_max])
# axes.set_ylim([y_min,y_max])


print plt.gca()
scat.set_offset_position('data')
print(scat.get_offsets()[0].dtype)

def update(frame_num):
    if len(files) > 0:
        cur_file = files.pop(0)
        print("showing: " + cur_file)
    
        x,y = read_file(cur_file)
        fig.clear()
        scat = plt.scatter(x, y, c='blue',label='start', s=masses, edgecolors='none')
        plt.xlim(x_min - abs(0.1*x_min),x_max + abs(0.1*x_max))
        plt.ylim(y_min - abs(0.1*y_min),y_max + abs(0.1*y_max))
        plt.draw()
        # newpos = np.array([x,y],np.float64)
        # scat.set_offset_position('data')
        # scat.set_offsets(newpos)

animation = FuncAnimation(fig,update,interval=int(1000/speed))    

#time.sleep(1)
plt.show()
