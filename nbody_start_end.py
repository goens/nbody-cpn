import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np 
import glob
import time
import matplotlib.cm as cmx
import matplotlib.colors as colors
from sys import argv

#todo:
#http://matplotlib.org/examples/mplot3d/scatter3d_demo.html


mass_multiplicator = 1
N = int(argv[1])

def read_file(filename):
    
    x = []
    y = []
    masses = []
    
    for line in open(filename):
        values = line.split(',')
        x.append(float(values[0]))
        y.append(float(values[1]))
        masses.append(float(values[3]))
    return x,y

def get_masses(filename):
    
    masses = []
    
    for line in open(filename):
        values = line.split(',')
        masses.append(float(values[3]))
    return masses

#http://stackoverflow.com/questions/14720331/how-to-generate-random-colors-in-matplotlib
def get_cmap(N):
    '''Returns a function that maps each index in 0, 1, ... N-1 to a distinct 
    RGB color.'''
    color_norm  = colors.Normalize(vmin=0, vmax=N-1)
    scalar_map = cmx.ScalarMappable(norm=color_norm, cmap='hsv') 
    def map_index_to_rgb_color(index):
        return scalar_map.to_rgba(index)
    return map_index_to_rgb_color

file_pattern = "spec/output/nbody_sim*"
files = sorted(glob.glob(file_pattern))
if N > len(files) or N <= 0:
    N = 1
    
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
               

fig = plt.figure()
colors = get_cmap(N)

step = int(len(files) / N)
j = 0

for i in range(0,len(files),step):
    x,y = read_file(files[i])
    plt.scatter(x, y, c=colors(j),label=files[i], s=masses, edgecolors='none')
    print(str(colors(j)))
    j = j+1


# axes = plt.gca()
# axes.set_xlim([x_min,x_max])
# axes.set_ylim([y_min,y_max])


plt.show()
