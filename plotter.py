import numpy
import sys
import matplotlib.pyplot as plt
import matplotlib.animation as animation

NUM_DIMENSIONS = 2

if len(sys.argv) != 4:
	print("Usage: python plotter.py <constants filename> <input filename (binary)> <output filename (pngs)>")
	sys.exit(1)

constants_file = open(sys.argv[1], "r")
input_file = open(sys.argv[2], "rb")
output_filename = sys.argv[3]

# Read constants
constants = {}
for line in constants_file:
	line = line.strip().split("=")
	constants[line[0]] = line[1]

num_particles = int(constants["num_particles"])
num_steps = int(constants["num_steps"])
box_size = float(constants["box_size"])
write_interval = int(constants["write_interval"])

num_frames = num_steps // write_interval

# Read input
# The input file is a binary file with the following format:
# num_steps sets of NUM_DIMENSIONS sets of num_particles doubles

# Read all the data
data = numpy.fromfile(input_file, dtype=numpy.float64)

# Reshape the data
data = data.reshape(num_frames, NUM_DIMENSIONS, num_particles).tolist()

# Plot each step and save to different file
for frame in range(num_frames):
	step = frame * write_interval
	x_arr = data[frame][0]
	y_arr = data[frame][1]
	plt.scatter(x_arr, y_arr)
	# plt.xlim(0, box_size)
	# plt.ylim(0, box_size)
	plt.savefig("plots/%s_%03d.png" % (output_filename, step))
	plt.clf()