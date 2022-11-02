import numpy
import sys
import numpy as np
from PIL import Image
import subprocess

NUM_DIMENSIONS = 2
IMAGE_SIZE = 300
DIR_NAME = "output_media"
VIDEO_DURATION = 10


def run_command(command: str, use_shell = False) -> tuple[str, str]:
	process = subprocess.Popen(command.split(" "), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=use_shell)
	stdout, stderr = process.communicate()
	# decode bytes to string
	stdout = stdout.decode("utf-8")
	stderr = stderr.decode("utf-8")
	if stderr:
		print(f"---ERROR from subprocess---\n{command}\n---------------------------\n\n{stderr}\n\n---End of ERROR from subprocess---\nExiting...")
		exit()
	return stdout, stderr

STAR_COLOURS = [
	[155, 176, 255],
	[170, 191, 255],
	[202, 215, 255],
	[248, 247, 255],
	[255, 244, 234],
	[255, 210, 161],
	[255, 204, 111]
]

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


# Create output folder
run_command(f"rm -r {DIR_NAME}")
run_command(f"mkdir {DIR_NAME}")

# Write each step and save to different file by coloring each image darker
for frame in range(num_frames):
	step = frame * write_interval
	x_arr = data[frame][0]
	y_arr = data[frame][1]
	# Find the centre of mass of the system
	x_centre = sum(x_arr) / num_particles
	y_centre = sum(y_arr) / num_particles
	# Move the centre of mass to (0, 0)
	x_arr = [x - x_centre for x in x_arr]
	y_arr = [y - y_centre for y in y_arr]
	# Find the maximum distance from the centre of mass (0, 0)
	max_dist = max([np.sqrt(x ** 2 + y ** 2) for x, y in zip(x_arr, y_arr)])
	# Bring all particles closer to the centre of mass so that the maximum distance is box_size / 2
	x_arr = [x * box_size / 2 / max_dist for x in x_arr]
	y_arr = [y * box_size / 2 / max_dist for y in y_arr]
	# Move the centre of mass to the centre of the image
	offset = box_size / 2
	x_arr = [x + offset for x in x_arr]
	y_arr = [y + offset for y in y_arr]

	# Create the image
	# Initialise the image and conver it to a numpy array
	image = Image.new("RGB", (IMAGE_SIZE, IMAGE_SIZE))
	image = np.array(image)
	image.fill(0)
	# Draw the particles
	for i in range(num_particles):
		x = x_arr[i]
		y = y_arr[i]
		# Convert to pixel coordinates
		x = int((x / box_size) * IMAGE_SIZE)
		y = int((y / box_size) * IMAGE_SIZE)
		# Check if the particle is in the image
		if x >= IMAGE_SIZE or y >= IMAGE_SIZE or x < 0 or y < 0:
			continue
		# Draw the particle
		image[x, y] = STAR_COLOURS[i % len(STAR_COLOURS)]
	# Convert the numpy array back to an image
	image = Image.fromarray(image)
	# Save the image
	image.save(f"{DIR_NAME}/{output_filename}_{step}.png")

# Turn images into mp4
run_command(f"ffmpeg -framerate {num_frames / VIDEO_DURATION} -pattern_type glob -i '{DIR_NAME}/{output_filename}_*.png' -c:v libx264 -pix_fmt yuv420p {DIR_NAME}/{output_filename}.mp4")
# Remove images
run_command("rm output/figures/" + output_filename + "*.png")

# Close the files
constants_file.close()
input_file.close()



