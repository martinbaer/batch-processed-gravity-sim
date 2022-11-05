import numpy
import sys
import pygame
import subprocess
import numpy as np

NUM_DIMENSIONS = 2
IMAGE_SIZE = 500
DIR_NAME = "output_media"
VIDEO_DURATION = 30
OUTPUT_NAME = "simulation"

def run_command(command: str, use_shell = False, ignore_error = False) -> tuple[str, str]:
	"""
	Runs a command and returns the output and error
	"""
	process = subprocess.Popen(command.split(" "), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=use_shell)
	stdout, stderr = process.communicate()
	# decode bytes to string
	stdout = stdout.decode("utf-8")
	stderr = stderr.decode("utf-8")
	if stderr and not ignore_error:
		print(f"---ERROR from subprocess---\n{command}\n----------STDERR-----------\n{stderr}\n----------STDOUT-----------\n{stdout}\n---End of ERROR from subprocess---\nExiting...")
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

if len(sys.argv) != 2:
	print("Usage: python plotter.py <constants filename>")
	sys.exit(1)

constants_file = open(sys.argv[1], "r")

# Read constants
constants = {}
for line in constants_file:
	line = line.strip().split("=")
	constants[line[0]] = line[1]

num_particles = int(constants["num_particles"])
num_steps = int(constants["num_steps"])
box_size = float(constants["box_size"])
write_interval = int(constants["write_interval"])
input_file = open(constants["output_filename"])

num_frames = num_steps // write_interval

# Read input
# The input file is a binary file with the following format:
# num_steps sets of NUM_DIMENSIONS sets of num_particles doubles

# Read all the data
data = numpy.fromfile(input_file, dtype=numpy.float64)

# Reshape the data
data = data.reshape(num_frames, NUM_DIMENSIONS, num_particles).tolist()


# Create output folder
run_command(f"rm -r {DIR_NAME}", ignore_error=True)
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
	# Find the 99% confidence interval from the centre of mass (0, 0)
	distances = [np.sqrt(x ** 2 + y ** 2) for x, y in zip(x_arr, y_arr)]
	distances.sort()
	max_dist = distances[int(num_particles * 0.99)]

	# Bring all particles closer to the centre of mass so that the maximum distance is box_size / 2
	x_arr = [x * box_size / 2 / max_dist for x in x_arr]
	y_arr = [y * box_size / 2 / max_dist for y in y_arr]
	# Move the centre of mass to the centre of the image
	offset = box_size / 2
	x_arr = [x + offset for x in x_arr]
	y_arr = [y + offset for y in y_arr]

	# Create a black pygame image
	surface = pygame.Surface((IMAGE_SIZE, IMAGE_SIZE))
	surface.fill((0, 0, 0))

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
		surface.set_at((x, y), STAR_COLOURS[i % len(STAR_COLOURS)])
	# Save the image
	pygame.image.save(surface, f"{DIR_NAME}/{OUTPUT_NAME}_{str(step).zfill(len(str(num_frames * write_interval)))}.png")


# Turn images into mp4
# run_command(f"ffmpeg -framerate {num_frames / VIDEO_DURATION} -pattern_type glob -i '{DIR_NAME}/{output_filename}_*.png' -c:v libx264 -pix_fmt yuv420p {DIR_NAME}/{output_filename}.mp4", use_shell=True)
run_command(f"sh .utils/ffmpeg.sh {num_frames // VIDEO_DURATION} {DIR_NAME}/{OUTPUT_NAME}.mp4", ignore_error=True)

# Remove images
for frame in range(num_frames):
	step = frame * write_interval
	run_command(f'rm {DIR_NAME}/{OUTPUT_NAME}_{str(step).zfill(len(str(num_frames * write_interval)))}.png')

# Close the files
constants_file.close()
input_file.close()



