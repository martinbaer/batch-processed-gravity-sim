import numpy
import sys
import pygame
import subprocess
import numpy as np
import os

NUM_DIMENSIONS = 2
IMAGE_SIZE = 500
VIDEO_DURATION = 30
OUTPUT_DIR = "output"
TEMP_DIR = "temp_pngs"




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
write_interval = int(constants["write_interval"])
input_filename = constants["output_filename"]
input_file = open(constants["output_filename"])

num_frames = num_steps // write_interval

# Read input
# The input file is a binary file with the following format:
# num_steps sets of NUM_DIMENSIONS sets of num_particles doubles

# Read all the data
data = numpy.fromfile(input_file, dtype=numpy.float64)

# Find the centre of mass of the entire simulation (every particle at every step)
centre_of_mass_x = 0
centre_of_mass_y = 0
for frame in range(num_frames):
	for i in range(num_particles):
		centre_of_mass_x += data[frame * 2 * num_particles + i]
		centre_of_mass_y += data[frame * 2 * num_particles + num_particles + i]
centre_of_mass_x /= len(data) / 2
centre_of_mass_y /= len(data) / 2

# Move the centre of mass to the centre of the screen
for frame in range(num_frames):
	for i in range(num_particles):
		data[frame * 2 * num_particles + i] -= centre_of_mass_x
		data[frame * 2 * num_particles + num_particles + i] -= centre_of_mass_y

# Find the distances from the centre of mass
distances = []
# for element in data:
# 	distances.append(np.sqrt((element - centre_of_mass_x)**2 + (element - centre_of_mass_y)**2))
for frame in range(num_frames):
	for i in range(num_particles):
		distances.append(np.sqrt((data[frame * 2 * num_particles + i] - centre_of_mass_x)**2 + (data[frame * 2 * num_particles + num_particles + i] - centre_of_mass_y)**2))

# Find the 99% percentile interval from the centre of mass (0, 0) and use this to scale the image
view_radius = np.percentile(distances, 95)
scale_factor = IMAGE_SIZE / (2 * view_radius)

# Scale the data to fit the screen and offset it to the centre
for frame in range(num_frames):
	for i in range(num_particles):
		data[frame * 2 * num_particles + i] *= scale_factor
		data[frame * 2 * num_particles + i] += IMAGE_SIZE / 2
		data[frame * 2 * num_particles + num_particles + i] *= scale_factor
		data[frame * 2 * num_particles + num_particles + i] += IMAGE_SIZE / 2

# Reshape the data
data = data.reshape(num_frames, NUM_DIMENSIONS, num_particles).tolist()


# Create output folder
run_command(f"rm -r {TEMP_DIR}", ignore_error=True)
run_command(f"mkdir {TEMP_DIR}")

output_filename = f"{os.path.splitext(input_filename)[0]}.mp4"
run_command(f"rm {output_filename}", ignore_error=True)

digits_needed = len(str(num_frames * write_interval))

# Write each step and save to different file by coloring each image darker
for frame in range(num_frames):
	step = frame * write_interval
	x_arr = data[frame][0]
	y_arr = data[frame][1]

	# Create a black pygame image
	surface = pygame.Surface((IMAGE_SIZE, IMAGE_SIZE))
	surface.fill((0, 0, 0))

	# Draw the particles
	for i in range(num_particles):
		x = x_arr[i]
		y = y_arr[i]
		# Convert to pixel coordinates
		x = int(x)
		y = int(y)
		# Check if the particle is in the image
		if x >= IMAGE_SIZE or y >= IMAGE_SIZE or x < 0 or y < 0:
			continue
		# Draw the particle
		surface.set_at((x, y), STAR_COLOURS[i % len(STAR_COLOURS)])
	# Save the image
	pygame.image.save(surface, f"{TEMP_DIR}/{str(step).zfill(digits_needed)}.png")

# Turn images into mp4
# run_command(f"ffmpeg -framerate {num_frames / VIDEO_DURATION} -pattern_type glob -i '{DIR_NAME}/{output_filename}_*.png' -c:v libx264 -pix_fmt yuv420p {DIR_NAME}/{output_filename}.mp4", use_shell=True)

run_command(f"sh ./shell_scripts/ffmpeg.sh {num_frames // VIDEO_DURATION} {output_filename}", ignore_error=True)

# Remove images
run_command(f"rm -r {TEMP_DIR}", ignore_error=True)

# Close the files
constants_file.close()
input_file.close()



