import os
import numpy as np
import matplotlib.pyplot as plt
import gudhi as gd
import yaml
from tqdm import tqdm

def read_particle_data(file_path, num_particles, number_of_time_steps, include_apc=True, skip_init_frame=False):
  """
  Read the particle data from a target file.

  Args:
      file_path: Res particle data target file
      num_particles: Number of particles in the simulation including the apc
      number_of_time_steps: Number of time steps in a simulation typically given by t_final / write_every
      include_apc: Default True. Will add the apc to each time snapshot

  Returns:
      3D array of shape (number_of_time_steps (+ 1), num_particles, NUM_DATA_PER_PARTICLE)
  """
  NUM_DATA_PER_PARTICLE = 6
  time_step_array = np.empty((number_of_time_steps + 1, num_particles, NUM_DATA_PER_PARTICLE))
    # each entry holds a num_particles by NUM_DATA_PER_PARTICLE shaped array
    # which in turn holds all the particle information
  if not include_apc:
    time_step_array = np.empty((number_of_time_steps + 1, num_particles - 1, NUM_DATA_PER_PARTICLE))
  time_step_labels = []
  apc_coord_and_radius = np.empty(3)
  apc_cell_recorded = False

  # Array to hold the coordinate, velocity, and radius information

  index_iter = 0

  if not include_apc:
    with open(file_path, 'r') as file:
      current_particle_index = 0
      for line in file:
        line_split = line.strip().split()

        # Start processing line
        if line[0] == 't':
          if current_particle_index != 0:
            index_iter += 1
          current_particle_index = 0
          time_step_labels.append(float(line_split[1]))
        elif not apc_cell_recorded and line[0] == '1':
          apc_cell_recorded = True
          apc_coord_and_radius[0] = float(line_split[1]) # x
          apc_coord_and_radius[1] = float(line_split[2]) # y
          apc_coord_and_radius[2] = float(line_split[5]) # radius
        elif line_split[0] == '1':
          continue
        else:
          time_step_array[index_iter][current_particle_index][0] = float(line_split[0])
          time_step_array[index_iter][current_particle_index][1] = float(line_split[1])
          time_step_array[index_iter][current_particle_index][2] = float(line_split[2])
          time_step_array[index_iter][current_particle_index][3] = float(line_split[3])
          time_step_array[index_iter][current_particle_index][4] = float(line_split[4])
          time_step_array[index_iter][current_particle_index][5] = float(line_split[5])
          current_particle_index += 1
  else : # include_apc = True
    with open(file_path, 'r') as file:
      current_particle_index = 0
      for line in file:
        line_split = line.strip().split()

        # Start processing line
        if line[0] == 't':
          if current_particle_index != 0:
            index_iter += 1
          current_particle_index = 0
          time_step_labels.append(float(line_split[1]))
        else:
          time_step_array[index_iter][current_particle_index][0] = float(line_split[0])
          time_step_array[index_iter][current_particle_index][1] = float(line_split[1])
          time_step_array[index_iter][current_particle_index][2] = float(line_split[2])
          time_step_array[index_iter][current_particle_index][3] = float(line_split[3])
          time_step_array[index_iter][current_particle_index][4] = float(line_split[4])
          time_step_array[index_iter][current_particle_index][5] = float(line_split[5])
          current_particle_index += 1
  return time_step_array

def calculate_alpha_complex_pd(coordinates):
  """
  Calculate the birth death pairs of dimension 0 and dimension 1 alpha complexes.

  Args:
      coordinates: An array where each element is an array containing an x and y coordinate.

  Returns:
      The birth, death pairs of dimension 0 and dimension 1
      dimension0 [birth, death], dimension1 [birth, death]
  """
  simplex_tree = gd.AlphaComplex(coordinates).create_simplex_tree()
  diagram = simplex_tree.persistence()

  dim0 = []
  dim1 = []
  for pair in diagram:
    (dim, (birth, death)) = pair
    if (death == float("inf")):
      continue
    if dim == 0:
      dim0.append(np.array([birth, death]))
    elif dim == 1:
      dim1.append(np.array([birth, death]))
  return np.array(dim0), np.array(dim1)

# Input one coordinate or list of coordinates
def dist_to_yx(coordinate):
  return abs(coordinate[0] - coordinate[1]) / np.sqrt(2)

# Averaged distance for each dimensional measure
def dist_measure_one_frame(dim0, dim1):
  dim0_measure = np.mean(dist_to_yx(dim0))
  dim1_measure = np.mean(dist_to_yx(dim1))
  return dim0_measure, dim1_measure


# One frame/timestamp via cell points
def one_frame_cell_points(data_for_one_frame):
  coord_data = data_for_one_frame[:, 1:3]
  dim0, dim1 = calculate_alpha_complex_pd(coord_data)
  # print(dim0.shape)
  # print(dim1.shape)
  dim0_dist_measure, dim1_dist_measure = dist_measure_one_frame(dim0, dim1)
  return dim0_dist_measure, dim1_dist_measure

# Average measure of one realization
def one_realization(data_with_all_frames):
  dim0_dist_measures = []
  dim1_dist_measures = []
  debug_count = 0
  for i in range(len(data_with_all_frames)):
    dim0_dist_measure, dim1_dist_measure = one_frame_cell_points(data_with_all_frames[i])
    dim0_dist_measures.append(dim0_dist_measure)
    dim1_dist_measures.append(dim1_dist_measure)
    debug_count += 1


  return np.mean(dim0_dist_measures), np.mean(dim1_dist_measures)

def run():
  NUM_FREQUENCIES = 14
  NUM_TIME_STEPS = 30   # t_final / (dt * write_every)
  NUM_PARTICLES = 300   # Particles including APC
  NUM_PARTICLE_DATA = 6
  NUM_REALIZATIONS = 10
  INCLUDE_APC = True
  SKIP_INIT_FRAME = True

  # Save data into memory (not used currently)
  configs = np.empty(NUM_FREQUENCIES, dtype='object')

  # Save list of frequencies for graphing
  frequency_list = np.empty(NUM_FREQUENCIES)

  # Measures
  dist_measure_dim0_list = np.empty(NUM_FREQUENCIES)
  dist_measure_dim1_list = np.empty(NUM_FREQUENCIES)

  # Runtime progress bar
  total_realizations = NUM_FREQUENCIES * NUM_REALIZATIONS
  rbar = tqdm(total=total_realizations)
  
  # Get yaml configs
  dirname = os.path.dirname(__file__)
  data_dir = ""
  with open(os.path.join(dirname, 'analysis_config.yaml')) as f:
    analysis_config = yaml.safe_load(f)
    data_dir = analysis_config['experiment_directory']
  data_dir = os.path.abspath(os.path.join(dirname, data_dir))

  # Iterator
  frequency_iterator = 0
  with os.scandir(data_dir) as experiments:
    for dir in experiments:

      # Accumulate means across realizations
      means_dim0_dist = np.empty(NUM_REALIZATIONS)
      means_dim1_dist = np.empty(NUM_REALIZATIONS)
      realization_iterator = 0
      is_yaml_retrieved = False
      with os.scandir(dir) as realizations:
        for sample in realizations:

          if not is_yaml_retrieved:
            is_yaml_retrieved = True
            with open(os.path.join(data_dir, dir.name, sample.name, 'config.yaml')) as f:
              data = yaml.safe_load(f)
              frequency_list[frequency_iterator] = data['frequency']


          # positions[index_iterator] = read_particle_data(os.path.join(data_dir, dir.name, 'Results.txt'), NUM_PARTICLES, NUM_TIME_STEPS, include_apc=True)
          realization_data = read_particle_data(os.path.join(data_dir, dir.name, sample.name, 'Results.txt'), NUM_PARTICLES, NUM_TIME_STEPS, include_apc=True)
          dim0_dist_measure, dim1_dist_measure = one_realization(realization_data)
          means_dim0_dist[realization_iterator] = dim0_dist_measure
          means_dim1_dist[realization_iterator] = dim1_dist_measure

          realization_iterator += 1

          rbar.update(1)

      dist_measure_dim0_list[frequency_iterator] = np.mean(means_dim0_dist)
      dist_measure_dim1_list[frequency_iterator] = np.mean(means_dim1_dist)

      frequency_iterator+=1


  rbar.close()

  index_sort = np.argsort(frequency_list)

  plt.scatter(frequency_list[index_sort], dist_measure_dim0_list[index_sort])
  plt.xlabel('Frequency')
  plt.ylabel('Distance Measure')
  plt.title('Distance Measure Dimension 0 vs Frequency')
  plt.savefig(os.path.join(data_dir, 'dist_dim0_vs_freq.png'))

  plt.figure()
  plt.scatter(frequency_list[index_sort], dist_measure_dim1_list[index_sort])
  plt.xlabel('Frequency')
  plt.ylabel('Distance Measure')
  plt.title('Distance Measure Dimension 1 vs Frequency')
  plt.savefig(os.path.join(data_dir, 'dist_dim1_vs_freq.png'))