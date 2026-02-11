import numpy as np

def read_particle_data(file_path, num_particles, number_of_timestamps, include_apc=True, skip_init_frame=False):
  """
  Read the particle data from a target file.

  Args:
      file_path: Res particle data target file.
      num_particles: Number of particles in the simulation including the apc.
      number_of_timestamps: Number of time steps in a simulation typically given by t_final / write_every.
      include_apc: Default True. Will add the apc to each time snapshot.

  Returns:
      3D array of shape (number_of_timestamps (+ 1), num_particles, NUM_DATA_PER_PARTICLE).
  """
  NUM_DATA_PER_PARTICLE = 6
  timestamp_array = np.empty()
  
  true_num_particles = num_particles
  true_num_timestamps = number_of_timestamps
  
  # Set size of array based on whether to include apc and whether to skip initial frame
  if not include_apc:
    true_num_particles = num_particles - 1
  if not skip_init_frame:
    true_num_timestamps = number_of_timestamps + 1

  timestamp_array = np.empty((true_num_timestamps, true_num_particles, NUM_DATA_PER_PARTICLE))
  
  apc_coord_and_radius = np.empty(3)
  apc_cell_recorded = False

  # Array to hold the coordinate, velocity, and radius information

  timestamp_index = 0
  
  is_first_frame_skipped = not skip_init_frame

  if not include_apc:
    with open(file_path, 'r') as file:
      current_particle_index = 0
      for line in file:
        line_split = line.strip().split()

        # Start processing line
        if line[0] == 't':
          if not is_first_frame_skipped and not int(line_split[1]) == 0:
            is_first_frame_skipped = True
          if current_particle_index != 0:
            timestamp_index += 1
          current_particle_index = 0
        elif not apc_cell_recorded and line[0] == '1' and is_first_frame_skipped:
          apc_cell_recorded = True
          apc_coord_and_radius[0] = float(line_split[1]) # x
          apc_coord_and_radius[1] = float(line_split[2]) # y
          apc_coord_and_radius[2] = float(line_split[5]) # radius
        elif line_split[0] == '1':
          continue
        elif is_first_frame_skipped:
          timestamp_array[timestamp_index][current_particle_index][0] = float(line_split[0])
          timestamp_array[timestamp_index][current_particle_index][1] = float(line_split[1])
          timestamp_array[timestamp_index][current_particle_index][2] = float(line_split[2])
          timestamp_array[timestamp_index][current_particle_index][3] = float(line_split[3])
          timestamp_array[timestamp_index][current_particle_index][4] = float(line_split[4])
          timestamp_array[timestamp_index][current_particle_index][5] = float(line_split[5])
          current_particle_index += 1
  else : # include_apc = True
    with open(file_path, 'r') as file:
      current_particle_index = 0
      for line in file:
        line_split = line.strip().split()

        # Start processing line
        if line[0] == 't':
          if not is_first_frame_skipped and not int(line_split[1]) == 0:
            is_first_frame_skipped = True
            continue
          if current_particle_index != 0:
            timestamp_index += 1
          current_particle_index = 0
        elif is_first_frame_skipped:
          timestamp_array[timestamp_index][current_particle_index][0] = float(line_split[0])
          timestamp_array[timestamp_index][current_particle_index][1] = float(line_split[1])
          timestamp_array[timestamp_index][current_particle_index][2] = float(line_split[2])
          timestamp_array[timestamp_index][current_particle_index][3] = float(line_split[3])
          timestamp_array[timestamp_index][current_particle_index][4] = float(line_split[4])
          timestamp_array[timestamp_index][current_particle_index][5] = float(line_split[5])
          current_particle_index += 1
  return timestamp_array