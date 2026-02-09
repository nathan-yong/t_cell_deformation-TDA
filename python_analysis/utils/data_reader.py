import numpy as np

def read_particle_data(file_path, num_particles, number_of_time_steps, include_apc=True, skip_init_frame=False):
  """
  Read the particle data from a target file.

  Args:
      file_path: Res particle data target file.
      num_particles: Number of particles in the simulation including the apc.
      number_of_time_steps: Number of time steps in a simulation typically given by t_final / write_every.
      include_apc: Default True. Will add the apc to each time snapshot.

  Returns:
      3D array of shape (number_of_time_steps (+ 1), num_particles, NUM_DATA_PER_PARTICLE).
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