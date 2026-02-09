import os
import numpy as np
from tqdm import tqdm
import yaml

# Local imports
from utils.data_reader import read_particle_data
from tda.alpha_complexes import one_realization

def main():
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
    



if __name__ == "__main__":
    main()