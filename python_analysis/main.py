import os
import numpy as np
from tqdm import tqdm
import yaml

# Local imports
from utils.data_reader import read_particle_data
from tda.alpha_complexes import one_realization

def main():
  NUM_EXPERIMENT_VARIABLES = 0
  NUM_TIMESTAMPS = 0   # t_final / (dt * write_every)
  NUM_PARTICLES = 0   # Particles including APC
  NUM_REALIZATIONS = 0
  INCLUDE_APC = True
  SKIP_INIT_FRAME = True
  
  # Get yaml configs
  analysis_config = {}
  dirname = os.path.dirname(__file__)
  data_dir = ""
  with open(os.path.join(dirname, 'analysis_config.yaml')) as f:
    analysis_config = yaml.safe_load(f)
    try:
      NUM_EXPERIMENT_VARIABLES = analysis_config['number_of_experiment_variables']
      NUM_TIMESTAMPS = analysis_config['number_of_timestampes']
      NUM_PARTICLES = analysis_config['number_of_particles']
      NUM_REALIZATIONS = analysis_config['number_of_realizations']
      INCLUDE_APC = analysis_config['include_apc']
      SKIP_INIT_FRAME = analysis_config['skip_init_frame']
      data_dir = analysis_config['target_experiment_directory']
    except KeyError as e:
      print(f"Missing key in analysis_config.yaml: {e}")
      return
  data_dir = os.path.abspath(os.path.join(dirname, data_dir))

  # Save list of frequencies for graphing
  frequency_list = np.empty(NUM_EXPERIMENT_VARIABLES)

  # Measures
  dist_measure_dim0_list = np.empty(NUM_EXPERIMENT_VARIABLES)
  dist_measure_dim1_list = np.empty(NUM_EXPERIMENT_VARIABLES)

  # Runtime progress bar
  total_realizations = NUM_EXPERIMENT_VARIABLES * NUM_REALIZATIONS
  rbar = tqdm(total=total_realizations)
  
  

  # Iterator
  frequency_iterator = 0
  
  experiment_analyses_dir = os.path.join(data_dir, 'analyses')
  if not os.path.exists(experiment_analyses_dir):
      os.makedirs(experiment_analyses_dir)
  
  with os.scandir(data_dir) as experiment:
    for experiment_sample in experiment:

      # Accumulate means across realizations
      means_dim0_dist = np.empty(NUM_REALIZATIONS)
      means_dim1_dist = np.empty(NUM_REALIZATIONS)
      realization_iterator = 0
      is_yaml_retrieved = False
      
      experiment_sample_analyses_dir = os.path.join(data_dir, experiment_sample.name, 'analyses')
      if not os.path.exists(experiment_sample_analyses_dir):
          os.makedirs(experiment_sample_analyses_dir)
          
      
      with os.scandir(experiment_sample) as realizations:
        for realization_sample in realizations:
          # Create analyses directory if it doesn't exist
          realization_analyses_dir = os.path.join(data_dir, experiment_sample.name, realization_sample.name, 'analyses')
          if not os.path.exists(realization_analyses_dir):
              os.makedirs(realization_analyses_dir)
          
          if not is_yaml_retrieved:
            is_yaml_retrieved = True
            with open(os.path.join(data_dir, experiment_sample.name, realization_sample.name, 'config.yaml')) as f:
              data = yaml.safe_load(f)
              frequency_list[frequency_iterator] = data['frequency']

          realization_data = read_particle_data(os.path.join(data_dir, experiment_sample.name, realization_sample.name, 'Results.txt'), 
                                                NUM_PARTICLES, NUM_TIMESTAMPS, include_apc=INCLUDE_APC, 
                                                SKIP_INIT_FRAME=SKIP_INIT_FRAME)
          if analysis_config['alpha_complexes_with_particle_coords']:
            # Function should return one np array of size (NUM_TIMESTAMPS, 2) for each dimension.
            dim0_dist_measure, dim1_dist_measure = one_realization(realization_data) # TODO: convert function
            # print data into file
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