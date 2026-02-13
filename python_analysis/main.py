import os
import numpy as np
import yaml
import re

# Local imports
import utils.progress_bar as progress_bar
import utils.data_reader as data_reader
import utils.data_printer as data_printer
from tda.alpha_complexes import alpha_complexes_with_particle_coords


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
      NUM_TIMESTAMPS = analysis_config['number_of_timestamps']
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

  num_analyses = 0
  if analysis_config['alpha_complexes_with_particle_coords']:
    num_analyses += 1

  true_frame_count = NUM_TIMESTAMPS
  if not SKIP_INIT_FRAME:
    true_frame_count += 1

  # Runtime progress bar
  total = NUM_EXPERIMENT_VARIABLES * NUM_REALIZATIONS * true_frame_count * num_analyses
  progress_bar.start_progress_bar(total)

  # Iterator
  frequency_iterator = 0
  
  experiment_analyses_dir = os.path.join(data_dir, 'analyses')
  if not os.path.exists(experiment_analyses_dir):
      os.makedirs(experiment_analyses_dir)
  
  with os.scandir(data_dir) as experiment:
    
    is_first_experiment_sample_print = True
    # LOOP: experiment sample
    for experiment_sample in experiment:
      if experiment_sample.name == 'analyses':
        continue

      # Accumulate means across realizations
      realization_iterator = 0
      is_yaml_retrieved = False
      experiment_sample_analyses_dir = os.path.join(data_dir, experiment_sample.name, 'analyses')
      
      if not os.path.exists(experiment_sample_analyses_dir):
          os.makedirs(experiment_sample_analyses_dir)
          
      
      with os.scandir(experiment_sample) as realizations:
        if experiment_sample.name == 'analyses':
          continue
        
        is_first_realization_print = True
        # LOOP: realizations
        for realization_sample in realizations:
          if realization_sample.name == 'analyses':
            continue 
          # Create analyses directory if it doesn't exist
          realization_analyses_dir = os.path.join(data_dir, experiment_sample.name, realization_sample.name, 'analyses')
          if not os.path.exists(realization_analyses_dir):
              os.makedirs(realization_analyses_dir)

          realization_data = data_reader.read_particle_data(os.path.join(data_dir, experiment_sample.name, realization_sample.name, 'Results.txt'), 
                                                            NUM_PARTICLES, 
                                                            NUM_TIMESTAMPS, 
                                                            include_apc=INCLUDE_APC, 
                                                            skip_init_frame=SKIP_INIT_FRAME)
          
          # Deploy all methods selected in the analysis_config
          if analysis_config['alpha_complexes_with_particle_coords']:
            alpha_complexes_with_particle_coords(realization_data, realization_analyses_dir, analysis_config['timestamp_difference'], SKIP_INIT_FRAME)
            # print to experiment_sample_analyses_dir. Copy the files created by function.
            data_printer.print_experiment_sample(os.path.join(realization_analyses_dir, "alpha_complexes_with_particle_coords_dim0.txt"), 
                                                experiment_sample_analyses_dir, 
                                                "freq" + str(experiment_sample.name) + "_alpha_complexes_with_particle_coords_dim0.txt",
                                                realization_sample.name, is_first_realization_print)
            data_printer.print_experiment_sample(os.path.join(realization_analyses_dir, "alpha_complexes_with_particle_coords_dim1.txt"), 
                                                experiment_sample_analyses_dir, 
                                                "freq" + str(experiment_sample.name) + "_alpha_complexes_with_particle_coords_dim1.txt",
                                                realization_sample.name, is_first_realization_print)
            is_first_realization_print = False
          realization_iterator += 1

      # print to experiment_analyses_dir. Copy realization file.
      data_printer.print_experiment(os.path.join(experiment_sample_analyses_dir, f"freq{experiment_sample.name}_alpha_complexes_with_particle_coords_dim0.txt"), 
                                    experiment_analyses_dir, 
                                    f"{analysis_config['experiment_name']}_alpha_complexes_with_particle_coords_dim0.txt",
                                    re.findall(r'\d+', experiment_sample.name)[0], is_first_experiment_sample_print)
      data_printer.print_experiment(os.path.join(experiment_sample_analyses_dir, f"freq{experiment_sample.name}_alpha_complexes_with_particle_coords_dim1.txt"), 
                                    experiment_analyses_dir, 
                                    f"{analysis_config['experiment_name']}_alpha_complexes_with_particle_coords_dim1.txt",
                                    re.findall(r'\d+', experiment_sample.name)[0], is_first_experiment_sample_print)

      is_first_experiment_sample_print = False
      frequency_iterator+=1
  
  progress_bar.stop_progress_bar()



if __name__ == "__main__":
    main()