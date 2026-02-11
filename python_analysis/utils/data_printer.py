import os

def print_realizaton(filename, measure_data, realization_analyses_dir, timestamp_difference, is_first_frame_skipped):
  """
  Print the distance measures for one dimension for each timestamp in a realization to a file. Each timestamp will be separated by a line of dashes.

  Args:
      filename: Name of the file to print the data to.
      measure_data: A list of lists. Each inner list contains the distance measures for one timestamp.
      realization_analyses_dir: The directory to print the file to.
      timestamp_difference: The time difference between each timestamp in the simulation.
      is_first_frame_skipped: Whether the initial frame of the simulation was skipped when reading in the data. If so, the first timestamp will be 0 + timestamp_difference instead of 0.
  """
  with open(os.path.join(realization_analyses_dir, filename), 'w') as f:
    for i in range(len(measure_data)):
      if is_first_frame_skipped:
        f.write(f'Timestamp: {(i+1) * timestamp_difference}\n')
      else:
        f.write(f'Timestamp: {i * timestamp_difference}\n')
      for measure in measure_data[i]:
        f.write(f'{measure}\n')