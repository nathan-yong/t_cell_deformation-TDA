import numpy as np
import gudhi as gd

def calculate_alpha_complex_pd(coordinates):
  """
  Calculate the birth death pairs of dimension 0 and dimension 1 alpha complexes.

  Args:
      coordinates: An array where each element is an array containing an x and y coordinate.

  Returns: 
      dim0: An array of birth death pairs for dimension 0.
      dim1: An array of birth death pairs for dimension 1.
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