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
      dim0.append([birth, death])
    elif dim == 1:
      dim1.append([birth, death])
  return dim0, dim1

# Input one coordinate or list of coordinates
def dist_to_yx(coordinate):
  return abs(coordinate[0] - coordinate[1]) / np.sqrt(2)

# One frame/timestamp via cell points
# Get dim0 and dim1 birth, death pairs
# Then calculate distance measures for each pair
def one_frame_dist_measures(coord_data):
  dim0, dim1 = calculate_alpha_complex_pd(coord_data)
  dim0_dist_measures = dist_to_yx(dim0)
  dim1_dist_measures = dist_to_yx(dim1)
  return dim0_dist_measures, dim1_dist_measures

# Distance measures for all frames/timestamps in one realization
def alpha_complexes_with_particle_coords(data_with_all_frames):
  for i in range(len(data_with_all_frames)):
    dim0_dist_measure, dim1_dist_measure = one_frame_dist_measures(data_with_all_frames[i][:, 1:3])
    # print out data into file to separate timestamps
