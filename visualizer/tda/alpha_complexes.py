import numpy as np
import gudhi as gd
import plotly.graph_objects as go


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

# One frame/timestamp
# Get dim0 and dim1 birth, death pairs
# Then calculate distance measures for each pair
def one_frame_dist_measures(coord_data):
  dim0, dim1 = calculate_alpha_complex_pd(coord_data)
  
  dim0_dist_measures = []
  for pair in dim0:
    dim0_dist_measures.append(dist_to_yx(pair))
  
  dim1_dist_measures = []
  for pair in dim1:
    dim1_dist_measures.append(dist_to_yx(pair))
  return dim0_dist_measures, dim1_dist_measures

# Distance measures for all frames/timestamps in one realization
def alpha_complexes_with_particle_coords(data_with_all_frames, realization_analyses_dir, timestamp_difference, is_first_frame_skipped):
  dim0_dist_measure_list = []
  dim1_dist_measure_list = []
  for i in range(len(data_with_all_frames)):
    dim0_dist_measure, dim1_dist_measure = one_frame_dist_measures(data_with_all_frames[i][:, 1:3])
    dim0_dist_measure_list.append(dim0_dist_measure)
    dim1_dist_measure_list.append(dim1_dist_measure)
  
def delunauy_plotly_visualization(coordinates):
    fig = go.Figure()

    # Add traces, one for each slider step
    for step in np.arange(0, 5, 0.1):
        fig.add_trace(
            go.Scatter(
                visible=False,
                line=dict(color="#00CED1", width=6),
                name="𝜈 = " + str(step),
                x=np.arange(0, 10, 0.01),
                y=np.sin(step * np.arange(0, 10, 0.01))))

    # Make 10th trace visible
    fig.data[10].visible = True

    # Create and add slider
    steps = []
    for i in range(len(fig.data)):
        step = dict(
            method="update",
            args=[{"visible": [False] * len(fig.data)},
                {"title": "Slider switched to step: " + str(i)}],  # layout attribute
        )
        step["args"][0]["visible"][i] = True  # Toggle i'th trace to "visible"
        steps.append(step)

    sliders = [dict(
        active=10,
        currentvalue={"prefix": "Frequency: "},
        pad={"t": 50},
        steps=steps
    )]

    fig.update_layout(
        sliders=sliders
    )

    return fig