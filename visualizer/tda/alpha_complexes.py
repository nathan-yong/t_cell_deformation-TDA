import numpy as np
import gudhi as gd
import plotly.graph_objects as go
import networkx as nx


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
  return np.array(dim0), np.array(dim1)

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
  
def delaunay_plotly_visualization(coordinates):
    # 1. Compute Alpha Complex
    alpha_complex = gd.AlphaComplex(points=coordinates)
    simplex_tree = alpha_complex.create_simplex_tree()

    # Compute persistence
    persistence = simplex_tree.persistence()

    # Extract birth and death times for dimension 1
    # persistence is a list of (dim, (birth, death))
    dim1_birth_times = [p[1][0] for p in persistence if p[0] == 1]
    dim1_death_times = [p[1][1] for p in persistence if p[0] == 1 and p[1][1] != float('inf')]

    # 2. Extract Simplices
    edges = []
    triangles = []
    filtrations = set()

    for simplex, filtration in simplex_tree.get_filtration():
        filtrations.add(filtration)
        if len(simplex) == 2:
            edges.append((simplex, filtration))
        elif len(simplex) == 3:
            triangles.append((simplex, filtration))

    sorted_filtrations = sorted(list(filtrations))
    if 0.0 not in sorted_filtrations:
        sorted_filtrations.insert(0, 0.0)

    alphas = [np.sqrt(f) for f in sorted_filtrations]

    # 3. Create Figure
    fig = go.Figure()

    # Trace 0: Points (always visible)
    fig.add_trace(go.Scatter(
        x=coordinates[:, 0],
        y=coordinates[:, 1],
        mode='markers',
        marker=dict(size=5, color='black'),
        name='Points'
    ))

    # Trace 1: Alpha Complex Edges
    fig.add_trace(go.Scatter(
        x=[],
        y=[],
        mode='lines',
        line=dict(color='blue', width=1),
        name='Alpha Complex'
    ))

    # Trace 2: Alpha Radius Circles
    fig.add_trace(go.Scatter(
        x=[],
        y=[],
        mode='lines',
        line=dict(color='rgba(255, 0, 0, 0.3)', width=1),
        name='Alpha Radius'
    ))

    # Trace 3: Birth Cycles (New Feature)
    fig.add_trace(go.Scatter(
        x=[],
        y=[],
        mode='lines',
        line=dict(color='red', width=4),
        name='New Cycle (Birth)'
    ))

    # Trace 4: Death Triangles (Feature Death)
    fig.add_trace(go.Scatter(
        x=[],
        y=[],
        mode='lines',
        line=dict(color='green', width=4),
        name='Filled Hole (Death)'
    ))

    # Precompute circle unit coordinates
    theta = np.linspace(0, 2 * np.pi, 40)
    unit_circle_x = np.cos(theta)
    unit_circle_y = np.sin(theta)

    # 4. Create Frames
    frames = []

    # Initialize graph for pathfinding (used for Births)
    G = nx.Graph()
    G.add_nodes_from(range(len(coordinates)))

    # Helper to check if a value is in a list with tolerance
    def is_in_list(val, target_list):
        return any(np.isclose(val, b) for b in target_list)

    # Process edges incrementally
    edges.sort(key=lambda x: x[1])
    current_edge_idx = 0

    for i, (alpha, filt_val) in enumerate(zip(alphas, sorted_filtrations)):

        # --- 1. Process Edges & Births ---
        step_birth_x = []
        step_birth_y = []

        while current_edge_idx < len(edges) and edges[current_edge_idx][1] <= filt_val:
            edge_simplex, edge_filt = edges[current_edge_idx]
            u, v = edge_simplex

            # Check for Birth: Is this edge creating a cycle at a birth time?
            if is_in_list(edge_filt, dim1_birth_times):
                if nx.has_path(G, u, v):
                    path = nx.shortest_path(G, u, v)
                    cycle_nodes = path + [u]
                    cx = coordinates[cycle_nodes, 0]
                    cy = coordinates[cycle_nodes, 1]
                    step_birth_x.extend(cx)
                    step_birth_x.append(None)
                    step_birth_y.extend(cy)
                    step_birth_y.append(None)

            G.add_edge(u, v)
            current_edge_idx += 1

        # --- 2. Process Deaths (Triangles) ---
        step_death_x = []
        step_death_y = []

        # Check triangles that appear exactly at this filtration step
        # (We could optimize by sorting triangles, but N is small enough)
        for tri_simplex, tri_filt in triangles:
            if np.isclose(tri_filt, filt_val):
                # Is this triangle killing a feature?
                if is_in_list(tri_filt, dim1_death_times):
                    # Highlight the triangle (3 edges)
                    u, v, w = tri_simplex
                    # Draw u->v->w->u
                    tri_nodes = [u, v, w, u]
                    tx = coordinates[tri_nodes, 0]
                    ty = coordinates[tri_nodes, 1]
                    step_death_x.extend(tx)
                    step_death_x.append(None)
                    step_death_y.extend(ty)
                    step_death_y.append(None)

        # --- 3. Construct Frame Data ---

        # A. Wireframe (Edges)
        vis_edge_x = []
        vis_edge_y = []
        # Re-scan edges currently in graph for visualization
        # (Scanning up to current_edge_idx is safe because we sorted them)
        for k in range(current_edge_idx):
            p0 = coordinates[edges[k][0][0]]
            p1 = coordinates[edges[k][0][1]]
            vis_edge_x.extend([p0[0], p1[0], None])
            vis_edge_y.extend([p0[1], p1[1], None])

        # B. Circles
        circ_x = []
        circ_y = []
        for point in coordinates:
            cx = point[0] + alpha * unit_circle_x
            cy = point[1] + alpha * unit_circle_y
            circ_x.extend(cx)
            circ_x.append(None)
            circ_y.extend(cy)
            circ_y.append(None)

        frames.append(go.Frame(
            data=[
                go.Scatter(x=vis_edge_x, y=vis_edge_y),
                go.Scatter(x=circ_x, y=circ_y),
                go.Scatter(x=step_birth_x, y=step_birth_y),
                go.Scatter(x=step_death_x, y=step_death_y)
            ],
            traces=[1, 2, 3, 4],
            name=str(i)
        ))

    fig.frames = frames

    # 5. Configure Slider
    sliders = [dict(
        active=0,
        currentvalue={"prefix": "Alpha (Radius): "},
        pad={"t": 50},
        steps=[dict(
            method="animate",
            args=[[str(i)], dict(mode="immediate", frame=dict(duration=0, redraw=True), transition=dict(duration=0))],
            label=f"{a:.3f}"
        ) for i, a in enumerate(alphas)]
    )]

    fig.update_layout(
        sliders=sliders,
        title="Interactive Alpha Complex: Births (Red) & Deaths (Green)",
        xaxis_title="X Position",
        yaxis_title="Y Position",
        showlegend=True,
        autosize=True,
        width=None,
        height=None,
        margin=dict(l=20, r=20, t=40, b=80),
        xaxis=dict(scaleanchor="y", scaleratio=1),
        yaxis=dict(constrain='domain')
    )

    return fig