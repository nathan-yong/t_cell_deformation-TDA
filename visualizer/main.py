from nicegui import ui, events
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go

import utils.data_reader as data_reader
import tda.alpha_complexes as alpha_complexes

# --------------------------------------------------------------------------
# SECTION A: NICEGUI GLOBALS
# --------------------------------------------------------------------------
# State to store our circles
# [{
# 'id': 0,
# 'x': 21,
# 'y': 21
# 'radius': 1.0
# 'color': "blue"
# }, ...]
circles = []
next_id = 1
# Track which circle is being dragged
dragging = None
simulation_size_L = 43.0
# Manual
num_particles_to_spawn = 10
num_particles = 0
# File upload
usingFirstContactTimes = False
particle_data = None
contact_times_string = ""
use_contact_times = False
num_frames = 0
selected_frame = 1


def handle_mouse(e):
    global dragging

    # Mouse Down: Check if we clicked near an existing circle
    if e.type == "mousedown":
        for c in circles:
            # Simple distance check (radius is roughly 20 units)
            if ((e.image_x - c["x"]) ** 2 + (e.image_y - c["y"]) ** 2) ** 0.5 < c[
                "radius"
            ]:
                dragging = c
                return

    # Mouse Move: Update coordinates if dragging
    if e.type == "mousemove" and dragging:
        dragging["x"] = e.image_x
        dragging["y"] = e.image_y
        update_image()

    # Mouse Up: Stop dragging
    if e.type == "mouseup":
        dragging = None


def spawn_particles():
    global next_id, circles, num_particles, num_particles_to_spawn
    next_id = 1
    circles.clear()
    num_particles_to_spawn = int(num_particles_to_spawn)

    average_radius = 1.0
    n = int(np.sqrt(num_particles_to_spawn))
    start_location = simulation_size_L / 2 - (n * average_radius)
    x = start_location
    y = start_location
    coord_step_size = average_radius * 2
    row_iterator = 0
    for i in range(num_particles_to_spawn):
        circles.append(
            {"id": next_id, "x": x, "y": y, "radius": average_radius, "color": "blue"}
        )
        x += coord_step_size
        row_iterator += 1
        if row_iterator == n:
            row_iterator = 0
            x = start_location
            y += coord_step_size
        next_id += 1

    num_particles = len(circles)

    update_image()


def add_particle():
    global next_id, circles, num_particles

    circles.append(
        {
            "id": next_id,
            "x": simulation_size_L / 2,
            "y": simulation_size_L / 2,
            "radius": 1,
            "color": "blue",
        }
    )
    next_id += 1
    num_particles = len(circles)

    update_image()

def update_image():
    # Generate SVG content based on the circles list
    content = ""
    for c in circles:
        content += f'<circle cx="{c["x"]}" cy="{c["y"]}" r="{c["radius"]}" fill="{c["color"]}" fill-opacity="0.4" stroke="{c["color"]}" stroke-width="0.1" />'
    ii.content = content

async def simulation_results_file_upload(e: events.UploadEventArguments):
    global num_frames, particle_data
    particle_data_string = await e.file.text()
    particle_data = data_reader.read_particle_data_from_string(
        particle_data_string, contact_times_string, use_contact_times
    )
    num_frames = len(particle_data)
    simulation_load_number.max = num_frames
    simulation_load_slider._props['max'] = num_frames
    simulation_load_slider.update()
    display_loaded_frame()

async def simulation_contacts_file_upload(e):
    global contact_times_string, use_contact_times
    contact_times_string = await e.file.text()
    use_contact_times = True



def display_loaded_frame():
    global selected_frame, particle_data, next_id, circles, num_particles
    selected_frame = int(selected_frame)
    circles.clear()
    if particle_data is None:
        return
    for p in particle_data[selected_frame - 1]:
        color = "blue"
        if p[0] == 1:
            color = "green"
        elif p[6]:
            color = "red"
        circles.append(
            {
                "id": p[0],
                "x": p[1],
                "y": p[2],
                "radius": p[5],
                "color": color,
            }
        )
    next_id = max([c["id"] for c in circles]) + 1
    num_particles = len(circles)
    update_image()

def run_alpha_complexes_analysis():
    if not circles:
        return
    coords = np.array([[c["x"], c["y"]] for c in circles])
    fig = alpha_complexes.delaunay_plotly_visualization(coords)
    ac_visualizer.update_figure(fig)
    dim0_pd_birth_death, dim1_pd_birth_death = alpha_complexes.calculate_alpha_complex_pd(coords)
    dim0_pd_x_list = dim0_pd_birth_death[:, 0]
    dim0_pd_y_list = dim0_pd_birth_death[:, 1]
    dim1_pd_x_list = dim1_pd_birth_death[:, 0]
    dim1_pd_y_list = dim1_pd_birth_death[:, 1]
    dim0_pd_ax.clear()
    dim0_pd_ax.plot(dim0_pd_birth_death[:, 0], dim0_pd_birth_death[:, 1], 'o')
    dim1_pd_ax.clear()
    dim1_pd_ax.plot(dim1_pd_birth_death[:, 0], dim1_pd_birth_death[:, 1], 'o')

    dim0_pd_plot.update()
    dim1_pd_plot.update()




# UI Layout
with ui.row().classes("w-full h-[95vh] no-wrap items-stretch bg-slate-50 p-4"):
    with ui.column().classes("flex-[7] h-full items-center justify-start p-5"):
        with ui.tabs().classes("w-full border-b") as tabs:
            ui.tab("Simulation", icon="play_arrow")
            ui.tab("Analysis Visualization", icon="analytics")
            ui.tab("Persistence Diagrams", icon="insights")
        
        with ui.tab_panels(tabs, value="Simulation").classes("w-full h-full bg-transparent"):
            with ui.tab_panel("Simulation").classes("w-full h-full items-center justify-center"):
                with ui.row():
                    ui.label("Number of Particles: ")
                    ui.label().bind_text_from(globals(), "num_particles")

                # Interactive Image Canvas
                ii = ui.interactive_image(
                    size=(simulation_size_L, simulation_size_L),
                    on_mouse=handle_mouse,
                    events=["mousedown", "mousemove", "mouseup"],
                    cross=False,
                ).classes("w-max-full h-full aspect-square border-4 border-black")
                ii.style("background-color: #ddd;")

            with ui.tab_panel("Analysis Visualization").classes("w-full h-full items-center justify-center"):
                ac_visualizer = ui.plotly(alpha_complexes.delaunay_plotly_visualization(np.array([[1,1], [2,1], [1,2], [2,4]]))).classes('w-max-full h-full aspect-square')

            with ui.tab_panel("Persistence Diagrams").classes("w-full h-full items-center justify-center"):
                with ui.row().classes("w-full gap-4 no-wrap items-center justify-center"):
                    with ui.card().classes("w-auto h-auto shadow-lg items-center p-4"):
                        with ui.pyplot(close=False) as dim0_pd_plot:
                            # Example plot (replace with actual alpha complex analysis results)
                            dim0_pd_ax = dim0_pd_plot.fig.add_subplot()
                            dim0_pd_ax.set_xlabel("Birth")
                            dim0_pd_ax.set_ylabel("Death")
                            dim0_pd_ax.set_title("Persistence Diagram Betti 0")
                            dim0_pd_x_list = [0, 1, 2]
                            dim0_pd_y_list = [0, 1, 0] 
                            dim0_pd_ax.plot(dim0_pd_x_list, dim0_pd_y_list, 'o')
                    with ui.card().classes("w-auto h-auto shadow-lg  items-center p-4"):
                        with ui.pyplot(close=False) as dim1_pd_plot:
                            # Example plot (replace with actual alpha complex analysis results)
                            dim1_pd_ax = dim1_pd_plot.fig.add_subplot()
                            dim1_pd_ax.set_xlabel("Birth")
                            dim1_pd_ax.set_ylabel("Death")
                            dim1_pd_ax.set_title("Persistence Diagram Betti 1")
                            dim1_pd_x_list = [0, 1, 2]
                            dim1_pd_y_list = [0, 1, 0] 
                            dim1_pd_ax.plot(dim1_pd_x_list, dim1_pd_y_list, 'o')
            
    with ui.column().classes("flex=[3] h-full items-center justify-start p-4"):
        with ui.card().classes("w-full h-auto shadow-lg p-0"):

            with ui.tabs().classes("w-full border-b") as tabs:
                ui.tab("Manual", icon="settings")
                ui.tab("Upload", icon="file_upload")
                ui.tab("Analysis", icon="analytics")

            with ui.tab_panels(tabs, value="Manual").classes("w-full bg-transparent"):
                with ui.tab_panel("Manual"):
                    with ui.column().classes("gap-4"):
                        ui.label("Configuration").classes("text-lg font-bold")

                        with ui.card().classes("w-full"):
                            ui.label("Number of Particles").classes("font-bold")
                            ui.number(
                                value=10, min=0, max=50, step=1, format="%i"
                            ).bind_value(globals(), "num_particles_to_spawn")
                            ui.slider(min=0, max=50, value=10).bind_value(
                                globals(), "num_particles_to_spawn"
                            ).classes("w-full")

                        ui.button("Spawn Particles", on_click=spawn_particles).classes(
                            "w-full py-2"
                        ).props("elevated color=primary")

                        ui.button("Add Particles", on_click=add_particle).classes(
                            "w-full py-2"
                        ).props("elevated color=primary")

                with ui.tab_panel("Upload"):
                    with ui.column().classes("gap-4 items-center"):
                        ui.label("Import Simulation").classes(
                            "text-lg font-bold self-start"
                        )

                        ui.upload(
                            on_upload=simulation_results_file_upload,
                            label="Upload Results.txt",
                        ).classes("w-full").props("accept=.txt")

                        firstContactTimesCheckbox = ui.checkbox(
                            "Use First Contact Times"
                        ).bind_value(globals(), "usingFirstContactTimes")

                        ui.upload(
                            on_upload=simulation_contacts_file_upload,
                            label="Upload First_contact_times.txt",
                        ).bind_visibility_from(
                            firstContactTimesCheckbox, "value"
                        ).classes(
                            "w-full"
                        ).props(
                            "accept=.txt"
                        )

                        with ui.card().classes("w-full"):
                            ui.label("Frame").classes("font-bold")
                            simulation_load_number = ui.number(
                                value=1, min=1, max=1, step=1, format="%i"
                            ).bind_value(globals(), "selected_frame").on(
                                "update:model-value", display_loaded_frame
                            )
                            simulation_load_slider = ui.slider(min=1, max=1, value=1).bind_value(
                                globals(), "selected_frame"
                            ).on("update:model-value", display_loaded_frame).classes(
                                "w-full"
                            )
                            
                with ui.tab_panel("Analysis"):
                    with ui.column().classes("gap-4 items-center"):
                        ui.label("Select Analysis").classes(
                            "text-lg font-bold self-start"
                        )
                        analysis_selection = ["Alpha Complexes", "Euclidean Distance Transform"]
                        analysis_selector = ui.select(analysis_selection, value=analysis_selection[0]).classes("w-full")

                        with ui.card().classes("w-full").bind_visibility_from(analysis_selector, "value", value="Alpha Complexes"):
                            ui.label("Alpha Complexes Analysis").classes("font-bold")
                            ui.button("Run Alpha Complexes", on_click=run_alpha_complexes_analysis).classes(
                                "w-full py-2"
                            ).props("elevated color=primary")


ui.run()
