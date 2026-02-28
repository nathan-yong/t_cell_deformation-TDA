from nicegui import ui
import numpy as np

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
number_of_particles = 10


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


def spawn_circles():
    global next_id, circles, number_of_particles
    number_of_particles = int(number_of_particles)
    circles.clear()

    average_radius = 1.0
    n = int(np.sqrt(number_of_particles))
    start_location = simulation_size_L / 2 - (n * average_radius)
    x = start_location
    y = start_location
    coord_step_size = average_radius * 2
    row_iterator = 0
    for i in range(number_of_particles):
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

    update_image()


def update_image():
    # Generate SVG content based on the circles list
    content = ""
    for c in circles:
        content += f'<circle cx="{c["x"]}" cy="{c["y"]}" r="{c["radius"]}" fill="rgba(0, 0, 255, 0.7)" stroke="{c["color"]}" stroke-width="0.1" />'
    ii.content = content


# UI Layout
with ui.row().classes("w-full h-[95vh] no-wrap items-stretch bg-slate-50 p-4"):
    with ui.column().classes("flex-[7] h-full items-center justify-center p-5"):
        # The interactive image
        # Note: 'crosshair' cursor helps with precision
        ii = ui.interactive_image(
            size=(simulation_size_L, simulation_size_L),
            on_mouse=handle_mouse,
            events=["mousedown", "mousemove", "mouseup"],
            cross=False,
        ).classes("w-max-full h-full aspect-square border-4 border-black")
        ii.style("background-color: #ddd;")

    with ui.column().classes("flex=[3] h-full items-center justify-start p-4"):
        with ui.card().classes("w-full h-auto shadow-lg p-0"):

            with ui.tabs().classes("w-full border-b") as tabs:
                ui.tab("Manual", icon="settings")
                ui.tab("Upload", icon="file_upload")

            with ui.tab_panels(tabs, value="Manual").classes("w-full bg-transparent"):
                with ui.tab_panel("Manual"):
                    with ui.column().classes("gap-4"):
                        ui.label("Configuration").classes("text-lg font-bold")

                        with ui.card().classes("w-full"):
                            ui.label("Number of Particles").classes("font-bold")
                            ui.number(
                                value=10, min=0, max=50, step=1, format="%i"
                            ).bind_value(globals(), "number_of_particles")
                            ui.slider(min=0, max=50, value=10).bind_value(
                                globals(), "number_of_particles"
                            ).classes("w-full")

                        ui.button("Spawn Circles", on_click=spawn_circles).classes(
                            "w-full py-2"
                        ).props("elevated color=primary")

                with ui.tab_panel("Upload"):
                    with ui.column().classes("gap-4 items-center"):
                        ui.label("Import Simulation").classes(
                            "text-lg font-bold self-start"
                        )

                        ui.upload(
                            on_upload=lambda e: ui.notify(f"Uploaded {e.file.name}"),
                            label="Upload Results.txt",
                        ).classes("w-full")

ui.run()
