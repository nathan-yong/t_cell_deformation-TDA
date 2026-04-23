import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Circle

# --- Configuration ---
INPUT_FILE = 'w3adh075rep0/Results.txt'
OUTPUT_FILE = 'w3adh075rep0.mp4'
FPS = 30
# Update LIMITS based on your simulation box size
LIMITS = (0, 43) 

def load_frames(filename):
    """Parses Results.txt into a list of frames."""
    frames = []
    current_particles = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            # Lines starting with 't' separate the frames
            if line.startswith('t'):
                if current_particles:
                    frames.append(current_particles)
                    current_particles = []
            else:
                parts = line.split()
                if len(parts) == 6:
                    try:
                        # Format: id, x, y, xvel, yvel, radius
                        current_particles.append({
                            'id': int(parts[0]),
                            'x': float(parts[1]),
                            'y': float(parts[2]),
                            'r': float(parts[5])
                        })
                    except ValueError:
                        continue
        
        # Add the final frame in the file
        if current_particles:
            frames.append(current_particles)
    return frames

# 1. Load the data
print(f"Loading data from {INPUT_FILE}...")
simulation_frames = load_frames(INPUT_FILE)
total_frames = len(simulation_frames)
print(f"Found {total_frames} frames.")

# 2. Setup the Plot
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_xlim(LIMITS)
ax.set_ylim(LIMITS)
ax.set_aspect('equal')
ax.set_title("Particle Simulation")

# Add a text object to display the frame count in the video
frame_text = ax.text(0.02, 0.95, '', transform=ax.transAxes, 
                     weight='bold', fontsize=12, color='black')

def update(frame_idx):
    """Function called for each frame of the animation."""
    # Print progress to the command line
    print(f"Working on frame: {frame_idx}/{total_frames - 1}", end='\r')
    
    # Remove previous circles from the plot
    for patch in ax.patches:
        patch.remove()
        
    # Update the frame counter text in the video
    frame_text.set_text(f'Frame: {frame_idx}')
    
    current_data = simulation_frames[frame_idx]
    
    for p in current_data:
        # ID 1 is the center cell -> Green; Others -> Blue
        color = 'green' if p['id'] == 1 else 'blue'
        
        # Create circle with position (x, y) and radius (r)
        c = Circle((p['x'], p['y']), p['r'], color=color, alpha=0.6)
        ax.add_patch(c)
        
    return ax.patches + [frame_text]

# 3. Create and Save the Animation
# 'blit=True' makes the rendering faster by only re-drawing changed parts
ani = animation.FuncAnimation(
    fig, update, frames=total_frames, 
    blit=True, interval=1000/FPS
)

print(f"\nStarting video encoding...")
try:
    # Requires FFmpeg installed on your machine
    writer = animation.FFMpegWriter(fps=FPS, metadata=dict(artist='Me'), bitrate=3800)
    ani.save(OUTPUT_FILE, writer=writer)
    print(f"\nSuccessfully saved animation to {OUTPUT_FILE}")
except Exception as e:
    print(f"\nError: {e}")
    print("Tip: Make sure 'ffmpeg' is installed and added to your system PATH.")

plt.close()