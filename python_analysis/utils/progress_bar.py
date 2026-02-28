from tqdm import tqdm

progress_bar = None

def start_progress_bar(total_frames):
  global progress_bar
  progress_bar = tqdm(total=total_frames)
  
def increment_progress_bar(increment):
  global progress_bar
  progress_bar.update(increment)
  
def stop_progress_bar():
  global progress_bar
  progress_bar.close()