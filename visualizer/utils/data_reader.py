from io import StringIO

import numpy as np
import copy


def read_particle_data_from_string(particle_string, contact_times_string, use_contact_times=False):
    num_frames = 0
    contact_times_id = []
    contact_times_time = []
    for line in contact_times_string.splitlines():
        line_split = line.strip().split()
        contact_times_id.append(int(line_split[0]))
        contact_times_time.append(float(line_split[1]))
    


    frame_array = []
    frame_data = []
    current_time = 0
    for line in particle_string.splitlines():
        line_split = line.strip().split()

        # Start processing line
        if line[0] == "t":
            current_time = float(line_split[1])
            if num_frames > 0:
                frame_array.append(copy.deepcopy(frame_data))
            frame_data = []
            num_frames += 1
        else:
            if use_contact_times and int(line_split[0]) in contact_times_id and current_time >= contact_times_time[contact_times_id.index(int(line_split[0]))]: 
                frame_data.append(
                    [
                        float(line_split[0]),
                        float(line_split[1]),
                        float(line_split[2]),
                        float(line_split[3]),
                        float(line_split[4]),
                        float(line_split[5]),
                        True
                    ]
                )
            else:
                frame_data.append(
                    [
                        float(line_split[0]),
                        float(line_split[1]),
                        float(line_split[2]),
                        float(line_split[3]),
                        float(line_split[4]),
                        float(line_split[5]),
                        False
                    ]
                )
    return frame_array
