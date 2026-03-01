import numpy as np
import copy


def read_particle_data_from_string(particle_string):
    num_frames = 0

    frame_array = []
    frame_data = []
    for line in particle_string.splitlines():
        line_split = line.strip().split()

        # Start processing line
        if line[0] == "t":
            if num_frames > 0:
                frame_array.append(copy.deepcopy(frame_data))
            frame_data = []
            num_frames += 1
        else:
            frame_data.append(
                [
                    float(line_split[0]),
                    float(line_split[1]),
                    float(line_split[2]),
                    float(line_split[3]),
                    float(line_split[4]),
                    float(line_split[5]),
                ]
            )
    return frame_array
