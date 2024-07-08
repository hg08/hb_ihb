#!/usr/bin/env python3
import os
import json

# List of trajectories
trajectories = [
    {
        "name": "128w-pos-1.xyz",
        "dt": 0.5, # fs 
        "sizeX": 15.6404, # Angstrom 
        "sizeY": 15.6404,
        "sizeZ": 31.2808,
        "numAtoms": 384 
    }
]

# Write metadata to json file
for traj in trajectories:
    metaFile = traj["name"].replace(".xyz", ".json") 
    metaFile = os.path.join('m2_traj', metaFile)
    with open(metaFile, "w") as outfile:
        json.dump(traj, outfile, indent=4)
    print(f"Metadata written to {metaFile}")
