# sets up the root folder for PULSE pipeline

import os

# PREPROCESS SETUP
def setup_for_preprocess():
    paths_to_create = [r'./output', r'./output/preprocess']
    for path in paths_to_create:
        if not os.path.exists(path):
            os.makedirs(path)

if __name__ == "__main__":
    setup_for_preprocess()