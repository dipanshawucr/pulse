import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier


training_examples = 145
total_data = 15929

cutoff = 0.6

features = pd.read_table('/home/wonjunetai/Documents/KimLab/pulse/input/param_files/ML_Training_data.txt', sep=',')


