import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier


training_examples = 145
total_data = 15929

num_iterations = 5

cutoff = 0.6

features = pd.read_table('/home/wonjunetai/Documents/KimLab/pulse/input/param_files/ML_Training_data.txt', sep=',')

Y = np.zeros(total_data)
N = Y.shape[0]

