import scipy.io
import pandas as pd
import numpy as np
from datetime import datetime, timedelta

# Load the MATLAB file
data = scipy.io.loadmat('data/metadata/LSEO24h_temp.mat')

# Extract the data
T = pd.DataFrame(data['T'])  # Temperature matrix (9390 x 60)
t = pd.DataFrame(data['t'].flatten())  # Time array (9390,) - MATLAB datenum format
dep = data['dep'].flatten()  # Depth array (60,)
SN = data['SN'].flatten()  # Serial numbers (60,)
dep_sn = {"sn" : SN, "dep" : dep}
dep_sn = pd.DataFrame(dep_sn)
# Convert MATLAB datenum to Python datetime
# MATLAB datenum starts from January 0, 0000 (which is December 31, 1 BC in proleptic Gregorian)
# Python datetime starts from January 1, 1 AD
# Difference is 366 days (year 0 in MATLAB)

# Create long format dataframe
# Create DataFrame

# Save dep_sn to CSV
output_path = 'data/metadata/dep_sn.csv'
dep_sn.to_csv(output_path, index=False)
# Save t to CSV
output_path = 'data/metadata/t.csv'
t.to_csv(output_path, index=False)
# Save T to CSV
output_path = 'data/metadata/T.csv'
T.to_csv(output_path, index=False)

