import scipy.io
import pandas as pd


# Load the MATLAB file
data = scipy.io.loadmat('data/metadata/LSEO24h_temp.mat')

# Extract the data
T = pd.DataFrame(data['T'])  # Temperature matrix (9390 x 60)
t = data['t'].flatten()  # Time array (9390,) - MATLAB datenum format
dep = data['dep'].flatten()  # Depth array (60,)
SN = data['SN'].flatten()  # Serial numbers (60,)
dep_sn = {"sn" : SN, "dep" : dep}
dep_sn = pd.DataFrame(dep_sn)
# Possible script to convert to other date time format

from datetime import datetime, timedelta

def matlab_datenum_to_python_datetime(matlab_datenum):
    """
    Convert a MATLAB serial date number to a Python datetime object.
    
    MATLAB serial date 1 is Jan 1, 0000.
    Python ordinal date 1 is Jan 1, 0001.
    There is a one-year difference, plus the non-existent year 0 in Python's ordinal system.
    The most reliable base difference is an offset of 366 days.
    """
    days = matlab_datenum % 1
    # datetime.fromordinal(int(matlab_datenum)) converts the integer part
    # timedelta(days=days) adds the fractional part (time of day)
    # - timedelta(days=366) adjusts for the epoch difference (Jan 1, 0000 vs Jan 1, 0001)
    python_datetime = datetime.fromordinal(int(matlab_datenum)) + timedelta(days=days) - timedelta(days=366)
    return python_datetime

# Example Usage:
# MATLAB serial date 731885.75 corresponds to '31-Oct-2003, 6:00 PM'
dt = for date in t : matlab_datenum_to_python_datetime(t[date])
print(t)
python_dt = matlab_datenum_to_python_datetime(t)

print(f"MATLAB datenum: {matlab_date_number}")
print(f"Python datetime: {python_dt}")

# Create DataFrames

# Save dep_sn to CSV
output_path = 'data/metadata/dep_sn.csv'
dep_sn.to_csv(output_path, index=False)
# Save t to CSV
output_path = 'data/metadata/t.csv'
t.to_csv(output_path, index=False)
# Save T to CSV
output_path = 'data/metadata/Time.csv'
T.to_csv(output_path, index=False)

