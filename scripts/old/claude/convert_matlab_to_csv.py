import scipy.io
import pandas as pd
import numpy as np

# Load the MATLAB file
data = scipy.io.loadmat('data/metadata/LSEO24h_temp.mat')

# Extract the data
T = data['T']  # Temperature matrix (9390 x 60)
t = data['t'].flatten()  # Time array (9390,)
dep = data['dep'].flatten()  # Depth array (60,)
SN = data['SN'].flatten()  # Serial numbers (60,)

# Create a DataFrame with time as index
# Each column represents a depth, with temperature values
df = pd.DataFrame(T, columns=[f'depth_{d:.1f}m' for d in dep])

# Add the time column
df.insert(0, 'time', t)

# Save to CSV
output_path = 'output/LSEO24h_temp.csv'
df.to_csv(output_path, index=False)

print(f"✓ Converted successfully!")
print(f"✓ Shape: {df.shape} (rows × columns)")
print(f"✓ Saved to: {output_path}")
print(f"\nFirst few rows:")
print(df.head())
print(f"\nColumn names:")
print(df.columns.tolist()[:10], "... (60 depth columns total)")
