import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Load data from CSV file
df = pd.read_csv('output.csv')

# Extract coordinates from data frame
x = df['X']
y = df['Y']
z = df['Z']

# Create figure and axes
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Create scatter plot of satellite locations
ax.scatter(x, y, z)

# Set labels and title
ax.set_xlabel('X coordinate')
ax.set_ylabel('Y coordinate')
ax.set_zlabel('Z coordinate')
ax.set_title('Satellite Locations')

# Add satellite ID and time error number as text annotations
for i, row in df.iterrows():
    ax.text(row['X'], row['Y'], row['Z'],
            f"{row['PRN']} ({row['dts_Li']})")

# Show plot
plt.show()