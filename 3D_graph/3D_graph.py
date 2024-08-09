



import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import make_interp_spline

lengths = {
    '0.01': ([100, 45.33, 58.79, 80.75, 84.95, 93.05, 94.73, 97.53, 98.57, 99.80], np.arange(0, 10)),
    '0.1': ([100, 44.93, 61.11, 81.98, 84.27, 93.05, 93.69, 97.70, 98.87, 98.22], np.arange(0, 10)),
    '0.2': ([100, 45.23, 57.53, 78.63, 80.24, 90.00, 89.94, 96.09, 97.83, 98.30, 99.23, 98.80], np.arange(0, 12)),
    '0.4': ([100, 49.02, 59.58, 79.18, 80.91, 89.22, 88.39, 88.12, 90.53, 91.16, 95.77, 98.10, 98.37, 99.33], np.arange(0, 14)),
    '0.6': ([100, 49.35, 58.62, 76.24, 79.61, 86.83, 88.35, 87.10, 87.51, 88.22, 93.36, 95.81, 97.52, 98.04, 99.26],
            np.arange(0, 15)),
    '1.0': ([100, 62.17, 61.29, 76.72, 84.27, 90.91, 92.65, 90.30, 84.20, 86.20, 89.99, 94.01, 96.78, 97.68, 98.03],
            np.arange(0, 15)),
    '1.5': ([100, 66.34, 56.75, 72.83, 82.28, 89.28, 92.12, 87.80, 89.26, 80.51, 88.85, 94.18, 96.73, 98.41, 98.90],
            np.arange(0, 15)),
    '2.0': ([100, 72.57, 57.42, 74.68, 85.26, 90.61, 93.35, 91.14, 83.49, 80.54, 89.94, 93.06, 96.10, 98.00, 99.73],
            np.arange(0, 15)),
    '2.5': (
    [99.87, 76.17, 57.24, 72.70, 85.47, 91.28, 93.42, 90.41, 89.82, 78.90, 85.95, 92.15, 94.24, 96.76, 98.42, 92.50],
    np.arange(0, 16)),
    '3.0': (
    [98.42, 73.28, 56.65, 73.10, 84.91, 90.93, 93.29, 91.58, 92.81, 77.79, 78.21, 86.76, 91.43, 94.09, 97.41, 94.40],
    np.arange(0, 16))
}


thicknesses = {
    '0.01': [100, 100, 100, 100, 100, 100, 99, 98, 95, 1],
    '0.1': [100, 100, 100, 100, 100, 100, 100, 98, 94, 5],
    '0.2': [100, 100, 100, 100, 100, 97, 96, 93, 82, 19, 3, 1],
    '0.4': [100, 100, 100, 100, 100, 95, 91, 86, 72, 46, 33, 23, 10, 3],
    '0.6': [100, 100, 100, 100, 100, 98, 97, 93, 78, 69, 56, 44, 25, 17, 5],
    '1.0': [100, 100, 100, 100, 100, 100, 95, 92, 48, 47, 40, 34, 26, 18, 7],
    '1.5': [100, 100, 100, 100, 100, 100, 100, 99, 45, 42, 38, 36, 24, 20, 10],
    '2.0': [100, 100, 100, 100, 100, 100, 100, 98, 32, 32, 28, 26, 21, 18, 9],
    '2.5': [100, 100, 100, 100, 100, 100, 100, 99, 26, 26, 25, 24, 17, 16, 12, 1],
    '3.0': [100, 100, 100, 100, 100, 100, 99, 91, 27, 21, 17, 16, 12, 12, 8, 1]
}


# gradient color map
cmap = plt.cm.viridis  # Choose a colormap
norm = plt.Normalize(vmin=0, vmax=len(lengths)-1)
colors = [cmap(norm(i)) for i in range(len(lengths))]

branch_lengths = list(lengths.keys())

fig = plt.figure(figsize=(14, 10))
ax = fig.add_subplot(111, projection='3d')

for i in range(len(branch_lengths) - 1, -1, -1):
    length = branch_lengths[i]
    values, x = lengths[length]


    spl = make_interp_spline(x, values, k=3)
    x_smooth = np.linspace(x.min(), x.max(), 500)
    values_smooth = spl(x_smooth)


    for j in range(len(x_smooth) - 1):
        ax.plot([x_smooth[j], x_smooth[j + 1]], [float(length), float(length)], [values_smooth[j], values_smooth[j + 1]],
                color=colors[i], linewidth=max(thicknesses[length][int(x_smooth[j])] / 14.0, 2.0))  # Minimum linewidth of 2.0

ax.set_xlabel('Order of Adjacent Branch', fontsize=20, labelpad=20)
ax.set_ylabel('Branch Length', fontsize=20, labelpad=20)
ax.set_zlabel('UFBoot Value', fontsize=20, labelpad=20)
ax.set_title('UFBoot Values in Clade t0', fontsize=20)

ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
ax.tick_params(axis='z', labelsize=16)

#thickness
legend_lines = []
legend_labels = []
for length, color in zip(branch_lengths, colors):
    legend_lines.append(plt.Line2D([0], [0], color=color, linewidth=4))
    legend_labels.append(f'Length {length}')

legend = ax.legend(legend_lines, legend_labels, fontsize=18, title='JB Length', title_fontsize=18,
                   bbox_to_anchor=(1.05, 1), loc='upper left')

# for revoming white space
fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)

plt.show()
