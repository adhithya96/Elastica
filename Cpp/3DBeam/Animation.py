
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
import os

x1 = []
y1 = []
z1 = []
x2 = []
y2 = []
z2 = []

FFMpegWriter = manimation.writers['ffmpeg']
metadata = dict(title='Movie Test', artist='Matplotlib',
                comment='a red circle following a blue sine wave')
writer = FFMpegWriter(fps=2, metadata=metadata)

fig = plt.figure()
ax2 = fig.subplots()

ax2.set_xlim([0, 10])
ax2.set_ylim([0, 0.001])
ax2.set_title("Beam Deformed configuration")
ax2.set_xlabel('x-coordinate')
ax2.set_ylabel('y-coordinate')
plot1, = ax2.plot([], [], marker='o', color='g')
plot2, = ax2.plot([], [], marker='o', color='r')

print(__file__)  # __vsc_ipynb__file

CSD = os.path.dirname(__file__)
data1_name = "EBBE3D1_"
data2_name = "EBBE3D2_"

with writer.saving(fig, "writer_test.mp4", 100):
    for i in range(1, 10):
        data1_path = os.path.join(CSD, data1_name + str((i - 1) * 10) + ".txt")
        data2_path = os.path.join(CSD, data2_name + str((i - 1) * 10) + ".txt")
        f = open(data1_path)
        x1 = []
        y1 = []
        z1 = []
        x2 = []
        y2 = []
        z2 = []

        for row in f:
            row = row.split(' ')
            x1.append(float(row[0]))
            y1.append(float(row[1]))
            z1.append(float(row[2]))

        g = open(data2_path)
        for row in g:
            row = row.split(' ')
            x2.append(float(row[0]))
            y2.append(float(row[1]))
            z2.append(float(row[2]))

        plot1.set_data(x1, y1)
        plot2.set_data(x2, y2)

        writer.grab_frame()
