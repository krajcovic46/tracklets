import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.io import fits

folder = "../misc/data/Oct2017/NEA/PR25_R_7/"

results_list = list()
with open("../ODResults/myresults.txt", 'r') as f:
    results_list = f.readlines()

print(results_list)
first_line = results_list.pop(0).split("\t"); results_list.pop(0)

slope = first_line[0]
intercept = first_line[1][:-1]

first_object = results_list.pop(0).split("\t")
last_object = results_list.pop(-3).split("\t")

x1 = first_object[-2]; y1 = first_object[-1][:-1]
x2 = last_object[-2]; y2 = last_object[-1][:-1]

fits_list = []
for filename in os.listdir(folder):
    if filename.endswith(".fit"):
        fits_list.append(folder + filename)

image_concat = []
for image in fits_list:
    image_concat.append(fits.getdata(image))

final_image = np.sum(image_concat, axis=0)

fig = plt.imshow(final_image, cmap="gray", vmin=10000, vmax=25000)
plt.plot([float(x1), float(x2)], [float(y1), float(y2)], "r-")
plt.savefig("result/fig.png")
plt.colorbar()
plt.show()
