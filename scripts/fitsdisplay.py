import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.io import fits

# TODO - intelligently set which folder to use; draw a line using first and last point (get slope and intercept from
# regression in java)

folder = "../misc/data/Oct2017/NEA/2017_PR25_R_7/"

fits_list = []
for filename in os.listdir(folder):
    if filename.endswith(".fit"):
        fits_list.append(folder + filename)

image_concat = []
for image in fits_list:
    image_concat.append(fits.getdata(image))

final_image = np.amax(image_concat, axis=0)

fig = plt.imshow(final_image, cmap="gray", vmin=1000, vmax=6000)
plt.savefig("result/fig.png")
plt.colorbar()
plt.show()
