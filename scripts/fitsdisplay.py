import numpy as np
import matplotlib
import os
# matplotlib.rc_file("../../templates/matplotlibrc")
import matplotlib.pyplot as plt
from astropy.io import fits

fits_list = []
for filename in os.listdir("frames/20120305/GEO_tracking/Object_1"):
    if filename.startswith("FIELD"):
        fits_list.append("frames/20120305/GEO_tracking/Object_1/" + filename)

image_concat = []
for image in fits_list:
    image_concat.append(fits.getdata(image))

final_image = np.amax(image_concat, axis=0)

# hdu_list = fits.open("frames.tar/frames/20120305/GEO_tracking/Object_1/FIELD01_R.091")
# final_image = hdu_list[0].data
# hdu_list.close()

plt.imshow(final_image, cmap="gray", vmin=4200, vmax=5000)
plt.colorbar()
plt.show()
#
NBINS = 1000
histogram = plt.hist(final_image.flatten(), NBINS)
plt.show()
