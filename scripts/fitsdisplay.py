import numpy as np
import matplotlib
# matplotlib.rc_file("../../templates/matplotlibrc")
import matplotlib.pyplot as plt

from astropy.io import fits
hdu_list = fits.open("../misc/FIELD01_R.091")
image_data = hdu_list[0].data
hdu_list.close()

plt.imshow(image_data, cmap="gray")

plt.colorbar()
