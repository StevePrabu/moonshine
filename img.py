from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.wcs import WCS

obsid = 1058715072
ch = 1

hdu = fits.open('processing/{0}/{0}-img-narrowband-{1}-image-pb.fits'.format(obsid, str(ch).zfill(4)))
data = hdu[0].data[0,0,:,:]

hdu2 = fits.open('processing/{0}/{0}-img-3C44-{1}-image-pb.fits'.format(obsid, str(ch).zfill(4)))
data2 = hdu2[0].data[0,0,:,:]


plt.figure(figsize=(10,3))
plt.subplot(131)
plt.imshow(data, origin='lower', vmax=1, vmin=-1)
plt.colorbar()


plt.subplot(132)
plt.imshow(data2, origin='lower', vmax=1, vmin=-1)
plt.colorbar()


plt.subplot(133)
plt.imshow(data-data2, origin='lower', vmax=1, vmin=-1)
plt.colorbar()



plt.savefig('img.png')