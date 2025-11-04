from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord, EarthLocation, get_moon, get_body
import astropy.units as u
from astropy.time import Time
from datetime import datetime

def radec2pixel(wcs, ra, dec):
    skycoord = SkyCoord(ra*u.deg, dec*u.deg, frame='icrs')
    x, y = wcs.world_to_pixel(skycoord)
    return x, y

def getMoonPixel(wcs, utc):
    # ra = 333.6071
    # dec = -17.0267
    t = Time(utc, scale='utc')
    #location = EarthLocation(lon=116.67*u.deg, lat=-26.7*u.deg,
    #    height=377*u.m)
    moon_icrs = get_moon(t)
    return radec2pixel(wcs, moon_icrs.ra.deg, moon_icrs.dec.deg)

def main(args):

    ## open fits file
    hdu = fits.open(args.inputFits)
    wcs = WCS(hdu[0].header, naxis=2)
    data = hdu[0].data[0,0,:,:] 
    utc = datetime.strptime(hdu[0].header['DATE-OBS'][:-2],
                    '%Y-%m-%dT%H:%M:%S')
    h,w = data.shape
    cal_col, cal_row = getMoonPixel(wcs, utc)
    
    mask = np.zeros((h,w))
    yy, xx = np.ogrid[:h, :w]
    dist = np.sqrt((yy-cal_col)**2 + (xx-cal_row)**2)

    mask[dist < 36] = 1
    
    hdun = fits.PrimaryHDU(mask, hdu[0].header)
    hdun.writeto('mask.fits')

if __name__ == "__main__":
    parser = ArgumentParser('masker')
    parser.add_argument('--inputFits', required=True, help='Input fits file')
    args = parser.parse_args()

    main(args)