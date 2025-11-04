from astropy.io import fits
import numpy as np
from astropy.wcs import WCS
from astropy.coordinates import EarthLocation, get_moon, get_body, AltAz
from astropy.time import Time
from datetime import datetime
import astropy.units as u
import math
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from scipy import ndimage
from astropy.nddata import Cutout2D
from matplotlib.patches import Circle
from scipy.ndimage import gaussian_filter
from scipy.interpolate import griddata

obs_list = [1058715072,
1063478696,
1084650488,
1084650608,
1094130472,
1094130592,
1094130712,
1094130832,
1103630400,
1103630520,
1103630640,
1108983696,
1108983816,
1108983936,
1108984776,
1108984896,
1108985016,
1108985976,
1108986096,
1108987296,
1108987416,
1108988736,
1108988856,
1108988976,
1108990176,
1108990296,
1108990416,
1113085752,
1113085928,
1117824136,
1117824256,
1117824376,
1117824496,
1117824616,
1120152064,
1120152184,
1120152304,
1120152424,
1120152544,
1122477280,
1122477408,
1122477536,
1148511160,
1148511336,
1150829488,
1150842928,
1157899096,
1157899216,
1157906656,
1157906776,
1157906896,
1157907016,
1157907136,
1160302768,
1160302944,
1160303128,
1160303304,
1160303488,
1160316136,
1162638200,
1165057048,
1165057224,
1165057408,
1165057584,
1165057768,
1176847384,
1176847568,
1176847744,
1176847928,
1176848104,
1179181048,
1179181224,
1179181584,
1179181768,
1181515120,
1181601000,
1181601544,
1193308160,
1193308336,
1193321536,
1193321656,
1193321776,
1193321896,
1193394600,
1193394776,
1193394960,
1193395136,
1193401952,
1193402072,
1195642328,
1195642504,
1195728592,
1195728768,
1195728952,
1195729128,
1195729312,
1195733960,
1207519072,
1207519256,
1207519432,
1207519616,
1207519792,
1209847496,
1209847616,
1209847736,
1209847856,
1209847976,
1209851776,
1209851960,
1209852136,
1209852312,
1209852496,
1209852720,
1209852896,
1209853080,
1209853256,
1209853440,
1209933896,
1209934016,
1209934136,
1212271776,
1212271960,
1212272496,
1212272768,
1212272952,
1212273128,
1212273312,
1212273488,
1214597832,
1214597952,
1214598072,
1214598192,
1214598312,
1214673656,
1214673776,
1214673896,
1216919664,
1217005824,
1221659952,
1223979784,
1223979960,
1223980144,
1223980320,
1223980504,
1223986384,
1224066392,
1224066576,
1224066752,
1224066936,
1226313672,
1226313848,
1226314032,
1226314208,
1226315336,
1226400120,
1226400296,
1226400480,
1226400656,
1226400840,
1226401216,
1226401344,
1226401464,
1228734112,
1228734288,
1228734472,
1228734648,
1228734832,
1228820552,
1228820728,
1228820912,
1228821088,
1240524392,
1240524568,
1240524752,
1240524928,
1240525112,
1242858056,
1242858232,
1242858416,
1242858592,
1242944488,
1242944664,
1245191488,
1245191664,
1250003680,
1250004624,
1252333040,
1254652200,
1254652376,
1254652552,
1254652728,
1254652904,
1254659480,
1254739152,
1254739328,
1254745672,
1256985968,
1256986144,
1256986320,
1256986496,
1256986672,
1256987984,
1257072584,
1257072760,
1257072936,
1257073112,
1257074424,
1259319984,
1259320160,
1259321880,
1259406432,
1259406608,
1259406784,
1259406960,
1259407136,
1259408320,
1271109488,
1271109664,
1271109840,
1271110016,
1271195736,
1271195912,
1271196088,
1271196264,
1271196440,
1273529400,
1273529576,
1273529752,
1273530104,
1275860736,
1275860912,
1275861088,
1275861264,
1275861440,
1275863296,
1275865232,
1275949400,
1275949576,
1275949752,
1275949928,
1275950104,
1278179072,
1287657568,
1287657744,
1287657920,
1287658096,
1287658272,
1287669912,
1287744008,
1287744184,
1287744360,
1287744536,
1287744712,
1289992056,
1289992232,
1290077968,
1290078144,
1290078320,
1290078496,
1290078672,
1318338384,
1318338416,
1318338448,
1323083496,
1323083672,
1323083848,
1323084024,
1323084200,
1323085392,
1323086592,
1323170112,
1323170288,
1323170464,
1323170640,
1339627416,
1339627472,
1339627528,
1339627584,
1339627640,
1341944496,
1341944616,
1341944736,
1344347896,
1344355416,
1346767872,
1346768168,
1346768464,
1346768760,
1346769352,
1349087312,
1349087368,
1349087424,
1349087480,
1349087536,
1351421048,
1351421104,
1351421160,
1351421216,
1351421272,
1351507600,
1351507656,
1351507712,
1353754088,
1353754208,
1353754328,
1353841560,
1353841616,
1353841672,
1353841728,
1365631952,
1365632008,
1365632064,
1365632120,
1365632176,
1367965608,
1367965664,
1367965720,
1367965776,
1367965832,
1370286472,
1370299208,
1370299264,
1370299320,
1370299376,
1370299432,
1370372608,
1370385632,
1370385688,
1370385744,
1370385800,
1370385856,
1372698736,
1375112416,
1379778776,
1379780576,
1379866736,
1382092824,
1382092880,
1382092944,
1382179088,
1382179144,
1382179200,
1382179256,
1382179312,
1382179376,
1382179496,
1382179616,
1382179736,
1382179856,
1382191008,
1384513032,
1384513088,
1384513144,
1384513200,
1384513256,
1384517440,
1412850176,
1412850736,
1412850792,
1412850848,
1412850904,
1412850960,
1412854360,
1412937168,
1412937224,
1412937280,
1412937336,
1412937392,
1415271008,
1415271064,
1415271120,
1415271176,
1415271232,
1417518584,
1417605024,
1417605080,
1417605136,
1417605192,
1417605248,
1436464104,
1436471144]

crab_obs = [1108983696,1108983816,1108983936,1108984776,1108984896,1108985016,1108985976,1108986096,
1108987296,1108987416,1108988736,1108988856,1108988976,1108990176,1108990296,1108990416]

def getMoonRADEC(UTC):
    """
    Gets the Moon's RA and DEC 
    """
    # Observer location
    mwa = EarthLocation(lon=116.67083333*u.deg, 
                        lat=-26.70331941*u.deg,
                        height=377.827*u.m)
    
    # Observation time
    t = Time(UTC, scale='utc')
    
    # Get Moon position in ICRS (RA/Dec)
    moon = get_moon(t, mwa)
    
    return moon.ra.deg, moon.dec.deg, moon


def radec_to_pixel(wcs, ra, dec):
    """
    Convert RA, Dec to image pixel coordinates
    """
    skycoord = SkyCoord(ra*u.deg, dec*u.deg, frame="icrs")
    x, y = wcs.world_to_pixel(skycoord)
    return x, y

def getNoise(data):
    """
    Calculates noise in image through recursive
    flagging of hot pixels
    """
    tmp = np.copy(data)
    tmp[np.abs(tmp) > 3*np.std(tmp)] = 0
    tmp[np.abs(tmp) > 3*np.std(tmp)] = 0
    return np.std(tmp)

def datetime_to_decimal_year(dt: datetime) -> float:
    year_start = datetime(dt.year, 1, 1)
    year_end = datetime(dt.year + 1, 1, 1)

    year_length = (year_end - year_start).total_seconds()
    seconds_into_year = (dt - year_start).total_seconds()

    return dt.year + seconds_into_year / year_length

def getCalPixel(cal, wcs):

    if cal == '3C444':
        ra = 333.6071
        dec = -17.0267

    elif cal == 'Crab':
        ra = 83.6333
        dec = 22.0144
    else:
        print(cal)

    return radec_to_pixel(wcs, ra, dec)


# def flux_3C444(frequency_mhz):
#     S1, nu1 = 146, 74     # Jy, MHz
#     S2, nu2 = 79.6, 160   # Jy, MHz
#     alpha = (math.log(S2/S1) / math.log(nu2/nu1))
#     return S1 * (frequency_mhz / nu1) ** alpha

def flux_3C444(freq):
    s_array = np.array([39.5384187, 22.0061198, 4.5379081,
               3.4024946, 3.0192216, 2.7268868,
               2.3310504, 2.2057688, 0.9764044,
               0.4799112])
    a_array = np.array([-0.78, -0.78, -0.68,
                        -0.62, -0.77, -0.78,
                        -0.58, -0.56, -0.77,
                        -0.39])
    s = s_array*(freq/154)**a_array
    return np.sum(s)

def flagData(data):
    """
    Return
    """
    # remove inf, nans, threshold
    noise_at_pc = getNoise(data[1000-75:1000+75,
                                1000-75:1000+75])
    
    if not np.isfinite(noise_at_pc):
        return True
    
    if noise_at_pc > 5:
        return True

def sourceFindV2(img, hdu):

    ## get noise
    noise = getNoise(img)

    if noise > 5:
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

    if np.max(np.abs(img)) < 5:
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

    beamVolume = 1.1331*hdu[0].header['BMIN']*hdu[0].header['BMAJ']
    pix2deg = np.abs(hdu[0].header['CDELT1'])

    h, w = img.shape
    yc, xc = h // 2, w // 2  # image centre
    yy, xx = np.ogrid[:h, :w]
    dist = np.sqrt((yy - yc)**2 + (xx - xc)**2)

    img[dist > 18] = np.nan

    n_pix = np.sum(dist <= 18)
    n_beams = n_pix*pix2deg**2/beamVolume
    s_int = np.nansum(img)
    s_err = np.sqrt(n_beams)*noise

    return noise, _, s_int, np.sqrt(s_err**2+noise**2), _, _




def sourceFind(img, hdu):
    """
    Function takes in img cuttout and does souce finding
    """
    ## local_noise, s_peak, s_int, s_err, type, dist

    ## get noise
    noise = getNoise(img)

    ## create SNR
    snr = img/noise

    beamVolume = 1.1331*hdu[0].header['BMIN']*hdu[0].header['BMAJ']
    pix2deg = np.abs(hdu[0].header['CDELT1'])

    mask = (snr>=3)
    labels, num = ndimage.label(mask)



    ### search for positive island
    for i in range(1, num + 1):
        region = (labels == i)

        ## peak intensity
        s_peak = np.nanmax(region*img)
        s_sum = np.nansum(region*img)
        n_pix = np.sum(region)
        n_beams = n_pix*pix2deg**2/beamVolume

        s_int = s_sum*pix2deg**2/beamVolume
        s_err = np.sqrt(n_beams)*noise

        ## calculate centroid (RA, DEC)
        row, col = np.where(region*img > 0)
        row_mean = np.mean(row)
        col_mean = np.mean(col)

        dist_from_pc = np.sqrt((row_mean-100)**2 + (col_mean-100)**2)

        if dist_from_pc <= 20 and s_int >= 3*noise:

            ## local_noise, s_peak, s_int, s_err, type, dist
            #return noise, s_peak, s_int, s_err, 'positive', dist_from_pc
            return noise, s_peak, s_int, np.sqrt(s_err**2 + noise**2), 'positive', dist_from_pc
        
        
    ## search for negative moon
    img*=-1
    snr = img/noise
    mask = (snr>=3)
    labels, num = ndimage.label(mask)

    ### search for positive island
    for i in range(1, num + 1):
        region = (labels == i)

        ## peak intensity
        s_peak = np.nanmax(region*img)
        s_sum = np.nansum(region*img)
        n_pix = np.sum(region)
        n_beams = n_pix*pix2deg**2/beamVolume

        s_int = s_sum*pix2deg**2/beamVolume
        s_err = np.sqrt(n_beams)*noise

        ## calculate centroid (RA, DEC)
        row, col = np.where(region*img > 0)
        row_mean = np.mean(row)
        col_mean = np.mean(col)

        dist_from_pc = np.sqrt((row_mean-100)**2 + (col_mean-100)**2)

        if dist_from_pc <= 20 and s_int >= 3*noise:
            ## local_noise, s_peak, s_int, s_err, type, dist
            #return noise, -s_peak, -s_int, s_err, 'negative', dist_from_pc
            return noise, -s_peak, -s_int, np.sqrt(s_err**2 + noise**2), 'negative', dist_from_pc
        
        
    ## return noise limit
    return noise, np.nan, np.nan, np.nan, 'limit', np.nan


def getCalFlux(img, hdu):

    ## get noise
    noise = getNoise(img)

    ## create SNR
    snr = img/noise

    beamVolume = 1.1331*hdu[0].header['BMIN']*hdu[0].header['BMAJ']
    pix2deg = np.abs(hdu[0].header['CDELT1'])

    mask = (snr>=3)
    labels, num = ndimage.label(mask)


    s_int_array = []
    s_err_array = []

    ### search for positive island
    for i in range(1, num + 1):
        region = (labels == i)

        ## peak intensity
        s_peak = np.nanmax(region*img)
        s_sum = np.nansum(region*img)
        n_pix = np.sum(region)
        n_beams = n_pix*pix2deg**2/beamVolume
        s_err = np.sqrt(n_beams)*noise
        s_int = s_sum*pix2deg**2/beamVolume
        
        s_int_array.append(s_int)
        s_err_array.append(s_err)

    s_int_array=np.array(s_int_array)
    s_err_array = np.array(s_err_array)
    
    cal_flux = np.nanmax(s_int_array)
    mask = (s_int_array == cal_flux)
    return cal_flux, float(s_err_array[mask])


def moon_az_el(utc_time_str):
    # Parse time (assumes UTC)
    t = Time(utc_time_str, scale='utc')

    # Observer location
    location = EarthLocation(lon=116.67083333*u.deg, lat=-26.70331941*u.deg, height=377.827*u.m)

    # AltAz frame for that time/location, include pressure/temperature for refraction
    altaz = AltAz(obstime=t,
                  location=location)  

    # Get Moon coord and transform to local AltAz
    moon_icrs = get_moon(t)            # geocentric ICRS position of the Moon
    moon_altaz = moon_icrs.transform_to(altaz)

    az_deg = moon_altaz.az.to(u.deg).value
    el_deg = moon_altaz.alt.to(u.deg).value

    return az_deg, el_deg, moon_icrs.ra.deg, moon_icrs.ra.deg

## initialise output values
obs_ary = []
frq_ary = []
chn_ary = []
utc_ary = []
noise_loc_ary = []
s_peak_ary = []
s_int_ary = []
s_err_ary = []
bmj_ary = []
scale_ary = [] 
dist2pc_ary = []
type_ary = []
az_ary = []
el_ary = []

c = 0
for obs in obs_list:

    print('progrees {}/{}'.format(c, len(obs_list)))
    c += 1
    if obs in  crab_obs:
        continue

    ## load observation data
    hdu = fits.open('processing/{0}/{0}-img-narrowband-0000-image-pb.fits'.format(obs))
    utc = datetime.strptime(hdu[0].header['DATE-OBS'][:-2], '%Y-%m-%dT%H:%M:%S')
    wcs = WCS(hdu[0].header, naxis=2)
    pc_ra = hdu[0].header['CRVAL1'] +360
    pc_dec = hdu[0].header['CRVAL2']
    moon_ra, moon_dec, _ = getMoonRADEC(utc)
    moon_col, moon_row = radec_to_pixel(wcs, moon_ra, moon_dec)
    row = int(moon_row)
    col = int(moon_col)
    ## calculate distance of moon from pc
    dist = np.sqrt((pc_ra - moon_ra)**2 + (pc_dec - moon_dec)**2)
    
    print('pc ra {} dec {} moon ra {} dec {}\ndist {}'.format(pc_ra,pc_dec, moon_ra, moon_dec, dist))
    
    
    
    ## get moon az el
    az_deg, el_deg, ra_deg, dec_deg = \
        moon_az_el(utc)

    cal_col, cal_row = getCalPixel('3C444', wcs)
    cal_col = int(cal_col)
    cal_row = int(cal_row)

    ## go through the channels and get the flux density
    for channel in range(24):
        hdu = fits.open('processing/{0}/{0}-img-narrowband-{1}-image-pb.fits'.format(obs, str(channel).zfill(4)))
        data = hdu[0].data[0,0,:,:]
    
        if flagData(data):
            continue

        ### amp scaling (using model of calibrator)
        # cal_flux = np.nanmax(data[cal_row-5:cal_row+5,
        #                               cal_col-5:cal_col+5])
        cutout = Cutout2D(data, (cal_col, cal_row), (300, 300), wcs)
        cal_flux, cal_flux_err = getCalFlux(cutout.data,hdu)
        scaling = flux_3C444(hdu[0].header['CRVAL3']/1e6)/cal_flux
        data *= scaling
        
        ## make img cuttout
        cutout = Cutout2D(data, (col, row), (200, 200), wcs)
        img = cutout.data

        local_noise, s_peak, s_int, s_err, mtype, _ = sourceFind(img, hdu)
        s_err = np.sqrt(s_err**2 + cal_flux_err**2)

        ## append to array
        if np.isfinite(s_int):

            obs_ary.append(obs)
            frq_ary.append(hdu[0].header['CRVAL3']/1e6)
            chn_ary.append(channel)
            utc_ary.append(datetime_to_decimal_year(utc))
            noise_loc_ary.append(local_noise)
            s_peak_ary.append(s_peak)
            s_int_ary.append(s_int)
            s_err_ary.append(s_err)
            bmj_ary.append(hdu[0].header['BMAJ'])
            scale_ary.append(scaling)
            dist2pc_ary.append(dist)
            type_ary.append(mtype)
            az_ary.append(az_deg)
            el_ary.append(el_deg)

        
obs_ary=np.array(obs_ary)
frq_ary=np.array(frq_ary)
chn_ary=np.array(chn_ary)
utc_ary=np.array(utc_ary)
noise_loc_ary=np.array(noise_loc_ary)
s_peak_ary=np.array(s_peak_ary)
s_int_ary=np.array(s_int_ary)
s_err_ary=np.array(s_err_ary)
bmj_ary=np.array(bmj_ary)
scale_ary=np.array(scale_ary)
dist2pc_ary=np.array(dist2pc_ary)
type_ary=np.array(type_ary)
az_ary=np.array(az_ary)
el_ary=np.array(el_ary)

## stack and save to disk
data = np.column_stack((obs_ary, frq_ary, chn_ary,
    utc_ary, noise_loc_ary, s_peak_ary, s_int_ary,
    s_err_ary, bmj_ary, scale_ary, dist2pc_ary,
    type_ary, az_ary, el_ary))

np.save('group1data.npy', data)


