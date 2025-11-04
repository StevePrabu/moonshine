import numpy as np
import matplotlib.pyplot as plt

def binning(t,y,dy,width):
    """
    simple function to bin lightcurves with errors.
    
    Parameters
    ----------
    t : array, time 
    
    y : array, rate/flux
    
    dy : array, rate/flux uncertainty
    
    width : number, width of each bin in time (in the same unit as the input t array)
    
    Returns
    -------
    binned_t : binned time array
    
    binned_y : binned rate/flux array
    
    binned_dy : binned rate/flux uncertainty array
    
    """
    width = float(width)
    nbins = int(np.ceil((t[-1]-t[0])/width))+1
    binned_y = np.zeros(nbins)
    binned_y_weights = np.zeros(nbins)
    
    for i in range(len(y)):
        bin_no = int((t[i]-t[0])/width)
        weight = dy[i]**-2.0
        binned_y[bin_no] += y[i] * weight
        binned_y_weights[bin_no] += weight
        
    binned_y /= binned_y_weights
    binned_dy = np.sqrt(1.0/(binned_y_weights)) 
    binned_t = np.arange(start=t[0]+(width/2.0),stop=t[0]+nbins*width,step=width)
    return binned_t, binned_y, binned_dy

def getMoonValues(freq):

    t_moon = np.random.uniform(150, 230, 1)+25.3*(freq/150)**-2.5
    t_sky_60 = np.random.uniform(2340, 3206, 1)
    sky_alpha = np.random.uniform(-2.364, -2.9, 1)

    h = 6.62607015e-34
    k = 1.380649e-23
    c = 299792458
    omega = 6.4*1e-5
    t_sky = np.random.uniform(3206, 2340,1)*(freq/60)**np.random.uniform(-2.364, -2.9, 1)
    s = 2*k*(t_moon-t_sky)*omega*(freq*1e6/c)**2*1*1e26
    return s

def getMaxDist(f):
    c = 299792458
    wavelenght = c/(f*1e6)   

    fwhm = 1.22*np.degrees(wavelenght/5.5)
    return fwhm/2

data = np.load('group1data.npy')

dist = np.array(data[:,10], dtype=np.float32)
freq = np.array(data[:,1], dtype=np.float32)
s_int = np.array(data[:,6], dtype=np.float32)
s_err = np.array(data[:,7], dtype=np.float32)

filt_freq = []
filt_s = []
filt_s_err = []


for d,f, s, sr in zip(dist, freq, s_int, s_err):
    if d < getMaxDist(f):
        filt_freq.append(f)
        filt_s.append(s)
        filt_s_err.append(sr)
    else:
        pass

filt_freq = np.array(filt_freq)
filt_s = np.array(filt_s)
filt_s_err = np.array(filt_s_err)

idx = np.argsort(filt_freq)
freq_sort = filt_freq[idx]
s_sort = filt_s[idx]
s_err_sort = filt_s_err[idx]

bin_freq, bin_s, bin_s_err = binning(freq_sort, s_sort, s_err_sort,1)

plt.errorbar(bin_freq, bin_s, yerr=bin_s_err, fmt='.k', alpha=1)
#plt.errorbar(filt_freq, filt_s, yerr=filt_s_err, fmt='.', color='grey', alpha=0.2)
plt.hlines(0, 70, 250, color='black')
plt.grid(linestyle='dashed')
plt.xlim(70, 250)
plt.xscale('log')

for i in range(100):
    plt.plot(freq, getMoonValues(freq), color='blue', alpha=0.2, zorder=-1)

plt.savefig('test.png')