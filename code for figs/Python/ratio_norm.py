import numpy as np
from scipy.special import erf

def phi(x):
    return 0.5*(1.+erf(x/np.sqrt(2.)))
    
# PDF for ratio of two normal random variables N(mean1,sdev1^2) / N(mean2,sdev2^2)
def pdf(xpoints, (mean1, sdev1), (mean2, sdev2)):
    xpoints = np.array(xpoints)
    a = np.sqrt( (xpoints**2)/(sdev1**2) + 1./(sdev2**2))
    b = mean1/(sdev1**2)*xpoints + mean2/(sdev2**2)
    c = (mean1**2)/(sdev1**2) + (mean2**2)/(sdev2**2)
    d = np.exp( 0.5*(b**2-c*(a**2)) / (a**2) )
    p = b*d/(a**3) * 1./(np.sqrt(2.*np.pi)*sdev1*sdev2) * (phi(b/a) - phi(-b/a)) + 1./(np.pi*sdev1*sdev2*(a**2))*np.exp(-0.5*c)
    
    return p
