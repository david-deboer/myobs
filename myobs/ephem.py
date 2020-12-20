from astropy import units as u
from astropy.coordinates import Angle
import numpy as np
from . import dateutil
from astropy.convolution import convolve, Gaussian1DKernel


C0 = 299792458.0  # m/s


def to_Angle(val, unit='degree'):
    if isinstance(val, (u.quantity.Quantity, Angle)):
        return val
    return Angle(val, unit)


def to_Time(times, toffset=None):
    """Convert to Time."""
    return dateutil.get_astropytime(times, toffset)


class BaseEphem:
    param = ['times', 'ra', 'dec', 'az', 'el', 'x', 'y', 'z', 'D',
             'dt', 'dra', 'ddec', 'daz', 'del', 'dx', 'dy', 'dz', 'dD']

    def __init__(self):
        """
        Provides an init'd base class for ephemerides per:
            times :  Time
            dt :  np.array
            ra, dec, az, el:  Angle
            dra, ddec, daz, delev:  deg/sec (Quantity)
            x, y, z:  m (Quantity)
            dx, dy, dz:  m/s (Quantity)
        """
        self.initall()

    def initall(self):
        for par in self.param:
            setattr(self, par, [])

    def to_Time(self, times='times', toffset=None):
        """
        Convert the array to Time and set as self.times.

        If str, then convert class 'times' variable.
        """
        if isinstance(times, str) and times == 'times':
            times = self.times
        self.times = to_Time(times, toffset)

    def to_Angle(self, ang, angle_unit=u.degree, angv=None):
        """
        Convert attribute ang to class array.

        """
        if angv is None:
            angv = getattr(self, ang)
        setattr(self, ang, to_Angle(angv, angle_unit))

    def calc_dt(self):
        """
        Compute self.dt in sec as ndarray Quantity.
        """
        self.dt = [0.0]  # in seconds
        for i, _t in enumerate(self.times[1:]):
            self.dt.append((_t - self.times[i]).value * 3600.0 * 24.0)
        self.dt[0] = self.dt[1]
        self.dt = np.array(self.dt) * u.second

    def _smrt(self, arr, smooth):
        try:
            anit = arr.unit
        except AttributeError:
            anit = None
        if smooth is not None:
            if anit is not None:
                return u.Quantity(convolve(arr, Gaussian1DKernel(smooth), boundary='extend'), anit)
            else:
                return convolve(arr, Gaussian1DKernel(smooth), boundary='extend')
        else:
            return arr

    def vis(self, arr, val=0.0, horizon=0.0):
        """
        Get filter for visible above horizon.  Those below the horizon are set to val.

        Usage is to call e.g. visible_doppler = self.vis(self.doppler)
        """
        varr = np.array(arr)
        varr[np.where(self.el < horizon)] = val
        return varr

    def dbydt(self, par, smooth=None, unwrap=False):
        if self.dt is None or len(self.dt) != len(getattr(self, par)):
            self.calc_dt()
        deriv = f"d{par}"
        setattr(self, deriv, [0.0])
        if unwrap:
            _param = np.unwrap(getattr(self, par))
        else:
            _param = getattr(self, par)
        _param = self._smrt(_param, smooth)
        for i, _pp in enumerate(_param[1:]):
            getattr(self, deriv).append((_pp - _param[i])/self.dt[i+1])
        getattr(self, deriv)[0] = getattr(self, deriv)[1]
        setattr(self, deriv, u.Quantity(getattr(self, deriv)))
