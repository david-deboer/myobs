"""
Read cut-pasted jpl horizons files
"""
import numpy as np
from astropy.time import Time
import astropy.units as u
from . import ephem


AU = 149597870700.0  # m
mon = {'Jan': '01', 'Feb': '02', 'Mar': '03', 'Apr': '04', 'May': '05', 'Jun': '06',
       'Jul': '07', 'Aug': '08', 'Sep': '09', 'Oct': '10', 'Nov': '11', 'Dec': '12'}


class Horizons(ephem.BaseEphem):
    def __init__(self, fname):
        self.initall()
        self.f = ephem.BaseEphem()  # as read from file
        self.read(fname)
        self.at()

    def interp(self, par, times):
        """
        Interpolate attribute par onto times.

        par : str
            string of the attribute to use
        times : Time
            times to be interpolated onto.
        """
        clpar = getattr(self.f, f"{par}")
        if not len(clpar):
            return clpar
        try:
            parunit = clpar.unit
        except AttributeError:
            return np.interp(times.jd, self.f.times.jd, clpar)
        return u.Quantity(np.interp(times.jd, self.f.times.jd, clpar), parunit)

    def at(self, times=None):
        """
        Put all data at "times".  None to data read from "fname".
        """
        if times is None:
            self.times = self.f.times
            for par in self.param:
                if par != 'times':
                    setattr(self, par, getattr(self.f, f"{par}"))
        else:
            self.to_Time(times)
            for par in self.param:
                if par != 'times':
                    setattr(self, par, self.interp(par, times))

    def read(self, fname):
        self.fname = fname
        self.f.initall()
        with open(fname, 'r') as fp:
            for line in fp:
                if len(line) < 50:
                    break
                data = line.strip().split()
                # Get date/time
                edate = data[0].split('-')
                edate[1] = mon[edate[1]]
                self.f.times.append(Time(f"{'-'.join(edate)} {data[1]}"))
                # Get RA
                hr, mn, sc = [float(_x) for _x in data[2:5]]
                self.f.ra.append(hr + mn/60 + sc/3600)
                # Get Dec
                sn = 1.0 if data[5][0] == '+' else -1.0
                dgr = float(data[5][1:])
                amn, asc = [float(_x) for _x in data[6:8]]
                self.f.dec.append(sn*(dgr + amn/60 + asc/3600))
                # Get range and velocity
                self.f.D.append(AU * float(data[10]))
                self.f.dD.append(1000.0 * float(data[11]))
        self.f.times = Time(self.f.times)
        self.f.to_Angle('ra', angle_unit='hourangle')
        self.f.to_Angle('dec', angle_unit='degree')
        self.f.D, self.f.dD = np.array(self.f.D)*u.m, np.array(self.f.dD)*(u.m/u.s)
        self.f.x = self.f.D * np.cos(self.f.ra) * np.cos(self.f.dec)
        self.f.y = self.f.D * np.sin(self.f.ra) * np.sin(self.f.dec)
        self.f.z = self.f.D * np.sin(self.f.dec)
        smoothfactor = int(len(self.f.times) / 20.0)
        self.f.dbydt('ra', smoothfactor, True)
        self.f.dbydt('dec', smoothfactor, False)
        self.f.dx = (self.f.D * (-np.sin(self.f.ra)*np.cos(self.f.dec)*self.f.dra.to('radian/s') -
                     np.cos(self.f.ra)*np.sin(self.f.dec)*self.f.ddec.to('radian/s')))
        self.f.dy = (self.f.D * (np.cos(self.f.ra)*np.cos(self.f.dec)*self.f.dra.to('radian/s') -
                     np.sin(self.f.ra)*np.sin(self.f.dec)*self.f.ddec.to('radian/s')))
        self.f.dz = self.f.D * np.cos(self.f.dec)*self.f.ddec.to('radian/s')
        # Cast to m/s (as opposed to m rad / s)
        self.f.dx = self.f.dx.value * u.m / u.second
        self.f.dy = self.f.dy.value * u.m / u.second
        self.f.dz = self.f.dz.value * u.m / u.second
