"""
Read cut-pasted jpl horizons files
"""
import numpy as np
from astropy.time import Time
import astropy.units as u
from . import ephem


mon = {'Jan': '01', 'Feb': '02', 'Mar': '03', 'Apr': '04', 'May': '05', 'Jun': '06',
       'Jul': '07', 'Aug': '08', 'Sep': '09', 'Oct': '10', 'Nov': '11', 'Dec': '12'}


class Horizons(ephem.BaseEphem):
    def __init__(self, fname):
        self.initall()
        self.read(fname)

    def read(self, fname):
        self.fname = fname
        self.initall()
        with open(fname, 'r') as fp:
            for line in fp:
                if len(line) < 50:
                    break
                data = line.strip().split()
                # Get date/time
                edate = data[0].split('-')
                edate[1] = mon[edate[1]]
                self.times.append(Time(f"{'-'.join(edate)} {data[1]}"))
                # Get RA
                hr, mn, sc = [float(_x) for _x in data[2:5]]
                self.ra.append(hr + mn/60 + sc/3600)
                # Get Dec
                sn = 1.0 if data[5][0] == '+' else -1.0
                dgr = float(data[5][1:])
                amn, asc = [float(_x) for _x in data[6:8]]
                self.dec.append(sn*(dgr + amn/60 + asc/3600))
                # Get range and velocity
                self.D.append(ephem.C.au * float(data[10]))
                self.Ddot.append(1000.0 * float(data[11]))
        self.times = Time(self.times)
        self.to_Angle('ra', angle_unit='hourangle')
        self.to_Angle('dec', angle_unit='degree')
        self.D, self.Ddot = np.array(self.D)*u.m, np.array(self.Ddot)*(u.m/u.s)
        self.x = self.D * np.cos(self.ra) * np.cos(self.dec)
        self.y = self.D * np.sin(self.ra) * np.sin(self.dec)
        self.z = self.D * np.sin(self.dec)
        smoothfactor = int(len(self.times) / 20.0)
        self.dbydt('ra', smoothfactor, True)
        self.dbydt('dec', smoothfactor, False)
        self.xdot = (self.D * (-np.sin(self.ra)*np.cos(self.dec)*self.radot.to('radian/s') -
                     np.cos(self.ra)*np.sin(self.dec)*self.decdot.to('radian/s')))
        self.ydot = (self.D * (np.cos(self.ra)*np.cos(self.dec)*self.radot.to('radian/s') -
                     np.sin(self.ra)*np.sin(self.dec)*self.decdot.to('radian/s')))
        self.zdot = self.D * np.cos(self.dec)*self.decdot.to('radian/s')
        # Cast to m/s (as opposed to m rad / s)
        self.xdot = self.xdot.value * u.m / u.second
        self.ydot = self.ydot.value * u.m / u.second
        self.zdot = self.zdot.value * u.m / u.second
