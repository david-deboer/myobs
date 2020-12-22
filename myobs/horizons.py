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
    std_hdr = {'Date__(UT)__HR:MN': 0,
               'R.A._____(ICRF)_____DEC': 23,
               'APmag': None,
               'S-brt': None,
               'delta': 63,
               'deldot': 80
               }

    def __init__(self, fname, skip_verify, hdr=None):
        super().__init__()
        self.read(fname=fname, skip_verify=skip_verify, hdr=hdr)

    def read(self, fname, skip_verify=False, hdr=None):
        self.fname = fname
        if hdr is None:
            hdr = self.std_hdr
        self.hdr = hdr
        self.initall()
        SOE = False
        header_verified = skip_verify
        with open(fname, 'r') as fp:
            for line in fp:
                # HEADER
                if line.strip().startswith('Center-site'):
                    print(f"Observers location: {line.split(':')[1].strip()}")
                if not header_verified:
                    data = line.strip().split()
                    if len(data) >= len(self.hdr):
                        is_ok = True
                        for i, hk in enumerate(self.hdr.keys()):
                            if hk != data[i]:
                                is_ok = False
                        if is_ok:
                            header_verified
                if line.strip().startswith('$$SOE'):
                    SOE = True
                    continue
                if not SOE:
                    continue
                if not header_verified:
                    raise ValueError("Header not verified.")
                # HEADER
                if line.startswith('$$EOE') or len(line) < 80:
                    break
                # Get date/time
                data = line[self.hdr['Date__(UT)__HR:MN']:self.hdr['R.A._____(ICRF)_____DEC']].strip().split()  # noqa
                edate = data[0].split('-')
                edate[1] = mon[edate[1]]
                self.times.append(Time(f"{'-'.join(edate)} {data[1]}"))
                # Get RA
                data = line[self.hdr['R.A._____(ICRF)_____DEC']:self.hdr['APmag']].strip().split()
                hr, mn, sc = [float(_x) for _x in data[0:3]]
                self.ra.append(hr + mn/60 + sc/3600)
                # Get Dec
                sn = 1.0 if data[3][0] == '+' else -1.0
                dgr = float(data[3][1:])
                amn, asc = [float(_x) for _x in data[4:6]]
                self.dec.append(sn*(dgr + amn/60 + asc/3600))
                # Get range and velocity
                data = line[self.hdr['delta']:].strip().split()
                self.D.append(ephem.C.au * float(data[0]))
                self.Ddot.append(1000.0 * float(data[1]))
        self.times = Time(self.times)
        self.to_Angle('ra', angle_unit='hourangle')
        self.to_Angle('dec', angle_unit='degree')
        self.D = u.Quantity(self.D)
        self.Ddot = np.array(self.Ddot)*(u.m/u.s)
        self.x = self.D * np.cos(self.ra) * np.cos(self.dec)
        self.y = self.D * np.sin(self.ra) * np.cos(self.dec)
        self.z = self.D * np.sin(self.dec)
        smoothfactor = int(len(self.times) / 100.0)
        self.dbydt('ra', smoothfactor, True)
        self.dbydt('dec', smoothfactor, False)
        dra = self.radot.to('radian/s')
        ddec = self.decdot.to('radian/s')
        self.xdot = (self.D*(np.sin(self.ra)*np.cos(self.dec)*dra +
                             np.cos(self.ra)*np.sin(self.dec)*ddec))
        self.ydot = (self.D*(np.cos(self.ra)*np.cos(self.dec)*dra -
                             np.sin(self.ra)*np.sin(self.dec)*ddec))
        self.zdot = self.D * np.cos(self.dec)*ddec
        # Cast to m/s (as opposed to m rad / s) and include radial
        self.xdot = self.Ddot*np.cos(self.ra)*np.cos(self.dec) - self.xdot.value * u.m / u.second
        self.ydot = self.Ddot*np.sin(self.ra)*np.cos(self.dec) + self.ydot.value * u.m / u.second
        self.zdot = self.Ddot*np.sin(self.dec) + self.zdot.value * u.m / u.second
