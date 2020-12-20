import numpy as np
from argparse import Namespace
from datetime import datetime
from . import ephem, observer
import astropy.units as u


class Track:
    satpar = ['since', 'x', 'y', 'z', 'vx', 'vy', 'vz']

    def __init__(self, fname):
        self.fname = fname
        for sp in self.satpar:
            setattr(self, sp, [])
        self.dtime = []
        self.read_file()

    def read_file(self):
        """
        Read in the .out file generated by satpos.
        """
        self.headerlines = []
        with open(self.fname, 'r') as fp:
            for line in fp:
                if line[0] == '#':
                    self.headerlines.append(line.strip()[1:])
                    continue
                data = line.strip().split()
                df = [float(data[i]) for i in range(len(self.satpar))]
                for i, sp in enumerate(self.satpar):
                    getattr(self, sp).append(df[i])
                yr, mn, dy = int(data[14]), int(data[15]), int(data[16])
                tod = line[199:].strip().split(':')
                hr, mt = int(tod[0]), int(tod[1])
                sc, ms = int(tod[2].split('.')[0]), int(tod[2].split('.')[1])
                self.dtime.append(datetime(year=yr, month=mn, day=dy,
                                  hour=hr, minute=mt, second=sc, microsecond=ms))
        for sp in self.satpar:
            if sp == 'since':
                self.since = np.array(self.since) * 60.0  # convert to sec
            elif sp in ['x', 'y', 'z']:
                setattr(self, sp, np.array(getattr(self, sp))*1000.0*u.m)  # convert to m
            elif sp in ['vx', 'vy', 'vz']:
                setattr(self, sp, np.array(getattr(self, sp))*1000.0*u.m/u.s)  # convert to m
        self.parse_header()

    def parse_header(self, default_epoch='2019-04-29 '):
        self.headers = {'scname': Namespace(text='spacecraft name', type=str, nval=0),
                        'satnum': Namespace(text='satnum', type=int, nval=0),
                        'period': Namespace(text='period', type=float, nval=0),
                        'sublon': Namespace(text='starting lon lat h', type=float, nval=0),
                        'sublat': Namespace(text='starting lon lat h', type=float, nval=1),
                        'height': Namespace(text='starting lon lat h', type=float, nval=2),
                        'epoch': Namespace(text='starting epoch', type=str, nval=0),
                        'scenario': Namespace(text='scenario epoch', type=float, nval=0)}
        for hdrline in self.headerlines:
            for hdr, valns in self.headers.items():
                if valns.text in hdrline:
                    if hdr == 'scname' or hdr == 'epoch':
                        data = hdrline.split(':')[1].strip()
                        if data.startswith('0 '):
                            data = data[2:]
                        data = [data]
                    else:
                        try:
                            data = hdrline.split(':')[1].strip().split()
                        except IndexError:
                            continue
                    setattr(self, hdr, valns.type(data[valns.nval]))
                    self.headers[hdr].value = valns.type(data[valns.nval])

                    break

    def calc(self, loc, freq):
        """Pipeline for general processing."""
        self.view(loc)
        self.rates(freq)
        self.subsat()

    def location(self, name, lon=None, lat=None, alt=None):
        """Make a location instance."""
        self.location = observer.Pointing(name, lon, lat, alt)
        self.loc = Namespace(x=self.location.loc.x,
                             y=self.location.loc.y,
                             z=self.location.loc.z)

    def view(self, name=None, lon=None, lat=None, alt=None):
        """Compute distance, enu, az/el and if ever viewable."""
        if name is not None:
            self.location(name, lon, lat, alt)
        self.R = Namespace(x=(self.x-self.loc.x),
                           y=(self.y-self.loc.y),
                           z=(self.z-self.loc.z))
        self.D = np.sqrt(self.R.x**2 + self.R.y**2 + self.R.z**2)
        robs = np.sqrt(self.loc.x**2 + self.loc.y**2 + self.loc.z**2)
        cxyz = self.loc.x*self.R.x + self.loc.y*self.R.y + self.loc.z*self.R.z
        self.el = np.rad2deg(np.arcsin(cxyz / (robs*self.D)))
        self.viewable = len(np.where(self.el > 0.0)[0]) > 1
        self.enu = Namespace()
        slon = np.sin(np.deg2rad(self.location.loc.lon))
        clon = np.cos(np.deg2rad(self.location.loc.lon))
        slat = np.sin(np.deg2rad(self.location.loc.lat))
        clat = np.cos(np.deg2rad(self.location.loc.lat))
        self.enu.x = -slon*self.R.x + clon*self.R.y
        self.enu.y = -slat*clon*self.R.x - slat*slon*self.R.y + clat*self.R.z
        self.enu.z = clat*clon*self.R.x + clat*slon*self.R.y + slat*self.R.z
        alp = np.arctan2(self.enu.y, self.enu.x)
        self.az = 90.0*u.deg - alp
        self.az[np.where(alp > 90*u.deg)] = 360.0*u.deg + self.az[np.where(alp > 90*u.deg)]

    def vis(self, arr, val=0.0, horizon=0.0):
        """
        Get filter for visible above horizon.  Those below the horizon are set to val.

        Usage is to call e.g. visible_doppler = self.vis(self.doppler)
        """
        varr = np.array(arr)
        varr[np.where(self.el < horizon)] = val
        return varr

    def rates(self, f=982E6):
        self.freq = f
        self.V = (self.vx*self.R.x + self.vy*self.R.y + self.vx*self.R.z) / self.D
        self.doppler = (np.array(self.V) / ephem.C0) * f
        self.drift = [0.0]
        for i in range(1, len(self.doppler)):
            dt = self.since[i] - self.since[i-1]
            dd = self.doppler[i] - self.doppler[i-1]
            self.drift.append(dd/dt)
        self.drift[0] = self.drift[1]

    def subsat(self):
        """Calculate the subsatellite point.  self.lon/self.lat"""
        self.lon = np.rad2deg(np.arctan2(self.y, self.x))
        self.rsat = np.sqrt(self.x**2 + self.y**2 + self.z**2)
        r = self.z / self.rsat
        self.lat = np.rad2deg(np.arcsin(r))

    def footprint(self, i):
        """
        Return array circumscribing the footprint for that subsat point/distance.
        """
        return np.array(footprint(self.lon[i], self.lat[i], self.rsat[i]))

    def waterfall(self, pwr=1.0, Tsys=50.0, bw=1.0, show_below=False):
        """Class wrapper to get waterfall."""
        from my_utils import transmitters
        t = self.since - self.since[0]
        f = self.freq + self.doppler
        r = 1.0 * self.D
        if not show_below:
            r[np.where(self.el < 0.0)] = 1E38
        self.ch, self.wf = transmitters.waterfall(t, f, r, pwr, Tsys, bw)

    def get_pointing(self):
        times = ephem.to_Time(self.dtime)
        self.location.xyz2pointing(self.x, self.y, self.z, times, self.el)
        print("This will yield self.loc.az, self.loc.el, self.loc.ra, self.loc.dec")


def footprint(lngs, lats, rsat):
    """Actual implementation of footprint."""
    Rearth = observer.earth_radius(lats)
    lngs = np.deg2rad(lngs)
    lats = np.deg2rad(lats)
    cg = Rearth/rsat
    g = np.arccos(cg)
    footprint = []
    for latf in np.arange(lats-g, lats+g, 0.01):
        cll = cg/np.cos(latf)/np.cos(lats) - np.tan(latf)*np.tan(lats)
        lngf = lngs - np.arccos(cll)
        footprint.append([np.rad2deg(lngf), np.rad2deg(latf)])
    for latf in np.arange(lats+g, lats-g, -0.01):
        cll = cg/np.cos(latf)/np.cos(lats) - np.tan(latf)*np.tan(lats)
        lngf = lngs + np.arccos(cll)
        footprint.append([np.rad2deg(lngf), np.rad2deg(latf)])
    latf = lats-g
    cll = cg/np.cos(latf)/np.cos(lats) - np.tan(latf)*np.tan(lats)
    lngf = lngs - np.arccos(cll)
    footprint.append([np.rad2deg(lngf), np.rad2deg(latf)])
    return footprint
