from astropy import units as u
from astropy.coordinates import SkyCoord, AltAz, EarthLocation, Angle
import numpy as np
from . import dateutil, horizons, ephem
from argparse import Namespace


LOCATIONLIST = {
                'pks': {'name': 'Parkes NSW',
                        'lat': -32.995126867833136,
                        'lon': 148.26238459073372,
                        'alt': 480.1,
                        'key': 'pks'},
                'gbt': {'name': 'GBT WV',
                        'lat': 38.43286455266003,
                        'lon': -79.83971215724313,
                        'alt': 147.8,
                        'key': 'gbt'},
                '0-0': {'name': "Origin",
                        'lat': 0.0,
                        'lon': 0.0,
                        'alt': 0.0,
                        'key': '0-0'}
               }


def getloc(name):
    if name in LOCATIONLIST.keys():
        loc = LOCATIONLIST[name]
        loc['r'] = earth_radius(loc['lat']) + loc['alt']
        return loc
    print("Options:")
    for k, v in LOCATIONLIST.items():
        print(f"\t{k}:  {v}")
    return {}


def earth_radius(lat):
    """
    Return earth's radius in m given lat.
    """
    req, rpo = 6378137.0, 6356752.0
    lat = ephem.to_Angle(lat, 'degree').to('rad').value
    rec = req * np.cos(lat)
    rps = rpo * np.sin(lat)
    a = (req*rec)
    b = (rpo*rps)
    return np.sqrt((a**2 + b**2) / (rec**2 + rps**2))


def celestialxyz(ra, dec, distance, ra_unit=u.hr, dec_unit=u.degree):
    x = distance * np.cos(ephem.to_Angle(ra, ra_unit))*np.cos(ephem.to_Angle(dec, dec_unit))
    y = distance * np.sin(ephem.to_Angle(ra, ra_unit))*np.cos(ephem.to_Angle(dec, dec_unit))
    z = distance * np.sin(ephem.to_Angle(dec, dec_unit))
    return Namespace(x=x, y=y, z=z)


class Pointing(ephem.BaseEphem):
    def __init__(self, name, lat=None, lon=None, alt=None):
        """
        Pointing class that holds: EarthLocation, SkyCoord and AltAz.....

        Parameters
        ----------
        name : str
            Name of location
        lat : float
            latitude in deg
        lon : float
            longitude in deg
        alt : float
            altitude in m
        """
        super().__init__()
        alt_offset = 0.0
        if isinstance(name, dict):
            self.name, lon, lat, alt = name['name'], name['lon'], name['lat'], name['alt']
        elif isinstance(name, str):
            if alt is not None:
                alt_offset = 1.0 * alt
            if lon is None:
                if name in LOCATIONLIST.keys():
                    loc = LOCATIONLIST[name]
                    self.name, lon, lat, alt = loc['name'], loc['lon'], loc['lat'], loc['alt']
                else:
                    self.loc = EarthLocation.of_site(name)
                    return
            else:
                self.name = name
        self.loc = EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=(alt+alt_offset)*u.m)
        self.altaz = None

    def pointing(self, ra, dec, times=None, toffset=None, ra_unit=u.hourangle, dec_unit=u.deg):
        """
        Given ra,dec provide az,el.

        Parameters
        ----------
        ra : float, list, array
            Right Ascension of source
        dec : float, list, array
            Declination of source
        times :
            converted by get_astropytime
        toffset :
            offset to times per get_astropytime [hr]
        ra_unit : astropy unit
            unit of ra value
        dec_unit : astropy unit
            unit of dec value
        """
        if times is None:
            if self.altaz is None:
                self.altaz = AltAz(location=self.loc, obstime=self.times)
        else:
            self.times = dateutil.get_astropytime(times, toffset)
            self.altaz = AltAz(location=self.loc, obstime=self.times)
            self.calc_dt()
        if isinstance(ra, (int, float)):
            ra = [ra] * len(self.times)
        if isinstance(dec, (int, float)):
            dec = [dec] * len(self.times)
        self.ra = Angle(ra, ra_unit)
        self.dec = Angle(dec, dec_unit)
        ptg = SkyCoord(ra=self.ra, dec=self.dec, frame='icrs')
        self.az = ptg.transform_to(self.altaz).az
        self.el = ptg.transform_to(self.altaz).alt

    def geocentric(self, voffset=None, isolate_offset=False):
        self.gcrs = Namespace(r=None, v=None, Ddot=None, acc=None)
        if isolate_offset:
            if voffset is None:
                raise ValueError("Can't isolate 'None' offset for geocentric.")
            self.gcrs.Ddot = voffset
        else:
            self.gcrs.r, self.gcrs.v = self.loc.get_gcrs_posvel(self.times)
            self.gcrs.Ddot = (np.cos(self.ra)*np.cos(self.dec)*self.gcrs.v.x +
                              np.sin(self.ra)*np.cos(self.dec)*self.gcrs.v.y +
                              np.sin(self.dec)*self.gcrs.v.z)
            if voffset is not None:
                self.gcrs.Ddot += voffset
        self.gcrs.acc = np.diff(self.gcrs.Ddot.value) / self.dt[1:]
        self.gcrs.acc = np.insert(self.gcrs.acc, 0, self.gcrs.acc[0])

    def barycentric(self, baryfile='barycenter.dat'):
        # # Ideally do something like this:
        # from astropy.coordinates import solar_system_ephemeris
        # from astropy.coordinates import ICRS
        # with solar_system_ephemeris.set('jpl'):
        #     self.altaz.transform_to(ICRS(obstime=self.times))
        # self.icrs.r, self.icrs.v = self.loc.get_icrs_posvel(self.times)
        # # But instead for now...
        print("EPH235:  Be sure to check dates of barycenter file.")
        sun = horizons.Horizons(baryfile)
        self.icrs = Namespace(v=Namespace(x=None, y=None, z=None), Ddot=None, acc=[0.0])
        self.icrs.v.x = self.gcrs.v.x + sun.interp('xdot', self.times)
        self.icrs.v.y = self.gcrs.v.y + sun.interp('ydot', self.times)
        self.icrs.v.z = self.gcrs.v.z + sun.interp('zdot', self.times)
        self.icrs.Ddot = (np.cos(self.ra)*np.cos(self.dec)*self.icrs.v.x +
                          np.sin(self.ra)*np.cos(self.dec)*self.icrs.v.y +
                          np.sin(self.dec)*self.icrs.v.z)
        self.icrs.acc = np.diff(self.icrs.Ddot) / self.dt[1:]
        self.icrs.acc = np.insert(self.icrs.acc, 0, self.icrs.acc[0])

    def xyz2pointing(self, x, y, z, times=None, toffset=None, el=None):
        """
        Given x,y,z provide az,el,ra,dec.

        If el has already been computed, don't redo it.

        Parameters
        ----------
        x : array
            array of x values
        y : array
            array of y values
        z : array
            array of z values
        times :
            converted by get_astropytime
        toffset :
            offset to times per get_astropytime [hr]
        el : array or None
        """
        if times is None:
            if self.altaz is None:
                self.altaz = AltAz(location=self.loc, obstime=self.times)
        else:
            self.times = dateutil.get_astropytime(times, toffset)
            self.altaz = AltAz(location=self.loc, obstime=self.times)
        print("below doesn't work - but should be fairly straightforward")
        # ptg = SkyCoord(x=x*u.m, y=y*u.m, z=z*u.m, frame='icrs')
        # self.az = ptg.transform_to(self.altaz).az.value
        # self.el = ptg.transform_to(self.altaz).alt.value
        # self.ra = ptg.transform_to(self.altaz).ra.value
        # self.dec = ptg.transform_to(self.altaz).dec.value

    def pointing_at(self, az, el, times=None, toffset=None):
        """Given an az,el provide the ra,dec"""
        if times is None:
            if self.altaz is None:
                self.altaz = AltAz(location=self.loc, obstime=self.times)
        else:
            self.times = dateutil.get_astropytime(times, toffset)
            self.altaz = AltAz(location=self.loc, obstime=self.times)
        ptg = SkyCoord(az=az*u.deg, el=el*u.deg, frame='icrs')
        self.az = el
        self.el = az
        self.ra = ptg.transform_to(self.altaz).ra.value
        self.dec = ptg.transform_to(self.altaz).dec.value
