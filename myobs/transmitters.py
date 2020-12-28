from astropy.coordinates import EarthLocation
import astropy.units as u
import numpy as np
from . import ephem, observer
from argparse import Namespace


def noisedist(Tsys, bw, tau, N):
    p = ephem.Const.k_B.value * Tsys * bw
    s = np.sqrt(bw * tau)
    #  pn = np.random.rayleigh(p/s, N)
    pn = np.fabs(np.random.normal(p, p/s, N))
    return pn


def waterfall(t, Rxfreq, r, pwr=1.0/1000.0, Tsys=50.0, bw=1.0, minsmear=4.0, flo=None, fhi=None):
    return Waterfall(t=t, Rxfreq=Rxfreq, r=r, pwr=pwr, Tsys=Tsys, bw=bw, minsmear=minsmear,
                     flo=flo, fhi=fhi)


class Waterfall:
    def __init__(self, t, Rxfreq, r, pwr=1.0/1000.0, Tsys=50.0, bw=1.0, minsmear=4.0,
                 flo=None, fhi=None):
        """
        Parameters
        ----------
        t : list/array of float
            Times in sec
        Rxfreq : list/array of float
            Doppler shifted frequency [Hz]
        r : list/array of float
            Distance to transmitter [m]
        pwr : float
            EIRP [W]
        Tsys : float
            System temperature [K]
        bw : float
            Bandwidth [Hz]
        """
        if r is None:
            self.pwr = {0: [pwr] * len(Rxfreq)}
        else:
            self.pwr = {0: pwr / (4.0 * np.pi * r**2)}
        self.t = {0: t}
        self.f = {0: Rxfreq}
        self.auto_ctr = 1
        self.r = r
        self.bw = bw
        self.Tsys = Tsys
        self.minsmear = minsmear
        if flo is None:
            self.flo = 3.0*Rxfreq.min()/2.0 - Rxfreq.max()/2.0
        else:
            self.flo = flo
        if fhi is None:
            self.fhi = 3.0*Rxfreq.max()/2.0 - Rxfreq.min()/2.0
        else:
            self.fhi = fhi
        self.numch = int((self.fhi - self.flo)/bw)
        self.ch = np.linspace(self.flo, self.fhi, self.numch)
        self.wf = []
        self.tstart = t[0]
        self.tstop = t[-1]
        self.int_time = t[1] - t[0]
        print("<Waterfall>")
        print(f"Nch: {self.numch}, tau: {self.int_time}, N: {len(t)}")
        print(f"flo: {self.flo}, fhi: {self.fhi}")
        for i in range(len(Rxfreq)):
            self.wf.append(noisedist(Tsys, bw, self.int_time, self.numch))
        self.wf = np.array(self.wf)
        self.add_signal(key=0)

    def add_signal(self, t=None, f=None, p=None, key=None):
        if key is None:
            key = self.auto_ctr
        self.auto_ctr += 1
        if t is None:
            t = self.t[0]
        else:
            if np.fabs((t[1] - t[0]) - self.int_time) / self.int_time > 0.1:
                print("Timing doesn't match - skipping adding signal.")
                return
            self.t[key] = t
        if f is None:
            f = self.f[0]
        else:
            self.f[key] = f
        if p is None:
            p = self.pwr[0]
        else:
            if self.r is None:
                self.pwr[key] = [p] * len(self.f)
            else:
                self.pwr[key] = p / (4.0 * np.pi * self.r**2)
        chnum = (f - self.flo)/self.bw
        tbin = (t - self.tstart)/self.int_time
        for i, (_t, _f, _p) in enumerate(zip(t, f, p)):
            this_time = int(tbin[i])
            if this_time > len(f) - 2:
                this_time = len(f) - 2
            this_chan = int(chnum[i])
            smear = (f[this_time+1] - f[this_time])/self.bw
            if smear == 0.0:
                smear = self.minsmear
            elif abs(smear) < self.minsmear:
                smear = self.minsmear * np.sign(smear)
            for k in range(this_chan, this_chan+int(smear), int(np.sign(smear))):
                self.wf[this_time, k] += p[this_time]/abs(smear)


def _split(i, lon, lat, alt):
    lat1, lat2 = np.deg2rad(lat[i-1]), np.deg2rad(lat[i])
    dlat = lat2 - lat1
    lon1, lon2 = np.deg2rad(lon[i-1]), np.deg2rad(lon[i])
    dlon = lon2 - lon1
    dalt = alt[i] - alt[i-1]
    return lon1, lon2, dlon, lat1, lat2, dlat, dalt


class Moving(ephem.BaseEphem):
    def __init__(self, f=None, date='2020-12-11'):
        """
        Moving transmitter, ignore antenna movement as negligible.
        """
        self.freq = f
        self.date = date  # Not used, but could if wanted "real" datetimes.
        self.initall()

    def observatory(self, name, lon=None, lat=None, alt=None):
        self.obs = observer.Pointing(name=name, lat=lat, lon=lon, alt=alt)

    def calc_doppler_drift(self, smooth=10, drift_smooth=10):
        vel = self._smrt(self.Ddot, smooth)  # additional smoothing in velocity
        self.doppler = (vel/self.c0) * self.freq
        self.drift = [0.0]
        for i in range(1, len(self.doppler)):
            self.drift.append(self.doppler[i].value-self.doppler[i-1].value)
        self.drift[0] = self.drift[1]
        self.drift = self._smrt(np.array(self.drift), drift_smooth) / self.dt.value

    def trajectory(self, lon, lat, alt, utc=12.0, speed=None, total_time=None, dt=10.0, smooth=10):
        """
        Take a list of waypoints and compute trajectory every dt sec over a great circle.

        Parameters
        ----------
        lon : list of float
            longitude in degrees
        lat : list of float
            latitude in degrees
        alt : list of float
            altitude in m
        utc : list of float, float
            utc.  If speed is not None, float is starting UTC.
        speed : None, float
            speed in m/s or None.  If None just use utc list, else calculate
            where utc is starting time.  Either speed or total_time can be specified,
            but not both.
        total_time : None, float
            total_time in sec or None (see above).
            Either speed or total_time can be specified, but not both.
        dt : float
            time interval for trajectory in s
        smooth : int or None
            smoothing factor for velocity calc
        """
        self.initall()
        calc_utc = False
        if speed is not None and total_time is not None:
            raise ValueError("Can't specify speed _and_ total_time.")
        if speed is not None or total_time is not None:
            calc_utc = True
            if not isinstance(utc, (float, int)):
                print("utc should just be starting time if given speed/total_time "
                      "-- using first term")
                utc = [utc[0]]
            else:
                utc = [utc]
        # Calculate waypoints
        self.waypt = Namespace(lon=lon, lat=lat, alt=alt, distance=[0.0], speed=[0.0])
        self.total_distance = 0.0
        Rearth = [observer.earth_radius(lat[0])]
        ca = [0.0]
        for i in range(1, len(lon)):
            lon1, lon2, dlon, lat1, lat2, dlat, dalt = _split(i, lon, lat, alt)
            _a = np.sin(dlat/2.0)**2 + np.cos(lat1)*np.cos(lat2)*np.sin(dlon/2.0)**2
            _rth = observer.earth_radius((lat1+lat2)/2.0)
            _ca = 2.0 * np.arctan2(np.sqrt(_a), np.sqrt(1.0 - _a))
            _D = np.sqrt((_rth*_ca)**2 + dalt**2)
            Rearth.append(_rth)
            ca.append(_ca)
            self.waypt.distance.append(_D)
            self.total_distance += _D
        print("Total distance is {:.1f} mi".format(self.total_distance / 1609.344))
        if total_time is not None:
            speed = self.total_distance / total_time
        self.traj = Namespace(utc=np.array([]), lon=np.array([]),
                              lat=np.array([]), alt=np.array([]))
        self.waypoint_indices = []
        for i in range(1, len(lon)):
            if calc_utc:
                dt_waypt = self.waypt.distance[i] / speed
                utc.append(utc[i-1] + dt_waypt/3600.0)
            else:
                dt_waypt = (utc[i] - utc[i-1]) * 3600.0
                speed = self.waypt.distance[i] / dt_waypt
            self.waypt.speed.append(speed)
            # Calculate waypoints into trajectory at dt
            lon1, lon2, dlon, lat1, lat2, dlat, dalt = _split(i, lon, lat, alt)
            _altinc = self.waypt.distance[i] * np.sin(np.arctan2(dalt, Rearth[i]*ca[i]))
            N_time = int(np.floor(dt_waypt / dt))
            if np.fabs(N_time) < 1.0:
                continue
            _f = np.linspace(0.0, 1.0-1.0/N_time, N_time)
            _ta, _tb = np.sin((1.0 - _f)*ca[i]) / np.sin(ca[i]), np.sin(_f*ca[i]) / np.sin(ca[i])
            _x = _ta*np.cos(lat1)*np.cos(lon1) + _tb*np.cos(lat2)*np.cos(lon2)
            _y = _ta*np.cos(lat1)*np.sin(lon1) + _tb*np.cos(lat2)*np.sin(lon2)
            _z = _ta*np.sin(lat1) + _tb*np.sin(lat2)
            self.traj.utc = np.append(self.traj.utc, utc[i-1] +
                                      (dt/3600.0) * np.linspace(0, N_time-1, N_time))
            self.traj.lat = np.append(self.traj.lat, np.arctan2(_z, np.sqrt(_x**2+_y**2)))
            self.traj.lon = np.append(self.traj.lon, np.arctan2(_y, _x))
            self.traj.alt = np.append(self.traj.alt, alt[i-1] + _f * _altinc)
            self.waypoint_indices.append(len(self.traj.utc)-1)
        self.waypt.utc = utc
        self.traj.lon = np.rad2deg(self.traj.lon)
        self.traj.lat = np.rad2deg(self.traj.lat)
        self.traj.loc = EarthLocation(lon=self.traj.lon*u.deg, lat=self.traj.lat*u.deg,
                                      height=self.traj.alt*u.m)
        delta_x = self.traj.loc.x - self.obs.loc.x
        delta_y = self.traj.loc.y - self.obs.loc.y
        delta_z = self.traj.loc.z - self.obs.loc.z
        cosang = ((self.traj.loc.x*self.obs.loc.x +
                   self.traj.loc.y*self.obs.loc.y +
                   self.traj.loc.z*self.obs.loc.z) /
                  (np.sqrt(self.traj.loc.x**2 + self.traj.loc.y**2 + self.traj.loc.z**2) *
                   np.sqrt(self.obs.loc.x**2 + self.obs.loc.y**2 + self.obs.loc.z**2)))
        rtraj = np.sqrt(self.traj.loc.x**2 + self.traj.loc.y**2 + self.traj.loc.z**2)
        robs = np.sqrt(self.obs.loc.x**2 + self.obs.loc.y**2 + self.obs.loc.z**2)
        self.D = np.sqrt(delta_x**2 + delta_y**2 + delta_z**2)
        self.el = np.arcsin((rtraj * cosang - robs)/self.D)
        self.utc = self.traj.utc
        self.dt = np.array([dt] * len(self.traj.utc)) * u.second
        self.dbydt('D', smooth, False)


class Local(ephem.BaseEphem):
    def __init__(self, f, E, N, H=18.0, R=34.0):
        """
        Fixed transmitter and moving antenna.

        Parameters
        ----------
        E  m - distance E of the antenna
        N  m - distance N of the antenna
        H  m - alt of rotation point; Parkes~18
        R  m - length of rotation arm; Parkes~34
        """
        self.initall()
        self.freq = f  # Hz - Frequency
        self.E, self.N, self.H, self.R = E, N, H, R
        self.X2 = E*E + N*N + H*H + R*R

    def trajectory(self, cls, smooth=None):
        for par in self.param:
            setattr(self, par, getattr(cls, par))
        self.D = np.sqrt(self.X2 - 2.0*self.R*(self.E*np.cos(self.az)*np.cos(self.el) +
                                               self.N*np.sin(self.az)*np.cos(self.el) -
                                               self.H*np.sin(self.el)))
        self.dbydt('D', smooth)

    def calc_doppler_drift(self):
        self.doppler = (self.Ddot/self.c0) * self.freq
        self.drift = [0.0]
        for i in range(1, len(self.doppler)):
            self.drift.append(self.doppler[i].value-self.doppler[i-1].value)
        self.drift[0] = self.drift[1]
        self.drift = np.array(self.drift) / self.dt.value

    def calc_adDdt(self, dazdt=None, deldt=None, n=None):
        """
        Given D and UTC analytically calculate per above the rate of change in distance.
        """
        dratio = self.R / self.D
        if dazdt is not None and deldt is not None:
            dx = ((self.E*np.cos(self.az)*np.cos(self.el) -
                   self.N*np.sin(self.az)*np.cos(self.el))*np.deg2rad(dazdt) +
                  (self.H*np.cos(self.el) - self.E*np.sin(self.az)*np.sin(self.el) -
                   self.N*np.cos(self.az)*np.sin(self.el))*np.deg2rad(deldt)) * dratio
            self.adDdt = dx
