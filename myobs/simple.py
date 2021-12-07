import numpy as np
import matplotlib.pyplot as plt
from . import my_dateutil as dateutil


def simple(f=982e6, dec=-62, lat=-32):
    """
    Do the simple geometric calculations for sanity check.
    """
    print(f"Freq: {f/1e6:.1f} MHz")
    print(f"Lat: {lat:.2f} deg")
    print(f"Dec: {dec:.2f} deg")
    dec = np.deg2rad(dec)
    lat = np.deg2rad(lat)
    re = 6378000.0  # m
    w = (2.0 * np.pi) / (24 * 3600)  # rad/sec
    c0 = 3e8  # m/s

    print(f"Max rate: {re*w:.2f} m/s")
    print(f"Max accel:  {re*w*w:.5f} m/s/s")
    print(f"Mag df: {re*w*f/c0*np.cos(lat):.2f} Hz")
    print(f"Mag drift: {re*w*w*f/c0*np.cos(lat):.2f} Hz/s")

    # theta = lng + w(t - t0) - RA
    theta = np.arange(-np.pi, np.pi, np.pi/180.)
    t = theta / w
    dt = np.diff(t)
    ts = dateutil.get_astropytime('2019-04-29', t/3600.)

    r = re * np.cos(theta) * np.cos(lat) * np.cos(dec) + re * np.sin(lat) * np.sin(dec)
    dr = -w * re * np.sin(theta) * np.cos(lat) * np.cos(dec)
    ddr = -w * w * re * np.cos(theta) * np.cos(lat) * np.cos(dec)
    el = np.rad2deg(np.arcsin(r / re))

    df = (dr / c0) * f
    ddf = (ddr / c0) * f
    a_ddf = np.diff(df) / dt  # These should all be identical to ddf
    b_ddf = ((np.diff(dr) / dt) / c0) * f  # These should all be identical to ddf
    checksum = np.sum(np.fabs(a_ddf - ddf[1:]))
    if checksum > 0.1:
        print("Not identical: {}".format(checksum))
        print("...{}".format(np.sum(np.fabs(b_ddf-ddf[1:]))))

    plt.figure("drift")
    plt.plot(ts.datetime[np.where(el > 0.0)], ddf[np.where(el > 0.0)], label='calc')
    plt.legend()
    plt.figure("delta_f")
    plt.plot(t[np.where(el > 0.0)]/3600, df[np.where(el > 0.0)])

    plt.figure('el')
    plt.plot(t/3600, el)
