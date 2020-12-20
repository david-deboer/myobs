"""Miscellaneous collection of modules to count satellite/etc"""
from argparse import Namespace
from . import sattrack
from os.path import join


def filter_drift(f=982e6, rng=[0.025, 0.05], absvalue=True,
                 trackfilelist='viewable.csv', path='output', addtag='out'):
    """
    Uses the satellites written in find_viewable below and
    filters on drift range from the .out files.

    Writes outsats.csv with the list and returns a dict of the relevant portion
    of track.
    """
    outsats = open('outsats.csv', 'w')
    drifts = {}
    with open(trackfilelist, 'r') as satslist:
        loc = satslist.readline().strip()
        print("Location: ", loc)
        for line in satslist:
            data = line.strip().split(',')
            if data[0] == 'file':
                continue
            fname = data[0]
            if addtag is not None:
                fname = f"{fname}.{addtag}"
            fname = join(path, fname)
            s = sattrack.Track(fname)
            s.view(loc)
            s.rates(f)
            for i, drift in enumerate(s.vis(s.drift)):
                chkdrift = drift
                if absvalue:
                    chkdrift = abs(drift)
                if s.za[i] < 90.0 and chkdrift > rng[0] and chkdrift < rng[1]:
                    drifts.setdefault(s.satnum, Namespace(t=[], drift=[], za=[], period=[]))
                    drifts[s.satnum].t.append(s.dtime[i])
                    drifts[s.satnum].drift.append(drift)
                    drifts[s.satnum].za.append(s.za[i])
                    drifts[s.satnum].period.append(s.period)
                    print(f"{fname},{s.satnum},{s.since[i]},{drift},{s.x[i]},{s.y[i]},{s.z[i]},"
                          f"{s.za[i]},{s.period}", file=outsats)
    outsats.close()
    return drifts


def showdist():
    """
    Reads in viewable.csv and notviewable.csv and shows various things about them.
    """
    import numpy as np
    view = Namespace(iarr=[], period=[], zmin=[], zmax=[])
    notv = Namespace(iarr=[], period=[], zmin=[], zmax=[])
    with open('viewable.csv', 'r') as fp:
        loc = fp.readline().strip()
        print("Location ", loc)
        for line in fp:
            data = line.split(',')
            view.iarr.append(int(data[1].split()[0]))
            view.period.append(float(data[1].split('=')[1].split()[0]))
            view.zmin.append(float(data[2]))
            view.zmax.append(float(data[3]))
    view.iarr = np.array(view.iarr)
    view.period = np.array(view.period)
    view.zmin = np.array(view.zmin)
    view.zmax = np.array(view.zmax)

    with open('notviewable.csv', 'r') as fp:
        loc = fp.readline().strip()
        for line in fp:
            data = line.split(',')
            notv.iarr.append(int(data[1].split()[0]))
            notv.period.append(float(data[1].split('=')[1].split()[0]))
            try:
                notv.zmin.append(float(data[2]))
                notv.zmax.append(float(data[3]))
            except ValueError:
                notv.zmin.append(0.0)
                notv.zmax.append(0.0)
    notv.iarr = np.array(notv.iarr)
    notv.period = np.array(notv.period)
    notv.zmin = np.array(notv.zmin)
    notv.zmax = np.array(notv.zmax)


def find_viewable(loc, trackfilelist='spactive.list', path='output', verbose=False):
    """
    Run satpos over the desired TLEs ('satpos active 1', ...)
    and ls that list in the path.

    This writes viewable.csv and notviewable.csv of those tracks for loc.
    """
    viewable = open('viewable.csv', 'w')
    print(loc, file=viewable)
    notviewable = open('notviewable.csv', 'w')
    print(loc, file=notviewable)
    hdrstr = "file,scname,satnum,orbit,period,sublon,zamin,zamax"
    print(hdrstr, file=viewable)
    print(hdrstr, file=notviewable)
    count = Namespace(leo=0, meo=0, geo=0, deep=0, other=0, viewable=0, notviewable=0)
    with open(trackfilelist, 'r') as fp:
        i = 0
        for line in fp:
            fname = line.strip()
            if verbose:
                print(f"Reading {fname}")
            s = sattrack.Track(join(path, fname))
            s.view(loc)
            if s.period > 1500.0:
                count.deep += 1
                orbit = 'deep'
            elif s.period < 1450.0 and s.period > 1420.0:
                count.geo += 1
                orbit = 'geo'
            elif s.period < 1000.0 and s.period > 500.0:
                count.meo += 1
                orbit = 'meo'
            elif s.period < 200.0:
                count.leo += 1
                orbit = 'leo'
            else:
                count.other += 1
                orbit = 'other'
            try:
                zamin = s.za.min()
                zamax = s.za.max()
            except ValueError:
                zamin = '!'
                zamax = '!'
            fnp = fname.split('.')[0]
            pline = f"{fnp},{s.scname},{s.satnum},{orbit},{s.period},{s.sublon},{zamin},{zamax}"
            if s.viewable:
                print(pline, file=viewable)
                count.viewable += 1
            else:
                print(pline, file=notviewable)
                count.notviewable += 1
            i += 1
    viewable.close()
    notviewable.close()

    print(f"LEO: {count.leo}")
    print(f"MEO: {count.meo}")
    print(f"GEO: {count.geo}")
    print(f"DEEP: {count.deep}")
    print(f"OTHER: {count.other}")
    print(f"VIEWABLE: {count.viewable}")
    print(f"NOTVIEWABLE: {count.notviewable}")


def satpos_script(tlefile):
    """
    Write a bash script to check entry numbers within a tle file repetitiously via satpos.
    """
    if '.' not in tlefile:
        tlefile = f"{tlefile}.tle"
    tprename = tlefile.split('.')[0]
    tot = 0
    with open(tlefile, 'r') as fp:
        for line in fp:
            tot += 1
    tot = int(tot / 3.0)
    outfile = f'satpos_{tprename}.sh'
    print(f"Writing {tot} entries to {outfile}.")
    with open(outfile, 'w') as fp:
        for i in range(tot):
            print(f"satpos {tprename} {i+1}", file=fp)


def generate_complete_set(epoch=None, path='tle', fmname='master.dat'):
    """
    Goes through all tle files to find unique satellites and produce a completeset.
    [N.B. should use most recent epoch -- currently doesn't check that.]
    """
    import os.path
    satellites = {}
    sats_by_file = {}
    total_count = 0
    flistname = os.path.join(path, fmname)
    with open(flistname, 'r') as fp:
        for line in fp:
            fname = os.path.join(path, f"{line.strip().split(':')[0]}")
            sats_by_file[fname] = []
            with open(fname, 'r') as fptle:
                for tleline in fptle:
                    data = tleline.split()
                    if not len(data):
                        continue
                    if data[0] not in ['1', '2']:
                        scname = tleline.strip()
                        line0 = tleline.strip('\n')
                    elif data[0] == '1':
                        line1 = tleline.strip('\n')
                        this_epoch = float(line1[3])
                    elif data[0] == '2':
                        line2 = tleline.strip('\n')
                        key = data[1]
                        total_count += 1
                        satellites.setdefault(key, {'scname': scname, 'files': [], 'epoch': epoch})
                        if abs(this_epoch-epoch) <= abs(this_epoch-satellites[key]['epoch']):
                            satellites[key]['line0'] = line0
                            satellites[key]['line1'] = line1
                            satellites[key]['line2'] = line2
                            satellites[key]['epoch'] = this_epoch
                        satellites[key]['files'].append(fname)
    satlist = list(satellites.keys())

    print("Total satellites listed: {}".format(total_count))
    print("Total unique satellites:  {}".format(len(satlist)))

    # Still do, even though aren't using
    for i in range(len(satlist)):
        this_sat = satlist.pop()
        for fname in sats_by_file.keys():
            if fname in satellites[this_sat]['files']:
                satdes = '{}:{}'.format(this_sat, satellites[this_sat]['scname'])
                sats_by_file[fname].append(satdes)
                break
    with open('completeset.tle', 'w') as fp:
        for this_sat in satellites.keys():
            print(satellites[this_sat]['line0'], file=fp)
            print(satellites[this_sat]['line1'], file=fp)
            print(satellites[this_sat]['line2'], file=fp)
