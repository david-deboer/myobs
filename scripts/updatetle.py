#! /usr/local/anaconda3/bin/python
import requests
from os import path
import argparse

ap = argparse.ArgumentParser()
ap.add_argument('--base-url', dest='base_url', help="Base url for tles",
                default='http://www.celestrak.com/NORAD/elements/')
ap.add_argument('--base-path', dest='base_path', help="Base path for tles",
                default='.')
args = ap.parse_args()


def updatetle(base_path, base_url):
    master_file = requests.get(path.join(base_url, 'master.asp'))
    master = master_file.text.splitlines()
    tlefiles = {}
    base_path = path.expanduser(base_path)

    ignore = ['debris', 'cesium']

    print('Reading Celestrak master file')
    for line in master:
        data = line.split('"')
        for word in data:
            if '.txt' in word and not word.startswith('/'):
                description = line.split('>')[2].split('(')[0].strip()
                try:
                    num = line.split('[')[1].split(']')[0]
                except IndexError:
                    num = '-'
                tlefiles[word] = "{} [{}]".format(description, num)
                break

    with open(path.join(base_path, 'master.dat'), 'w') as master:
        for lll in tlefiles:
            useThis = True
            for ig in ignore:
                if ig in tlefiles[lll].lower():
                    useThis = False
            if useThis:
                a = lll.split('.')
                outfile = a[0]+'.tle'
                print('Reading %s:  %s' % (lll, tlefiles[lll]))
                sat = requests.get(path.join(base_url, lll)).text.splitlines()  # noqa
                with open(path.join(base_path, outfile), 'w') as fp:
                    for line in sat:
                        print(line, file=fp)
                print("{}:  {}".format(outfile, tlefiles[lll]), file=master)


if __name__ == '__main__':
    updatetle(base_path=args.base_path, base_url=args.base_url)
