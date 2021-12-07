#! /usr/bin/env python
import requests
from os import path


base_path = 'http://www.celestrak.com/NORAD/elements/'

if __name__ == '__main__':
    master_file = requests.get(path.join(base_path, 'master.asp'))
    master = master_file.text.splitlines()
    tlefiles = {}

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
                tlefiles[word] = f"{description} [{num}]"
                break

    with open('master.dat', 'w') as master:
        for lll in tlefiles:
            useThis = True
            for ig in ignore:
                if ig in tlefiles[lll].lower():
                    useThis = False
            if useThis:
                a = lll.split('.')
                outfile = a[0]+'.tle'
                print('Reading %s:  %s' % (lll, tlefiles[lll]))
                sat = requests.get(path.join(base_path, lll)).text.splitlines()  # noqa
                with open(outfile, 'w') as fp:
                    for line in sat:
                        print(line, file=fp)
                print("{}:  {}".format(outfile, tlefiles[lll]), file=master)
