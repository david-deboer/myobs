#! /usr/bin/env python
import requests

master_file = requests.get("http://www.celestrak.com/NORAD/elements/master.asp")
master = master_file.text.splitlines()
tlefiles = {}

ignore = ['debris']

print('Reading Celestrak master file')
for line in master:
    data = line.split('"')
    for word in data:
        if '.txt' in word:
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
            sat = requests.get(f"http://www.celestrak.com/NORAD/elements/{lll}").text.splitlines()
            with open(outfile, 'w') as fp:
                for line in sat:
                    print(line, file=fp)
            print("{}:  {}".format(outfile, tlefiles[lll]), file=master)
