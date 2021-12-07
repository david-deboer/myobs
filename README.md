satpos  implements the sgp4 orbital tracking code (get ref to where the actual c code came from) and some other modules to look at the data.

The TLE files are maintained /tle and may be updated using the script 'updatetle.py' while in that directory.  Note that when updated should git commit -am 'TLE update on YY-MM-DD' so that old ones may be found via git log.

--satpos--
Given TLE (from a file, location set in satpos.cfg) it writes a file:
    'sp{filename}{entry_number_in_file:04d}.out' which contains a track over the period specified in satpos.cfg
e.g. satpos active 123


--sattrack--
from satpos import sattrack
reads in the .out files from above and computes various parameters.

--satcensus--
from satpos import satcensus
random modules using sattrack (or not) to look at satellite data.
