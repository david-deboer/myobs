from datetime import datetime
from astropy.time import Time
from astropy.time import TimeDelta
from numpy import ndarray


def get_astropytime(dates, times=None):
    """
    Get an astropy time object based on provided date/time formats.

    Take in various incarnations of dates/times and return an astropy.Time object or None.
    No time zone is allowed.

    Returns:  either astropy Time or None

    Parameters
    ----------
    dates :  a date in various formats.  If a list, recursively goes through it.
                    astropy Time
                    datetime
                    int, long, float:  interpreted as gps_second or julian date
                                       depending on appropriate range
                    string:  '<' - PAST_DATE
                             '>' - future_date()
                             'now' or 'current'
                             'YYYY/M/D' or 'YYYY-M-D'
             the following just return None, ignoring times
                    string:  'none' return None
                    None/False:  return None
    times : offset from dates in hours
                float, int:  hours in decimal time
                string:  HH[:MM[:SS]] or hours in decimal time

    Returns
    -------
    astropy.Time or None

    """
    # Process dates
    if dates is None or dates is False:
        return None
    if isinstance(dates, list) or isinstance(times, (list, ndarray)):
        if isinstance(times, ndarray):
            times = list(times)
        if not isinstance(times, list):
            times = [times] * len(dates)
        if not isinstance(dates, list):
            dates = [dates] * len(times)
        if len(dates) != len(times):
            raise ValueError("dates/times list lengths must match.")
        return_Time = []
        if len(dates) > 1000:
            print("Converting {} time entries - could take a moment.".format(len(dates)))
        for _date, _time in zip(dates, times):
            return_Time.append(get_astropytime(_date, _time))
        return Time(return_Time)
    if isinstance(dates, str):
        if dates.lower() == 'none':
            return None
        if dates == '<':
            return Time('2000-01-01', scale='utc')
        if dates == '>':
            return Time.now() + TimeDelta(1000, format='jd')
        if dates.lower() == 'now' or dates.lower() == 'current':
            return Time.now()
    if isinstance(dates, Time):
        return_Time = dates
    elif isinstance(dates, datetime):
        return_Time = Time(dates, format='datetime')
    else:
        try:
            dates = float(dates)
            if dates > 1000000000.0:
                return_Time = Time(dates, format='gps')
            elif dates > 2400000.0 and dates < 2500000.0:
                return_Time = Time(dates, format='jd')
            else:
                raise ValueError(f'Invalid format:  date as a number should be gps time '
                                 f'or julian date, not {dates}.')
        except ValueError:
            dates = dates.replace('/', '-')
            try:
                return_Time = Time(dates, scale='utc')
            except ValueError:
                raise ValueError(
                    f'Invalid format:  YYYY[/-]M[/-]D [HH:MM:SS], not {dates}')
    # add on times
    if times is None or abs(times) < 1E-6:
        return return_Time
    try:
        times = float(times)
        return return_Time + TimeDelta(times * 3600.0, format='sec')
    except ValueError:
        pass
    sign_of_times = 1.0
    if times[0] == '-':
        sign_of_times = -1.0
        times = times[1:]
    add_time = 0.0
    for i, d in enumerate(times.split(':')):
        add_time += (float(d)) * 3600.0 / (60.0**i)
    add_time *= sign_of_times
    return return_Time + TimeDelta(add_time, format='sec')
