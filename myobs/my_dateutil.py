"""This is copied from my_utils."""

from datetime import datetime


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
    from astropy.time import Time
    from astropy.time import TimeDelta
    from numpy import ndarray

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


def sort_date_keys(keys, date_format=None):
    data = {}
    if isinstance(date_format, str) and '%' in date_format:
        date_format = [date_format]
    for dd in list(keys):
        if date_format is None:
            data[string_to_date(dd)] = dd
        else:
            data[string_to_date(dd, date_format)] = dd
    return [data[x] for x in sorted(list(data.keys()))]


def get_now(date, fmt='%Y%m%d'):
    return date_to_string(datetime.now(), fmt=fmt)


def date_to_string(date, fmt='%m/%d/%y'):
    return datetime.strftime(date, fmt)


def string_to_date(date, strings_to_try=['%m/%d/%y', '%m/%d/%Y', '%Y%m%d'], return_format=False):
    if isinstance(strings_to_try, str):
        strings_to_try = [strings_to_try]
    for fmt in strings_to_try:
        try:
            dt = datetime.strptime(date, fmt)
            if return_format:
                return fmt
            return dt
        except ValueError:
            continue


class datekey:
    def __init__(self, strings_to_try=['%m/%d/%y', '%m/%d/%Y', '%Y%m%d']):
        if isinstance(strings_to_try, str) and '%' in strings_to_try:
            self.default_date_format = strings_to_try
            self.strings_to_try = [strings_to_try]
        else:
            self.strings_to_try = strings_to_try
            self.default_date_format = None

    def set_format(self, datefmt):
        if '%' in datefmt:
            self.default_date_format = datefmt
            return
        self.default_date_format = string_to_date(datefmt, return_format=True,
                                                  strings_to_try=self.strings_to_try)

    def string_to_date(self, date, strings_to_try=None):
        if strings_to_try is None:
            strings_to_try = self.strings_to_try
        elif strings_to_try is str:
            strings_to_try = [strings_to_try]
        if self.default_date_format is None:
            self.set_format(date)
        return string_to_date(date, strings_to_try=strings_to_try)

    def date_to_string(self, date, datefmt=None):
        if datefmt is None:
            datefmt = self.default_date_format
        return date_to_string(date, datefmt)
