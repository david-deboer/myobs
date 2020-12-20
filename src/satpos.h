#ifndef _MODELS_H
#define _MODELS_H

/* If MAIN is defined, this header will generate definitions.  Otherwise,
   only declarations will be generated. */

#ifdef MAIN   /* generate DEFINITIONS */
	#define CLASS
	#define EQUALS(N) = N
struct satTime {
	int Year;
	int Month;
	int Day;
	int Hour;
	int Minute;
	double Second;
	double time;
};
struct observer {
	char Name[12];
	char Code[8];
    char TLEpath[80];
	double lat;
	double lng;
	double alt;
	double h;
	double range;
	double Horizon;
	int timeZone;
	int daylightSavings;
	struct satTime tstart;
	struct satTime tstop;
	double tstep;
};
#else       /* generate DECLARATIONS */
	#define CLASS extern
	#define EQUALS(N)
#endif

#ifndef PI
  #define PI (3.14159654)
#endif

#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define SQ(X)  ( (X)*(X) )
#undef VERBOSE

int readObserver(struct observer *obs);
int readLocation(struct observer *obs);
int invxyz(double x, double y, double z, double *lng, double *lat, double *h);
int teme2pef(struct observer obs, double jd, double *ro, struct observer *subsat);
void addTime(struct satTime *base_t, struct satTime diff_t);
void resetTime(struct satTime *t);
int getline(FILE *fp, char *s, int lim);

#endif
