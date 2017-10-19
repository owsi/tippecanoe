#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "projection.hpp"

struct projection projections[] = {
	{"EPSG:4326", lonlat2tile, tile2lonlat, "urn:ogc:def:crs:OGC:1.3:CRS84"},
	{"EPSG:3857", epsg3857totile, tiletoepsg3857, "urn:ogc:def:crs:EPSG::3857"},
	{"EPSG:3395", lonlat2epsg3395tile, epsg3395tiletolonlat, "urn:ogc:def:crs:EPSG::3395"},
	{NULL, NULL},
};

struct projection *projection = &projections[0];

// http://wiki.openstreetmap.org/wiki/Slippy_map_tilenames
void lonlat2tile(double lon, double lat, int zoom, long long *x, long long *y) {
	// Must limit latitude somewhere to prevent overflow.
	// 89.9 degrees latitude is 0.621 worlds beyond the edge of the flat earth,
	// hopefully far enough out that there are few expectations about the shape.
	if (lat < -89.9) {
		lat = -89.9;
	}
	if (lat > 89.9) {
		lat = 89.9;
	}

	if (lon < -360) {
		lon = -360;
	}
	if (lon > 360) {
		lon = 360;
	}

	double lat_rad = lat * M_PI / 180;
	unsigned long long n = 1LL << zoom;

	long long llx = n * ((lon + 180) / 360);
	long long lly = n * (1 - (log(tan(lat_rad) + 1 / cos(lat_rad)) / M_PI)) / 2;

	*x = llx;
	*y = lly;
}

// http://wiki.openstreetmap.org/wiki/Slippy_map_tilenames
void tile2lonlat(long long x, long long y, int zoom, double *lon, double *lat) {
	unsigned long long n = 1LL << zoom;
	*lon = 360.0 * x / n - 180.0;
	*lat = atan(sinh(M_PI * (1 - 2.0 * y / n))) * 180.0 / M_PI;
}

void lonlat2epsg3395tile(double lon, double lat, int zoom, long long *x, long long *y) {
	// Must limit latitude somewhere to prevent overflow.
	// 89.9 degrees latitude is 0.621 worlds beyond the edge of the flat earth,
	// hopefully far enough out that there are few expectations about the shape.
	if (lat < -89.9) {
		lat = -89.9;
	}
	if (lat > 89.9) {
		lat = 89.9;
	}

	if (lon < -360) {
		lon = -360;
	}
	if (lon > 360) {
		lon = 360;
	}
	double a = 6378137.0;
	double e = 0.081819190842622;
    unsigned long long n = 1LL << zoom;

	double lat_rad = lat*M_PI/180.0;
	double lon_rad = lon*M_PI/180.0;
	double dx = a * lon_rad;
	double d1 = tan(M_PI/4 + lat_rad/2);
	double sinLat = e * sin(lat_rad);
	double dt = 1.0 - sinLat;
	double dt1 = 1.0 + sinLat;
	double dy = a * log(d1 * pow(dt/dt1, e/2));

    double circumference = 2.0 * M_PI * 6378137.0;
	double res = circumference / 256.0;
	res = res / exp2((double)zoom);
    double originShift = circumference / 2.0;
	double px = (dx + originShift) / res;
	double py = (dy + originShift) / res;
    long long llx = ceil(px/256.0) - 1;
    long long lly = ceil(py/256.0) - 1;

    // Store as Google tile (flip y)
	*x = llx;
	*y = (n - 1) - lly;
}

void epsg3395tiletolonlat(long long x, long long y, int zoom, double *lon, double *lat) {
	unsigned long long n = 1LL << zoom;
	double a = 6378137.0;
	double e = 0.081819190842622;
    double circumference = 2.0 * M_PI * 6378137.0;
    // Convert back to TMS from Google tile (flip y)
    long long ny = (n - 1) - y;

	double res = circumference / 256.0;
	res = res / n;

	// First unmap back to pixels
	long long px = (x + 1) * 256;
	long long py = (ny + 1) * 256;

	// Now back to meters
    double originShift = circumference / 2.0;
 	double dx = px/res - originShift;
 	double dy = py/res - originShift;

 	// Now iterate to arrive at lat (radians)
 	double t = exp(-dy/a);
 	double tlat = M_PI/2 - atan(t);
 	double EPSILON = 0.00000001;
 	bool changing = true;
 	double term = e * sin(tlat);
 	double prev = M_PI/2 - (2 * atan(t * pow((1 - term)/(1 + term), e/2)));
 	while (changing) {
 	    term = e * sin(tlat);
 	    tlat = M_PI/2 - (2 * atan(t * pow((1 - term)/(1 + term), e/2)));

        double diff = tlat - prev;
        if (diff < 0.0) {
            diff = diff * -1.0;
        }
 	    if (diff <= EPSILON ) {
 	        changing = false;
 	     } else {
 	        prev = tlat;
 	     }

 	}

 	// Back to decimal degrees
    *lat = tlat * (180.0/M_PI);
 	*lon = (dx/a) * (180.0/M_PI);
}


void epsg3857totile(double ix, double iy, int zoom, long long *x, long long *y) {
	*x = ix * (1LL << 31) / 6378137.0 / M_PI + (1LL << 31);
	*y = ((1LL << 32) - 1) - (iy * (1LL << 31) / 6378137.0 / M_PI + (1LL << 31));

	if (zoom != 0) {
		*x >>= (32 - zoom);
		*y >>= (32 - zoom);
	}
}

void tiletoepsg3857(long long ix, long long iy, int zoom, double *ox, double *oy) {
	if (zoom != 0) {
		ix <<= (32 - zoom);
		iy <<= (32 - zoom);
	}

	*ox = (ix - (1LL << 31)) * M_PI * 6378137.0 / (1LL << 31);
	*oy = ((1LL << 32) - 1 - iy - (1LL << 31)) * M_PI * 6378137.0 / (1LL << 31);
}

unsigned long long encode(unsigned int wx, unsigned int wy) {
	unsigned long long out = 0;

	int i;
	for (i = 0; i < 32; i++) {
		unsigned long long v = ((wx >> (32 - (i + 1))) & 1) << 1;
		v |= (wy >> (32 - (i + 1))) & 1;
		v = v << (64 - 2 * (i + 1));

		out |= v;
	}

	return out;
}

static unsigned char decodex[256];
static unsigned char decodey[256];

void decode(unsigned long long index, unsigned *wx, unsigned *wy) {
	static int initialized = 0;
	if (!initialized) {
		for (size_t ix = 0; ix < 256; ix++) {
			size_t xx = 0, yy = 0;

			for (size_t i = 0; i < 32; i++) {
				xx |= ((ix >> (64 - 2 * (i + 1) + 1)) & 1) << (32 - (i + 1));
				yy |= ((ix >> (64 - 2 * (i + 1) + 0)) & 1) << (32 - (i + 1));
			}

			decodex[ix] = xx;
			decodey[ix] = yy;
		}

		initialized = 1;
	}

	*wx = *wy = 0;

	for (size_t i = 0; i < 8; i++) {
		*wx |= ((unsigned) decodex[(index >> (8 * i)) & 0xFF]) << (4 * i);
		*wy |= ((unsigned) decodey[(index >> (8 * i)) & 0xFF]) << (4 * i);
	}
}

void set_projection_or_exit(const char *optarg) {
	struct projection *p;
	for (p = projections; p->name != NULL; p++) {
		if (strcmp(p->name, optarg) == 0) {
			projection = p;
			break;
		}
		if (strcmp(p->alias, optarg) == 0) {
			projection = p;
			break;
		}
	}
	if (p->name == NULL) {
		fprintf(stderr, "Unknown projection (-s): %s\n", optarg);
		exit(EXIT_FAILURE);
	}
}
