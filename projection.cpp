#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "projection.hpp"

// input projections of source files (GeoJson)
// The default is 4326
struct projection projections[] = {
	{"EPSG:4326", lonlat2tile, tile2lonlat, "urn:ogc:def:crs:OGC:1.3:CRS84"},
	{"EPSG:3857", epsg3857totile, tiletoepsg3857, "urn:ogc:def:crs:EPSG::3857"},
	{NULL, NULL},
};

// 4326 is the lingua franca of tile projections
// everything is from-to 4326
// Note: 3857 maps perfectly back to 4326 (lossless)
//       3395 is accurate to 8 decimal places (but an approx when converting back to 4326)
// The default is 3857
struct tileprojection tileprojections[] = {
        {"EPSG:3857", epsg4326_to_3857, epsg3857_to_4326, "urn:ogc:def:crs:EPSG::3857"},
        {"EPSG:3395", epsg4326_to_3395, epsg3857_to_4326, "urn:ogc:def:crs:EPSG::3395"}, // using simpler 3857 for now
        {NULL, NULL},
};

struct projection *projection = &projections[0]; // projection of source (GeoJson)
struct tileprojection *tileprojection = &tileprojections[0]; // projections of tiles (vector tiles)

// http://wiki.openstreetmap.org/wiki/Slippy_map_tilenames
void lonlat2tile(double lon, double lat, int zoom, long long *x, long long *y) {
    tileprojection->project(lon, lat, zoom, x, y);
}

// http://wiki.openstreetmap.org/wiki/Slippy_map_tilenames
void tile2lonlat(long long x, long long y, int zoom, double *lon, double *lat) {
    tileprojection->unproject(x, y, zoom, lon, lat);
}

void epsg3857totile(double ix, double iy, int zoom, long long *x, long long *y) {
    double lon, lat;
    epsg3857_to_4326(ix, iy, zoom, &lon, &lat); // first get it as 4326
    tileprojection->project(lon, lat, zoom, x, y);
}

void tiletoepsg3857(long long ix, long long iy, int zoom, double *ox, double *oy) {
    double lon, lat;
    tileprojection->unproject(ix, iy, zoom, &lon, &lat);  // first get as 4326
    long long x, y;
    epsg4326_to_3857(lon, lat, 32, &x, &y); // project to 3857
    *ox = (double)x;
    *oy = (double)y;
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

// http://wiki.openstreetmap.org/wiki/Slippy_map_tilenames
void epsg4326_to_3857(double lon, double lat, int zoom, long long *x, long long *y) {
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

void epsg4326_to_3395(double lon, double lat, int zoom, long long *x, long long *y) {
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

// http://wiki.openstreetmap.org/wiki/Slippy_map_tilenames
// This is a lossless transformation (perfect mapping)
void epsg3857_to_4326(long long x, long long y, int zoom, double *lon, double *lat) {
    unsigned long long n = 1LL << zoom;
    *lon = 360.0 * x / n - 180.0;
    *lat = atan(sinh(M_PI * (1 - 2.0 * y / n))) * 180.0 / M_PI;
    fprintf(stdout, "\n\n\n***epsg3857_to_4326([%lld, %lld, %d]) yields [%.05f, %.05f] ***\n\n\n", x, y, zoom, *lon, *lat);
}

// This is an iterative transformation, valid to 8 decimal places
void epsg3395_to_4326(long long x, long long y, int zoom, double *lon, double *lat) {
    if (zoom < 0) {
        zoom = 32;
    }
    unsigned long long n = 1LL << zoom;
    *lon = 360.0 * x / n - 180.0;

    double a = 6378137.0;
    double e = 0.081819190842622;
    double circumference = 2.0 * M_PI * a;
    // Convert back to TMS from Google tile (flip y)
    long long ny = (n - 1) - y;

    double res = circumference / 256.0;
    res = res / n;

    // First unmap back to pixels
    long long py = (ny + 1) * 256;

    // Now back to meters
    double originShift = circumference / 2.0;
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
    fprintf(stdout, "\n\n\n***epsg3395_to_4326([%lld, %lld, %d]) yields [%.05f, %.05f] ***\n\n\n", x, y, zoom, *lon, *lat);
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

void set_output_projection_or_exit(const char *optarg) {
    struct tileprojection *p;
    for (p = tileprojections; p->name != NULL; p++) {
        if (strcmp(p->name, optarg) == 0) {
            tileprojection = p;
            break;
        }
        if (strcmp(p->alias, optarg) == 0) {
            tileprojection = p;
            break;
        }
    }
    if (p->name == NULL) {
        fprintf(stderr, "Unknown output (tile) projection (-O): %s\n", optarg);
        exit(EXIT_FAILURE);
    }
}
