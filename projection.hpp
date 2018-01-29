#ifndef PROJECTION_HPP
#define PROJECTION_HPP

void lonlat2tile(double lon, double lat, int zoom, long long *x, long long *y);
void epsg3857totile(double ix, double iy, int zoom, long long *x, long long *y);
void tile2lonlat(long long x, long long y, int zoom, double *lon, double *lat);
void tiletoepsg3857(long long x, long long y, int zoom, double *ox, double *oy);
unsigned long long encode(unsigned int wx, unsigned int wy);
void decode(unsigned long long index, unsigned *wx, unsigned *wy);
void epsg4326_to_3857(double lon, double lat, int zoom, long long *x, long long *y);
void epsg4326_to_3395(double lon, double lat, int zoom, long long *x, long long *y);
void epsg3857_to_4326(long long x, long long y, int zoom, double *lon, double *lat);
void epsg3395_to_4326(long long x, long long y, int zoom, double *lon, double *lat);

void set_projection_or_exit(const char *optarg);
void set_output_projection_or_exit(const char *optarg);

struct projection {
	const char *name;
	void (*project)(double ix, double iy, int zoom, long long *ox, long long *oy);
	void (*unproject)(long long ix, long long iy, int zoom, double *ox, double *oy);
	const char *alias;
};

struct tileprojection {
	const char *name;
	void (*project)(double ix, double iy, int zoom, long long *ox, long long *oy);
	void (*unproject)(long long ix, long long iy, int zoom, double *ox, double *oy);
	const char *alias;
};

extern struct projection *projection;
extern struct projection projections[];
extern struct tileprojection *outprojection;
extern struct tileprojection tileprojections[];

#endif
