#ifndef G2T_PROFILES_HPP
#define G2T_PROFILES_HPP

#define _USE_MATH_DEFINES

#include <cmath>
#include <array>
#include <string>
#include <vector>

#include "G2TOptions.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace g2t {

const int MAXZOOMLEVEL = 32;

/// TMS Global Mercator Profile (EPSG:3857), ported from gdal2tiles.py
class GlobalMercator {
public:
  explicit GlobalMercator(int tileSize = 256)
    : tile_size(tileSize),
      initialResolution(2.0 * M_PI * 6378137.0 / tileSize),
      originShift(2.0 * M_PI * 6378137.0 / 2.0) {}

  void LatLonToMeters(double lat, double lon, double &mx, double &my) const {
    mx = lon * originShift / 180.0;
    my = std::log(std::tan((90.0 + lat) * M_PI / 360.0)) / (M_PI / 180.0);
    my = my * originShift / 180.0;
  }

  void MetersToLatLon(double mx, double my, double &lat, double &lon) const {
    lon = (mx / originShift) * 180.0;
    lat = (my / originShift) * 180.0;
    lat = 180.0 / M_PI * (2.0 * std::atan(std::exp(lat * M_PI / 180.0)) - M_PI / 2.0);
  }

  double Resolution(int zoom) const {
    return initialResolution / std::pow(2.0, zoom);
  }

  void PixelsToMeters(double px, double py, int zoom, double &mx, double &my) const {
    double res = Resolution(zoom);
    mx = px * res - originShift;
    my = py * res - originShift;
  }

  void MetersToPixels(double mx, double my, int zoom, double &px, double &py) const {
    double res = Resolution(zoom);
    px = (mx + originShift) / res;
    py = (my + originShift) / res;
  }

  void PixelsToTile(double px, double py, int &tx, int &ty) const {
    tx = static_cast<int>(std::ceil(px / static_cast<double>(tile_size))) - 1;
    ty = static_cast<int>(std::ceil(py / static_cast<double>(tile_size))) - 1;
  }

  void MetersToTile(double mx, double my, int zoom, int &tx, int &ty) const {
    double px, py;
    MetersToPixels(mx, my, zoom, px, py);
    PixelsToTile(px, py, tx, ty);
  }

  void TileBounds(int tx, int ty, int zoom, double b[4]) const {
    PixelsToMeters(tx * tile_size, ty * tile_size, zoom, b[0], b[1]);
    PixelsToMeters((tx + 1) * tile_size, (ty + 1) * tile_size, zoom, b[2], b[3]);
  }

  void TileLatLonBounds(int tx, int ty, int zoom, double b[4]) const {
    double tb[4];
    TileBounds(tx, ty, zoom, tb);
    double minLat, minLon, maxLat, maxLon;
    MetersToLatLon(tb[0], tb[1], minLat, minLon);
    MetersToLatLon(tb[2], tb[3], maxLat, maxLon);
    b[0] = minLat; b[1] = minLon; b[2] = maxLat; b[3] = maxLon;
  }

  int ZoomForPixelSize(double pixelSize) const {
    for (int i = 0; i < MAXZOOMLEVEL; ++i) {
      if (pixelSize > Resolution(i))
        return std::max(0, i - 1);
    }
    return MAXZOOMLEVEL - 1;
  }

  int tile_size;

private:
  double initialResolution;
  double originShift;
};

/// TMS Global Geodetic Profile (EPSG:4326), ported from gdal2tiles.py
class GlobalGeodetic {
public:
  GlobalGeodetic(bool tmscompatible, int tileSize = 256)
    : tile_size(tileSize) {
    resFact = tmscompatible ? 180.0 / tileSize : 360.0 / tileSize;
  }

  void LonLatToPixels(double lon, double lat, int zoom, double &px, double &py) const {
    double res = resFact / std::pow(2.0, zoom);
    px = (180.0 + lon) / res;
    py = (90.0 + lat) / res;
  }

  void PixelsToTile(double px, double py, int &tx, int &ty) const {
    tx = static_cast<int>(std::ceil(px / static_cast<double>(tile_size))) - 1;
    ty = static_cast<int>(std::ceil(py / static_cast<double>(tile_size))) - 1;
  }

  void LonLatToTile(double lon, double lat, int zoom, int &tx, int &ty) const {
    double px, py;
    LonLatToPixels(lon, lat, zoom, px, py);
    PixelsToTile(px, py, tx, ty);
  }

  double Resolution(int zoom) const {
    return resFact / std::pow(2.0, zoom);
  }

  int ZoomForPixelSize(double pixelSize) const {
    for (int i = 0; i < MAXZOOMLEVEL; ++i) {
      if (pixelSize > Resolution(i))
        return std::max(0, i - 1);
    }
    return MAXZOOMLEVEL - 1;
  }

  void TileBounds(int tx, int ty, int zoom, double b[4]) const {
    double res = resFact / std::pow(2.0, zoom);
    b[0] = tx * tile_size * res - 180.0;
    b[1] = ty * tile_size * res - 90.0;
    b[2] = (tx + 1) * tile_size * res - 180.0;
    b[3] = (ty + 1) * tile_size * res - 90.0;
  }

  void TileLatLonBounds(int tx, int ty, int zoom, double b[4]) const {
    double tb[4];
    TileBounds(tx, ty, zoom, tb);
    b[0] = tb[1]; b[1] = tb[0]; b[2] = tb[3]; b[3] = tb[2];
  }

  int tile_size;

private:
  double resFact;
};

/// The tile currently being processed
struct TileDetail {
  int tx = 0;
  int ty = 0;
  int ty_tms = 0;
  int tz = 0;
  int rx = 0;
  int ry = 0;
  int rxsize = 0;
  int rysize = 0;
  int wx = 0;
  int wy = 0;
  int wxsize = 0;
  int wysize = 0;
  int querysize = 0;
};

/// Plain object holding tile job configuration for a dataset
struct TileJobInfo {
  std::string src_file;
  int nb_data_bands = 0;
  std::string output_file_path;
  std::string tile_extension;
  int tile_size = 0;
  std::string tile_driver;
  std::vector<std::array<int, 4>> tminmax;
  int tminz = 0;
  int tmaxz = 0;
  std::string in_srs_wkt;
  double out_geo_trans[6] = {0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  double ominy = 0.0;
  bool is_epsg_4326 = false;
  bool exclude_transparent = false;
  const Options *options = nullptr;
};

} // namespace g2t

#endif /* G2T_PROFILES_HPP */
