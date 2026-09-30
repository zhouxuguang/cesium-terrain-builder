#ifndef G2T_OPTIONS_HPP
#define G2T_OPTIONS_HPP

#include <string>

namespace g2t {

/// Command line options (mirrors gdal2tiles.py, without mpi/kml/webviewer/mapml)
struct Options {
  std::string profile = "mercator";
  std::string resampling = "average";
  std::string s_srs = "";
  std::string zoom_str = "";
  int tminz = -1;               ///< -1 means auto
  int tmaxz = -1;               ///< -1 means auto
  bool resume = false;
  std::string srcnodata = "";
  bool tmscompatible = false;
  bool xyz = false;
  bool verbose = false;
  bool exclude_transparent = false;
  bool quiet = false;
  int nb_processes = 1;
  int tilesize = 256;
  std::string tiledriver = "PNG";
  std::string excluded_values = "";
  double excluded_values_pct_threshold = 50.0;
  double nodata_values_pct_threshold = 100.0;
  int webp_quality = 75;
  bool webp_lossless = false;
  int jpeg_quality = 75;
  std::string title = "";
};

} // namespace g2t

#endif /* G2T_OPTIONS_HPP */
