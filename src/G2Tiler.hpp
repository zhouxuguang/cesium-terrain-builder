#ifndef G2T_TILER_HPP
#define G2T_TILER_HPP

#include <array>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "config.hpp"
#include "gdal_priv.h"
#include "ogr_spatialref.h"

#include "G2TOptions.hpp"
#include "G2TProfiles.hpp"

namespace g2t {

class G2TException;
class GDAL2Tiles;

} // namespace g2t

/// Exception thrown for all gdal2tiles errors
class CTB_DLL g2t::G2TException : public std::runtime_error {
public:
  explicit G2TException(const std::string &msg) : std::runtime_error(msg) {}
};

namespace g2t {

/// Scale down the query dataset to the tile dataset
CTB_DLL void scaleQueryToTile(GDALDataset *dsquery, GDALDataset *dstile,
                              const Options &options,
                              const std::string &tilefilename);

/// Create a base tile from the input raster
CTB_DLL void createBaseTile(const TileJobInfo &tile_job_info,
                            const TileDetail &tile_detail);

/// Create an overview tile from no more than 4 underlying base tiles
CTB_DLL void createOverviewTile(int base_tz,
                                const std::vector<std::pair<int, int>> &base_tiles,
                                const std::string &output_folder,
                                const TileJobInfo &tile_job_info,
                                const Options &options);

/// Group base tiles that belong to the same overview tile
CTB_DLL std::vector<std::vector<std::pair<int, int>>>
groupOverviewBaseTiles(int base_tz, const std::string &output_folder,
                       const TileJobInfo &tile_job_info);

/// Count the number of overview tiles
CTB_DLL int countOverviewTiles(const TileJobInfo &tile_job_info);

/// Post process the command line options (defaults, zoom parsing, validation)
CTB_DLL void postProcessOptions(Options &options, const std::string &input_file,
                                const std::string &output_folder);

} // namespace g2t

/// Ported GDAL2Tiles class
class CTB_DLL g2t::GDAL2Tiles {
public:
  GDAL2Tiles(const std::string &input_file, const std::string &output_folder,
             const Options &options);
  ~GDAL2Tiles();

  void openInput();
  void generateMetadata();
  void run();

  static int getYTile(int ty, int tz, const Options &options);

private:
  void generateBaseTiles(TileJobInfo &tile_job_info,
                         std::vector<TileDetail> &tile_details);
  void geoQuery(GDALDataset *ds, double ulx, double uly, double lrx, double lry,
                int querysize, int &rx, int &ry, int &rxsize, int &rysize,
                int &wx, int &wy, int &wxsize, int &wysize) const;

  std::string m_input_file;
  std::string m_output_folder;
  Options m_options;

  int m_tile_size;
  std::string m_tiledriver;
  std::string m_tileext;
  std::string m_tmp_vrt_filename;

  bool m_scaledquery;
  int m_querysize;
  bool m_overviewquery;

  int m_tminz;
  int m_tmaxz;

  GDALDataset *m_input_dataset;
  GDALDataset *m_warped_input_dataset;
  OGRSpatialReference *m_in_srs;
  OGRSpatialReference *m_out_srs;
  std::string m_in_srs_wkt;

  int m_dataBandsCount;
  double m_out_gt[6];
  double m_ominx;
  double m_omaxx;
  double m_omaxy;
  double m_ominy;

  GlobalMercator *m_mercator;
  GlobalGeodetic *m_geodetic;
  std::vector<std::array<int, 4>> m_tminmax;
  bool m_isepsg4326;

  int m_nativezoom;
  std::vector<int> m_tsize;
};

#endif /* G2T_TILER_HPP */
