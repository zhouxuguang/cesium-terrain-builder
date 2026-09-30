#define _USE_MATH_DEFINES

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <functional>
#include <iostream>
#include <map>
#include <mutex>
#include <sstream>
#include <thread>

#include "cpl_conv.h"
#include "cpl_error.h"
#include "cpl_minixml.h"
#include "cpl_string.h"
#include "cpl_vsi.h"
#include "gdal.h"
#include "gdal_alg.h"
#include "gdal_priv.h"
#include "gdal_utils.h"
#include "gdalwarper.h"
#include "ogr_spatialref.h"

#include "G2Tiler.hpp"

namespace g2t {

// ---------------------------------------------------------------------------
// Small helpers
// ---------------------------------------------------------------------------

static std::string
joinPath(const std::string &a, const std::string &b) {
  if (a.empty())
    return b;
  char last = a[a.size() - 1];
  if (last == '/' || last == '\\')
    return a + b;
  return a + "/" + b;
}

static void
makedirs(const std::string &path) {
  if (path.empty())
    return;

  VSIStatBufL stat;
  if (VSIStatL(path.c_str(), &stat) == 0)
    return;

  if (VSIMkdirRecursive(path.c_str(), 0755) != 0) {
    if (VSIStatL(path.c_str(), &stat) != 0)
      throw G2TException("Cannot create the output directory: " + path);
  }
}

static bool
fileExists(const std::string &path) {
  VSIStatBufL stat;
  if (VSIStatL(path.c_str(), &stat) != 0)
    return false;
  return VSI_ISREG(stat.st_mode) != 0;
}

static std::string
dirname(const std::string &path) {
  std::string::size_type pos = path.find_last_of("/\\");
  if (pos == std::string::npos)
    return ".";
  return path.substr(0, pos);
}

static std::string
proj4String(OGRSpatialReference *srs) {
  if (srs == nullptr)
    return "";
  char *pszProj4 = nullptr;
  srs->exportToProj4(&pszProj4);
  std::string result = pszProj4 ? pszProj4 : "";
  CPLFree(pszProj4);
  return result;
}

static std::string
wktString(OGRSpatialReference *srs) {
  if (srs == nullptr)
    return "";
  char *pszWkt = nullptr;
  srs->exportToWkt(&pszWkt);
  std::string result = pszWkt ? pszWkt : "";
  CPLFree(pszWkt);
  return result;
}

static bool
hasGeoreference(GDALDataset *ds) {
  double gt[6];
  ds->GetGeoTransform(gt);
  bool identity = gt[0] == 0.0 && gt[1] == 1.0 && gt[2] == 0.0 &&
                  gt[3] == 0.0 && gt[4] == 0.0 && gt[5] == 1.0;
  return !identity || ds->GetGCPCount() != 0;
}

static std::vector<std::string>
split(const std::string &s, char sep) {
  std::vector<std::string> parts;
  std::string::size_type start = 0;
  for (;;) {
    std::string::size_type pos = s.find(sep, start);
    if (pos == std::string::npos) {
      parts.push_back(s.substr(start));
      break;
    }
    parts.push_back(s.substr(start, pos - start));
    start = pos + 1;
  }
  return parts;
}

static int
nbDataBands(GDALDataset *ds) {
  GDALRasterBand *mask = ds->GetRasterBand(1)->GetMaskBand();
  int flags = mask->GetMaskFlags();
  if ((flags & GMF_ALPHA) || ds->GetRasterCount() == 4 || ds->GetRasterCount() == 2)
    return ds->GetRasterCount() - 1;
  return ds->GetRasterCount();
}

static std::vector<double>
setupNoDataValues(GDALDataset *ds, const Options &options) {
  std::vector<double> in_nodata;

  if (!options.srcnodata.empty()) {
    std::vector<std::string> tokens = split(options.srcnodata, ',');
    std::vector<double> nds;
    for (size_t i = 0; i < tokens.size(); ++i)
      nds.push_back(CPLAtof(tokens[i].c_str()));

    if (static_cast<int>(nds.size()) < ds->GetRasterCount()) {
      for (int i = 0; i < ds->GetRasterCount(); ++i)
        in_nodata.push_back(nds[i % nds.size()]);
    } else {
      in_nodata = nds;
    }
  } else {
    for (int i = 1; i <= ds->GetRasterCount(); ++i) {
      GDALRasterBand *band = ds->GetRasterBand(i);
      int hasNoData = 0;
      double value = band->GetNoDataValue(&hasNoData);
      if (hasNoData) {
        if (band->GetRasterDataType() == GDT_Byte &&
            (value != static_cast<double>(static_cast<int>(value)) || value < 0 || value > 255)) {
          in_nodata.clear();
          break;
        }
        in_nodata.push_back(value);
      }
    }
  }

  return in_nodata;
}

static OGRSpatialReference *
setupInputSrs(GDALDataset *ds, const Options &options, std::string &in_srs_wkt) {
  OGRSpatialReference *in_srs = nullptr;

  if (!options.s_srs.empty()) {
    in_srs = new OGRSpatialReference();
    if (in_srs->SetFromUserInput(options.s_srs.c_str()) != OGRERR_NONE) {
      delete in_srs;
      throw G2TException("Invalid value for --s_srs option");
    }
    in_srs_wkt = wktString(in_srs);
  } else {
    const char *projection = ds->GetProjectionRef();
    if ((projection == nullptr || projection[0] == '\0') && ds->GetGCPCount() != 0)
      projection = ds->GetGCPProjection();
    if (projection != nullptr && projection[0] != '\0') {
      in_srs_wkt = projection;
      in_srs = new OGRSpatialReference();
      if (in_srs->importFromWkt(projection) != OGRERR_NONE) {
        delete in_srs;
        in_srs = nullptr;
      }
    }
  }

  if (in_srs != nullptr)
    in_srs->SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

  return in_srs;
}

static OGRSpatialReference *
setupOutputSrs(OGRSpatialReference *in_srs, const Options &options) {
  OGRSpatialReference *out_srs = new OGRSpatialReference();

  if (options.profile == "mercator") {
    out_srs->importFromEPSG(3857);
  } else if (options.profile == "geodetic") {
    out_srs->importFromEPSG(4326);
  } else {
    if (in_srs != nullptr)
      out_srs->operator=(*in_srs);
  }

  out_srs->SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
  return out_srs;
}

// ---------------------------------------------------------------------------
// VRT XML manipulation
// ---------------------------------------------------------------------------

static void
setXmlAttr(CPLXMLNode *node, const char *name, const std::string &value) {
  CPLXMLNode *attr = CPLCreateXMLNode(node, CXT_Attribute, name);
  CPLCreateXMLNode(attr, CXT_Text, value.c_str());
}

static CPLXMLNode *
findXmlNode(CPLXMLNode *root, const char *name) {
  for (CPLXMLNode *node = root; node != nullptr; node = node->psNext) {
    if (node->eType == CXT_Element && strcmp(node->pszValue, name) == 0)
      return node;
    if (node->psChild != nullptr) {
      CPLXMLNode *found = findXmlNode(node->psChild, name);
      if (found != nullptr)
        return found;
    }
  }
  return nullptr;
}

static std::string
getVrtXml(GDALDataset *ds) {
  char **metadata = ds->GetMetadata("xml:VRT");
  if (metadata != nullptr && metadata[0] != nullptr)
    return metadata[0];
  return "";
}

static std::string
serializeXml(CPLXMLNode *root) {
  char *xml = CPLSerializeXMLTree(root);
  std::string result = xml ? xml : "";
  CPLFree(xml);
  CPLDestroyXMLNode(root);
  return result;
}

static std::string
addAlphaBandToVrt(const std::string &vrt_string) {
  CPLXMLNode *root = CPLParseXMLString(vrt_string.c_str());
  if (root == nullptr)
    return vrt_string;

  int nb_bands = 0;
  CPLXMLNode *last_band = nullptr;
  for (CPLXMLNode *node = root->psChild; node != nullptr; node = node->psNext) {
    if (node->eType == CXT_Element && strcmp(node->pszValue, "VRTRasterBand") == 0) {
      nb_bands++;
      CPLXMLNode *color = CPLGetXMLNode(node, "ColorInterp");
      if (color != nullptr && color->psChild != nullptr &&
          color->psChild->eType == CXT_Text &&
          strcmp(color->psChild->pszValue, "Alpha") == 0) {
        CPLDestroyXMLNode(root);
        throw G2TException("Alpha band already present");
      }
      last_band = node;
    } else if (nb_bands) {
      break;
    }
  }

  CPLXMLNode *band = CPLCreateXMLNode(nullptr, CXT_Element, "VRTRasterBand");
  setXmlAttr(band, "dataType", "Byte");
  setXmlAttr(band, "band", std::to_string(nb_bands + 1));
  setXmlAttr(band, "subClass", "VRTWarpedRasterBand");
  CPLXMLNode *color = CPLCreateXMLNode(band, CXT_Element, "ColorInterp");
  CPLCreateXMLNode(color, CXT_Text, "Alpha");

  if (last_band != nullptr) {
    band->psNext = last_band->psNext;
    last_band->psNext = band;
  } else {
    band->psNext = root->psChild;
    root->psChild = band;
  }

  CPLXMLNode *warp_options = findXmlNode(root->psChild, "GDALWarpOptions");
  if (warp_options != nullptr) {
    CPLXMLNode *dst_alpha = CPLCreateXMLNode(warp_options, CXT_Element, "DstAlphaBand");
    CPLCreateXMLNode(dst_alpha, CXT_Text, std::to_string(nb_bands + 1).c_str());
    CPLXMLNode *option = CPLCreateXMLNode(warp_options, CXT_Element, "Option");
    setXmlAttr(option, "name", "INIT_DEST");
    CPLCreateXMLNode(option, CXT_Text, "0");
  }

  return serializeXml(root);
}

static std::string
addWarpOptionsToVrt(const std::string &vrt_string,
                    const std::vector<std::pair<std::string, std::string>> &warp_options) {
  if (warp_options.empty())
    return vrt_string;

  CPLXMLNode *root = CPLParseXMLString(vrt_string.c_str());
  if (root == nullptr)
    return vrt_string;

  CPLXMLNode *options = findXmlNode(root->psChild, "GDALWarpOptions");
  if (options == nullptr) {
    CPLDestroyXMLNode(root);
    return vrt_string;
  }

  for (size_t i = 0; i < warp_options.size(); ++i) {
    CPLXMLNode *option = CPLCreateXMLNode(nullptr, CXT_Element, "Option");
    setXmlAttr(option, "name", warp_options[i].first);
    CPLCreateXMLNode(option, CXT_Text, warp_options[i].second.c_str());
    option->psNext = options->psChild;
    options->psChild = option;
  }

  return serializeXml(root);
}

static GDALDataset *
updateNoDataValues(GDALDataset *warped_vrt_dataset,
                   const std::vector<double> &nodata_values, const Options &options) {
  std::string vrt_string = getVrtXml(warped_vrt_dataset);

  std::vector<std::pair<std::string, std::string>> warp_options;
  warp_options.push_back(std::make_pair(std::string("INIT_DEST"), std::string("NO_DATA")));
  warp_options.push_back(std::make_pair(std::string("UNIFIED_SRC_NODATA"), std::string("YES")));
  vrt_string = addWarpOptionsToVrt(vrt_string, warp_options);

  GDALDataset *corrected = static_cast<GDALDataset *>(
      GDALOpen(vrt_string.c_str(), GA_ReadOnly));
  if (corrected == nullptr)
    throw G2TException("Cannot open the modified warped VRT");

  std::string joined;
  for (size_t i = 0; i < nodata_values.size(); ++i) {
    if (i)
      joined += " ";
    std::ostringstream os;
    os << nodata_values[i];
    joined += os.str();
  }
  corrected->SetMetadataItem("NODATA_VALUES", joined.c_str());
  return corrected;
}

static GDALDataset *
updateAlphaValueForNonAlphaInputs(GDALDataset *warped_vrt_dataset, const Options &options) {
  int count = warped_vrt_dataset->GetRasterCount();
  if (count == 1 || count == 3) {
    std::string vrt_string = getVrtXml(warped_vrt_dataset);
    vrt_string = addAlphaBandToVrt(vrt_string);
    GDALDataset *corrected = static_cast<GDALDataset *>(
        GDALOpen(vrt_string.c_str(), GA_ReadOnly));
    if (corrected != nullptr)
      return corrected;
  }
  return warped_vrt_dataset;
}

// ---------------------------------------------------------------------------
// Reprojection
// ---------------------------------------------------------------------------

static GDALDataset *
warpByGdalWarp(GDALDataset *src, const std::vector<std::string> &args) {
  std::vector<char *> argv;
  for (size_t i = 0; i < args.size(); ++i)
    argv.push_back(CPLStrdup(args[i].c_str()));
  argv.push_back(nullptr);

  GDALWarpAppOptions *warp_options = GDALWarpAppOptionsNew(argv.data(), nullptr);
  GDALDatasetH src_ds = static_cast<GDALDatasetH>(src);
  int usage_error = 0;
  GDALDatasetH result = GDALWarp("", nullptr, 1, &src_ds, warp_options, &usage_error);
  GDALWarpAppOptionsFree(warp_options);
  for (size_t i = 0; i + 1 < argv.size(); ++i)
    CPLFree(argv[i]);

  if (result == nullptr)
    throw G2TException("gdal.Warp failed while reprojecting the dataset");
  return static_cast<GDALDataset *>(result);
}

static GDALDataset *
reprojectDataset(GDALDataset *from_dataset, OGRSpatialReference *from_srs,
                 OGRSpatialReference *to_srs, const Options &options) {
  if (from_srs == nullptr || to_srs == nullptr)
    throw G2TException("from and to SRS must be defined to reproject the dataset");

  if (proj4String(from_srs) != proj4String(to_srs) || from_dataset->GetGCPCount() != 0) {
    const char *authority_name = to_srs->GetAuthorityName(nullptr);
    const char *authority_code = to_srs->GetAuthorityCode(nullptr);
    if (from_srs->IsGeographic() && authority_name != nullptr &&
        strcmp(authority_name, "EPSG") == 0 && authority_code != nullptr &&
        strcmp(authority_code, "3857") == 0) {
      double from_gt[6];
      if (from_dataset->GetGeoTransform(from_gt) == CE_None &&
          from_gt[2] == 0 && from_gt[4] == 0 && from_gt[5] < 0) {
        double minlon = from_gt[0];
        double maxlon = from_gt[0] + from_dataset->GetRasterXSize() * from_gt[1];
        double maxlat = from_gt[3];
        double minlat = from_gt[3] + from_dataset->GetRasterYSize() * from_gt[5];
        const double MAX_LAT = 85.0511287798066;
        bool adjust = false;
        if (minlon < -180.0) { minlon = -180.0; adjust = true; }
        if (maxlon > 180.0) { maxlon = 180.0; adjust = true; }
        if (maxlat > MAX_LAT) { maxlat = MAX_LAT; adjust = true; }
        if (minlat < -MAX_LAT) { minlat = -MAX_LAT; adjust = true; }

        if (adjust) {
          OGRCoordinateTransformation *ct =
              OGRCreateCoordinateTransformation(from_srs, to_srs);
          double west = minlon, south = minlat, east = maxlon, north = maxlat;
          ct->Transform(1, &west, &south);
          ct->Transform(1, &east, &north);
          OGRCoordinateTransformation::DestroyCT(ct);

          std::vector<std::string> args;
          args.push_back("-of");
          args.push_back("VRT");
          args.push_back("-te");
          args.push_back(std::to_string(west));
          args.push_back(std::to_string(south));
          args.push_back(std::to_string(east));
          args.push_back(std::to_string(north));
          args.push_back("-s_srs");
          args.push_back(wktString(from_srs));
          args.push_back("-t_srs");
          args.push_back("EPSG:3857");
          return warpByGdalWarp(from_dataset, args);
        }
      }
    }

    std::string src_wkt = wktString(from_srs);
    std::string dst_wkt = wktString(to_srs);
    GDALDataset *to_dataset = static_cast<GDALDataset *>(GDALAutoCreateWarpedVRT(
        static_cast<GDALDatasetH>(from_dataset), src_wkt.c_str(), dst_wkt.c_str(),
        GRA_NearestNeighbour, 0.0, nullptr));
    if (to_dataset == nullptr)
      throw G2TException("GDALAutoCreateWarpedVRT failed");
    return to_dataset;
  }

  return from_dataset;
}

// ---------------------------------------------------------------------------
// Scaling / tile writing
// ---------------------------------------------------------------------------

static std::vector<std::string>
getCreationOptions(const Options &options, const std::string &driver) {
  std::vector<std::string> copts;
  if (driver == "WEBP") {
    if (options.webp_lossless)
      copts.push_back("LOSSLESS=True");
    else
      copts.push_back("QUALITY=" + std::to_string(options.webp_quality));
  } else if (driver == "JPEG") {
    copts.push_back("QUALITY=" + std::to_string(options.jpeg_quality));
  }
  return copts;
}

static GDALResampleAlg
resampleAlgFromName(const std::string &name) {
  if (name == "near") return GRA_NearestNeighbour;
  if (name == "bilinear") return GRA_Bilinear;
  if (name == "cubic") return GRA_Cubic;
  if (name == "cubicspline") return GRA_CubicSpline;
  if (name == "lanczos" || name == "antialias") return GRA_Lanczos;
  if (name == "mode") return GRA_Mode;
  if (name == "max") return GRA_Max;
  if (name == "min") return GRA_Min;
  if (name == "med") return GRA_Med;
  if (name == "q1") return GRA_Q1;
  if (name == "q3") return GRA_Q3;
  return GRA_NearestNeighbour;
}

void
scaleQueryToTile(GDALDataset *dsquery, GDALDataset *dstile, const Options &options,
                 const std::string &tilefilename) {
  int querysize = dsquery->GetRasterXSize();
  int tile_size = dstile->GetRasterXSize();
  int tilebands = dstile->GetRasterCount();

  double gtquery[6] = {0.0, tile_size / static_cast<double>(querysize), 0.0,
                       0.0, 0.0, tile_size / static_cast<double>(querysize)};
  dsquery->SetGeoTransform(gtquery);
  double gttile[6] = {0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  dstile->SetGeoTransform(gttile);

  if (options.resampling == "average" &&
      (!options.excluded_values.empty() || options.nodata_values_pct_threshold < 100)) {
    std::vector<std::string> args;
    args.push_back("-r");
    args.push_back("average");
    args.push_back("-wo");
    args.push_back("NODATA_VALUES_PCT_THRESHOLD=" +
                   std::to_string(options.nodata_values_pct_threshold));
    if (!options.excluded_values.empty()) {
      args.push_back("-wo");
      args.push_back("EXCLUDED_VALUES=" + options.excluded_values);
      args.push_back("-wo");
      args.push_back("EXCLUDED_VALUES_PCT_THRESHOLD=" +
                     std::to_string(options.excluded_values_pct_threshold));
    }

    std::vector<char *> argv;
    for (size_t i = 0; i < args.size(); ++i)
      argv.push_back(CPLStrdup(args[i].c_str()));
    argv.push_back(nullptr);

    GDALWarpAppOptions *warp_options = GDALWarpAppOptionsNew(argv.data(), nullptr);
    GDALDatasetH src_ds = static_cast<GDALDatasetH>(dsquery);
    int usage_error = 0;
    GDALDatasetH result = GDALWarp("", static_cast<GDALDatasetH>(dstile), 1,
                                   &src_ds, warp_options, &usage_error);
    GDALWarpAppOptionsFree(warp_options);
    for (size_t i = 0; i + 1 < argv.size(); ++i)
      CPLFree(argv[i]);

    if (result == nullptr)
      throw G2TException("gdal.Warp failed on " + tilefilename + " (error " +
                         std::to_string(usage_error) + ")");
    if (result != static_cast<GDALDatasetH>(dstile))
      GDALClose(result);
  } else if (options.resampling == "average") {
    for (int i = 1; i <= tilebands; ++i) {
      GDALRasterBandH dst_band = static_cast<GDALRasterBandH>(dstile->GetRasterBand(i));
      CPLErr res = GDALRegenerateOverviews(
          static_cast<GDALRasterBandH>(dsquery->GetRasterBand(i)), 1, &dst_band,
          "average", nullptr, nullptr);
      if (res != CE_None)
        throw G2TException("RegenerateOverview() failed on " + tilefilename);
    }
  } else {
    GDALResampleAlg alg = resampleAlgFromName(options.resampling);
    CPLErr res = GDALReprojectImage(
        static_cast<GDALDatasetH>(dsquery), nullptr,
        static_cast<GDALDatasetH>(dstile), nullptr, alg, 0.0, 0.0, nullptr,
        nullptr, nullptr);
    if (res != CE_None)
      throw G2TException("ReprojectImage() failed on " + tilefilename);
  }
}

static GDALDataset *
removeAlphaBand(GDALDataset *src_ds) {
  if (src_ds->GetRasterBand(src_ds->GetRasterCount())->GetColorInterpretation() !=
      GCI_AlphaBand)
    return src_ds;

  int new_band_count = src_ds->GetRasterCount() - 1;
  GDALDriver *mem_driver = GetGDALDriverManager()->GetDriverByName("MEM");
  GDALDataset *dst_ds = mem_driver->Create(
      "", src_ds->GetRasterXSize(), src_ds->GetRasterYSize(), new_band_count,
      src_ds->GetRasterBand(1)->GetRasterDataType(), nullptr);

  double gt[6];
  if (src_ds->GetGeoTransform(gt) == CE_None)
    dst_ds->SetGeoTransform(gt);
  const char *projection = src_ds->GetProjectionRef();
  if (projection != nullptr && projection[0] != '\0')
    dst_ds->SetProjection(projection);

  int xsize = src_ds->GetRasterXSize();
  int ysize = src_ds->GetRasterYSize();
  std::vector<GByte> buffer(static_cast<size_t>(xsize) * ysize);
  for (int i = 1; i <= new_band_count; ++i) {
    src_ds->GetRasterBand(i)->RasterIO(GF_Read, 0, 0, xsize, ysize, buffer.data(),
                                       xsize, ysize, GDT_Byte, 0, 0, nullptr);
    dst_ds->GetRasterBand(i)->RasterIO(GF_Write, 0, 0, xsize, ysize, buffer.data(),
                                       xsize, ysize, GDT_Byte, 0, 0, nullptr);
  }

  return dst_ds;
}

static void
writeTile(GDALDataset *dstile, const TileJobInfo &tile_job_info,
          const Options &options, const std::string &tilefilename) {
  GDALDriver *out_driver = GetGDALDriverManager()->GetDriverByName(
      tile_job_info.tile_driver.c_str());
  GDALDataset *out_src = dstile;
  if (tile_job_info.tile_driver == "JPEG")
    out_src = removeAlphaBand(dstile);

  std::vector<std::string> copts = getCreationOptions(options, tile_job_info.tile_driver);
  CPLStringList creation_options;
  for (size_t i = 0; i < copts.size(); ++i)
    creation_options.AddString(copts[i].c_str());

  GDALDataset *created = out_driver->CreateCopy(tilefilename.c_str(), out_src, FALSE,
                         creation_options.List(), nullptr, nullptr);
  if (created != nullptr)
    GDALClose(created);
  if (out_src != dstile)
    GDALClose(out_src);

  std::string aux_xml = tilefilename + ".aux.xml";
  VSIStatBufL stat;
  if (VSIStatL(aux_xml.c_str(), &stat) == 0)
    VSIUnlink(aux_xml.c_str());
}

// Thread local cached source dataset
struct ThreadDataset {
  GDALDataset *ds = nullptr;
  std::string file;
  ~ThreadDataset() {
    if (ds != nullptr)
      GDALClose(ds);
  }
};

static thread_local ThreadDataset tls_dataset;

static GDALDataset *
getSourceDataset(const std::string &filename) {
  if (tls_dataset.ds != nullptr && tls_dataset.file == filename)
    return tls_dataset.ds;
  if (tls_dataset.ds != nullptr) {
    GDALClose(tls_dataset.ds);
    tls_dataset.ds = nullptr;
    tls_dataset.file.clear();
  }
  GDALDataset *ds = static_cast<GDALDataset *>(GDALOpen(filename.c_str(), GA_ReadOnly));
  if (ds == nullptr)
    throw G2TException("Cannot open the source dataset: " + filename);
  tls_dataset.ds = ds;
  tls_dataset.file = filename;
  return ds;
}

// ---------------------------------------------------------------------------
// Base tile creation
// ---------------------------------------------------------------------------

void
createBaseTile(const TileJobInfo &tile_job_info, const TileDetail &tile_detail) {
  const Options &options = *tile_job_info.options;
  int dataBandsCount = tile_job_info.nb_data_bands;
  int tilebands = dataBandsCount + 1;
  int tile_size = tile_job_info.tile_size;

  GDALDriver *mem_driver = GetGDALDriverManager()->GetDriverByName("MEM");
  if (mem_driver == nullptr)
    throw G2TException("The 'MEM' driver was not found");

  GDALDataset *ds = getSourceDataset(tile_job_info.src_file);
  GDALRasterBand *alphaband = ds->GetRasterBand(1)->GetMaskBand();

  std::string tilefilename = joinPath(
      joinPath(joinPath(tile_job_info.output_file_path,
                        std::to_string(tile_detail.tz)),
               std::to_string(tile_detail.tx)),
      std::to_string(tile_detail.ty) + "." + tile_job_info.tile_extension);
  makedirs(dirname(tilefilename));

  GDALDataset *dstile = mem_driver->Create("", tile_size, tile_size, tilebands,
                                           GDT_Byte, nullptr);
  if (dstile == nullptr)
    throw G2TException("Cannot create the in-memory tile");
  dstile->GetRasterBand(tilebands)->SetColorInterpretation(GCI_AlphaBand);

  int rx = tile_detail.rx, ry = tile_detail.ry;
  int rxsize = tile_detail.rxsize, rysize = tile_detail.rysize;
  int wx = tile_detail.wx, wy = tile_detail.wy;
  int wxsize = tile_detail.wxsize, wysize = tile_detail.wysize;

  std::vector<GByte> data;
  std::vector<GByte> alpha;

  if (rxsize != 0 && rysize != 0 && wxsize != 0 && wysize != 0) {
    alpha.resize(static_cast<size_t>(wxsize) * wysize);
    alphaband->RasterIO(GF_Read, rx, ry, rxsize, rysize, alpha.data(), wxsize,
                        wysize, GDT_Byte, 0, 0, nullptr);

    if (tile_job_info.exclude_transparent) {
      bool all_transparent = true;
      for (size_t i = 0; i < alpha.size(); ++i) {
        if (alpha[i] != 0) { all_transparent = false; break; }
      }
      if (all_transparent) {
        GDALClose(dstile);
        return;
      }
    }

    data.resize(static_cast<size_t>(wxsize) * wysize * dataBandsCount);
    std::vector<int> band_map(dataBandsCount);
    for (int i = 0; i < dataBandsCount; ++i)
      band_map[i] = i + 1;
    ds->RasterIO(GF_Read, rx, ry, rxsize, rysize, data.data(), wxsize, wysize,
                 GDT_Byte, dataBandsCount, band_map.data(), 0, 0, 0, nullptr);
  }

  if (!data.empty()) {
    std::vector<int> band_map(dataBandsCount);
    for (int i = 0; i < dataBandsCount; ++i)
      band_map[i] = i + 1;
    int alpha_band = tilebands;

    if (tile_size == tile_detail.querysize) {
      dstile->RasterIO(GF_Write, wx, wy, wxsize, wysize, data.data(), wxsize,
                       wysize, GDT_Byte, dataBandsCount, band_map.data(), 0, 0,
                       0, nullptr);
      dstile->RasterIO(GF_Write, wx, wy, wxsize, wysize, alpha.data(), wxsize,
                       wysize, GDT_Byte, 1, &alpha_band, 0, 0, 0, nullptr);
    } else {
      GDALDataset *dsquery = mem_driver->Create("", tile_detail.querysize,
                                                tile_detail.querysize, tilebands,
                                                GDT_Byte, nullptr);
      dsquery->GetRasterBand(tilebands)->SetColorInterpretation(GCI_AlphaBand);
      dsquery->RasterIO(GF_Write, wx, wy, wxsize, wysize, data.data(), wxsize,
                        wysize, GDT_Byte, dataBandsCount, band_map.data(), 0, 0,
                        0, nullptr);
      dsquery->RasterIO(GF_Write, wx, wy, wxsize, wysize, alpha.data(), wxsize,
                        wysize, GDT_Byte, 1, &alpha_band, 0, 0, 0, nullptr);
      scaleQueryToTile(dsquery, dstile, options, tilefilename);
      GDALClose(dsquery);
    }
  }

  writeTile(dstile, tile_job_info, options, tilefilename);
  GDALClose(dstile);
}

// ---------------------------------------------------------------------------
// Overview tile creation
// ---------------------------------------------------------------------------

int
countOverviewTiles(const TileJobInfo &tile_job_info) {
  int tile_number = 0;
  for (int tz = tile_job_info.tmaxz - 1; tz >= tile_job_info.tminz; --tz) {
    int tminx = tile_job_info.tminmax[tz][0];
    int tminy = tile_job_info.tminmax[tz][1];
    int tmaxx = tile_job_info.tminmax[tz][2];
    int tmaxy = tile_job_info.tminmax[tz][3];
    tile_number += (1 + std::abs(tmaxx - tminx)) * (1 + std::abs(tmaxy - tminy));
  }
  return tile_number;
}

std::vector<std::vector<std::pair<int, int>>>
groupOverviewBaseTiles(int base_tz, const std::string &output_folder,
                       const TileJobInfo &tile_job_info) {
  std::map<std::pair<int, int>, size_t> index;
  std::vector<std::vector<std::pair<int, int>>> groups;

  int tminx = tile_job_info.tminmax[base_tz][0];
  int tminy = tile_job_info.tminmax[base_tz][1];
  int tmaxx = tile_job_info.tminmax[base_tz][2];
  int tmaxy = tile_job_info.tminmax[base_tz][3];

  for (int ty = tmaxy; ty >= tminy; --ty) {
    int overview_ty = ty >> 1;
    for (int tx = tminx; tx <= tmaxx; ++tx) {
      int overview_tx = tx >> 1;
      std::pair<int, int> key(overview_tx, overview_ty);
      std::map<std::pair<int, int>, size_t>::iterator it = index.find(key);
      size_t group;
      if (it == index.end()) {
        group = groups.size();
        index[key] = group;
        groups.push_back(std::vector<std::pair<int, int>>());
      } else {
        group = it->second;
      }
      groups[group].push_back(std::make_pair(tx, ty));
    }
  }

  int overview_tz = base_tz - 1;
  for (int tx = tminx; tx <= tmaxx; ++tx) {
    int overview_tx = tx >> 1;
    makedirs(joinPath(joinPath(output_folder, std::to_string(overview_tz)),
                      std::to_string(overview_tx)));
  }

  return groups;
}

void
createOverviewTile(int base_tz, const std::vector<std::pair<int, int>> &base_tiles,
                   const std::string &output_folder,
                   const TileJobInfo &tile_job_info, const Options &options) {
  int overview_tz = base_tz - 1;
  int overview_tx = base_tiles[0].first >> 1;
  int overview_ty = base_tiles[0].second >> 1;
  int overview_ty_real = GDAL2Tiles::getYTile(overview_ty, overview_tz, options);

  std::string tilefilename = joinPath(
      joinPath(joinPath(output_folder, std::to_string(overview_tz)),
               std::to_string(overview_tx)),
      std::to_string(overview_ty_real) + "." + tile_job_info.tile_extension);
  if (options.resume && fileExists(tilefilename))
    return;

  GDALDriver *mem_driver = GetGDALDriverManager()->GetDriverByName("MEM");
  int tile_size = tile_job_info.tile_size;
  int tilebands = tile_job_info.nb_data_bands + 1;

  GDALDataset *dsquery = mem_driver->Create("", 2 * tile_size, 2 * tile_size,
                                            tilebands, GDT_Byte, nullptr);
  dsquery->GetRasterBand(tilebands)->SetColorInterpretation(GCI_AlphaBand);
  GDALDataset *dstile = mem_driver->Create("", tile_size, tile_size, tilebands,
                                           GDT_Byte, nullptr);
  dstile->GetRasterBand(tilebands)->SetColorInterpretation(GCI_AlphaBand);

  std::vector<int> band_map(tilebands);
  for (int i = 0; i < tilebands; ++i)
    band_map[i] = i + 1;
  std::vector<GByte> mask(static_cast<size_t>(tile_size) * tile_size, 255);
  std::vector<GByte> base_data(static_cast<size_t>(tile_size) * tile_size * tilebands);

  bool any_usable = false;

  for (size_t t = 0; t < base_tiles.size(); ++t) {
    int base_tx = base_tiles[t].first;
    int base_ty = base_tiles[t].second;
    int base_ty_real = GDAL2Tiles::getYTile(base_ty, base_tz, options);

    std::string base_tile_path = joinPath(
        joinPath(joinPath(output_folder, std::to_string(base_tz)),
                 std::to_string(base_tx)),
        std::to_string(base_ty_real) + "." + tile_job_info.tile_extension);
    if (!fileExists(base_tile_path))
      continue;

    GDALDataset *dsquerytile = static_cast<GDALDataset *>(
        GDALOpen(base_tile_path.c_str(), GA_ReadOnly));
    if (dsquerytile == nullptr)
      continue;

    int tileposx = (base_tx % 2 == 0) ? 0 : tile_size;
    int tileposy;
    if (options.xyz && options.profile == "raster")
      tileposy = (base_ty % 2 == 0) ? 0 : tile_size;
    else
      tileposy = (base_ty % 2 == 0) ? tile_size : 0;

    GDALDataset *owned_tmp = nullptr;
    if (tile_job_info.tile_driver == "JPEG" && dsquerytile->GetRasterCount() == 3 &&
        tilebands == 2) {
      owned_tmp = mem_driver->Create("", dsquerytile->GetRasterXSize(),
                                     dsquerytile->GetRasterYSize(), 2, GDT_Byte,
                                     nullptr);
      std::vector<GByte> band1(static_cast<size_t>(tile_size) * tile_size);
      dsquerytile->GetRasterBand(1)->RasterIO(GF_Read, 0, 0, tile_size, tile_size,
                                              band1.data(), tile_size, tile_size,
                                              GDT_Byte, 0, 0, nullptr);
      owned_tmp->GetRasterBand(1)->RasterIO(GF_Write, 0, 0, tile_size, tile_size,
                                            band1.data(), tile_size, tile_size,
                                            GDT_Byte, 0, 0, nullptr);
      owned_tmp->GetRasterBand(2)->RasterIO(GF_Write, 0, 0, tile_size, tile_size,
                                            mask.data(), tile_size, tile_size,
                                            GDT_Byte, 0, 0, nullptr);
      owned_tmp->GetRasterBand(2)->SetColorInterpretation(GCI_AlphaBand);
      dsquerytile = owned_tmp;
    } else if (dsquerytile->GetRasterCount() == tilebands - 1) {
      owned_tmp = mem_driver->CreateCopy("", dsquerytile, FALSE, nullptr, nullptr,
                                         nullptr);
      owned_tmp->AddBand(GDT_Byte, nullptr);
      int alpha_band = tilebands;
      owned_tmp->RasterIO(GF_Write, 0, 0, tile_size, tile_size, mask.data(),
                          tile_size, tile_size, GDT_Byte, 1, &alpha_band, 0, 0, 0,
                          nullptr);
      dsquerytile = owned_tmp;
    } else if (dsquerytile->GetRasterCount() != tilebands) {
      int got = dsquerytile->GetRasterCount();
      GDALClose(dsquerytile);
      throw G2TException("Unexpected number of bands in base tile. Got " +
                         std::to_string(got) + ", expected " +
                         std::to_string(tilebands));
    }

    dsquerytile->RasterIO(GF_Read, 0, 0, tile_size, tile_size, base_data.data(),
                          tile_size, tile_size, GDT_Byte, tilebands,
                          band_map.data(), 0, 0, 0, nullptr);
    dsquery->RasterIO(GF_Write, tileposx, tileposy, tile_size, tile_size,
                      base_data.data(), tile_size, tile_size, GDT_Byte, tilebands,
                      band_map.data(), 0, 0, 0, nullptr);

    GDALClose(dsquerytile);
    any_usable = true;
  }

  if (!any_usable) {
    GDALClose(dsquery);
    GDALClose(dstile);
    return;
  }

  scaleQueryToTile(dsquery, dstile, options, tilefilename);
  writeTile(dstile, tile_job_info, options, tilefilename);

  GDALClose(dsquery);
  GDALClose(dstile);
}

// ---------------------------------------------------------------------------
// GDAL2Tiles
// ---------------------------------------------------------------------------

GDAL2Tiles::GDAL2Tiles(const std::string &input_file,
                       const std::string &output_folder, const Options &options)
    : m_input_file(input_file), m_output_folder(output_folder), m_options(options),
      m_input_dataset(nullptr), m_warped_input_dataset(nullptr), m_in_srs(nullptr),
      m_out_srs(nullptr), m_mercator(nullptr), m_geodetic(nullptr) {
  m_tile_size = m_options.tilesize ? m_options.tilesize : 256;
  m_tiledriver = m_options.tiledriver;
  if (m_tiledriver == "PNG")
    m_tileext = "png";
  else if (m_tiledriver == "WEBP")
    m_tileext = "webp";
  else
    m_tileext = "jpg";

  m_tmp_vrt_filename = CPLGenerateTempFilenameSafe("gdal2tiles") + ".vrt";

  m_scaledquery = true;
  m_querysize = 4 * m_tile_size;
  m_overviewquery = false;
  if (m_options.resampling == "near")
    m_querysize = m_tile_size;
  else if (m_options.resampling == "bilinear")
    m_querysize = 2 * m_tile_size;

  m_tminz = m_options.tminz;
  m_tmaxz = m_options.tmaxz;

  m_dataBandsCount = 0;
  m_nativezoom = 0;
  m_out_gt[0] = 0.0; m_out_gt[1] = 1.0; m_out_gt[2] = 0.0;
  m_out_gt[3] = 0.0; m_out_gt[4] = 0.0; m_out_gt[5] = 1.0;
  m_ominx = m_omaxx = m_omaxy = m_ominy = 0.0;
  m_isepsg4326 = false;
}

GDAL2Tiles::~GDAL2Tiles() {
  if (m_warped_input_dataset != nullptr && m_warped_input_dataset != m_input_dataset)
    GDALClose(m_warped_input_dataset);
  if (m_input_dataset != nullptr)
    GDALClose(m_input_dataset);
  delete m_in_srs;
  delete m_out_srs;
  delete m_mercator;
  delete m_geodetic;

  if (!m_tmp_vrt_filename.empty()) {
    VSIUnlink(m_tmp_vrt_filename.c_str());
    std::string aux = m_tmp_vrt_filename + ".aux.xml";
    VSIStatBufL stat;
    if (VSIStatL(aux.c_str(), &stat) == 0)
      VSIUnlink(aux.c_str());
  }
}

int
GDAL2Tiles::getYTile(int ty, int tz, const Options &options) {
  if (options.xyz && options.profile != "raster") {
    if (options.profile == "mercator" || options.profile == "geodetic")
      return static_cast<int>(std::pow(2.0, tz)) - 1 - ty;
    return ty;
  }
  return ty;
}

void
GDAL2Tiles::openInput() {
  GDALAllRegister();
  GDALDriverManager *driver_manager = GetGDALDriverManager();
  GDALDriver *out_drv = driver_manager->GetDriverByName(m_tiledriver.c_str());
  GDALDriver *mem_drv = driver_manager->GetDriverByName("MEM");

  if (out_drv == nullptr)
    throw G2TException("The '" + m_tiledriver +
                       "' driver was not found, is it available in this GDAL build?");
  if (mem_drv == nullptr)
    throw G2TException("The 'MEM' driver was not found, is it available in this GDAL build?");

  m_input_dataset = static_cast<GDALDataset *>(
      GDALOpen(m_input_file.c_str(), GA_ReadOnly));
  if (m_input_dataset == nullptr)
    throw G2TException("It is not possible to open the input file '" + m_input_file + "'.");
  if (m_input_dataset->GetRasterCount() == 0)
    throw G2TException("Input file '" + m_input_file + "' has no raster band");
  if (m_input_dataset->GetRasterBand(1)->GetColorTable() != nullptr)
    throw G2TException("Please convert this file to RGB/RGBA and run gdal2tiles on the result.");
  if (m_input_dataset->GetRasterBand(1)->GetRasterDataType() != GDT_Byte)
    throw G2TException("Please convert this file to 8-bit and run gdal2tiles on the result.");

  std::vector<double> in_nodata = setupNoDataValues(m_input_dataset, m_options);

  m_in_srs = setupInputSrs(m_input_dataset, m_options, m_in_srs_wkt);
  m_out_srs = setupOutputSrs(m_in_srs, m_options);

  m_warped_input_dataset = nullptr;
  if (m_options.profile != "raster") {
    if (m_in_srs == nullptr)
      throw G2TException("Input file has unknown SRS. Use --s_srs EPSG:xyz to provide source reference system.");
    if (!hasGeoreference(m_input_dataset))
      throw G2TException("There is no georeference - neither affine transformation (worldfile) nor GCPs.");

    if (proj4String(m_in_srs) != proj4String(m_out_srs) ||
        m_input_dataset->GetGCPCount() != 0) {
      m_warped_input_dataset = reprojectDataset(m_input_dataset, m_in_srs, m_out_srs, m_options);
      if (in_nodata.empty())
        m_warped_input_dataset =
            updateAlphaValueForNonAlphaInputs(m_warped_input_dataset, m_options);
    } else {
      GDALDriver *vrt_drv = driver_manager->GetDriverByName("VRT");
      m_warped_input_dataset = vrt_drv->CreateCopy("", m_input_dataset, FALSE,
                                                   nullptr, nullptr, nullptr);
    }

    if (!in_nodata.empty())
      m_warped_input_dataset =
          updateNoDataValues(m_warped_input_dataset, in_nodata, m_options);
  }

  if (m_warped_input_dataset == nullptr)
    m_warped_input_dataset = m_input_dataset;

  GDALDriver *vrt_drv = driver_manager->GetDriverByName("VRT");
  GDALDataset *tmp_ds = vrt_drv->CreateCopy(m_tmp_vrt_filename.c_str(),
                                            m_warped_input_dataset, FALSE, nullptr,
                                            nullptr, nullptr);
  if (tmp_ds == nullptr)
    throw G2TException("Cannot create the temporary VRT file: " + m_tmp_vrt_filename);
  GDALClose(tmp_ds);

  m_dataBandsCount = nbDataBands(m_warped_input_dataset);

  OGRSpatialReference srs4326;
  srs4326.importFromEPSG(4326);
  srs4326.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
  m_isepsg4326 = m_out_srs != nullptr && proj4String(&srs4326) == proj4String(m_out_srs);

  m_warped_input_dataset->GetGeoTransform(m_out_gt);
  if (m_out_gt[2] != 0.0 || m_out_gt[4] != 0.0)
    throw G2TException("Georeference of the raster contains rotation or skew. Such raster is not supported.");

  m_ominx = m_out_gt[0];
  m_omaxx = m_out_gt[0] + m_warped_input_dataset->GetRasterXSize() * m_out_gt[1];
  m_omaxy = m_out_gt[3];
  m_ominy = m_out_gt[3] - m_warped_input_dataset->GetRasterYSize() * m_out_gt[1];

  if (m_options.profile == "mercator") {
    m_mercator = new GlobalMercator(m_tile_size);
    m_tminmax.resize(MAXZOOMLEVEL);
    for (int tz = 0; tz < MAXZOOMLEVEL; ++tz) {
      int tminx, tminy, tmaxx, tmaxy;
      m_mercator->MetersToTile(m_ominx, m_ominy, tz, tminx, tminy);
      m_mercator->MetersToTile(m_omaxx, m_omaxy, tz, tmaxx, tmaxy);
      tminx = std::max(0, tminx); tminy = std::max(0, tminy);
      tmaxx = std::min(static_cast<int>(std::pow(2.0, tz)) - 1, tmaxx);
      tmaxy = std::min(static_cast<int>(std::pow(2.0, tz)) - 1, tmaxy);
      std::array<int, 4> bounds = {{tminx, tminy, tmaxx, tmaxy}};
      m_tminmax[tz] = bounds;
    }

    if (m_tminz < 0)
      m_tminz = m_mercator->ZoomForPixelSize(
          m_out_gt[1] * std::max(m_warped_input_dataset->GetRasterXSize(),
                                 m_warped_input_dataset->GetRasterYSize()) /
          static_cast<double>(m_tile_size));
    if (m_tmaxz < 0) {
      m_tmaxz = m_mercator->ZoomForPixelSize(m_out_gt[1]);
      m_tmaxz = std::max(m_tminz, m_tmaxz);
    }
    m_tminz = std::min(m_tminz, m_tmaxz);

  } else if (m_options.profile == "geodetic") {
    m_geodetic = new GlobalGeodetic(m_options.tmscompatible, m_tile_size);
    m_tminmax.resize(MAXZOOMLEVEL);
    for (int tz = 0; tz < MAXZOOMLEVEL; ++tz) {
      int tminx, tminy, tmaxx, tmaxy;
      m_geodetic->LonLatToTile(m_ominx, m_ominy, tz, tminx, tminy);
      m_geodetic->LonLatToTile(m_omaxx, m_omaxy, tz, tmaxx, tmaxy);
      tminx = std::max(0, tminx); tminy = std::max(0, tminy);
      tmaxx = std::min(static_cast<int>(std::pow(2.0, tz + 1)) - 1, tmaxx);
      tmaxy = std::min(static_cast<int>(std::pow(2.0, tz)) - 1, tmaxy);
      std::array<int, 4> bounds = {{tminx, tminy, tmaxx, tmaxy}};
      m_tminmax[tz] = bounds;
    }

    if (m_tminz < 0)
      m_tminz = m_geodetic->ZoomForPixelSize(
          m_out_gt[1] * std::max(m_warped_input_dataset->GetRasterXSize(),
                                 m_warped_input_dataset->GetRasterYSize()) /
          static_cast<double>(m_tile_size));
    if (m_tmaxz < 0) {
      m_tmaxz = m_geodetic->ZoomForPixelSize(m_out_gt[1]);
      m_tmaxz = std::max(m_tminz, m_tmaxz);
    }
    m_tminz = std::min(m_tminz, m_tmaxz);

  } else {
    double native = std::max(
        std::ceil(std::log10(m_warped_input_dataset->GetRasterXSize() /
                             static_cast<double>(m_tile_size)) /
                  std::log10(2.0)),
        std::ceil(std::log10(m_warped_input_dataset->GetRasterYSize() /
                             static_cast<double>(m_tile_size)) /
                  std::log10(2.0)));
    m_nativezoom = std::max(0, static_cast<int>(native));

    if (m_tminz < 0)
      m_tminz = 0;
    if (m_tmaxz < 0) {
      m_tmaxz = m_nativezoom;
      m_tmaxz = std::max(m_tminz, m_tmaxz);
    } else if (m_tmaxz > m_nativezoom) {
      int oversample_factor = 1 << (m_tmaxz - m_nativezoom);
      std::string resample_alg = m_options.resampling;
      if (!(resample_alg == "near" || resample_alg == "bilinear" ||
            resample_alg == "average" || resample_alg == "cubic" ||
            resample_alg == "cubicspline" || resample_alg == "lanczos" ||
            resample_alg == "mode"))
        resample_alg = "bilinear";

      std::vector<std::string> args;
      args.push_back("-of");
      args.push_back("VRT");
      args.push_back("-outsize");
      args.push_back(std::to_string(m_warped_input_dataset->GetRasterXSize() * oversample_factor));
      args.push_back(std::to_string(m_warped_input_dataset->GetRasterYSize() * oversample_factor));
      args.push_back("-r");
      args.push_back(resample_alg);

      std::vector<char *> argv;
      for (size_t i = 0; i < args.size(); ++i)
        argv.push_back(CPLStrdup(args[i].c_str()));
      argv.push_back(nullptr);
      GDALTranslateOptions *translate_options = GDALTranslateOptionsNew(argv.data(), nullptr);
      GDALDatasetH src_ds = static_cast<GDALDatasetH>(m_warped_input_dataset);
      int usage_error = 0;
      GDALDatasetH translated = GDALTranslate(m_tmp_vrt_filename.c_str(), src_ds,
                                              translate_options, &usage_error);
      GDALTranslateOptionsFree(translate_options);
      for (size_t i = 0; i + 1 < argv.size(); ++i)
        CPLFree(argv[i]);
      if (translated == nullptr)
        throw G2TException("gdal.Translate failed while oversampling");

      m_warped_input_dataset = static_cast<GDALDataset *>(translated);
      m_warped_input_dataset->GetGeoTransform(m_out_gt);
      m_nativezoom = m_tmaxz;
    }

    m_tminmax.resize(m_tmaxz + 1);
    m_tsize.resize(m_tmaxz + 1);
    for (int tz = 0; tz <= m_tmaxz; ++tz) {
      double tsize = std::pow(2.0, m_nativezoom - tz) * m_tile_size;
      int tmaxx = static_cast<int>(std::ceil(
          m_warped_input_dataset->GetRasterXSize() / tsize)) - 1;
      int tmaxy = static_cast<int>(std::ceil(
          m_warped_input_dataset->GetRasterYSize() / tsize)) - 1;
      std::array<int, 4> bounds = {{0, 0, tmaxx, tmaxy}};
      m_tminmax[tz] = bounds;
      m_tsize[tz] = static_cast<int>(std::ceil(tsize));
    }
  }
}

void
GDAL2Tiles::geoQuery(GDALDataset *ds, double ulx, double uly, double lrx,
                     double lry, int querysize, int &rx, int &ry, int &rxsize,
                     int &rysize, int &wx, int &wy, int &wxsize, int &wysize) const {
  double geotran[6];
  ds->GetGeoTransform(geotran);

  rx = static_cast<int>((ulx - geotran[0]) / geotran[1] + 0.001);
  ry = static_cast<int>((uly - geotran[3]) / geotran[5] + 0.001);
  rxsize = std::max(1, static_cast<int>((lrx - ulx) / geotran[1] + 0.5));
  rysize = std::max(1, static_cast<int>((lry - uly) / geotran[5] + 0.5));

  if (!querysize) {
    wxsize = rxsize;
    wysize = rysize;
  } else {
    wxsize = querysize;
    wysize = querysize;
  }

  wx = 0;
  if (rx < 0) {
    int rxshift = std::abs(rx);
    wx = static_cast<int>(wxsize * (static_cast<double>(rxshift) / rxsize));
    wxsize = wxsize - wx;
    rxsize = rxsize - static_cast<int>(rxsize * (static_cast<double>(rxshift) / rxsize));
    rx = 0;
  }
  if (rx + rxsize > ds->GetRasterXSize()) {
    wxsize = static_cast<int>(wxsize *
                              (static_cast<double>(ds->GetRasterXSize() - rx) / rxsize));
    rxsize = ds->GetRasterXSize() - rx;
  }

  wy = 0;
  if (ry < 0) {
    int ryshift = std::abs(ry);
    wy = static_cast<int>(wysize * (static_cast<double>(ryshift) / rysize));
    wysize = wysize - wy;
    rysize = rysize - static_cast<int>(rysize * (static_cast<double>(ryshift) / rysize));
    ry = 0;
  }
  if (ry + rysize > ds->GetRasterYSize()) {
    wysize = static_cast<int>(wysize *
                              (static_cast<double>(ds->GetRasterYSize() - ry) / rysize));
    rysize = ds->GetRasterYSize() - ry;
  }
}

void
GDAL2Tiles::generateBaseTiles(TileJobInfo &tile_job_info,
                              std::vector<TileDetail> &tile_details) {
  int tminx = m_tminmax[m_tmaxz][0];
  int tminy = m_tminmax[m_tmaxz][1];
  int tmaxx = m_tminmax[m_tmaxz][2];
  int tmaxy = m_tminmax[m_tmaxz][3];

  GDALDataset *ds = m_warped_input_dataset;
  int querysize = m_querysize;
  int tz = m_tmaxz;

  for (int tx = tminx; tx <= tmaxx; ++tx)
    makedirs(joinPath(joinPath(m_output_folder, std::to_string(tz)),
                      std::to_string(tx)));

  for (int ty = tmaxy; ty >= tminy; --ty) {
    for (int tx = tminx; tx <= tmaxx; ++tx) {
      int ytile = getYTile(ty, tz, m_options);
      std::string tilefilename = joinPath(
          joinPath(joinPath(m_output_folder, std::to_string(tz)),
                   std::to_string(tx)),
          std::to_string(ytile) + "." + m_tileext);
      if (m_options.resume && fileExists(tilefilename))
        continue;

      int rx = 0, ry = 0, rxsize = 0, rysize = 0;
      int wx = 0, wy = 0, wxsize = 0, wysize = 0;

      if (m_options.profile != "raster") {
        double b[4];
        if (m_options.profile == "mercator")
          m_mercator->TileBounds(tx, ty, tz, b);
        else
          m_geodetic->TileBounds(tx, ty, tz, b);

        geoQuery(ds, b[0], b[3], b[2], b[1], 0, rx, ry, rxsize, rysize, wx, wy,
                 wxsize, wysize);
        geoQuery(ds, b[0], b[3], b[2], b[1], querysize, rx, ry, rxsize, rysize,
                 wx, wy, wxsize, wysize);
      } else {
        int tsize = m_tsize[tz];
        int xsize = ds->GetRasterXSize();
        int ysize = ds->GetRasterYSize();
        querysize = m_tile_size;

        rx = tx * tsize;
        rxsize = 0;
        if (tx == tmaxx)
          rxsize = xsize % tsize;
        if (rxsize == 0)
          rxsize = tsize;

        ry = ty * tsize;
        rysize = 0;
        if (ty == tmaxy)
          rysize = ysize % tsize;
        if (rysize == 0)
          rysize = tsize;

        wx = 0; wy = 0;
        wxsize = static_cast<int>(rxsize / static_cast<double>(tsize) * m_tile_size);
        wysize = static_cast<int>(rysize / static_cast<double>(tsize) * m_tile_size);

        if (!m_options.xyz) {
          ry = ysize - (ty * tsize) - rysize;
          if (wysize != m_tile_size)
            wy = m_tile_size - wysize;
        }
      }

      if (rxsize == 0 || rysize == 0 || wxsize == 0 || wysize == 0)
        continue;

      TileDetail tile_detail;
      tile_detail.tx = tx;
      tile_detail.ty_tms = ty;
      tile_detail.ty = ytile;
      tile_detail.tz = tz;
      tile_detail.rx = rx;
      tile_detail.ry = ry;
      tile_detail.rxsize = rxsize;
      tile_detail.rysize = rysize;
      tile_detail.wx = wx;
      tile_detail.wy = wy;
      tile_detail.wxsize = wxsize;
      tile_detail.wysize = wysize;
      tile_detail.querysize = querysize;
      tile_details.push_back(tile_detail);
    }
  }

  tile_job_info.src_file = m_tmp_vrt_filename;
  tile_job_info.nb_data_bands = m_dataBandsCount;
  tile_job_info.output_file_path = m_output_folder;
  tile_job_info.tile_extension = m_tileext;
  tile_job_info.tile_driver = m_tiledriver;
  tile_job_info.tile_size = m_tile_size;
  tile_job_info.tminmax = m_tminmax;
  tile_job_info.tminz = m_tminz;
  tile_job_info.tmaxz = m_tmaxz;
  tile_job_info.in_srs_wkt = m_in_srs_wkt;
  for (int i = 0; i < 6; ++i)
    tile_job_info.out_geo_trans[i] = m_out_gt[i];
  tile_job_info.ominy = m_ominy;
  tile_job_info.is_epsg_4326 = m_isepsg4326;
  tile_job_info.exclude_transparent = m_options.exclude_transparent;
  tile_job_info.options = &m_options;
}

void
GDAL2Tiles::generateMetadata() {
  makedirs(m_output_folder);

  if (m_options.xyz)
    return;

  double south, west, north, east;
  if (m_options.profile == "mercator") {
    double minLat, minLon, maxLat, maxLon;
    m_mercator->MetersToLatLon(m_ominx, m_ominy, minLat, minLon);
    m_mercator->MetersToLatLon(m_omaxx, m_omaxy, maxLat, maxLon);
    south = std::max(-85.05112878, minLat);
    west = std::max(-180.0, minLon);
    north = std::min(85.05112878, maxLat);
    east = std::min(180.0, maxLon);
  } else {
    south = std::max(-90.0, m_ominy);
    west = std::max(-180.0, m_ominx);
    north = std::min(90.0, m_omaxy);
    east = std::min(180.0, m_omaxx);
  }

  std::string srs;
  if (m_options.profile == "mercator")
    srs = "EPSG:3857";
  else if (m_options.profile == "geodetic")
    srs = "EPSG:4326";
  else if (!m_options.s_srs.empty())
    srs = m_options.s_srs;
  else
    srs = wktString(m_out_srs);

  std::string title = m_options.title;

  std::ostringstream os;
  os << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
     << "<TileMap version=\"1.0.0\" tilemapservice=\"http://tms.osgeo.org/1.0.0\">\n"
     << "  <Title>" << title << "</Title>\n"
     << "  <Abstract></Abstract>\n"
     << "  <SRS>" << srs << "</SRS>\n";
  os << "  <BoundingBox minx=\"" << west << "\" miny=\"" << south
     << "\" maxx=\"" << east << "\" maxy=\"" << north << "\"/>\n";
  os << "  <Origin x=\"" << west << "\" y=\"" << south << "\"/>\n";
  os << "  <TileFormat width=\"" << m_tile_size << "\" height=\"" << m_tile_size
     << "\" mime-type=\"image/" << m_tileext << "\" extension=\"" << m_tileext
     << "\"/>\n";
  os << "  <TileSets profile=\"" << m_options.profile << "\">\n";
  for (int z = m_tminz; z <= m_tmaxz; ++z) {
    double upp;
    if (m_options.profile == "raster")
      upp = std::pow(2.0, m_nativezoom - z) * m_out_gt[1];
    else if (m_options.profile == "mercator")
      upp = 156543.0339 / std::pow(2.0, z);
    else
      upp = 0.703125 / std::pow(2.0, z);
    os << "        <TileSet href=\"" << z << "\" units-per-pixel=\"" << upp
       << "\" order=\"" << z << "\"/>\n";
  }
  os << "      </TileSets>\n</TileMap>\n";

  std::string filename = joinPath(m_output_folder, "tilemapresource.xml");
  VSILFILE *fp = VSIFOpenL(filename.c_str(), "wb");
  if (fp != nullptr) {
    std::string content = os.str();
    VSIFWriteL(content.c_str(), 1, content.size(), fp);
    VSIFCloseL(fp);
  }
}

// ---------------------------------------------------------------------------
// Progress reporting
// ---------------------------------------------------------------------------

static std::string
formatDuration(double seconds) {
  if (!(seconds >= 0.0))
    return "--:--:--";
  long total = static_cast<long>(seconds + 0.5);
  char buffer[32];
  std::snprintf(buffer, sizeof(buffer), "%02ld:%02ld:%02ld",
                total / 3600, (total % 3600) / 60, total % 60);
  return std::string(buffer);
}

/// Thread safe progress reporter (percentage, elapsed time and ETA)
class ProgressReporter {
public:
  ProgressReporter(const std::string &label, long total, bool enabled)
      : m_label(label), m_total(total > 0 ? total : 1), m_done(0), m_enabled(enabled),
        m_start(std::chrono::steady_clock::now()),
        m_step(total > 200 ? total / 200 : 1) {
    if (m_enabled)
      print(0);
  }

  void inc(long n = 1) {
    long previous = m_done.fetch_add(n);
    long done = previous + n;
    if (!m_enabled)
      return;
    if (done >= m_total || previous / m_step != done / m_step)
      print(done);
  }

  void finish() {
    if (!m_enabled)
      return;
    print(m_done.load());
    std::fprintf(stderr, "\n");
    std::fflush(stderr);
    m_enabled = false;
  }

private:
  void print(long done) {
    if (done > m_total)
      done = m_total;
    double elapsed = std::chrono::duration<double>(
                         std::chrono::steady_clock::now() - m_start).count();
    double ratio = static_cast<double>(done) / static_cast<double>(m_total);
    double eta = done > 0 ? elapsed * (m_total - done) / done : 0.0;

    std::lock_guard<std::mutex> lock(m_mutex);
    std::fprintf(stderr, "\r%s: %6.2f%% (%ld/%ld) elapsed %s ETA %s     ",
                 m_label.c_str(), ratio * 100.0, done, m_total,
                 formatDuration(elapsed).c_str(), formatDuration(eta).c_str());
    std::fflush(stderr);
  }

  std::string m_label;
  long m_total;
  std::atomic<long> m_done;
  bool m_enabled;
  std::chrono::steady_clock::time_point m_start;
  long m_step;
  std::mutex m_mutex;
};

// ---------------------------------------------------------------------------
// Parallel execution
// ---------------------------------------------------------------------------

static void
parallelFor(size_t count, int nb_threads,
            const std::function<void(size_t)> &func) {
  if (nb_threads <= 1 || count <= 1) {
    for (size_t i = 0; i < count; ++i)
      func(i);
    return;
  }

  std::atomic<size_t> next(0);
  std::mutex error_mutex;
  std::string error_message;

  int threads_to_use = static_cast<int>(std::min<size_t>(nb_threads, count));
  std::vector<std::thread> threads;
  for (int t = 0; t < threads_to_use; ++t) {
    threads.push_back(std::thread([&]() {
      for (;;) {
        size_t i = next.fetch_add(1);
        if (i >= count)
          break;
        try {
          func(i);
        } catch (const std::exception &e) {
          std::lock_guard<std::mutex> lock(error_mutex);
          if (error_message.empty())
            error_message = e.what();
        } catch (...) {
          std::lock_guard<std::mutex> lock(error_mutex);
          if (error_message.empty())
            error_message = "unknown error";
        }
      }
    }));
  }
  for (size_t i = 0; i < threads.size(); ++i)
    threads[i].join();

  if (!error_message.empty())
    throw G2TException(error_message);
}

void
GDAL2Tiles::run() {
  openInput();
  generateMetadata();

  TileJobInfo tile_job_info;
  std::vector<TileDetail> tile_details;
  generateBaseTiles(tile_job_info, tile_details);

  int nb_threads = m_options.nb_processes < 1 ? 1 : m_options.nb_processes;
  bool show_progress = !m_options.quiet;

  {
    ProgressReporter progress("Generating Base Tiles",
                              static_cast<long>(tile_details.size()), show_progress);
    if (nb_threads == 1) {
      for (size_t i = 0; i < tile_details.size(); ++i) {
        createBaseTile(tile_job_info, tile_details[i]);
        progress.inc();
      }
    } else {
      parallelFor(tile_details.size(), nb_threads, [&](size_t i) {
        createBaseTile(tile_job_info, tile_details[i]);
        progress.inc();
      });
    }
    progress.finish();
  }

  int overview_count = countOverviewTiles(tile_job_info);
  if (overview_count > 0) {
    ProgressReporter progress("Generating Overview Tiles", overview_count,
                              show_progress);
    for (int base_tz = tile_job_info.tmaxz; base_tz > tile_job_info.tminz; --base_tz) {
      std::vector<std::vector<std::pair<int, int>>> groups =
          groupOverviewBaseTiles(base_tz, m_output_folder, tile_job_info);
      if (nb_threads == 1) {
        for (size_t i = 0; i < groups.size(); ++i) {
          createOverviewTile(base_tz, groups[i], m_output_folder, tile_job_info, m_options);
          progress.inc();
        }
      } else {
        parallelFor(groups.size(), nb_threads, [&](size_t i) {
          createOverviewTile(base_tz, groups[i], m_output_folder, tile_job_info, m_options);
          progress.inc();
        });
      }
    }
    progress.finish();
  }
}

// ---------------------------------------------------------------------------
// Options post processing
// ---------------------------------------------------------------------------

void
postProcessOptions(Options &options, const std::string &input_file,
                   const std::string &output_folder) {
  (void)output_folder;

  GDALAllRegister();

  if (options.title.empty()) {
    std::string::size_type pos = input_file.find_last_of("/\\");
    options.title = (pos == std::string::npos) ? input_file : input_file.substr(pos + 1);
  }

  options.tminz = -1;
  options.tmaxz = -1;
  if (!options.zoom_str.empty()) {
    std::string::size_type dash = options.zoom_str.find('-');
    if (dash == std::string::npos) {
      options.tminz = CPLAtof(options.zoom_str.c_str());
      options.tmaxz = options.tminz;
    } else {
      std::string min_str = options.zoom_str.substr(0, dash);
      std::string max_str = options.zoom_str.substr(dash + 1);
      options.tminz = CPLAtof(min_str.c_str());
      if (!max_str.empty()) {
        options.tmaxz = CPLAtof(max_str.c_str());
        if (options.tmaxz < options.tminz)
          throw G2TException("max zoom less than min zoom");
      }
    }
  }

  if (options.resampling == "antialias")
    options.resampling = "lanczos";

  if (options.tiledriver == "WEBP") {
    if (GetGDALDriverManager()->GetDriverByName("WEBP") == nullptr)
      throw G2TException("WEBP driver is not available");
    if (!options.webp_lossless &&
        (options.webp_quality <= 0 || options.webp_quality > 100))
      throw G2TException("webp_quality should be in the range [1-100]");
  } else if (options.tiledriver == "JPEG") {
    if (GetGDALDriverManager()->GetDriverByName("JPEG") == nullptr)
      throw G2TException("JPEG driver is not available");
    if (options.jpeg_quality <= 0 || options.jpeg_quality > 100)
      throw G2TException("jpeg_quality should be in the range [1-100]");
  }
}

} // namespace g2t
