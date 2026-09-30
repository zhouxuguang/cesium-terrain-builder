/*******************************************************************************
 * A C++ port of the GDAL gdal2tiles.py tool.
 *
 * Tiles are generated following the same logic and command line interface as
 * the original Python utility, without the --mpi, KML, Web viewer and MapML
 * options.
 *******************************************************************************/

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "cpl_string.h"
#include "gdal.h"
#include "gdal_version.h"

#include "G2Tiler.hpp"

#ifdef _WIN32
#include <windows.h>
#include <shellapi.h>
#endif

using namespace std;

static const char *RESAMPLING_LIST[] = {
  "average", "near", "bilinear", "cubic", "cubicspline", "lanczos",
  "antialias", "mode", "max", "min", "med", "q1", "q3"
};

static const char *PROFILE_LIST[] = { "mercator", "geodetic", "raster" };

static bool
inList(const string &value, const char * const *list, size_t count) {
  for (size_t i = 0; i < count; ++i)
    if (value == list[i])
      return true;
  return false;
}

static void
printUsage(const char *name) {
  cout <<
    "Usage: " << name << " [options] input_file [output]\n"
    "\n"
    "Options:\n"
    "  -p, --profile <profile>       Tile cutting profile (mercator,geodetic,raster)\n"
    "                                - default 'mercator' (Google Maps compatible)\n"
    "  -r, --resampling <method>     Resampling method (average,near,bilinear,cubic,\n"
    "                                cubicspline,lanczos,antialias,mode,max,min,med,q1,q3)\n"
    "                                - default 'average'\n"
    "  -s, --s_srs <SRS>             The spatial reference system used for the source data\n"
    "  -z, --zoom <zooms>            Zoom levels to render (format:'2-5', '10-' or '10').\n"
    "  -e, --resume                  Resume mode. Generate only missing files.\n"
    "  -a, --srcnodata <value>       Value in the input dataset considered as transparent\n"
    "  -d, --tmscompatible           When using the geodetic profile, specifies the base\n"
    "                                resolution as 0.703125 or 2 tiles at zoom level 0.\n"
    "      --xyz                     Use XYZ tile numbering (OSM Slippy Map tiles)\n"
    "                                instead of TMS\n"
    "  -v, --verbose                 Print status messages to stdout\n"
    "  -x, --exclude                 Exclude transparent tiles from result tileset\n"
    "  -q, --quiet                   Disable messages and status to stdout\n"
    "      --processes <n>           Number of processes to use for tiling\n"
    "      --tilesize <pixels>       Width and height in pixel of a tile\n"
    "      --tiledriver <driver>     Which tile driver to use for the tiles (PNG,WEBP,JPEG)\n"
    "      --excluded-values <v>     Tuples of values that must be ignored as contributing\n"
    "                                source pixels during resampling\n"
    "      --excluded-values-pct-threshold <pct>\n"
    "                                Minimum percentage of source pixels that must be set at\n"
    "                                one of the --excluded-values (default 50)\n"
    "      --nodata-values-pct-threshold <pct>\n"
    "                                Minimum percentage of source pixels that must be at\n"
    "                                nodata to cause transparent output (default 100)\n"
    "\n"
    "WEBP options:\n"
    "      --webp-quality <q>        quality of webp image, integer between 1 and 100\n"
    "                                (default 75)\n"
    "      --webp-lossless           use lossless compression for the webp image\n"
    "\n"
    "JPEG options:\n"
    "      --jpeg-quality <q>        quality of jpeg image, integer between 1 and 100\n"
    "                                (default 75)\n"
    "\n"
    "  -h, --help                    Show this help message and exit\n"
    "      --version                 Show the program version and exit\n";
}

/// Fetch the value of an option, supporting both '--opt value' and '--opt=value'
static string
optionValue(int argc, char *argv[], int &i, const char *arg, const char *opt_name) {
  string value;
  const char *equals = strchr(arg, '=');
  if (equals != nullptr) {
    value = equals + 1;
  } else {
    if (i + 1 >= argc) {
      cerr << "  Error: option '" << opt_name << "' requires an argument" << endl;
      exit(1);
    }
    value = argv[++i];
  }
  return value;
}

static int
runMain(vector<string> args) {
  vector<char *> argv_storage;
  for (size_t i = 0; i < args.size(); ++i)
    argv_storage.push_back(&args[i][0]);
  int argc = static_cast<int>(argv_storage.size());
  char **argv = argv_storage.empty() ? nullptr : argv_storage.data();

  g2t::Options options;
  vector<string> positional;

  for (int i = 1; i < argc; ++i) {
    string arg = argv[i];

    if (arg == "-h" || arg == "--help") {
      printUsage(argv[0]);
      return 0;
    } else if (arg == "--version") {
      cout << argv[0] << " " << GDALVersionInfo("RELEASE_NAME") << endl;
      return 0;
    } else if (arg == "-p" || arg.compare(0, 10, "--profile=") == 0 || arg == "--profile") {
      options.profile = optionValue(argc, argv, i, arg.c_str(), "--profile");
    } else if (arg == "-r" || arg.compare(0, 13, "--resampling=") == 0 || arg == "--resampling") {
      options.resampling = optionValue(argc, argv, i, arg.c_str(), "--resampling");
    } else if (arg == "-s" || arg.compare(0, 8, "--s_srs=") == 0 || arg == "--s_srs") {
      options.s_srs = optionValue(argc, argv, i, arg.c_str(), "--s_srs");
    } else if (arg == "-z" || arg.compare(0, 7, "--zoom=") == 0 || arg == "--zoom") {
      options.zoom_str = optionValue(argc, argv, i, arg.c_str(), "--zoom");
    } else if (arg == "-a" || arg.compare(0, 12, "--srcnodata=") == 0 || arg == "--srcnodata") {
      options.srcnodata = optionValue(argc, argv, i, arg.c_str(), "--srcnodata");
    } else if (arg == "--processes" || arg.compare(0, 12, "--processes=") == 0) {
      options.nb_processes = atoi(optionValue(argc, argv, i, arg.c_str(), "--processes").c_str());
    } else if (arg == "--tilesize" || arg.compare(0, 11, "--tilesize=") == 0) {
      options.tilesize = atoi(optionValue(argc, argv, i, arg.c_str(), "--tilesize").c_str());
    } else if (arg == "--tiledriver" || arg.compare(0, 13, "--tiledriver=") == 0) {
      options.tiledriver = optionValue(argc, argv, i, arg.c_str(), "--tiledriver");
    } else if (arg == "--excluded-values" || arg.compare(0, 18, "--excluded-values=") == 0) {
      options.excluded_values = optionValue(argc, argv, i, arg.c_str(), "--excluded-values");
    } else if (arg == "--excluded-values-pct-threshold" ||
               arg.compare(0, 32, "--excluded-values-pct-threshold=") == 0) {
      options.excluded_values_pct_threshold = atof(
          optionValue(argc, argv, i, arg.c_str(), "--excluded-values-pct-threshold").c_str());
    } else if (arg == "--nodata-values-pct-threshold" ||
               arg.compare(0, 30, "--nodata-values-pct-threshold=") == 0) {
      options.nodata_values_pct_threshold = atof(
          optionValue(argc, argv, i, arg.c_str(), "--nodata-values-pct-threshold").c_str());
    } else if (arg == "--webp-quality" || arg.compare(0, 15, "--webp-quality=") == 0) {
      options.webp_quality = atoi(
          optionValue(argc, argv, i, arg.c_str(), "--webp-quality").c_str());
    } else if (arg == "--webp-lossless") {
      options.webp_lossless = true;
    } else if (arg == "--jpeg-quality" || arg.compare(0, 15, "--jpeg-quality=") == 0) {
      options.jpeg_quality = atoi(
          optionValue(argc, argv, i, arg.c_str(), "--jpeg-quality").c_str());
    } else if (arg == "-e" || arg == "--resume") {
      options.resume = true;
    } else if (arg == "-d" || arg == "--tmscompatible") {
      options.tmscompatible = true;
    } else if (arg == "--xyz") {
      options.xyz = true;
    } else if (arg == "-v" || arg == "--verbose") {
      options.verbose = true;
    } else if (arg == "-x" || arg == "--exclude") {
      options.exclude_transparent = true;
    } else if (arg == "-q" || arg == "--quiet") {
      options.quiet = true;
    } else if (!arg.empty() && arg[0] == '-' && arg != "-") {
      cerr << "  Error: unknown option '" << arg << "'" << endl;
      printUsage(argv[0]);
      return 1;
    } else {
      positional.push_back(arg);
    }
  }

  if (positional.empty()) {
    cerr << "  Error: You need to specify at least an input file as argument" << endl;
    printUsage(argv[0]);
    return 1;
  }
  if (positional.size() > 2) {
    cerr << "  Error: Processing of several input files is not supported." << endl;
    return 1;
  }

  string input_file = positional[0];
  string output_folder;
  if (positional.size() == 2) {
    output_folder = positional[1];
  } else {
    string::size_type slash = input_file.find_last_of("/\\");
    string base = (slash == string::npos) ? input_file : input_file.substr(slash + 1);
    string::size_type dot = base.find_last_of('.');
    output_folder = (dot == string::npos) ? base : base.substr(0, dot);
  }

  try {
    GDALAllRegister();

    if (!inList(options.profile, PROFILE_LIST,
                sizeof(PROFILE_LIST) / sizeof(PROFILE_LIST[0]))) {
      cerr << "  Error: Invalid value for --profile option: " << options.profile << endl;
      return 1;
    }
    if (!inList(options.resampling, RESAMPLING_LIST,
                sizeof(RESAMPLING_LIST) / sizeof(RESAMPLING_LIST[0]))) {
      cerr << "  Error: Invalid value for --resampling option: " << options.resampling << endl;
      return 1;
    }
    if (options.tiledriver != "PNG" && options.tiledriver != "WEBP" &&
        options.tiledriver != "JPEG") {
      cerr << "  Error: Invalid value for --tiledriver option: " << options.tiledriver << endl;
      return 1;
    }

    g2t::postProcessOptions(options, input_file, output_folder);

    g2t::GDAL2Tiles tiler(input_file, output_folder, options);
    tiler.run();
  } catch (const std::exception &e) {
    cerr << "ERROR: " << e.what() << endl;
    return 1;
  }

  return 0;
}

#ifdef _WIN32
static string
wideToUtf8(const wchar_t *wide) {
  if (wide == nullptr)
    return "";
  int size = WideCharToMultiByte(CP_UTF8, 0, wide, -1, nullptr, 0, nullptr, nullptr);
  if (size <= 0)
    return "";
  string result(size - 1, '\0');
  WideCharToMultiByte(CP_UTF8, 0, wide, -1, &result[0], size, nullptr, nullptr);
  return result;
}

int
main() {
  int wide_argc = 0;
  LPWSTR *wide_argv = CommandLineToArgvW(GetCommandLineW(), &wide_argc);
  vector<string> args;
  if (wide_argv != nullptr) {
    for (int i = 0; i < wide_argc; ++i)
      args.push_back(wideToUtf8(wide_argv[i]));
    LocalFree(wide_argv);
  }
  if (args.empty())
    args.push_back("gdal2tiles");
  return runMain(args);
}
#else
int
main(int argc, char *argv[]) {
  vector<string> args;
  for (int i = 0; i < argc; ++i)
    args.push_back(argv[i]);
  return runMain(args);
}
#endif
