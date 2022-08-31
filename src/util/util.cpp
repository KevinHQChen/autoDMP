#include "util/util.hpp"

// sample commands:
// For SYSID:   .\robodrop.exe --vidSrc 2 --templateSource 0 --chanSource 0
// --tmThres 0.8 --ctrl 2 --dryRun 1 For control: .\robodrop.exe --vidSrc 2
// --templateSource 0 --chanSource 0 --tmThres 0.8 --ctrl 1 --dryRun 1 (replace
// with REPL if we get the chance)
config::config() {
  cli.set_config("--config", "config.toml", "Read config file", true);

  cli.add_option("-s,--videoSource", videoSource,
                 "set video source (0:Webcam, 1:From file, 2:Andor)");
  // app.set_help_all_flag("--help-all", "Expand all help");
  // app.add_flag("--version", "Get version");

  // cameraCli = cli.add_subcommand("cam", "Camera control");
  // cameraCli->add_option("-s,--videoSource", videoSource,
  //                       "set video source (0:Webcam, 1:From file, 2:Andor)");

  // cli.

  //     CLI::App *cameraApp =
  //     app.add_subcommand("cam", "Configure the app camera");
  // cameraApp->require_subcommand(0, 1); // 0 (default) or 1 camera

  // std::string mvcamera_config_file = "mvcamera_config.json";
  // CLI::App *mvcameraApp = cameraApp->add_subcommand(
  //     "mvcamera", "MatrixVision Camera Configuration");
  // mvcameraApp
  //     ->add_option("-c,--config", mvcamera_config_file, "Config filename")
  //     ->capture_default_str()
  //     ->check(CLI::ExistingFile);

  // std::string mock_camera_path;
  // CLI::App *mockcameraApp =
  //     cameraApp->add_subcommand("mock", "Mock Camera Configuration");
  // mockcameraApp->add_option("-p,--path", mock_camera_path, "Path")
  //     ->required()
  //     ->check(CLI::ExistingPath);

  // bool show_version = false;
  // cli.add_flag("--version", show_version, "Show version information");

  // cameraCli->require_subcommand(1);

  // cameraCli->add_option("-t,--templateSource", templateSource, "set
  // template source (0:Webcam, 1:From file, 2:Andor)");
  // cameraCli->add_option("-c,--chanSource", chanSource, "set channel source
  // (0:Webcam, 1:From file, 2:Andor)"); cameraCli->add_option("-T,--tmThres",
  // tmThres, "set template matching threshold");
  // cameraCli->add_option("-C,--ctrl", ctrl, "set control mode (0:SYSID,
  // 1:control)"); cameraCli->add_option("-D,--dryRun", dryRun, "set dry run
  // mode (0:off, 1:on)"); cameraCli->add_option("-V,--vidSrc", vidSrc, "set
  // video source (0:Webcam, 1:From file, 2:Andor)");
  // cameraCli->add_option("-O,--outFile", outFile, "set output file");
  // cameraCli->add_option("-F,--frameRate", frameRate, "set frame rate");

  // namespace po = boost::program_options;
  // try {
  //     po::options_description desc("Allowed options");
  //     desc.add_options()
  //         ("help,H", "produce help message")
  //         ("vidSrc", po::value<int>(&videoSource)->default_value(2),
  //             "set video source (0:WebCam, 1:From File, 2:Andor)")
  //         ("vidSrcFn",
  //         po::value<std::string>(&vidSrcFn)->default_value("C:/Users/khqc/RoboDrop/rawVideo2.avi"),
  //             "set video source filename")
  //         ("templateSource",
  //         po::value<int>(&templateSource)->default_value(0),
  //             "set template source (0:From File, 1:From Video Source)")
  //         ("chanSource", po::value<int>(&chanSource)->default_value(0),
  //             "set channel source (0:From File, 1:From Video Source)")
  //         ("tmThres", po::value<double>(&tmThres)->default_value(0.75),
  //             "set template matching threshold")
  //         ("saveRaw", po::value<int>(&saveRaw)->default_value(1),
  //             "save raw video (0:No, 1:Yes)")
  //         ("rawFPS", po::value<double>(&rawFPS)->default_value(40),
  //             "set raw FPS")
  //         ("saveProc", po::value<int>(&saveProc)->default_value(1),
  //             "save proc video (0:No, 1:Yes)")
  //         ("procFPS", po::value<double>(&procFPS)->default_value(40),
  //             "set proc FPS")
  //         ("imProc", po::value<int>(&imProc)->default_value(1),
  //             "enable image processing (0:No, 1:Yes)")
  //         ("ctrl", po::value<int>(&ctrl)->default_value(0),
  //             "engage controller (0:No, 1:Yes, 2:perform open loop sysID)")
  //         ("useGPU", po::value<int>(&useGPU)->default_value(0),
  //             "use GPU (0:No, 1:Yes)")
  //         ("numChannels", po::value<int>(&numChannels)->default_value(3),
  //         "set number of observable channels")
  //         ("dryRun", po::value<int>(&dryRun)->default_value(0), "dry Run
  //         (0:No, 1:Yes)")
  //     ;
  //     po::variables_map config_map;
  //     po::store(po::parse_command_line(argc, argv, desc), config_map);
  //     po::notify(config_map);

  //     if (config_map.count("help")) {
  //         std::cout << desc << "\n";
  //         paramsValid = 0;
  //     } else
  //         paramsValid = 1;
  // } catch(std::exception& e) {
  //     std::cout << "error: " << e.what() << "\n";
  //     paramsValid = 0;
  // } catch(...) {
  //     std::cout << "Exception of unknown type!\n";
  //     paramsValid = 0;
  // }
}

int config::parse(int argc, const char **argv) { CLI11_PARSE(cli, argc, argv); }

ffstream::ffstream(std::string filename) {
  // logFile.open("log1.txt");
  // add time stamp + put in folder

  // QDateTime now = QDateTime::currentDateTime();

  // In Windows (and MSVC) wchar_t is preferred over char, so it's a good
  // idea to use wide strings everywhere std::wstring filename =
  // L"ueva_log.txt";
  // filename.append(now.toString("yyyy_MM_dd_HH_mm_ss"));
  // filename.append(".txt");
  fileStream.open(filename);
  fileStream << filename << "\n";
}

ffstream::~ffstream() { fileStream.close(); }

mstream::mstream(std::string filename) {
  // logFile.open("log1.txt");
  // add time stamp + put in folder

  // QDateTime now = QDateTime::currentDateTime();

  multiStream.open(filename);

  multiStream << "====================================\n";
  multiStream << filename << "\n";
  multiStream << "====================================\n";
}

mstream::~mstream() { multiStream.close(); }

void saveData(std::string fileName, Eigen::MatrixXd matrix) {
  // https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html
  const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");

  std::ofstream file(fileName);
  if (file.is_open()) {
    file << matrix.format(CSVFormat);
    file.close();
  }
}

Eigen::MatrixXd openData(std::string fileToOpen) {

  // the inspiration for creating this function was drawn from here (I did
  // NOT copy and paste the code)
  // https://stackoverflow.com/questions/34247057/how-to-read-csv-file-and-assign-to-eigen-matrix

  // the input is the file: "fileToOpen.csv":
  // a,b,c
  // d,e,f
  // This function converts input file data into the Eigen matrix format

  // the matrix entries are stored in this variable row-wise. For example
  // if we have the matrix: M=[a b c
  //    d e f]
  // the entries are stored as matrixEntries=[a,b,c,d,e,f], that is the
  // variable "matrixEntries" is a row vector later on, this vector is
  // mapped into the Eigen matrix format
  std::vector<double> matrixEntries;

  // in this object we store the data from the matrix
  std::ifstream matrixDataFile(fileToOpen);

  // this variable is used to store the row of the matrix that contains
  // commas
  std::string matrixRowString;

  // this variable is used to store the matrix entry;
  std::string matrixEntry;

  // this variable is used to track the number of rows
  int matrixRowNumber = 0;

  while (std::getline(matrixDataFile,
                      matrixRowString)) // here we read a row by row of
                                        // matrixDataFile and store every line into
                                        // the string variable matrixRowString
  {
    std::stringstream matrixRowStringStream(matrixRowString); // convert matrixRowString that is a
                                                              // string to a stream variable.

    while (std::getline(matrixRowStringStream, matrixEntry,
                        ',')) // here we read pieces of the stream
                              // matrixRowStringStream until every comma, and store
                              // the resulting character into the matrixEntry
    {
      matrixEntries.push_back(std::stod(matrixEntry)); // here we convert the string to double
                                                       // and fill in the row vector storing
                                                       // all the matrix entries
    }
    matrixRowNumber++; // update the column numbers
  }

  // here we convet the vector variable into the matrix and return the
  // resulting object, note that matrixEntries.data() is the pointer to
  // the first memory location at which the entries of the vector
  // matrixEntries are stored;
  return Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
      matrixEntries.data(), matrixRowNumber, matrixEntries.size() / matrixRowNumber);
}
