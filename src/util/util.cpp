#include "util/util.hpp"

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
