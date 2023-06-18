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
  if (!fileStream.is_open()) {
    error("Could not open file");
  }
  fileStream << filename << "\n";
  fileStream.flush(); // "commits" current buffer to file (in case close() is not called)
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
