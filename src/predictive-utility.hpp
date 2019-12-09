#ifndef __PREDICTIVE_UTILITY_HPP__
#define __PREDICTIVE_UTILITY_HPP__

#include <fstream>
#include <random>
#include <chrono>

#include <iostream>
using std::cout;
using std::endl;

#include <vector>
using std::vector;
using std::pair;
using std::make_pair;

#include <cmath>

#include "ArgParse.hpp"

template<typename S, typename T> inline bool print_to_csv(const vector<pair<S, T> >& data, const string& fileName) {
  if (data.empty()) return true;
  std::ofstream fout(fileName);
  if (fout.fail()) {
    cout << "Could not open file [" << fileName << "].\n";
    return false;
  }
  for (const auto& pr : data) fout << pr.first << "," << pr.second << "\n";
  fout.close();
  return true;
}

template<typename T> inline bool print_to_csv(const vector<T>& data, const string& fileName) {
  if (data.empty()) return true;
  std::ofstream fout(fileName);
  if (fout.fail()) {
    cout << "Could not open file [" << fileName << "].\n";
    return false;
  }
  for (int i=0; i<data.size(); ++i) {
    fout << data[i];
    if (i!=data.size()-1) fout << ",";
  }
  fout.close();
  return true;
}

template<typename T> inline bool print_to_csv(const vector<vector<T> >& data, const string& fileName) {
  if (data.empty()) return true;
  std::ofstream fout(fileName);
  if (fout.fail()) {
    cout << "Could not open file [" << fileName << "].\n";
    return false;
  }
  for (const auto& entry : data) {
    for (int i=0; i<entry.size(); ++i) {
      fout << entry[i];
      if (i!=entry.size()-1) fout << ",";
    }
    fout << "\n";
  }
  fout.close();
  return true;
}

#endif // __PREDICTIVE_UTILITY_HPP__