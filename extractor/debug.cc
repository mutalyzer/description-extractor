// *******************************************************************
//   (C) Copyright 2014 Leiden Institute of Advanced Computer Science
//   Universiteit Leiden
//   All Rights Reserved
// *******************************************************************
// Extractor (library)
// *******************************************************************
// FILE INFORMATION:
//   File:     debug.cc (depends on extractor.h)
//   Author:   Jonathan K. Vis
//   Revision: 1.05a
//   Date:     2014/02/05
// *******************************************************************
// DESCRIPTION:
//   This source can be used to debug the Extractor library within
//   C/C++.
// *******************************************************************

#include "extractor.h"
using namespace mutalyzer;

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;

int main(int argc, char* argv[])
{
  static_cast<void>(argc);

  char* reference;
  char* sample;
  ifstream file(argv[1]);

  fprintf(stderr, "debug.cc --- loading files\n  reference: %s\n  sample:    %s\n", argv[1], argv[2]);

  file.seekg(0, file.end);
  size_t const reference_length = static_cast<size_t>(file.tellg()) - 1;
  file.seekg(0, file.beg);
  reference = new char[reference_length];
  file.read(reference, reference_length);
  file.close();
  file.open(argv[2]);
  file.seekg(0, file.end);
  size_t const sample_length = static_cast<size_t>(file.tellg()) - 1;
  file.seekg(0, file.beg);
  sample = new char[sample_length];
  file.read(sample, sample_length);
  file.close();

  fprintf(stderr, "debug.cc --- files loaded\n  reference: %ld\n  sample:    %ld\n", reference_length, sample_length);

  vector<Variant> result = extract(reference, reference_length, sample, sample_length);

  for (size_t i = 0; i < result.size(); ++i)
  {
    cout << result[i].reference_start << "--" << result[i].reference_end << ", "
         << result[i].sample_start << "--" << result[i].sample_end << ", " << result[i].type << endl;
  } // for


  delete[] reference;
  delete[] sample;

  return 0;
} // main

