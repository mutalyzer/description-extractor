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
  size_t reference_length;
  size_t sample_length;
  ifstream file(argv[1]);

  fprintf(stderr, "debug.cc --- loading files\n  reference: %s\n  sample:    %s\n", argv[1], argv[2]);

  file.seekg(0, file.end);
  reference_length = file.tellg();
  file.seekg(0, file.beg);
  reference = new char[reference_length];
  file.read(reference, reference_length);
  file.close();
  file.open(argv[2]);
  file.seekg(0, file.end);
  sample_length = file.tellg();
  file.seekg(0, file.beg);
  sample = new char[sample_length];
  file.read(sample, sample_length);
  file.close();

  fprintf(stderr, "debug.cc --- files loaded\n  reference: %ld\n  sample:    %ld\n", reference_length, sample_length);

  vector<Variant> result = extract(reference, reference_length, sample, sample_length);

  delete[] reference;
  delete[] sample;

  return 0;
} // main

