// *******************************************************************
//   (C) Copyright 2015 Leiden Institute of Advanced Computer Science
//   Universiteit Leiden
//   All Rights Reserved
// *******************************************************************
// Extractor (library)
// *******************************************************************
// FILE INFORMATION:
//   File:     debug.cc
//   Author:   Jonathan K. Vis
//   Revision: 2.3.0
//   Date:     2015/07/31
// *******************************************************************
// DESCRIPTION:
//   This source can be used to debug the Extractor library within
//   C/C++. It opens two files given as arguments and perform the
//   description extraction. Supply the -D__debug__ flag in the
//   Makefile for tracing.
// *******************************************************************

#include "extractor.h"
using namespace mutalyzer;

#include <cstdio>


#include <iostream>
#include <string>
using namespace std;


// Entry point.
int main(int argc, char* argv[])
{
/*
  if (argc < 3)
  {
    fprintf(stderr, "usage: %s reference sample\n", argv[0]);
    return 1;
  } // if
  fprintf(stderr, "HGVS description extractor\n");


  // Opening files.
  FILE* file = fopen(argv[1], "r");
  if (file == 0)
  {
    fprintf(stderr, "ERROR: could not open file `%s'\n", argv[1]);
    return 1;
  } // if
  fseek(file, 0, SEEK_END);
  size_t const reference_length = ftell(file);
  rewind(file);
  char_t* reference = new char_t[reference_length];
  size_t const ref_length = fread(reference, sizeof(char_t), reference_length, file);
  static_cast<void>(ref_length);
  fclose(file);

  file = fopen(argv[2], "r");
  if (file == 0)
  {
    fprintf(stderr, "ERROR: could not open file `%s'\n", argv[2]);
    delete[] reference;
    return 1;
  } // if
  fseek(file, 0, SEEK_END);
  size_t const sample_length = ftell(file);
  rewind(file);
  char_t* sample = new char_t[sample_length];
  size_t const alt_length = fread(sample, sizeof(char_t), sample_length, file);
  static_cast<void>(alt_length);
  fclose(file);
*/

  string header[4305];
  string protein[4305];

  for (int i = 0; i < 4305; ++i)
  {
    getline(cin, header[i]);
    getline(cin, protein[i]);
  } // for

  for (int i = 0; i < 4305; ++i)
  {
    cerr << i << endl;
    double best = 1.f;
    for (int j = i + 1; j < 4305; ++j)
    {
      vector<Variant> variant;
      extract(variant, protein[i].c_str(), protein[i].length(), protein[j].c_str(), protein[j].length(), TYPE_PROTEIN, "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF");

      for (std::vector<Variant>::iterator it = variant.begin(); it != variant.end(); ++it)
      {
        if (it->type >= FRAME_SHIFT && best > it->probability)
        {
          best = it->probability;
          fprintf(stdout, "%.9s %.9s %ld--%ld, %ld--%ld, %d, %.10e\n", header[i].c_str(), header[j].c_str(), it->reference_start, it->reference_end, it->sample_start, it->sample_end, it->type, 1.f - it->probability);
        } // if
      } // for
    } // for
    fprintf(stdout, "\n");
  } // for

/*
  // The actual extraction.
  std::vector<Variant> variant;
  size_t const weight = extract(variant, reference, reference_length - 1, sample, sample_length - 1, TYPE_PROTEIN, "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF");


  // Printing the variants.
  fprintf(stdout, "Variants (%ld / %ld):\n", variant.size(), weight);
  for (std::vector<Variant>::iterator it = variant.begin(); it != variant.end(); ++it)
  {
    if (it->type >= FRAME_SHIFT)
    {
      fprintf(stdout, "%ld--%ld, %ld--%ld, %d, %lf, %ld--%ld\n", it->reference_start, it->reference_end, it->sample_start, it->sample_end, it->type, 1.f - it->probability, it->transposition_start, it->transposition_end);
    } // if
    else
    {
      fprintf(stdout, "%ld--%ld, %ld--%ld, %d, %ld, %ld--%ld\n", it->reference_start, it->reference_end, it->sample_start, it->sample_end, it->type, it->weight, it->transposition_start, it->transposition_end);
    } // else
  } // for


  // Cleaning up.
  delete[] reference;
  delete[] sample;
*/
  return 0;
} // main

