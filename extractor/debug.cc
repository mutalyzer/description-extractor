// *******************************************************************
//   (C) Copyright 2014 Leiden Institute of Advanced Computer Science
//   Universiteit Leiden
//   All Rights Reserved
// *******************************************************************
// Extractor (library)
// *******************************************************************
// FILE INFORMATION:
//   File:     debug.cc
//   Author:   Jonathan K. Vis
//   Revision: 2.1.0
//   Date:     2014/08/06
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

// Entry point.
int main(int argc, char* argv[])
{

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
  fread(reference, sizeof(char_t), reference_length, file);
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
  fread(sample, sizeof(char_t), sample_length, file);
  fclose(file);


  // The actual extraction.
  std::vector<Variant> variant;
  size_t const weight = extract(variant, reference, reference_length, sample, sample_length);


  // Printing the variants.
  fprintf(stdout, "\nVariants (%ld / %ld):\n", variant.size(), weight);
  for (std::vector<Variant>::iterator it = variant.begin() ; it != variant.end(); ++it)
  {
    //if (it->type != IDENTITY)
    {
      fprintf(stdout, "%ld--%ld, %ld--%ld, %d, %ld, %ld--%ld\n", it->reference_start, it->reference_end, it->sample_start, it->sample_end, it->type, it->weight, it->transposition_start, it->transposition_end);
    } // if
  } // for


  // Cleaning up.
  delete[] reference;
  delete[] sample;

  return 0;
} // main

