// *******************************************************************
// Extractor (library)
// *******************************************************************
// FILE INFORMATION:
//   File:     debug.cc
//   Author:   Jonathan K. Vis
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

/*
  int const N = 359;

  string header[N];
  string protein[N];

  fprintf(stdout, "\t");
  for (int i = 0; i < N; ++i)
  {
    cin >> header[i] >> protein[i];
    fprintf(stdout, "%s\t", header[i].c_str());
  } // for
  fprintf(stdout, "\n");

  for (int i = 0; i < N; ++i)
  {
    cerr << i << endl;
    fprintf(stdout, "%s\t", header[i].c_str());
    for (int j = 0; j < N; ++j)
    {
      if (i == j)
      {
        fprintf(stdout, "%.10e\t", 0.);
        continue;
      } // if
      vector<Variant> variants;
      extract(variants, protein[i].c_str(), protein[i].length(), protein[j].c_str(), protein[j].length(), TYPE_PROTEIN, "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF");

      size_t best = 0;
      for (std::vector<Variant>::iterator it = variants.begin(); it != variants.end(); ++it)
      {
        if (it->type >= FRAME_SHIFT && it->reference_end - it->reference_start > best)
        {
          best = it->reference_end - it->reference_start;
        } // if
      } // for
      size_t const length = min(protein[i].length(), protein[j].length());
      fprintf(stdout, "%.10e\t", 1. - static_cast<double>(best) / static_cast<double>(length));
    } // for
    fprintf(stdout, "\n");
  } // for
*/

  // The actual extraction.
  std::vector<Variant> variant;
  size_t const weight = extract(variant, reference, reference_length - 1, sample, sample_length - 1, TYPE_DNA);


  // Printing the variants.
  fprintf(stdout, "Variants (%ld / %ld):\n", variant.size(), weight);
  for (std::vector<Variant>::iterator it = variant.begin(); it != variant.end(); ++it)
  {
    if (it->type >= FRAME_SHIFT)
    {
      fprintf(stdout, "%ld--%ld, %ld--%ld, %d, %lf, %ld--%ld\n", it->reference_start, it->reference_end, it->sample_start, it->sample_end, it->type, 1.f - it->probability, it->transposition_start, it->transposition_end);
      char_t ref_DNA[(it->reference_end - it->reference_start) * 3];
      char_t alt_DNA[(it->reference_end - it->reference_start) * 3];
      backtranslation(ref_DNA, alt_DNA, reference, it->reference_start, sample, it->sample_start, it->reference_end - it->reference_start, it->type);
      fprintf(stdout, "ref_DNA: ");
      fwrite(ref_DNA, sizeof(char_t), (it->reference_end - it->reference_start) * 3, stdout);
      fprintf(stdout, "\nref_pro: ");
      fwrite(reference + it->reference_start, sizeof(char_t), (it->reference_end - it->reference_start), stdout);
      fprintf(stdout , "\nalt_DNA: ");
      fwrite(alt_DNA, sizeof(char_t), (it->reference_end - it->reference_start) * 3, stdout);
      fprintf(stdout, "\nalt_pro: ");
      fwrite(sample + it->sample_start, sizeof(char_t), (it->reference_end - it->reference_start), stdout);
      fprintf(stdout , "\n");
    } // if
    else
    {
      fprintf(stdout, "%ld--%ld, %ld--%ld, %d, %ld, %ld--%ld\n", it->reference_start, it->reference_end, it->sample_start, it->sample_end, it->type, it->weight, it->transposition_start, it->transposition_end);
    } // else
  } // for


  // Cleaning up.
  delete[] reference;
  delete[] sample;

  return 0;
} // main

