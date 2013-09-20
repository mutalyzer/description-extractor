// *******************************************************************
//   (C) Copyright 2013 Leiden Institute of Advanced Computer Science
//   Universiteit Leiden
//   All Rights Reserved
// *******************************************************************
// Extractor (library)
// *******************************************************************
// FILE INFORMATION:
//   File:     debug.cc (depends on extractor.h)
//   Author:   Jonathan K. Vis
//   Revision: 1.04b
//   Date:     2013/09/20
// *******************************************************************
// DESCRIPTION:
//   This source can be used to debug the Extractor library within
//   C/C++.
// *******************************************************************

#include <cstdio>
#include <cstdlib>
#include <vector>
using namespace std;

#include "extractor.h"
using namespace mutalyzer;

static size_t const LENGTH = 16;

static char const* const REFERENCE = "AAAAAAAAAAAAAAAA";

static char const* const SAMPLE[] =
{
  "AAAAAAAAAAAAAAAA",
  "GGGGGGGGGGGGGGGG",
  "CAAAAAAAAAAAAAAA",
  "AAAAAAAAAAAAAAAC",
  "CAAAAAAAAAAAAAAC",
  "TTTTTTTTTTTTTTTT",
  "CTTTTTTTTTTTTTTT",
  "TTTTTTTTTTTTTTTC",
  "CTTTTTTTTTTTTTTC",
  "AAAAAAACCAAAAAAA",
  "AAACCAAAAAACCAAA",
  "AAAAAAATTAAAAAAA",
  "TTAAAAAAAAAAAATT",
  

}; // SAMPLE


int main(int argc, char* argv[])
{
  char const* const reference = REFERENCE;
  char const* const sample = argc < 2 ? SAMPLE[0] : SAMPLE[atoi(argv[1])];

  vector<Variant> result = extract(reference, LENGTH, sample, LENGTH);

  // simple HGVS description (for illustration/debug only)
  for (size_t i = 0; i < result.size(); ++i)
  {
    if (result[i].type == VARIANT_IDENTITY)
    {
      printf("%ld_%ldidem%ld_%ld\n", result[i].reference_start + 1, result[i].reference_end, result[i].sample_start + 1, result[i].sample_end);
    } // if
    else if (result[i].type == VARIANT_REVERSE_COMPLEMENT)
    {
      printf("%ld_%ldinv\n", result[i].reference_start + 1, result[i].reference_end);
    } // if
    else if (result[i].type == VARIANT_TRANSPOSITION || result[i].type == VARIANT_REVERSE_COMPLEMENT_TRANSPOSITION)
    {
      printf("%ld_%ldins%ld_%ld\n", result[i].reference_start + 1, result[i].reference_end, result[i].sample_start, result[i].sample_end);
    } // if
    else
    {
      size_t const reference_length = result[i].reference_end - result[i].reference_start;
      size_t const sample_length = result[i].sample_end - result[i].sample_start;
      if (reference_length == 0)
      {
        fprintf(stdout, "%ld_%ldins", result[i].reference_start, result[i].reference_start + 1);
        fwrite(sample + result[i].sample_start, sizeof(char), sample_length, stdout);
        fputc('\n', stdout);
      } // if
      else if (sample_length == 0)
      {
        if (reference_length == 1)
        {
          fprintf(stdout, "%lddel\n", result[i].reference_start + 1);
        } // if
        else
        {
          fprintf(stdout, "%ld_%lddel\n", result[i].reference_start + 1, result[i].reference_end);
        } // else
      } // if
      else if (reference_length == 1 && sample_length == 1)
      {
        fprintf(stdout, "%ld%c>%c\n", result[i].reference_start + 1, reference[result[i].reference_start], sample[result[i].sample_start]);
      } // if
      else
      {
        if (reference_length == 1)
        {
          fprintf(stdout, "%lddelins", result[i].reference_start + 1);
        } // if
        else
        {
          fprintf(stdout, "%ld_%lddelins", result[i].reference_start + 1, result[i].reference_end);
        } // else
        fwrite(sample + result[i].sample_start, sizeof(char), sample_length, stdout);
        fputc('\n', stdout);
      } // else
    } // else
  } // for
  return 0;
} // main

