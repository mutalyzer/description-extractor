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
//   Revision: 1.03a
//   Date:     2013/07/25
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

int main(int argc, char* argv[])
{
  static_cast<void>(argc);
  static_cast<void>(argv);
                        // 000000000111111111122
                        // 123456789012345678901
  char const* reference = "ATAGATGATAGATAGATAGAT";
  char const* sample    = "ATAGATGATTGATAGATAGAT";

  vector<Variant> result = extract(reference, 21, sample, 21);

  // simple HGVS description (for illustration only)
  for (size_t i = 0; i < result.size(); ++i)
  {
    if (result[i].reverse_complement)
    {
      printf("%ld_%ldinv\n", result[i].reference_start + 1, result[i].reference_end);
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

