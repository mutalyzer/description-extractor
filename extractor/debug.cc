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

  int count_idem = 0;
  int count_inv = 0;
  int count_trans = 0;
  int count_trans_inv = 0;
  int count_ins = 0;
  int count_del = 0;
  int count_sub = 0;
  int count_indel = 0;

  int count_length_del = 0;
  int count_length_ins = 0;

  // simple HGVS description (for illustration/debug only)
  for (size_t i = 0; i < result.size(); ++i)
  {
    if (result[i].type == VARIANT_IDENTITY)
    {
      ++count_idem;
      //printf("%ld_%ldidem%ld_%ld\n", result[i].reference_start + 1, result[i].reference_end, result[i].sample_start + 1, result[i].sample_end);
    } // if
    else if (result[i].type == VARIANT_REVERSE_COMPLEMENT)
    {
      ++count_inv;
      printf("%ld_%ldinv\n", result[i].reference_start + 1, result[i].reference_end);
    } // if
    else if (result[i].type == VARIANT_TRANSPOSITION)
    {
      ++count_trans;
      printf("%ld_%ldins%ld_%ld\n", result[i].reference_start + 1, result[i].reference_end, result[i].sample_start, result[i].sample_end);
    } // if
    else if (result[i].type == VARIANT_REVERSE_COMPLEMENT_TRANSPOSITION)
    {
      ++count_trans_inv;
      printf("%ld_%ldinso%ld_%ld\n", result[i].reference_start + 1, result[i].reference_end, result[i].sample_start, result[i].sample_end);
    } // if
    else
    {
      size_t const reference_length = result[i].reference_end - result[i].reference_start;
      size_t const sample_length = result[i].sample_end - result[i].sample_start;
      if (reference_length == 0)
      {
        count_length_ins += sample_length;
        ++count_ins;
        fprintf(stdout, "%ld_%ldins", result[i].reference_start, result[i].reference_start + 1);
        fwrite(sample + result[i].sample_start, sizeof(char), sample_length, stdout);
        fputc('\n', stdout);
      } // if
      else if (sample_length == 0)
      {
        count_length_del += reference_length;
        ++count_del;
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
        ++count_sub;
        fprintf(stdout, "%ld%c>%c\n", result[i].reference_start + 1, reference[result[i].reference_start], sample[result[i].sample_start]);
      } // if
      else
      {
        count_length_ins += sample_length;
        count_length_del += reference_length;
        ++count_indel;
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

  cout << endl;
  cout << "idem: " << count_idem << endl;
  cout << "inv: " << count_inv << endl;
  cout << "trans: " << count_trans << endl;
  cout << "trans_inv: " << count_trans_inv << endl;
  cout << "ins: " << count_ins << endl;
  cout << "del: " << count_del << endl;
  cout << "sub: " << count_sub << endl;
  cout << "indel: " << count_indel << endl;
  cout << endl;
  cout << "length_ins: " << count_length_ins << endl;
  cout << "length_del: " << count_length_del << endl;

  delete[] reference;
  delete[] sample;

  return 0;
} // main

