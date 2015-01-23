// *******************************************************************
//   (C) Copyright 2015 Leiden Institute of Advanced Computer Science
//   Universiteit Leiden
//   All Rights Reserved
// *******************************************************************
// Sandbox (testing playground)
// *******************************************************************
// FILE INFORMATION:
//   File:     sandbox.cc
//   Author:   Jonathan K. Vis
//   Revision: 1.1.2
//   Date:     2015/01/23
// *******************************************************************
// DESCRIPTION:
//  General testing playground: automatic tandem repeat detection
//  based on run length encoding with variable run lengths.
// *******************************************************************

#include <cstdio>
#include <cstdlib>

typedef char char_t;

bool string_match(char_t const* const string_1,
                  char_t const* const string_2,
                  size_t const        length)
{
  for (size_t i = 0; i < length; ++i)
  {
    if (string_1[i] != string_2[i])
    {
      return false;
    } // if
  } // for
  return true;
} // string_match

void tandem_repeat_annotation(char_t const* const string,
                              size_t const        start,
                              size_t const        end)
{
  size_t const length = end - start;
  size_t i = 0;
  while (i < length)
  {
    size_t max_count = 0;
    size_t max_k = 1;
    for (size_t k = 1; k < length / 2 + 1; ++k)
    {
      size_t count = 0;
      for (size_t j = i + k; j < length - k + 1; j += k)
      {
        if (!string_match(string + start + i, string + start + j, k))
        {
          break;
        } // if
        ++count;
      } // for
      if (count > 0 && count >= max_count)
      {
        max_count = count;
        max_k = k;
      } // if
    } // for
    for (size_t j = 0; j < max_k; ++j)
    {
      printf("%c", string[start + i + j]);
    } // for
    if (max_count > 0)
    {
      printf("%ld", max_count + 1);
    } // if
    printf(";");
    i += max_k * (max_count + 1);
  } // while
  printf("\n");
  return;
} // tandem_repeat_annotation

int main(int, char* [])
{
  char_t const* const string = "AACAAC";
  // length = 220
  //char_t const* const string = "CATGCTGGCCATATTCACTTGCCCACTTCTGCCCAGGGATCTATTTTTCTGTGGTGTGTATTCCCTGTGCCTTTGGGGGCATCTCTTATACTCATGAAATCAACAGAGGCTTGCATGTATCTATCTGTCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATGAGACAGGGTCTTGCTCTGTCACCCAGATTGGACTGCAGT";
  // length = 229
  //char_t const* const string = "ATATGTGAGTCAATTCCCCAAGTGAATTGCCTTCTATCTATCTATCTATCTATCTGTCTGTCTGTCTGTCTGTCTGTCTATCTATCTATATCTATCTATCATCTATCTATCCATATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATATCTATCGTCTATCTATCCAGTCTATCTACCTCCTATTAGTCTGTCTCTGGAGAACATTGACTAATACA";
  // length = 204
  //char_t const* const string = "GGCGACTGAGCAAGACTCAGTCTCAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAATTGTAAGGAGTTTTCTCAATTAATAACCCAAATAAGAGAATTCTTTCCATGTATCAATCATGATACTAAGCACTTTACACACATGTATGTTATGTAATCATTATATCATGCATGCAAGGTAATGAGT";

  tandem_repeat_annotation(string, 0, 6);

  return 0;
} // main

