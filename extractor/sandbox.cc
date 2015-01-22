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
//   Revision: 1.1.0
//   Date:     2015/01/22
// *******************************************************************
// DESCRIPTION:
//  General testing playground: automatic (small) repeat annotation
//  and protein frame shift construction code (in comments).
// *******************************************************************

#include <cstdio>
#include <cstdlib>

typedef char char_t;

static char_t const IUPAC_BASE[4] =
{
  'A',
  'C',
  'G',
  'T'
}; // IUPAC_BASE

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

int main(int, char* [])
{
  size_t const length = 220;
  //char_t const* const string = "CCATAGATAGATAGCC";

  // length = 220
  char_t const* const string = "CATGCTGGCCATATTCACTTGCCCACTTCTGCCCAGGGATCTATTTTTCTGTGGTGTGTATTCCCTGTGCCTTTGGGGGCATCTCTTATACTCATGAAATCAACAGAGGCTTGCATGTATCTATCTGTCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATGAGACAGGGTCTTGCTCTGTCACCCAGATTGGACTGCAGT";
  // length = 229
  //char_t const* const string = "ATATGTGAGTCAATTCCCCAAGTGAATTGCCTTCTATCTATCTATCTATCTATCTGTCTGTCTGTCTGTCTGTCTGTCTATCTATCTATATCTATCTATCATCTATCTATCCATATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATATCTATCGTCTATCTATCCAGTCTATCTACCTCCTATTAGTCTGTCTCTGGAGAACATTGACTAATACA";
  // length = 204
  //char_t const* const string = "GGCGACTGAGCAAGACTCAGTCTCAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAATTGTAAGGAGTTTTCTCAATTAATAACCCAAATAAGAGAATTCTTTCCATGTATCAATCATGATACTAAGCACTTTACACACATGTATGTTATGTAATCATTATATCATGCATGCAAGGTAATGAGT";

  for (size_t i = 0; i < length; ++i)
  {
    size_t count;
    //printf("%ld %c\n", i, string[i]);
    for (size_t k = 1; k < length / 2 + 1; ++k)
    {
      count = 0;
      //printf("  k = %ld\n", k);
      for (size_t j = i + k; j < length - k + 1; j += k)
      {
        /*
        printf("    %ld: ", j);
        for (size_t c = 0; c < k; ++c)
        {
          printf("%c", string[i + c]);
        } // for
        printf(" --- ");
        for (size_t c = 0; c < k; ++c)
        {
          printf("%c", string[j + c]);
        } // for
        printf("\n");
        */
        if (!string_match(string + i, string + j, k))
        {
          break;
        } // if
        ++count;
      } // for
      if (count > 0)
      {
        for (size_t j = 0; j < k; ++j)
        {
          printf("%c", string[i + j]);
        } // for
        printf("%ld;", count + 1);
        i += (count + 1) * k - 1;
        break;
      } // if
    } // for
    if (count == 0)
    {
      printf("%c;", string[i]);
    } // if
  } // for
  printf("\n");

  return 0;
} // main

