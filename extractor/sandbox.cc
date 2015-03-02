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
//   Revision: 1.0.0
//   Date:     2015/03/02
// *******************************************************************
// DESCRIPTION:
//  General testing playground: protein frame shift construction code.
// *******************************************************************

#include <cstdio>
#include <cstdlib>

typedef char char_t;

#include <map>

static unsigned int const FRAME_SHIFT_0            = 0x00;
static unsigned int const FRAME_SHIFT_1            = 0x01;
static unsigned int const FRAME_SHIFT_2            = 0x02;
static unsigned int const FRAME_SHIFT_REVERSE_IDEM = 0x04;
static unsigned int const FRAME_SHIFT_REVERSE_1    = 0x08;
static unsigned int const FRAME_SHIFT_REVERSE_2    = 0x10;

static char_t const* const codon_string = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";

static char_t const IUPAC_BASE[4] =
{
  'A',
  'C',
  'G',
  'T'
}; // IUPAC_BASE

char_t codon_table[64] = {'\0'};
std::multimap<char_t const, unsigned int const> codon_map;

unsigned int frame_shift(char_t const reference, char_t const sample_1, char_t const sample_2)
{
  std::pair<std::multimap<char_t const, unsigned int const>::const_iterator, std::multimap<char_t const, unsigned int const>::const_iterator> const range_1 = codon_map.equal_range(sample_1);
  std::pair<std::multimap<char_t const, unsigned int const>::const_iterator, std::multimap<char_t const, unsigned int const>::const_iterator> const range_2 = codon_map.equal_range(sample_2);

  unsigned int shift = FRAME_SHIFT_0;
  for (std::multimap<char_t const, unsigned int const>::const_iterator it_1 = range_1.first; it_1 != range_1.second; ++it_1)
  {
    for (std::multimap<char_t const, unsigned int const>::const_iterator it_2 = range_2.first; it_2 != range_2.second; ++it_2)
    {
      const unsigned int codon_1 = ((it_1->second << 0x2) | (it_2->second >> 0x4)) & 0x3f;
      const unsigned int codon_2 = ((it_1->second << 0x4) | (it_2->second >> 0x2)) & 0x3f;
      if (reference == codon_table[codon_1])
      {
        shift |= FRAME_SHIFT_1;
      } // if
      if (reference == codon_table[codon_2])
      {
        shift |= FRAME_SHIFT_2;
      } // if
      if (shift == (FRAME_SHIFT_1 | FRAME_SHIFT_2))
      {
        return shift;
      } // if
      printf("%d %d -> %d %d\n", it_1->second, it_2->second, codon_1, codon_2);
    } // for
  } // for
  return shift;
} // frame_shift

void print_codon(unsigned int const index)
{
  printf("%c%c%c", IUPAC_BASE[index >> 0x4], IUPAC_BASE[(index >> 0x2) & 0x3], IUPAC_BASE[index & 0x3]);
} // print_codon

unsigned int reverse_complement(unsigned int const index)
{
  return (~((index >> 0x4) | (((index >> 0x2) & 0x3) << 2) | ((index & 0x3) << 0x4)) & 0x3f);
} // reverse_complement

unsigned int frame_shift_reverse(char_t const reference, char_t const sample_1)
{
  std::pair<std::multimap<char_t const, unsigned int const>::const_iterator, std::multimap<char_t const, unsigned int const>::const_iterator> const range_1 = codon_map.equal_range(sample_1);

  unsigned int shift = FRAME_SHIFT_0;
  for (std::multimap<char_t const, unsigned int const>::const_iterator it_1 = range_1.first; it_1 != range_1.second; ++it_1)
  {
    if (reference == codon_table[reverse_complement(it_1->second)])
    {
      shift |= FRAME_SHIFT_REVERSE_IDEM;
    } // if
  } // for
  return shift;
} // frame_shift_reverse


int main(int, char* [])
{
  for (unsigned int i = 0; i < 64; ++i)
  {
    codon_table[i] = codon_string[i];
    codon_map.insert(std::pair<char_t const, unsigned int const>(codon_table[i], i));
  } // for

  unsigned int const shift = frame_shift_reverse('D', 'V');
  printf("%d\n", shift);
  return 0;
} // main

