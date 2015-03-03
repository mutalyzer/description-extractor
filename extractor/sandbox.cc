#include <cstdio>

typedef char char_t;
typedef unsigned long long uint64_t;

static unsigned int const FRAME_SHIFT_NONE      = 0x00;
static unsigned int const FRAME_SHIFT_1         = 0x01;
static unsigned int const FRAME_SHIFT_2         = 0x02;
static unsigned int const FRAME_SHIFT_REVERSE   = 0x04;
static unsigned int const FRAME_SHIFT_REVERSE_1 = 0x08;
static unsigned int const FRAME_SHIFT_REVERSE_2 = 0x10;


// default codon string
char_t const* const CODON_STRING = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";

char_t const IUPAC_BASE[4] =
{
  'A',
  'C',
  'G',
  'T'
}; // IUPAC_BASE

// bitstrings with ASCII index: for each acid (as ASCII value) which
// codons are coding for it (64 bits: 1 is coding)
static uint64_t acid_map[128]   = {0x0};
// for each codon what is the corresponding acid (as ASCII value)
static   char_t codon_table[64] = {'\0'};

void print_codon(FILE*stream, unsigned int const index)
{
  fprintf(stream, "%c%c%c", IUPAC_BASE[index >> 0x4],
                            IUPAC_BASE[(index >> 0x2) & 0x3],
                            IUPAC_BASE[index & 0x3]);
  return;
} // print_codon

void print_mapping(FILE* stream)
{
  for (unsigned int i = 0; i < 128; ++i)
  {
    if (acid_map[i] != 0x0)
    {
      fprintf(stream, "%c: ", static_cast<char_t>(i));
      for (unsigned int j = 0; j < 64; ++j)
      {
        if (((acid_map[i] >> j) & 0x1) == 0x1)
        {
          print_codon(stream, j);
          fprintf(stream, " (%d) ", j);
        } // if
      } // for
      fprintf(stream, "\n");
    } // if
  } // for
  for (unsigned int i = 0; i < 64; ++i)
  {
    print_codon(stream, i);
    fprintf(stream, " (%d): %c\n", i, codon_table[i]);
  } // for
  return;
} // print_mapping

// calculate all possible forwards frame shifts: find the sample codon
// in the combination of the two reference condons
unsigned int frame_shift_forward(char_t const reference_1,
                                 char_t const reference_2,
                                 char_t const sample)
{
  unsigned int shift = FRAME_SHIFT_NONE;
  for (unsigned int i = 0; i < 64; ++i)
  {
    // i-th bit set: so i-th codon is coding for this
    // acid (reference_1)
    if (((acid_map[static_cast<unsigned int>(reference_1)] >> i) & 0x1) == 0x1)
    {
      for (unsigned int j = 0; j < 64; ++j)
      {
        // j-th bit set: so j-th codon is coding for this
        // acid (reference_2)
        if (((acid_map[static_cast<unsigned int>(reference_2)] >> j) & 0x1) == 0x1)
        {
          // calculate the two frame shifted codons with the reference
          // codons.
          // codon_1 last two bases of reference_1 concat with first
          // base of reference_2
          unsigned int const codon_1 = ((i & 0xf) << 0x2) | (j >> 0x4);
          // codon_2 last base of reference_1 concat with first two
          // bases of reference_2
          unsigned int const codon_2 = ((i & 0x3) << 0x4) | ((j & 0x3c) >> 0x2);

          for (unsigned int k = 0; k < 64; ++k)
          {
            // k-th bit set: so k-th codon is coding for this
            // acid (sample)
            if (((acid_map[static_cast<unsigned int>(sample)] >> k) & 0x1) == 0x1)
            {
              if (codon_1 == k)
              {
                shift |= FRAME_SHIFT_2;
              } // if
              if (codon_2 == k)
              {
                shift |= FRAME_SHIFT_1;
              } // if
              if (shift == (FRAME_SHIFT_1 | FRAME_SHIFT_2))
              {
                //return shift;
              } // if
            } // if
          } // for
        } // if
      } // for
    } // if
  } // for
  return shift;
} // frame_shift_forward

// calculate all possible reverse frame shifts: find the sample codon
// in the combination of the (two) reference condon(s)
unsigned int frame_shift_reverse(char_t const reference_1,
                                 char_t const reference_2,
                                 char_t const sample)
{
  unsigned int shift = FRAME_SHIFT_NONE;
  for (unsigned int i = 0; i < 64; ++i)
  {
    // i-th bit set: so i-th codon is coding for this
    // acid (reference_1)
    if (((acid_map[static_cast<unsigned int>(reference_1)] >> i) & 0x1) == 0x1)
    {
      // reverse complement codon
      unsigned int const codon_0 = ((i >> 0x4) | (i & 0xc) | ((i & 0x3) << 0x4)) ^ 0x3f;
      for (unsigned int j = 0; j < 64; ++j)
      {
        // j-th bit set: so j-th codon is coding for this
        // acid (reference_2)
        if (((acid_map[static_cast<unsigned int>(reference_2)] >> j) & 0x1) == 0x1)
        {
          unsigned int const codon_1 = (((i & 0xc) >> 0x2) | ((i & 0x3) << 0x2) | (j & 0x30)) ^ 0x3f;
          unsigned int const codon_2 = ((i & 0x3) | ((j & 0x30) >> 0x2) | ((j & 0xc) << 0x2)) ^ 0x3f;
          for (unsigned int k = 0; k < 64; ++k)
          {
            // k-th bit set: so k-th codon is coding for this
            // acid (sample)
            if (((acid_map[static_cast<unsigned int>(sample)] >> k) & 0x1) == 0x1)
            {
              if (codon_0 == k)
              {
                shift |= FRAME_SHIFT_REVERSE;
              } // if
              if (codon_1 == k)
              {
                shift |= FRAME_SHIFT_REVERSE_2;
              } // if
              if (codon_2 == k)
              {
                shift |= FRAME_SHIFT_REVERSE_1;
              } // if
              if (shift == (FRAME_SHIFT_REVERSE | FRAME_SHIFT_REVERSE_1 | FRAME_SHIFT_REVERSE_2))
              {
                //return shift;
              } // if
            } // if
          } // for
        } // if
      } // for
    } // if
  } // for
  return shift;
} // frame_shift_reverse

// entry point
int main(int, char* [])
{
  // initialize acid_map and codon_map: assume CODON_STRING holds at
  // least 64 ASCII characters (from the lower 127)
  for (unsigned int i = 0; i < 64; ++i)
  {
    // update bitstring: set i-th bit
    acid_map[CODON_STRING[i] & 0x7f] |= (0x1ll << i);
    // update codon table
    codon_table[i] = CODON_STRING[i];
  } // for

//  printf("frame shift (forward) = %d\n", frame_shift_forward('D', 'Y', 'L'));
  printf("frame shift (reverse) = %d\n", frame_shift_reverse('D', 'Y', 'S'));

  return 0;
} // main

