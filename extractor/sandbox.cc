#include <cstddef>
#include <cstdio>

typedef char char_t;

typedef unsigned char       uint8_t;
typedef unsigned long long uint64_t;

static uint8_t const FRAME_SHIFT_NONE      = 0x00;
static uint8_t const FRAME_SHIFT_1         = 0x01;
static uint8_t const FRAME_SHIFT_2         = 0x02;
static uint8_t const FRAME_SHIFT_REVERSE   = 0x04;
static uint8_t const FRAME_SHIFT_REVERSE_1 = 0x08;
static uint8_t const FRAME_SHIFT_REVERSE_2 = 0x10;

char_t const* const CODON_STRING = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";

static uint8_t frame_shift_map[128][128][128] = {{{FRAME_SHIFT_NONE}}};

uint8_t calculate_frame_shift(uint64_t const acid_map[], size_t const reference_1, size_t const reference_2, size_t const sample)
{
  uint8_t shift = FRAME_SHIFT_NONE;
  for (size_t i = 0; i < 64; ++i)
  {
    if (((acid_map[reference_1] >> i) & 0x1) == 0x1)
    {
      size_t const codon_reverse = ((i >> 0x4) | (i & 0xc) | ((i & 0x3) << 0x4)) ^ 0x3f;
      for (size_t j = 0; j < 64; ++j)
      {
        if (((acid_map[reference_2] >> j) & 0x1) == 0x1)
        {
          size_t const codon_1 = ((i & 0x3) << 0x4) | ((j & 0x3c) >> 0x2);
          size_t const codon_2 = ((i & 0xf) << 0x2) | (j >> 0x4);
          size_t const codon_reverse_1 = (((i & 0xc) >> 0x2) | ((i & 0x3) << 0x2) | (j & 0x30)) ^ 0x3f;
          size_t const codon_reverse_2 = ((i & 0x3) | ((j & 0x30) >> 0x2) | ((j & 0xc) << 0x2)) ^ 0x3f;
          for (size_t k = 0; k < 64; ++k)
          {
            if (((acid_map[sample] >> k) & 0x1) == 0x1)
            {
              if (codon_1 == k)
              {
                shift |= FRAME_SHIFT_1;
              } // if
              if (codon_2 == k)
              {
                shift |= FRAME_SHIFT_2;
              } // if
              if (codon_reverse == k)
              {
                shift |= FRAME_SHIFT_REVERSE;
              } // if
              if (codon_reverse_1 == k)
              {
                shift |= FRAME_SHIFT_REVERSE_1;
              } // if
              if (codon_reverse_2 == k)
              {
                shift |= FRAME_SHIFT_REVERSE_2;
              } // if
            } // if
          } // for
        } // if
      } // for
    } // if
  } // for
  return shift;
} // calculate_frame_shift

void initialize_frame_shift_map(char_t const* const codon_string)
{
  uint64_t acid_map[128] = {0x0};
  for (size_t i = 0; i < 64; ++i)
  {
    acid_map[codon_string[i] & 0x7f] |= (0x1ll << i);
  } // for
  for (size_t i = 0; i < 128; ++i)
  {
    if (acid_map[i] != 0x0)
    {
      for (size_t j = 0; j < 128; ++j)
      {
        if (acid_map[j] != 0x0)
        {
          for (size_t k = 0; k < 128; ++k)
          {
            if (acid_map[k] != 0x0)
            {
              frame_shift_map[i][j][k] = calculate_frame_shift(acid_map, i, j, k);
            } // if
          } // if
        } // if
      } // for
    } // if
  } // for
  return;
} // initialize_frame_shift_map

uint8_t frame_shift(char_t const reference_1, char_t const reference_2, char_t const sample)
{
  return frame_shift_map[reference_1 & 0x7f][reference_2 & 0x7f][sample & 0x7f];
} // frame_shift

int main(int, char* [])
{
  initialize_frame_shift_map(CODON_STRING);

  char_t const* const reference = "MDYSL";
  //char_t const* const sample    = "MALFP";
  char_t const* const sample      = "MLFPW";
  size_t const start = 1;
  size_t const end   = 5;

  
  for (size_t i = start; i < end - 1; ++i)
  {
    printf("%c %c vs %c: %d\n", reference[i], reference[i + 1], sample[i], frame_shift(reference[i], reference[i + 1], sample[i]));
  } // for


  return 0;
} // main

