#include <cstddef>
#include <cstdio>

#include <vector>

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

struct FS_Substring
{
  size_t  reference_index;
  size_t  sample_index;
  size_t  length;
  uint8_t type;

  inline FS_Substring(size_t const  reference_index,
                      size_t const  sample_index,
                      size_t const  length,
                      uint8_t const type):
         reference_index(reference_index),
         sample_index(sample_index),
         length(length),
         type(type) { }

  inline FS_Substring(void): length(0) { }
}; // FS_Substring


uint8_t calculate_frame_shift(uint64_t const acid_map[],
                              size_t const   reference_1,
                              size_t const   reference_2,
                              size_t const   sample)
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
          } // for
        } // if
      } // for
    } // if
  } // for
  return;
} // initialize_frame_shift_map

uint8_t frame_shift(char_t const reference_1,
                    char_t const reference_2,
                    char_t const sample)
{
  return frame_shift_map[reference_1 & 0x7f][reference_2 & 0x7f][sample & 0x7f];
} // frame_shift

size_t lcs_frame_shift(std::vector<FS_Substring> &substring,
                       char_t const* const        reference,
                       size_t const               reference_start,
                       size_t const               reference_end,
                       char_t const* const        sample,
                       size_t const               sample_start,
                       size_t const               sample_end)
{
  size_t const reference_length = reference_end - reference_start;
  size_t const sample_length = sample_end - sample_start;

  size_t lcs[2][reference_length][5];
  for (size_t i = 0; i < reference_length; ++i)
  {
    lcs[0][i][0] = 0;
    lcs[0][i][1] = 0;
    lcs[0][i][2] = 0;
    lcs[0][i][3] = 0;
    lcs[0][i][4] = 0;
    lcs[1][i][0] = 0;
    lcs[1][i][1] = 0;
    lcs[1][i][2] = 0;
    lcs[1][i][3] = 0;
    lcs[1][i][4] = 0;
  } // for

  FS_Substring fs_substring[5];
  for (size_t i = 0; i < sample_length; ++i)
  {
    uint8_t const shift_reverse = frame_shift(reference[reference_end - 1], reference[reference_end - 2], sample[sample_start + i]);
    if ((shift_reverse & FRAME_SHIFT_REVERSE) == FRAME_SHIFT_REVERSE)
    {
      lcs[i % 2][0][2] = 1;
      fs_substring[2] = FS_Substring(reference_start, sample_start + i, lcs[i % 2][0][2], FRAME_SHIFT_REVERSE);
    } // if
    for (size_t j = 1; j < reference_length; ++j)
    {
      uint8_t const shift_forward = frame_shift(reference[reference_start + j - 1], reference[reference_start + j], sample[sample_start + i]);
      uint8_t const shift_reverse = frame_shift(reference[reference_end - j - 1], reference[reference_end - j], sample[sample_start + i]);
      if ((shift_forward & FRAME_SHIFT_1) == FRAME_SHIFT_1)
      {
        lcs[i % 2][j][0] = lcs[(i + 1) % 2][j - 1][0] + 1;
      } // if
      if ((shift_forward & FRAME_SHIFT_2) == FRAME_SHIFT_2)
      {
        lcs[i % 2][j][1] = lcs[(i + 1) % 2][j - 1][1] + 1;
      } // if
      if ((shift_reverse & FRAME_SHIFT_REVERSE) == FRAME_SHIFT_REVERSE)
      {
        lcs[i % 2][j][2] = lcs[(i + 1) % 2][j - 1][2] + 1;
      } // if
      if ((shift_reverse & FRAME_SHIFT_REVERSE_1) == FRAME_SHIFT_REVERSE_1)
      {
        lcs[i % 2][j][3] = lcs[(i + 1) % 2][j - 1][3] + 1;
      } // if
      if ((shift_reverse & FRAME_SHIFT_REVERSE_2) == FRAME_SHIFT_REVERSE_2)
      {
        lcs[i % 2][j][4] = lcs[(i + 1) % 2][j - 1][4] + 1;
      } // if
      if (lcs[i % 2][j][0] > fs_substring[0].length)
      {
        fs_substring[0] = FS_Substring(reference_start + j - lcs[i % 2][j][0], sample_start + i - lcs[i % 2][j][0] + 1, lcs[i % 2][j][0], FRAME_SHIFT_1);
      } // if
      if (lcs[i % 2][j][1] > fs_substring[1].length)
      {
        fs_substring[1] = FS_Substring(reference_start + j - lcs[i % 2][j][1], sample_start + i - lcs[i % 2][j][1] + 1, lcs[i % 2][j][1], FRAME_SHIFT_2);
      } // if
      if (lcs[i % 2][j][2] > fs_substring[2].length)
      {
        fs_substring[2] = FS_Substring(reference_start + j - lcs[i % 2][j][2], sample_start + i - lcs[i % 2][j][2], lcs[i % 2][j][2], FRAME_SHIFT_REVERSE);
      } // if
      if (lcs[i % 2][j][3] > fs_substring[3].length)
      {
        fs_substring[3] = FS_Substring(reference_start + j - lcs[i % 2][j][3], sample_start + i - lcs[i % 2][j][3] + 1, lcs[i % 2][j][3], FRAME_SHIFT_REVERSE_1);
      } // if
      if (lcs[i % 2][j][4] > fs_substring[4].length)
      {
        fs_substring[4] = FS_Substring(reference_start + j - lcs[i % 2][j][4], sample_start + i - lcs[i % 2][j][4] + 1, lcs[i % 2][j][4], FRAME_SHIFT_REVERSE_2);
      } // if
    } // for
  } // for
  substring = std::vector<FS_Substring>(1, fs_substring[0]);
  substring.push_back(fs_substring[1]);
  substring.push_back(fs_substring[2]);
  substring.push_back(fs_substring[3]);
  substring.push_back(fs_substring[4]);
  return 0;
} // lcs_frame_shift

FS_Substring maximize_frame_shift(std::vector<FS_Substring> const &substring)
{
  FS_Substring fs_substring;
  size_t max_length = 0;
  for (size_t i = 0; i < 5; ++i)
  {
    if (substring[i].length > max_length)
    {
      max_length = substring[i].length;
      fs_substring = substring[i];
    } // if
    else if (substring[i].length == max_length)
    {
      if ((fs_substring.type == FRAME_SHIFT_1 && substring[i].type == FRAME_SHIFT_2) ||
          (fs_substring.type == FRAME_SHIFT_2 && substring[i].type == FRAME_SHIFT_1) ||
          (fs_substring.type == FRAME_SHIFT_REVERSE_1 && substring[i].type == FRAME_SHIFT_REVERSE_2) ||
          (fs_substring.type == FRAME_SHIFT_REVERSE_2 && substring[i].type == FRAME_SHIFT_REVERSE_1))
      {
        fs_substring.type |= substring[i].type;
      } // if
    } // if
  } // for
  return fs_substring;
} // maximize_frame_shift

int main(int, char* [])
{
  initialize_frame_shift_map(CODON_STRING);

  char_t const* const reference = "MDYSL";
  size_t const start = 1;
  size_t const end   = 5;

  std::vector<FS_Substring> substring;
  FS_Substring fs_lcs;

  lcs_frame_shift(substring, reference, start, end, "MALFP", start, end);
  fs_lcs = maximize_frame_shift(substring);
  printf("fs_lcs: %d @ %ld, %ld length = %ld\n", fs_lcs.type, fs_lcs.reference_index, fs_lcs.sample_index, fs_lcs.length);
  lcs_frame_shift(substring, reference, start, end, "MATIP", start, end);
  fs_lcs = maximize_frame_shift(substring);
  printf("fs_lcs: %d @ %ld, %ld length = %ld\n", fs_lcs.type, fs_lcs.reference_index, fs_lcs.sample_index, fs_lcs.length);
  lcs_frame_shift(substring, reference, start, end, "MTIPW", start, end);
  fs_lcs = maximize_frame_shift(substring);
  printf("fs_lcs: %d @ %ld, %ld length = %ld\n", fs_lcs.type, fs_lcs.reference_index, fs_lcs.sample_index, fs_lcs.length);
  lcs_frame_shift(substring, reference, start, end, "MLFPG", start, end);
  fs_lcs = maximize_frame_shift(substring);
  printf("fs_lcs: %d @ %ld, %ld length = %ld\n", fs_lcs.type, fs_lcs.reference_index, fs_lcs.sample_index, fs_lcs.length);
  lcs_frame_shift(substring, reference, start, end, "MQGIV", start, end);
  fs_lcs = maximize_frame_shift(substring);
  printf("fs_lcs: %d @ %ld, %ld length = %ld\n", fs_lcs.type, fs_lcs.reference_index, fs_lcs.sample_index, fs_lcs.length);
  lcs_frame_shift(substring, reference, start, end, "MAGNS", start, end);
  fs_lcs = maximize_frame_shift(substring);
  printf("fs_lcs: %d @ %ld, %ld length = %ld\n", fs_lcs.type, fs_lcs.reference_index, fs_lcs.sample_index, fs_lcs.length);
  lcs_frame_shift(substring, reference, start, end, "MDRE*", start, end);
  fs_lcs = maximize_frame_shift(substring);
  printf("fs_lcs: %d @ %ld, %ld length = %ld\n", fs_lcs.type, fs_lcs.reference_index, fs_lcs.sample_index, fs_lcs.length);
  lcs_frame_shift(substring, reference, start, end, "MRE*S", start, end);
  fs_lcs = maximize_frame_shift(substring);
  printf("fs_lcs: %d @ %ld, %ld length = %ld\n", fs_lcs.type, fs_lcs.reference_index, fs_lcs.sample_index, fs_lcs.length);
  lcs_frame_shift(substring, reference, start, end, "MGNSP", start, end);
  fs_lcs = maximize_frame_shift(substring);
  printf("fs_lcs: %d @ %ld, %ld length = %ld\n", fs_lcs.type, fs_lcs.reference_index, fs_lcs.sample_index, fs_lcs.length);

  return 0;
} // main

