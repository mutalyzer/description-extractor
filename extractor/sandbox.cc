#include <cstddef>
#include <cstdio>

#include <vector>

typedef char char_t;

typedef unsigned char       uint8_t;
typedef unsigned long long uint64_t;

static unsigned int const IDENTITY            = 0x01;
static unsigned int const FRAME_SHIFT         = 0x20;

static uint8_t const FRAME_SHIFT_NONE      = 0x00;
static uint8_t const FRAME_SHIFT_1         = 0x01;
static uint8_t const FRAME_SHIFT_2         = 0x02;
static uint8_t const FRAME_SHIFT_REVERSE   = 0x04;
static uint8_t const FRAME_SHIFT_REVERSE_1 = 0x08;
static uint8_t const FRAME_SHIFT_REVERSE_2 = 0x10;

char_t const* const CODON_STRING = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";

static uint8_t frame_shift_map[128][128][128] = {{{FRAME_SHIFT_NONE}}};
static uint8_t frame_shift_count[128][128][5] = {{{0}}};

struct Variant
{
  size_t       reference_start;
  size_t       reference_end;
  size_t       sample_start;
  size_t       sample_end;
  unsigned int type;
  size_t       weight;
  size_t       transposition_start;
  size_t       transposition_end;

  inline Variant(size_t const       reference_start,
                 size_t const       reference_end,
                 size_t const       sample_start,
                 size_t const       sample_end,
                 unsigned int const type                = IDENTITY,
                 size_t const       weight              = 0,
                 size_t const       transposition_start = 0,
                 size_t const       transposition_end   = 0):
         reference_start(reference_start),
         reference_end(reference_end),
         sample_start(sample_start),
         sample_end(sample_end),
         type(type),
         weight(weight),
         transposition_start(transposition_start),
         transposition_end(transposition_end) { }

  inline Variant(void) { }
}; // Variant

struct Substring
{
  size_t  reference_index;
  size_t  sample_index;
  size_t  length;
  union
  {
    bool    reverse_complement;
    uint8_t type;
  }; // union

  inline Substring(size_t const  reference_index,
                   size_t const  sample_index,
                   size_t const  length,
                   bool const    reverse_complement = false):
         reference_index(reference_index),
         sample_index(sample_index),
         length(length),
         reverse_complement(reverse_complement) { }

  inline Substring(size_t const  reference_index,
                   size_t const  sample_index,
                   size_t const  length,
                   uint8_t const type):
         reference_index(reference_index),
         sample_index(sample_index),
         length(length),
         type(type) { }

  inline Substring(void): length(0) { }
}; // Substring

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
              uint8_t const shift = calculate_frame_shift(acid_map, i, j, k);
              frame_shift_map[i][j][k] = shift;
              if ((shift & FRAME_SHIFT_1) == FRAME_SHIFT_1)
              {
                ++frame_shift_count[i][j][0];
              } // if
              if ((shift & FRAME_SHIFT_2) == FRAME_SHIFT_2)
              {
                ++frame_shift_count[i][j][1];
              } // if
              if ((shift & FRAME_SHIFT_REVERSE) == FRAME_SHIFT_REVERSE)
              {
                ++frame_shift_count[i][j][2];
              } // if
              if ((shift & FRAME_SHIFT_REVERSE_1) == FRAME_SHIFT_REVERSE_1)
              {
                ++frame_shift_count[i][j][3];
              } // if
              if ((shift & FRAME_SHIFT_REVERSE_2) == FRAME_SHIFT_REVERSE_2)
              {
                ++frame_shift_count[i][j][4];
              } // if
            } // if
          } // for
          printf("%c%c: %d %d %d %d %d\n", static_cast<char_t>(i), static_cast<char_t>(j), frame_shift_count[i][j][0], frame_shift_count[i][j][1], frame_shift_count[i][j][2], frame_shift_count[i][j][3], frame_shift_count[i][j][4]);
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

void LCS_frame_shift(std::vector<Substring> &substring,
                     char_t const* const     reference,
                     size_t const            reference_start,
                     size_t const            reference_end,
                     char_t const* const     sample,
                     size_t const            sample_start,
                     size_t const            sample_end)
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

  Substring fs_substring[5];
  for (size_t i = 0; i < sample_length; ++i)
  {
    uint8_t const shift_reverse = frame_shift(reference[reference_end - 1], reference[reference_end - 2], sample[sample_start + i]);
    if ((shift_reverse & FRAME_SHIFT_REVERSE) == FRAME_SHIFT_REVERSE)
    {
      lcs[i % 2][0][2] = 1;
      fs_substring[2] = Substring(reference_start, sample_start + i, lcs[i % 2][0][2], FRAME_SHIFT_REVERSE);
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
        fs_substring[0] = Substring(reference_start + j - lcs[i % 2][j][0], sample_start + i - lcs[i % 2][j][0] + 1, lcs[i % 2][j][0], FRAME_SHIFT_1);
      } // if
      if (lcs[i % 2][j][1] > fs_substring[1].length)
      {
        fs_substring[1] = Substring(reference_start + j - lcs[i % 2][j][1], sample_start + i - lcs[i % 2][j][1] + 1, lcs[i % 2][j][1], FRAME_SHIFT_2);
      } // if
      if (lcs[i % 2][j][2] > fs_substring[2].length)
      {
        fs_substring[2] = Substring(reference_start + j - lcs[i % 2][j][2], sample_start + i - lcs[i % 2][j][2], lcs[i % 2][j][2], FRAME_SHIFT_REVERSE);
      } // if
      if (lcs[i % 2][j][3] > fs_substring[3].length)
      {
        fs_substring[3] = Substring(reference_start + j - lcs[i % 2][j][3], sample_start + i - lcs[i % 2][j][3] + 1, lcs[i % 2][j][3], FRAME_SHIFT_REVERSE_1);
      } // if
      if (lcs[i % 2][j][4] > fs_substring[4].length)
      {
        fs_substring[4] = Substring(reference_start + j - lcs[i % 2][j][4], sample_start + i - lcs[i % 2][j][4] + 1, lcs[i % 2][j][4], FRAME_SHIFT_REVERSE_2);
      } // if
    } // for
  } // for
  substring = std::vector<Substring>(1, fs_substring[0]);
  substring.push_back(fs_substring[1]);
  substring.push_back(fs_substring[2]);
  substring.push_back(fs_substring[3]);
  substring.push_back(fs_substring[4]);
  return;
} // LCS_frame_shift

void extractor_frame_shift(std::vector<Variant> &annotation,
                           char_t const* const   reference,
                           size_t const          reference_start,
                           size_t const          reference_end,
                           char_t const* const   sample,
                           size_t const          sample_start,
                           size_t const          sample_end)
{
  size_t const reference_length = reference_end - reference_start;
  size_t const sample_length = sample_end - sample_start;

  // First the base cases to end the recursion.
  if (reference_length <= 0 && sample_length <= 0)
  {
    return;
  } // if


  // Calculate the frame shift LCS of the two strings.
  std::vector<Substring> substring;
  LCS_frame_shift(substring, reference, reference_start, reference_end, sample, sample_start, sample_end);


  // Pick the ``best fitting'' frame shift LCS, i.e., pushed as far to
  // the start of the reference string as possible. Also update in
  // case of compound frame shifts.
  Substring lcs(0, 0, 0, FRAME_SHIFT_NONE);
  for (size_t i = 0; i < 5; ++i)
  {
    if (substring[i].length > lcs.length ||
        (substring[i].length == lcs.length && substring[i].reference_index < lcs.reference_index))
    {
      lcs = substring[i];
    } // if
    // Update compound frame shifts, e.g.,
    // FRAME_SHIFT_1 | FRAME_SHIFT_2
    else if (substring[i].length == lcs.length &&
             substring[i].reference_index == lcs.reference_index &&
             substring[i].sample_index == lcs.sample_index)
    {
      lcs.type |= substring[i].type;
    } // if
  } // for


  // No LCS found: no frame shift annotation.
  if (lcs.length <= 0)
  {
    return;
  } // if

  fprintf(stderr, "  LCS type = %d\n", lcs.type);
  fprintf(stderr, "    %ld--%ld: ", lcs.reference_index, lcs.reference_index + lcs.length);
  fprintf(stderr, " (%ld)\n    %ld--%ld: ", lcs.length, lcs.sample_index, lcs.sample_index + lcs.length);
  fprintf(stderr, " (%ld)", lcs.length);
  fputs("\n", stderr);


  size_t weight = 1;
  for (size_t i = 0; i < lcs.length; ++i)
  {
    size_t weight_compound = 0;
    if ((lcs.type & FRAME_SHIFT_1) == FRAME_SHIFT_1)
    {
      weight_compound += frame_shift_count[reference[lcs.reference_index + i] & 0x7f][reference[lcs.reference_index + i + 1] & 0x7f][0];
    } // if
    if ((lcs.type & FRAME_SHIFT_2) == FRAME_SHIFT_2)
    {
      weight_compound += frame_shift_count[reference[lcs.reference_index + i] & 0x7f][reference[lcs.reference_index + i + 1] & 0x7f][1];
    } // if
    if ((lcs.type & FRAME_SHIFT_REVERSE) == FRAME_SHIFT_REVERSE)
    {
      weight_compound += frame_shift_count[reference[lcs.reference_index + i] & 0x7f][reference[lcs.reference_index + i + 1] & 0x7f][2];
    } // if
    if ((lcs.type & FRAME_SHIFT_REVERSE_1) == FRAME_SHIFT_REVERSE_1)
    {
      weight_compound += frame_shift_count[reference[lcs.reference_index + i] & 0x7f][reference[lcs.reference_index + i + 1] & 0x7f][3];
    } // if
    if ((lcs.type & FRAME_SHIFT_REVERSE_2) == FRAME_SHIFT_REVERSE_2)
    {
      weight_compound += frame_shift_count[reference[lcs.reference_index + i] & 0x7f][reference[lcs.reference_index + i + 1] & 0x7f][4];
    } // if
    weight *= weight_compound;
  } // for


  // Recursively apply this function to the prefixes of the strings.
  std::vector<Variant> prefix;
  extractor_frame_shift(prefix, reference, reference_start, lcs.reference_index, sample, sample_start, lcs.sample_index);


  // Recursively apply this function to the suffixes of the strings.
  std::vector<Variant> suffix;
  extractor_frame_shift(suffix, reference, lcs.reference_index + lcs.length, reference_end, sample, lcs.sample_index + lcs.length, sample_end);


  // Add all variants (in order) to the annotation vector.
  annotation.insert(annotation.end(), prefix.begin(), prefix.end());
  annotation.push_back(Variant(lcs.reference_index, lcs.reference_index + lcs.length, lcs.sample_index, lcs.sample_index + lcs.length, FRAME_SHIFT | lcs.type, weight));
  annotation.insert(annotation.end(), suffix.begin(), suffix.end());

  return;
} // extractor_frame_shift

int main(int, char* [])
{
  initialize_frame_shift_map(CODON_STRING);

  std::vector<Variant> annotation;
  extractor_frame_shift(annotation, "MDYSLAAALTLHGH", 1, 11, "MTIPWRSPHFHGH", 1, 10);

  // Printing the variants.
  fprintf(stdout, "Annotation (%ld):\n", annotation.size());
  for (std::vector<Variant>::iterator it = annotation.begin(); it != annotation.end(); ++it)
  {
    fprintf(stdout, "%ld--%ld, %ld--%ld, %d, %ld, %ld--%ld\n", it->reference_start, it->reference_end, it->sample_start, it->sample_end, it->type, it->weight, it->transposition_start, it->transposition_end);
  } // for

  return 0;
} // main

