#include <cstddef>
#include <cstdio>

#include <vector>

typedef char char_t;

typedef unsigned char       uint8_t;
typedef unsigned long long uint64_t;

static char_t const IUPAC_ALPHA[16] =
{
  'x',  // 0x00
  'A',  // 0x01
  'C',  // 0x02
  'M',  // 0x03  A | C
  'G',  // 0x04
  'R',  // 0x05  A | G
  'S',  // 0x06  C | G
  'V',  // 0x07  A | C | G
  'T',  // 0x08
  'W',  // 0x09  A | T
  'Y',  // 0x0a  C | T
  'H',  // 0x0b  A | C | T
  'K',  // 0x0c  G | T
  'D',  // 0x0d  A | G | T
  'B',  // 0x0e  C | G | T
  'N'   // 0x0f  A | C | G | T
}; // IUPAC_ALPHA

static char_t const IUPAC_BASE[4] =
{
  'A',
  'C',
  'G',
  'T'
}; // IUPAC_BASE

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

static uint64_t acid_map[128] = {0x0ull};

static double acid_frequency[128] = {.0f};
static double frame_shift_frequency[128][128][5] = {{{.0f}}};

struct Variant
{
  size_t       reference_start;
  size_t       reference_end;
  size_t       sample_start;
  size_t       sample_end;
  unsigned int type;
  union
  {
    size_t       weight;
    double       probability;
  }; // union
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

size_t Dprint_codon(size_t const index,
                    FILE*        stream = stderr)
{
  return fprintf(stream, "%c%c%c", IUPAC_BASE[index >> 0x4],
                                   IUPAC_BASE[(index >> 0x2) & 0x3],
                                   IUPAC_BASE[index & 0x3]);
} // Dprint_codon

size_t Dprint_truncated(char_t const* const string,
                        size_t const        start,
                        size_t const        end,
                        size_t const        length = 40,
                        FILE*               stream = stderr)
{
  if (end - start <= length - 3)
  {
    return fwrite(string + start, sizeof(char_t), end - start, stream);
  } // if
  return fwrite(string + start, sizeof(char_t), length / 2 - 1, stream) +
         fputs("...", stream) +
         fwrite(string + end - (length / 2) + 1, sizeof(char_t), length / 2 - 1, stream);
} // Dprint_truncated

uint8_t calculate_frame_shift(size_t const reference_1,
                              size_t const reference_2,
                              size_t const sample)
{
  uint8_t shift = FRAME_SHIFT_NONE;
  for (size_t i = 0; i < 64; ++i)
  {
    if (((acid_map[reference_1] >> i) & 0x1ull) == 0x1ull)
    {
      size_t const codon_reverse = ((i >> 0x4) | (i & 0xc) | ((i & 0x3) << 0x4)) ^ 0x3f;
      for (size_t j = 0; j < 64; ++j)
      {
        if (((acid_map[reference_2] >> j) & 0x1ull) == 0x1ull)
        {
          size_t const codon_1 = ((i & 0x3) << 0x4) | ((j & 0x3c) >> 0x2);
          size_t const codon_2 = ((i & 0xf) << 0x2) | (j >> 0x4);
          size_t const codon_reverse_1 = (((i & 0xc) >> 0x2) | ((i & 0x3) << 0x2) | (j & 0x30)) ^ 0x3f;
          size_t const codon_reverse_2 = ((i & 0x3) | ((j & 0x30) >> 0x2) | ((j & 0xc) << 0x2)) ^ 0x3f;
          for (size_t k = 0; k < 64; ++k)
          {
            if (((acid_map[sample] >> k) & 0x1ull) == 0x1ull)
            {

/*
              Dprint_codon(i);
              fprintf(stderr, " ");
              Dprint_codon(j);
              fprintf(stderr, " vs ");
              Dprint_codon(k);
              fprintf(stderr, ": ");
              Dprint_codon(codon_1);
              fprintf(stderr, " (fs_1) ");
              Dprint_codon(codon_2);
              fprintf(stderr, " (fs_2) ");
              Dprint_codon(codon_reverse);
              fprintf(stderr, " (fs_inv) ");
              Dprint_codon(codon_reverse_1);
              fprintf(stderr, " (fs1_inv) ");
              Dprint_codon(codon_reverse_2);
              fprintf(stderr, " (fs2_inv): ");
*/

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
              //fprintf(stderr, "0x%x\n", shift);
            } // if
          } // for
        } // if
      } // for
    } // if
  } // for
  return shift;
} // calculate_frame_shift

uint8_t frame_shift(char_t const reference_1,
                    char_t const reference_2,
                    char_t const sample)
{
  return frame_shift_map[reference_1 & 0x7f][reference_2 & 0x7f][sample & 0x7f];
} // frame_shift

void initialize_frame_shift_map(char_t const* const codon_string)
{
  for (size_t i = 0; i < 64; ++i)
  {
    acid_map[codon_string[i] & 0x7f] |= (0x1ull << i);
  } // for
  //return;
  for (size_t i = 0; i < 128; ++i)
  {
    if (acid_map[i] != 0x0ull)
    {
      for (size_t j = 0; j < 128; ++j)
      {
        if (acid_map[j] != 0x0ull)
        {
          for (size_t k = 0; k < 128; ++k)
          {
            if (acid_map[k] != 0x0ull)
            {
              uint8_t const shift = calculate_frame_shift(i, j, k);
              frame_shift_map[i][j][k] = shift;

              if ((shift & FRAME_SHIFT_1) == FRAME_SHIFT_1)
              {
                ++frame_shift_count[i][j][0];
                frame_shift_frequency[i][j][0] += acid_frequency[k];
              } // if
              if ((shift & FRAME_SHIFT_2) == FRAME_SHIFT_2)
              {
                ++frame_shift_count[i][j][1];
                frame_shift_frequency[i][j][1] += acid_frequency[k];
              } // if
              if ((shift & FRAME_SHIFT_REVERSE) == FRAME_SHIFT_REVERSE)
              {
                ++frame_shift_count[i][j][2];
                frame_shift_frequency[i][j][2] += acid_frequency[k];
              } // if
              if ((shift & FRAME_SHIFT_REVERSE_1) == FRAME_SHIFT_REVERSE_1)
              {
                ++frame_shift_count[i][j][3];
                frame_shift_frequency[i][j][3] += acid_frequency[k];
              } // if
              if ((shift & FRAME_SHIFT_REVERSE_2) == FRAME_SHIFT_REVERSE_2)
              {
                ++frame_shift_count[i][j][4];
                frame_shift_frequency[i][j][4] += acid_frequency[k];
              } // if

            } // if
          } // for
        } // if
      } // for
    } // if
  } // for
  return;
} // initialize_frame_shift_map

void initialize_acid_frequency(void)
{
  for (size_t i = 0; i < 128; ++i)
  {
    acid_frequency[i] = .05f;
  } // for
  acid_frequency['A'] = .09515673f;
  acid_frequency['C'] = .01157279f;
  acid_frequency['D'] = .05151007f;
  acid_frequency['E'] = .05762795f;
  acid_frequency['F'] = .03890338f;
  acid_frequency['G'] = .07374416f;
  acid_frequency['H'] = .02266328f;
  acid_frequency['I'] = .06010209f;
  acid_frequency['K'] = .04406110f;
  acid_frequency['L'] = .10672657f;
  acid_frequency['M'] = .02819341f;
  acid_frequency['N'] = .03945573f;
  acid_frequency['P'] = .04425210f;
  acid_frequency['Q'] = .04439959f;
  acid_frequency['R'] = .05510809f;
  acid_frequency['S'] = .05802322f;
  acid_frequency['T'] = .05398938f;
  acid_frequency['U'] = .00000221f;
  acid_frequency['V'] = .07073316f;
  acid_frequency['W'] = .01531018f;
  acid_frequency['X'] = .00001106f;
  acid_frequency['Y'] = .02845373f;
  return;
} // initialize_acid_frequency

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

/*
  fprintf(stderr, "  LCS table (forward):\n     ");
  for (size_t i = 1; i < reference_length; ++i)
  {
    fprintf(stderr, "%4c%c", reference[reference_start + i - 1], reference[reference_start + i]);
  } // for
  fprintf(stderr, "\n");
*/

/*
  fprintf(stderr, "  LCS table (reverse):\n      %c(%c)", reference[reference_end - 1], reference[reference_end - 1]);
  for (size_t i = 1; i < reference_length; ++i)
  {
    fprintf(stderr, "%4c%c", reference[reference_end - i - 1], reference[reference_end - i]);
  } // for
  fprintf(stderr, "\n");
*/

  Substring fs_substring[5];
  for (size_t i = 0; i < sample_length; ++i)
  {

    //fprintf(stderr, "    %c", sample[sample_start + i]);

    uint8_t const shift_reverse = frame_shift(reference[reference_end - 1], reference[reference_end - 1], sample[sample_start + i]);
    if ((shift_reverse & FRAME_SHIFT_REVERSE) == FRAME_SHIFT_REVERSE)
    {
      lcs[i % 2][0][2] = 1;
    } // if
    if (lcs[i % 2][0][2] > fs_substring[2].length)
    {
      fs_substring[2] = Substring(reference_start - lcs[i % 2][0][2] + 1, sample_start + i - lcs[i % 2][0][2] + 1, lcs[i % 2][0][2], FRAME_SHIFT_REVERSE);
    } // if

    //fprintf(stderr, "%5d", shift_reverse);

    for (size_t j = 1; j < reference_length; ++j)
    {
      uint8_t const shift_forward = frame_shift(reference[reference_start + j - 1], reference[reference_start + j], sample[sample_start + i]);
      uint8_t const shift_reverse = frame_shift(reference[reference_end - j - 1], reference[reference_end - j], sample[sample_start + i]);

      //fprintf(stderr, "%5d", shift_reverse);

      if ((shift_forward & FRAME_SHIFT_1) == FRAME_SHIFT_1)
      {
        lcs[i % 2][j][0] = lcs[(i + 1) % 2][j - 1][0] + 1;
      } // if
      else
      {
        lcs[i % 2][j][0] = 0;
      } // else
      if ((shift_forward & FRAME_SHIFT_2) == FRAME_SHIFT_2)
      {
        lcs[i % 2][j][1] = lcs[(i + 1) % 2][j - 1][1] + 1;
      } // if
      else
      {
        lcs[i % 2][j][1] = 0;
      } // else
      if ((shift_reverse & FRAME_SHIFT_REVERSE) == FRAME_SHIFT_REVERSE)
      {
        lcs[i % 2][j][2] = lcs[(i + 1) % 2][j - 1][2] + 1;
      } // if
      else
      {
        lcs[i % 2][j][2] = 0;
      } // else
      if ((shift_reverse & FRAME_SHIFT_REVERSE_1) == FRAME_SHIFT_REVERSE_1)
      {
        lcs[i % 2][j][3] = lcs[(i + 1) % 2][j - 1][3] + 1;
      } // if
      else
      {
        lcs[i % 2][j][3] = 0;
      } // else
      if ((shift_reverse & FRAME_SHIFT_REVERSE_2) == FRAME_SHIFT_REVERSE_2)
      {
        lcs[i % 2][j][4] = lcs[(i + 1) % 2][j - 1][4] + 1;
      } // if
      else
      {
        lcs[i % 2][j][4] = 0;
      } // else

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
        fs_substring[2] = Substring(reference_end - j - 1, sample_start + i - lcs[i % 2][j][2] + 1, lcs[i % 2][j][2], FRAME_SHIFT_REVERSE);
      } // if
      if (lcs[i % 2][j][3] > fs_substring[3].length)
      {
        fs_substring[3] = Substring(reference_end - j - 1, sample_start + i - lcs[i % 2][j][3] + 1, lcs[i % 2][j][3], FRAME_SHIFT_REVERSE_1);
      } // if
      if (lcs[i % 2][j][4] > fs_substring[4].length)
      {
        fs_substring[4] = Substring(reference_end - j - 1, sample_start + i - lcs[i % 2][j][4] + 1, lcs[i % 2][j][4], FRAME_SHIFT_REVERSE_2);
      } // if
    } // for
    //fprintf(stderr, "\n");
  } // for
  substring = std::vector<Substring>(1, fs_substring[0]);
  substring.push_back(fs_substring[1]);
  substring.push_back(fs_substring[2]);
  substring.push_back(fs_substring[3]);
  substring.push_back(fs_substring[4]);
  return;
} // LCS_frame_shift

void backtranslation(size_t              reference_DNA[],
                     size_t              sample_DNA[],
                     char_t const* const reference,
                     size_t const        reference_start,
                     char_t const* const sample,
                     size_t const        sample_start,
                     size_t const        length,
                     uint8_t const       type)
{
  for (size_t i = 0; i < 3 * length; ++ i)
  {
    reference_DNA[i] = 0;
    sample_DNA[i] = 0;
  } // for
  for (size_t p = 0; p < length; ++p)
  {
    for (size_t i = 0; i < 64; ++i)
    {
      if (((acid_map[reference[reference_start + p] & 0x7f] >> i) & 0x1ull) == 0x1ull)
      {
        size_t const codon_reverse = ((i >> 0x4) | (i & 0xc) | ((i & 0x3) << 0x4)) ^ 0x3f;
        for (size_t k = 0; k < 64; ++k)
        {
          if (((acid_map[sample[sample_start + length - p - 1] & 0x7f] >> k) & 0x1ull) == 0x1ull)
          {
            if ((type & FRAME_SHIFT_REVERSE) == FRAME_SHIFT_REVERSE && codon_reverse == k)
            {
              reference_DNA[p * 3] |= 0x1 << (i >> 4);
              reference_DNA[p * 3 + 1] |= 0x1 << ((i >> 2) & 0x3);
              reference_DNA[p * 3 + 2] |= 0x1 << (i & 0x3);
              sample_DNA[(length - p) * 3 - 3] |= 0x1 << (codon_reverse >> 4);
              sample_DNA[(length - p) * 3 - 2] |= 0x1 << ((codon_reverse >> 2) & 0x3);
              sample_DNA[(length - p) * 3 - 1] |= 0x1 << (codon_reverse & 0x3);
            } // if
          } // if
        } // for

        for (size_t j = 0; j < 64; ++j)
        {
          if (((acid_map[reference[reference_start + p + 1] & 0x7f] >> j) & 0x1ull) == 0x1ull)
          {
            size_t const codon_1 = ((i & 0x3) << 0x4) | ((j & 0x3c) >> 0x2);
            size_t const codon_2 = ((i & 0xf) << 0x2) | (j >> 0x4);
            size_t const codon_reverse_1 = (((i & 0xc) >> 0x2) | ((i & 0x3) << 0x2) | (j & 0x30)) ^ 0x3f;
            size_t const codon_reverse_2 = ((i & 0x3) | ((j & 0x30) >> 0x2) | ((j & 0xc) << 0x2)) ^ 0x3f;

            for (size_t k = 0; k < 64; ++k)
            {
              if (((acid_map[sample[sample_start + p] & 0x7f] >> k) & 0x1ull) == 0x1ull)
              {
                if ((type & FRAME_SHIFT_1) == FRAME_SHIFT_1 && codon_1 == k)
                {
                  reference_DNA[p * 3] |= 0x1 << (i >> 4);
                  reference_DNA[p * 3 + 1] |= 0x1 << ((i >> 2) & 0x3);
                  reference_DNA[p * 3 + 2] |= 0x1 << (i & 0x3);
                  sample_DNA[p * 3] |= 0x1 << (codon_1 >> 4);
                  sample_DNA[p * 3 + 1] |= 0x1 << ((codon_1 >> 2) & 0x3);
                  sample_DNA[p * 3 + 2] |= 0x1 << (codon_1 & 0x3);
                } // if
                if ((type & FRAME_SHIFT_2) == FRAME_SHIFT_2 && codon_2 == k)
                {
                  reference_DNA[p * 3] |= 0x1 << (i >> 4);
                  reference_DNA[p * 3 + 1] |= 0x1 << ((i >> 2) & 0x3);
                  reference_DNA[p * 3 + 2] |= 0x1 << (i & 0x3);
                  sample_DNA[p * 3] |= 0x1 << (codon_2 >> 4);
                  sample_DNA[p * 3 + 1] |= 0x1 << ((codon_2 >> 2) & 0x3);
                  sample_DNA[p * 3 + 2] |= 0x1 << (codon_2 & 0x3);
                } // if
              } // if

              if (((acid_map[sample[sample_start + length - p - 1] & 0x7f] >> k) & 0x1ull) == 0x1ull)
              {
                if ((type & FRAME_SHIFT_REVERSE_1) == FRAME_SHIFT_REVERSE_1 && codon_reverse_1 == k)
                {
                  reference_DNA[p * 3] |= 0x1 << (i >> 4);
                  reference_DNA[p * 3 + 1] |= 0x1 << ((i >> 2) & 0x3);
                  reference_DNA[p * 3 + 2] |= 0x1 << (i & 0x3);
                  sample_DNA[(length - p) * 3 - 3] |= 0x1 << (codon_reverse_1 >> 4);
                  sample_DNA[(length - p) * 3 - 2] |= 0x1 << ((codon_reverse_1 >> 2) & 0x3);
                  sample_DNA[(length - p) * 3 - 1] |= 0x1 << (codon_reverse_1 & 0x3);
                } // if
                if ((type & FRAME_SHIFT_REVERSE_2) == FRAME_SHIFT_REVERSE_2 && codon_reverse_2 == k)
                {
                  reference_DNA[p * 3] |= 0x1 << (i >> 4);
                  reference_DNA[p * 3 + 1] |= 0x1 << ((i >> 2) & 0x3);
                  reference_DNA[p * 3 + 2] |= 0x1 << (i & 0x3);
                  sample_DNA[(length - p) * 3 - 3] |= 0x1 << (codon_reverse_2 >> 4);
                  sample_DNA[(length - p) * 3 - 2] |= 0x1 << ((codon_reverse_2 >> 2) & 0x3);
                  sample_DNA[(length - p) * 3 - 1] |= 0x1 << (codon_reverse_2 & 0x3);
                } // if
              } // if
            } // for
          } // if
        } // for
      } // if
    } // for
  } // for
  return;
} // backtranslation

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

  fputs("Extraction (frame shift)\n", stderr);
  fprintf(stderr, "  reference %ld--%ld:  ", reference_start, reference_end);
  Dprint_truncated(reference, reference_start, reference_end);
  fprintf(stderr, " (%ld)\n", reference_length);
  fprintf(stderr, "  sample %ld--%ld:     ", sample_start, sample_end);
  Dprint_truncated(sample, sample_start, sample_end);
  fprintf(stderr, " (%ld)\n", sample_length);


  // First the base cases to end the recursion.
  if (reference_length <= 0 || sample_length <= 0)
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
  Dprint_truncated(reference, lcs.reference_index, lcs.reference_index + lcs.length);
  fprintf(stderr, " (%ld)\n    %ld--%ld: ", lcs.length, lcs.sample_index, lcs.sample_index + lcs.length);
  Dprint_truncated(sample, lcs.sample_index, lcs.sample_index + lcs.length);
  fprintf(stderr, " (%ld)", lcs.length);
  fputs("\n", stderr);


  fprintf(stderr, "  Probability:\n");
  double probability = 1.f;
  for (size_t i = 0; i < lcs.length; ++i)
  {
    double probability_compound = .0f;
    if ((lcs.type & FRAME_SHIFT_1) == FRAME_SHIFT_1)
    {
      probability_compound += frame_shift_frequency[reference[lcs.reference_index + i] & 0x7f][reference[lcs.reference_index + i + 1] & 0x7f][0];
    } // if
    if ((lcs.type & FRAME_SHIFT_2) == FRAME_SHIFT_2)
    {
      probability_compound += frame_shift_frequency[reference[lcs.reference_index + i] & 0x7f][reference[lcs.reference_index + i + 1] & 0x7f][1];
    } // if
    if ((lcs.type & FRAME_SHIFT_REVERSE) == FRAME_SHIFT_REVERSE)
    {
      probability_compound += frame_shift_frequency[reference[lcs.reference_index + i] & 0x7f][reference[lcs.reference_index + i] & 0x7f][2];
    } // if
    if ((lcs.type & FRAME_SHIFT_REVERSE_1) == FRAME_SHIFT_REVERSE_1)
    {
      probability_compound += frame_shift_frequency[reference[lcs.reference_index + i] & 0x7f][reference[lcs.reference_index + i + 1] & 0x7f][3];
    } // if
    if ((lcs.type & FRAME_SHIFT_REVERSE_2) == FRAME_SHIFT_REVERSE_2)
    {
      probability_compound += frame_shift_frequency[reference[lcs.reference_index + i] & 0x7f][reference[lcs.reference_index + i + 1] & 0x7f][4];
    } // if
    //fprintf(stderr, "    probability_compound: %lf\n", probability_compound);
    probability *= probability_compound;
  } // for
    fprintf(stderr, "    probability: %lf\n", probability);


  // DNA reconstruction
  fprintf(stderr, "  Backtranslation:\n    ");

  size_t reference_DNA[lcs.length * 3];
  size_t sample_DNA[lcs.length * 3];
  backtranslation(reference_DNA, sample_DNA, reference, lcs.reference_index, sample, lcs.sample_index, lcs.length, lcs.type);

  for (size_t i = 0; i < lcs.length; ++i)
  {
    fprintf(stderr, "%c%c%c ", IUPAC_ALPHA[reference_DNA[i * 3]], IUPAC_ALPHA[reference_DNA[i * 3 + 1]], IUPAC_ALPHA[reference_DNA[i * 3 + 2]]);
  } // for
  fprintf(stderr, "\n    ");
  for (size_t i = 0; i < lcs.length; ++i)
  {
    fprintf(stderr, "%c%c%c ", IUPAC_ALPHA[sample_DNA[i * 3]], IUPAC_ALPHA[sample_DNA[i * 3 + 1]], IUPAC_ALPHA[sample_DNA[i * 3 + 2]]);
  } // for
  fprintf(stderr, "\n");

  // Recursively apply this function to the prefixes of the strings.
  std::vector<Variant> prefix;
  extractor_frame_shift(prefix, reference, reference_start, lcs.reference_index, sample, sample_start, lcs.sample_index);


  // Recursively apply this function to the suffixes of the strings.
  std::vector<Variant> suffix;
  extractor_frame_shift(suffix, reference, lcs.reference_index + lcs.length, reference_end, sample, lcs.sample_index + lcs.length, sample_end);


  // Add all variants (in order) to the annotation vector.
  annotation.insert(annotation.end(), prefix.begin(), prefix.end());
  Variant variant(lcs.reference_index, lcs.reference_index + lcs.length, lcs.sample_index, lcs.sample_index + lcs.length, FRAME_SHIFT | lcs.type);
  variant.probability = probability;
  annotation.push_back(variant);
  annotation.insert(annotation.end(), suffix.begin(), suffix.end());

  return;
} // extractor_frame_shift

int main(int, char* [])
{
  initialize_acid_frequency();
  initialize_frame_shift_map(CODON_STRING);

/*
  std::vector<Variant> annotation;
  extractor_frame_shift(annotation, "GQPVEFTWQSDDGISL",
                                    0, 16,
                                    "ATLIPVLNGVHQAGLS",
                                    0, 16);

  // Printing the variants.
  fprintf(stdout, "Annotation (%ld):\n", annotation.size());
  for (std::vector<Variant>::iterator it = annotation.begin(); it != annotation.end(); ++it)
  {
    fprintf(stdout, "  %ld--%ld, %ld--%ld, %d, %lf, %ld--%ld\n", it->reference_start, it->reference_end, it->sample_start, it->sample_end, it->type, 1.f - it->probability, it->transposition_start, it->transposition_end);
  } // for


  return 1;
*/
  {
    fprintf(stdout, "\n%s %s\n", "MDYSL", "MALFP");
    std::vector<Variant> annotation;
    extractor_frame_shift(annotation, "MDYSL", 1, 5, "MALFP", 1, 5);

    // Printing the variants.
    fprintf(stdout, "Annotation (%ld):\n", annotation.size());
    for (std::vector<Variant>::iterator it = annotation.begin(); it != annotation.end(); ++it)
    {
      fprintf(stdout, "  %ld--%ld, %ld--%ld, %d, %lf, %ld--%ld\n", it->reference_start, it->reference_end, it->sample_start, it->sample_end, it->type, 1.f - it->probability, it->transposition_start, it->transposition_end);
    } // for
  }

  {
    fprintf(stdout, "\n%s %s\n", "MDYSL", "MATIP");
    std::vector<Variant> annotation;
    extractor_frame_shift(annotation, "MDYSL", 1, 5, "MATIP", 1, 5);

    // Printing the variants.
    fprintf(stdout, "Annotation (%ld):\n", annotation.size());
    for (std::vector<Variant>::iterator it = annotation.begin(); it != annotation.end(); ++it)
    {
      fprintf(stdout, "  %ld--%ld, %ld--%ld, %d, %lf, %ld--%ld\n", it->reference_start, it->reference_end, it->sample_start, it->sample_end, it->type, 1.f - it->probability, it->transposition_start, it->transposition_end);
    } // for
  }

  {
    fprintf(stdout, "\n%s %s\n", "MDYSL", "MTIPW");
    std::vector<Variant> annotation;
    extractor_frame_shift(annotation, "MDYSL", 1, 5, "MTIPW", 1, 5);

    // Printing the variants.
    fprintf(stdout, "Annotation (%ld):\n", annotation.size());
    for (std::vector<Variant>::iterator it = annotation.begin(); it != annotation.end(); ++it)
    {
      fprintf(stdout, "  %ld--%ld, %ld--%ld, %d, %lf, %ld--%ld\n", it->reference_start, it->reference_end, it->sample_start, it->sample_end, it->type, 1.f - it->probability, it->transposition_start, it->transposition_end);
    } // for
  }

  {
    fprintf(stdout, "\n%s %s\n", "MDYSL", "MLFPG");
    std::vector<Variant> annotation;
    extractor_frame_shift(annotation, "MDYSL", 1, 5, "MLFPG", 1, 5);

    // Printing the variants.
    fprintf(stdout, "Annotation (%ld):\n", annotation.size());
    for (std::vector<Variant>::iterator it = annotation.begin(); it != annotation.end(); ++it)
    {
      fprintf(stdout, "  %ld--%ld, %ld--%ld, %d, %lf, %ld--%ld\n", it->reference_start, it->reference_end, it->sample_start, it->sample_end, it->type, 1.f - it->probability, it->transposition_start, it->transposition_end);
    } // for
  }

  {
    fprintf(stdout, "\n%s %s\n", "MDYSL", "MQGIV");
    std::vector<Variant> annotation;
    extractor_frame_shift(annotation, "MDYSL", 1, 5, "MQGIV", 1, 5);

    // Printing the variants.
    fprintf(stdout, "Annotation (%ld):\n", annotation.size());
    for (std::vector<Variant>::iterator it = annotation.begin(); it != annotation.end(); ++it)
    {
      fprintf(stdout, "  %ld--%ld, %ld--%ld, %d, %lf, %ld--%ld\n", it->reference_start, it->reference_end, it->sample_start, it->sample_end, it->type, 1.f - it->probability, it->transposition_start, it->transposition_end);
    } // for
  }

  {
    fprintf(stdout, "\n%s %s\n", "MDYSL", "MAGNS");
    std::vector<Variant> annotation;
    extractor_frame_shift(annotation, "MDYSL", 1, 5, "MAGNS", 1, 5);

    // Printing the variants.
    fprintf(stdout, "Annotation (%ld):\n", annotation.size());
    for (std::vector<Variant>::iterator it = annotation.begin(); it != annotation.end(); ++it)
    {
      fprintf(stdout, "  %ld--%ld, %ld--%ld, %d, %lf, %ld--%ld\n", it->reference_start, it->reference_end, it->sample_start, it->sample_end, it->type, 1.f - it->probability, it->transposition_start, it->transposition_end);
    } // for
  }

  {
    fprintf(stdout, "\n%s %s\n", "MDYSL", "MDRE*");
    std::vector<Variant> annotation;
    extractor_frame_shift(annotation, "MDYSL", 1, 5, "MDRE*", 1, 5);

    // Printing the variants.
    fprintf(stdout, "Annotation (%ld):\n", annotation.size());
    for (std::vector<Variant>::iterator it = annotation.begin(); it != annotation.end(); ++it)
    {
      fprintf(stdout, "  %ld--%ld, %ld--%ld, %d, %lf, %ld--%ld\n", it->reference_start, it->reference_end, it->sample_start, it->sample_end, it->type, 1.f - it->probability, it->transposition_start, it->transposition_end);
    } // for
  }

  {
    fprintf(stdout, "\n%s %s\n", "MDYSL", "IRE*S");
    std::vector<Variant> annotation;
    extractor_frame_shift(annotation, "MDYSL", 0, 5, "IRE*S", 0, 5);

    // Printing the variants.
    fprintf(stdout, "Annotation (%ld):\n", annotation.size());
    for (std::vector<Variant>::iterator it = annotation.begin(); it != annotation.end(); ++it)
    {
      fprintf(stdout, "  %ld--%ld, %ld--%ld, %d, %lf, %ld--%ld\n", it->reference_start, it->reference_end, it->sample_start, it->sample_end, it->type, 1.f - it->probability, it->transposition_start, it->transposition_end);
    } // for
  }

  {
    fprintf(stdout, "\n%s %s\n", "MDYSL", "TGNSP");
    std::vector<Variant> annotation;
    extractor_frame_shift(annotation, "MDYSL", 0, 5, "TGNSP", 0, 5);

    // Printing the variants.
    fprintf(stdout, "Annotation (%ld):\n", annotation.size());
    for (std::vector<Variant>::iterator it = annotation.begin(); it != annotation.end(); ++it)
    {
      fprintf(stdout, "  %ld--%ld, %ld--%ld, %d, %lf, %ld--%ld\n", it->reference_start, it->reference_end, it->sample_start, it->sample_end, it->type, 1.f - it->probability, it->transposition_start, it->transposition_end);
    } // for
  }

  {
    fprintf(stdout, "\n%s %s\n", "MDYSLAAALTLHGH", "MTIPWRSPHFHGH");
    std::vector<Variant> annotation;
    extractor_frame_shift(annotation, "MDYSLAAALTLHGH", 1, 11, "MTIPWRSPHFHGH", 1, 10);

    // Printing the variants.
    fprintf(stdout, "Annotation (%ld):\n", annotation.size());
    for (std::vector<Variant>::iterator it = annotation.begin(); it != annotation.end(); ++it)
    {
      fprintf(stdout, "  %ld--%ld, %ld--%ld, %d, %lf, %ld--%ld\n", it->reference_start, it->reference_end, it->sample_start, it->sample_end, it->type, 1.f - it->probability, it->transposition_start, it->transposition_end);
    } // for
  }

  return 0;
} // main

