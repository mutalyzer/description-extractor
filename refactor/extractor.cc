#include "extractor.h"

namespace mutalyzer
{

size_t global_reference_length = 0;

size_t weight_position = 0;


size_t extract(std::vector<Variant> &variant,
               char_t const* const   reference,
               size_t const          reference_length,
               char_t const* const   sample,
               size_t const          sample_length,
               int const             type)
{
  global_reference_length = reference_length;
  weight_position = ceil(log10(reference_length));

  size_t const prefix = prefix_match(reference, reference_length, sample, sample_length);
  size_t const suffix = suffix_match(reference, reference_length, sample, sample_length, prefix);

#if defined(__debug__)
  fputs("reference: ", stderr);
  Dprint_truncated(reference, 0, reference_length);
  fprintf(stderr, " (%ld)\n", reference_length);
  fputs("sample:    ", stderr);
  Dprint_truncated(sample, 0, sample_length);
  fprintf(stderr, " (%ld)\n", sample_length);
  fputs("prefix: ", stderr);
  Dprint_truncated(reference, 0, prefix);
  fprintf(stderr, " (%ld)\n", prefix);
  fputs("suffix: ", stderr);
  Dprint_truncated(reference, reference_length - suffix, reference_length);
  fprintf(stderr, " (%ld)\n", suffix);
  fputs("Construction IUPAC complement\n", stderr);
#endif

  char_t const* complement = type == TYPE_DNA ? IUPAC_complement(reference, reference_length) : 0;

#if defined(__debug__)
  fputs("complement: ", stderr);
  Dprint_truncated(complement, 0, reference_length);
  fprintf(stderr, " (%ld)\n", reference_length);
  fprintf(stderr, "position weight: %ld\n", weight_position);
#endif

  if (prefix > 0)
  {
    variant.push_back(Variant(0, prefix, 0, prefix));
  } // if

  size_t const weight = extractor(variant, reference, complement, prefix, reference_length - suffix, sample, prefix, sample_length - suffix);

  if (suffix > 0)
  {
    variant.push_back(Variant(reference_length - suffix, reference_length, sample_length - suffix, sample_length));
  } // if

  delete[] complement;

  return weight;
} // extract

size_t extractor(std::vector<Variant> &variant,
                 char_t const* const   reference,
                 char_t const* const   complement,
                 size_t const          reference_start,
                 size_t const          reference_end,
                 char_t const* const   sample,
                 size_t const          sample_start,
                 size_t const          sample_end,
                 bool const            transposition)
{
  size_t const reference_length = reference_end - reference_start;
  size_t const sample_length = sample_end - sample_start;
  size_t const weight_trivial = weight_position + WEIGHT_DELETION_INSERTION + WEIGHT_BASE * sample_length + (reference_length > 1 ? weight_position + WEIGHT_SEPARATOR : 0);
  size_t weight = 0;

#if defined(__debug__)
  fputs("Extraction\n", stderr);
  fprintf(stderr, "  reference %ld--%ld:  ", reference_start, reference_end);
  Dprint_truncated(reference, reference_start, reference_end);
  fprintf(stderr, " (%ld)\n", reference_length);
  fprintf(stderr, "  complement %ld--%ld: ", reference_start, reference_end);
  Dprint_truncated(complement, reference_start, reference_end);
  fprintf(stderr, " (%ld)\n", reference_length);
  fprintf(stderr, "  sample %ld--%ld:     ", sample_start, sample_end);
  Dprint_truncated(sample, sample_start, sample_end);
  fprintf(stderr, " (%ld)\n", sample_length);
  fprintf(stderr, "  trivial weight: %ld\n", weight_trivial);
#endif

  // insertions
  if (reference_length <= 0)
  {
    if (sample_length > 0)
    {
      // only simple insertions TODO: transpositions
      weight = 2 * weight_position + WEIGHT_SEPARATOR + WEIGHT_INSERTION + WEIGHT_BASE * sample_length;
      variant.push_back(Variant(reference_start, reference_end, sample_start, sample_end, SUBSTITUTION, weight));
    } // if
    return weight;
  } // if

  // deletions
  if (sample_length <= 0)
  {
    weight = weight_position + WEIGHT_DELETION + (reference_length > 1 ? weight_position + WEIGHT_SEPARATOR : 0);
    variant.push_back(Variant(reference_start, reference_end, sample_start, sample_end, SUBSTITUTION, weight));
    return weight;
  } // if

  // simple substitutions
  if (reference_length == 1 && sample_length == 1)
  {
    weight = weight_position + 2 * WEIGHT_BASE + WEIGHT_SUBSTITUTION;
    variant.push_back(Variant(reference_start, reference_end, sample_start, sample_end, SUBSTITUTION, weight));
    return weight;
  } // if

  // these are always delins: 2 vs 1 or 1 vs 2
  if (reference_length < 3 && sample_length < 3)
  {
    variant.push_back(Variant(reference_start, reference_end, sample_start, sample_end, SUBSTITUTION, weight_trivial));
    return weight_trivial;
  } // if


  std::vector<Substring> substring;
  size_t const length = LCS(substring, reference, complement, reference_start, reference_end, sample, sample_start, sample_end);

  // deletion/insertion TODO: transpositions
  if (length <= 0 || substring.size() <= 0)
  {
    variant.push_back(Variant(reference_start, reference_end, sample_start, sample_end, SUBSTITUTION, weight_trivial));
    return weight_trivial;
  } // if


  size_t difference = reference_length + sample_length;
  std::vector<Substring>::iterator lcs = substring.begin();
  for (std::vector<Substring>::iterator it = substring.begin() ; it != substring.end(); ++it)
  {
    if (static_cast<size_t>(abs(it->reference_index - it->sample_index)) < difference)
    {
      difference = abs(it->reference_index - it->sample_index);
      lcs = it;
    } // if
  } // for
  if (lcs->reverse_complement)
  {
    weight = 2 * weight_position + WEIGHT_SEPARATOR + WEIGHT_INVERSE;
  } // if

#if defined(__debug__)
  fprintf(stderr, "  LCS (x%ld)\n", substring.size());
  for (std::vector<Substring>::iterator it = substring.begin() ; it != substring.end(); ++it)
  {
    if (!it->reverse_complement)
    {
      fprintf(stderr, "    %ld--%ld: ", it->reference_index, it->reference_index + length);
      Dprint_truncated(reference, it->reference_index, it->reference_index + length);
      fprintf(stderr, " (%ld)\n    %ld--%ld: ", length, it->sample_index, it->sample_index + length);
      Dprint_truncated(sample, it->sample_index, it->sample_index + length);
      fprintf(stderr, " (%ld)", length);
    } // if
    else
    {
      fprintf(stderr, "    %ld--%ld: ", it->reference_index, it->reference_index + length);
      Dprint_truncated(complement, it->reference_index, it->reference_index + length);
      fprintf(stderr, " (%ld), (reverse complement)\n    %ld--%ld: ", length, it->sample_index, it->sample_index + length);
      Dprint_truncated(sample, it->sample_index, it->sample_index + length);
      fprintf(stderr, " (%ld)", length);
    } // else
    if (it == lcs)
    {
      fputs(" selected", stderr);
    } // if
    fputs("\n", stderr);
  } // for
#endif

  std::vector<Variant> prefix;
  weight += extractor(prefix, reference, complement, reference_start, lcs->reference_index, sample, sample_start, lcs->sample_index, transposition);

  if (weight > weight_trivial)
  {
    variant.push_back(Variant(reference_start, reference_end, sample_start, sample_end, SUBSTITUTION, weight_trivial));
    return weight_trivial;
  } // if

  std::vector<Variant> suffix;
  weight += extractor(suffix, reference, complement, lcs->reference_index + length, reference_end, sample, lcs->sample_index + length, sample_end, transposition);

  if (weight > weight_trivial)
  {
    variant.push_back(Variant(reference_start, reference_end, sample_start, sample_end, SUBSTITUTION, weight_trivial));
    return weight_trivial;
  } // if


  for (std::vector<Variant>::iterator it = prefix.begin(); it != prefix.end(); ++it)
  {
    variant.push_back(*it);
  } // for
  if (!lcs->reverse_complement)
  {
    variant.push_back(Variant(lcs->reference_index, lcs->reference_index + length, lcs->sample_index, lcs->sample_index + length));
  } // if
  else
  {
    variant.push_back(Variant(lcs->reference_index, lcs->reference_index + length, lcs->sample_index, lcs->sample_index + length, REVERSE_COMPLEMENT, 2 * weight_position + WEIGHT_SEPARATOR + WEIGHT_INVERSE));
  } // else
  for (std::vector<Variant>::iterator it = suffix.begin(); it != suffix.end(); ++it)
  {
    variant.push_back(*it);
  } // for
  return weight;
} // extractor


size_t LCS(std::vector<Substring> &substring,
           char_t const* const     reference,
           char_t const* const     complement,
           size_t const            reference_start,
           size_t const            reference_end,
           char_t const* const     sample,
           size_t const            sample_start,
           size_t const            sample_end)
{
  return LCS_1(substring, reference, complement, reference_start, reference_end, sample, sample_start, sample_end);
} // LCS

size_t LCS_1(std::vector<Substring> &substring,
             char_t const* const     reference,
             char_t const* const     complement,
             size_t const            reference_start,
             size_t const            reference_end,
             char_t const* const     sample,
             size_t const            sample_start,
             size_t const            sample_end)
{
  size_t const reference_length = reference_end - reference_start;
  size_t const sample_length = sample_end - sample_start;

  // Just a fancy way of allocation a continuous 2D array in heap
  // space.
  typedef size_t array[2][reference_length];
  array &LCS_line = *(reinterpret_cast<array*>(new size_t[2 * reference_length]));
  array &LCS_line_rc = *(reinterpret_cast<array*>(new size_t[2 * reference_length]));

  size_t length = 0;

  // Filling the LCS matrix (actually only the current and the
  // previous row).
  for (size_t i = 0; i < sample_length; ++i)
  {
    for (size_t j = 0; j < reference_length; ++j)
    {
      // A match
      if (reference[reference_start + j] == sample[sample_start + i])
      {
        if (i == 0 || j == 0)
        {
          LCS_line[i % 2][j] = 1;
        } // if
        else
        {
          LCS_line[i % 2][j] = LCS_line[(i + 1) % 2][j - 1] + 1;
        } // else

        // Check for a new maximal length.
        if (LCS_line[i % 2][j] > length)
        {
          length = LCS_line[i % 2][j];
          substring = std::vector<Substring>(1, Substring(j - length + reference_start + 1, i - length + sample_start + 1, length));
        } // if
        // Found a LCS of the same maximal length.
        else if (LCS_line[i % 2][j] == length)
        {
          substring.push_back(Substring(j - length + reference_start + 1, i - length + sample_start + 1, length));
        } // if
      } // if
      else
      {
        LCS_line[i % 2][j] = 0;
      } // else

      // If applicable check for a LCS in reverse complement space.
      // The same code is used as before but the complement string is
      // travesed backwards (towards the start).
      if (complement != 0 && complement[reference_start + j] == sample[sample_start + i])
      {
        if (i == 0 || j == 0)
        {
          LCS_line_rc[i % 2][j] = 1;
        } // if
        else
        {
          LCS_line_rc[i % 2][j] = LCS_line_rc[(i + 1) % 2][j - 1] + 1;
        } // else
        if (LCS_line_rc[i % 2][j] > length)
        {
          length = LCS_line_rc[i % 2][j];
          substring = std::vector<Substring>(1, Substring(j - length + reference_start + 1, i - length + sample_start + 1, length, true));
        } // if
        else if (LCS_line_rc[i % 2][j] == length)
        {
          substring.push_back(Substring(j - length + reference_start + 1, i - length + sample_start + 1, length, true));
        } // if
      } // if
      else
      {
        LCS_line_rc[i % 2][j] = 0;
      } // else
    } // for
  } // for

  delete[] &LCS_line;
  delete[] &LCS_line_rc;

  return length;
} // LCS_1

/*
size_t LCS_k(std::vector<Substring> &substring,
             char_t const* const     reference,
             char_t const* const     complement,
             size_t const            reference_start,
             size_t const            reference_end,
             char_t const* const     sample,
             size_t const            sample_start,
             size_t const            sample_end,
             size_t const            k)
{
  return 0;
} // LCS_k
*/


size_t prefix_match(char_t const* const reference,
                    size_t const        reference_length,
                    char_t const* const sample,
                    size_t const        sample_length)
{
  size_t i = 0;
  while (i < reference_length && i < sample_length && reference[i] == sample[i])
  {
    ++i;
  } // while
  return i;
} // prefix_match

size_t suffix_match(char_t const* const reference,
                    size_t const        reference_length,
                    char_t const* const sample,
                    size_t const        sample_length,
                    size_t const        prefix)

{
  size_t i = 0;
  while (i < reference_length - prefix && i < sample_length - prefix && reference[reference_length - i - 1] == sample[sample_length - i - 1])
  {
    ++i;
  } // while
  return i;
} // suffix_match

char_t IUPAC_base_complement(char_t const base)
{
  switch (base)
  {
    case 'A':
      return 'T';
    case 'C':
      return 'G';
    case 'G':
      return 'C';
    case 'T':
    case 'U':
      return 'A';
  } // switch
  return base;
} // IUPAC_base_complement

char_t const* IUPAC_complement(char_t const* const string,
                               size_t const        length)
{
  char_t* complement = new char_t[length];
  for (size_t i = 0; i < length; ++i)
  {
    complement[i] = IUPAC_base_complement(string[i]);
  } // for
  return complement;
} // IUPAC_complement


#if defined(__debug__)
size_t Dprint_truncated(char_t const* const string,
                        size_t const        start,
                        size_t const        end,
                        size_t const        length,
                        FILE*               stream)
{
  if (end - start <= length - 3)
  {
    return fwrite(string + start, sizeof(char_t), end - start, stream);
  } // if
  return fwrite(string + start, sizeof(char_t), length / 2 - 1, stream) +
         fputs("...", stream) +
         fwrite(string + end - (length / 2) + 1, sizeof(char_t), length / 2 - 1, stream);
} // Dprint_truncated
#endif

} // mutalyzer

