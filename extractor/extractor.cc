// *******************************************************************
//   (C) Copyright 2014 Leiden Institute of Advanced Computer Science
//   Universiteit Leiden
//   All Rights Reserved
// *******************************************************************
// Extractor (library)
// *******************************************************************
// FILE INFORMATION:
//   File:     extractor.cc (depends on extractor.h)
//   Author:   Jonathan K. Vis
//   Revision: 2.01a
//   Date:     2014/07/28
// *******************************************************************
// DESCRIPTION:
//   This library can be used to generete HGVS variant descriptions as
//   accepted by the Mutalyzer Name Checker.
// *******************************************************************

#include "extractor.h"

namespace mutalyzer
{

size_t global_reference_length = 0;

size_t weight_position = 0;

std::vector<Variant> extract(char_t const* const reference,
                             size_t const        reference_length,
                             char_t const* const sample,
                             size_t const        sample_length,
                             int const           type)
{
  std::vector<Variant> variant;
  extract(variant, reference, reference_length, sample, sample_length, type);
  return variant;
} // extract

size_t extract(std::vector<Variant> &variant,
               char_t const* const   reference,
               size_t const          reference_length,
               char_t const* const   sample,
               size_t const          sample_length,
               int const             type)
{
  global_reference_length = reference_length;
  weight_position = ceil(log10(reference_length / 4));

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
                 size_t const          sample_end)
{
  size_t const reference_length = reference_end - reference_start;
  size_t const sample_length = sample_end - sample_start;
  size_t const weight_trivial = weight_position + WEIGHT_DELETION_INSERTION + WEIGHT_BASE * sample_length + (reference_length != 1 ? weight_position + WEIGHT_SEPARATOR : 0);
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
      weight = 2 * weight_position + WEIGHT_SEPARATOR + WEIGHT_INSERTION + WEIGHT_BASE * sample_length;

      std::vector<Variant> transposition;
      size_t const weight_transposition = extractor_transposition(transposition, reference, complement, reference_start, reference_end, sample, sample_start, sample_end, weight) + 2 * weight_position + 3 * WEIGHT_SEPARATOR + WEIGHT_INSERTION;


#if defined(__debug__)
  fprintf(stderr, "Transpositions: %ld (trivial: %ld)\n", weight_transposition, weight);
  for (std::vector<Variant>::iterator it = transposition.begin(); it != transposition.end(); ++it)
  {
    fprintf(stderr, "  %ld--%ld, %ld--%ld, %d, %ld, %ld--%ld\n", it->reference_start, it->reference_end, it->sample_start, it->sample_end, it->type, it->weight, it->transposition_start, it->transposition_end);
  } // for
#endif


      if (weight > weight_transposition && transposition.size() > 0 && !(transposition.size() == 1 && transposition.front().type == SUBSTITUTION))
      {
        transposition.front().type |= TRANSPOSITION_OPEN;
        transposition.back().type |= TRANSPOSITION_CLOSE;
        variant.insert(variant.end(), transposition.begin(), transposition.end());
        return weight_transposition;
      } // if

      // simple insertion
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


  std::vector<Substring> substring;
  size_t const length = LCS(substring, reference, complement, reference_start, reference_end, sample, sample_start, sample_end);

  // deletion/insertion
  if (length <= 0 || substring.size() <= 0)
  {
    weight = weight_trivial;

    std::vector<Variant> transposition;
    size_t const weight_transposition = extractor_transposition(transposition, reference, complement, reference_start, reference_end, sample, sample_start, sample_end, weight)  + 2 * weight_position + 3 * WEIGHT_SEPARATOR + WEIGHT_INSERTION;


#if defined(__debug__)
  fprintf(stderr, "Transpositions: %ld (trivial: %ld)\n", weight_transposition, weight);
  for (std::vector<Variant>::iterator it = transposition.begin(); it != transposition.end(); ++it)
  {
    fprintf(stderr, "  %ld--%ld, %ld--%ld, %d, %ld, %ld--%ld\n", it->reference_start, it->reference_end, it->sample_start, it->sample_end, it->type, it->weight, it->transposition_start, it->transposition_end);
  } // for
#endif


    if (weight > weight_transposition && transposition.size() > 0 && !(transposition.size() == 1 && transposition.front().type == SUBSTITUTION))
    {
      transposition.front().type |= TRANSPOSITION_OPEN;
      transposition.back().type |= TRANSPOSITION_CLOSE;
      variant.insert(variant.end(), transposition.begin(), transposition.end());
      return weight_transposition;
    } // if

    variant.push_back(Variant(reference_start, reference_end, sample_start, sample_end, SUBSTITUTION, weight_trivial));
    return weight_trivial;
  } // if


  size_t diff = (reference_end - reference_start) + (sample_end - sample_start);
  std::vector<Substring>::iterator lcs = substring.begin();
  for (std::vector<Substring>::iterator it = substring.begin(); it != substring.end(); ++it)
  {
    size_t const prefix_diff = abs((it->reference_index - reference_start) - (it->sample_index - sample_start));
    size_t const suffix_diff = abs((reference_end - (it->reference_index + it->length)) - (sample_end - (it->sample_index + it->length)));
    if (prefix_diff + suffix_diff < diff)
    {
      diff = prefix_diff + suffix_diff;
      lcs = it;
    } // if
  } // for

  if (lcs->reverse_complement)
  {
    weight = 2 * weight_position + WEIGHT_SEPARATOR + WEIGHT_INVERSION;
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
    fputs("\n", stderr);
  } // for
#endif


  std::vector<Variant> prefix;
  weight += extractor(prefix, reference, complement, reference_start, lcs->reference_index, sample, sample_start, lcs->sample_index);

  if (weight > weight_trivial)
  {
    variant.push_back(Variant(reference_start, reference_end, sample_start, sample_end, SUBSTITUTION, weight_trivial));
    return weight_trivial;
  } // if

  std::vector<Variant> suffix;
  weight += extractor(suffix, reference, complement, lcs->reference_index + length, reference_end, sample, lcs->sample_index + length, sample_end);

  if (weight > weight_trivial)
  {
    variant.push_back(Variant(reference_start, reference_end, sample_start, sample_end, SUBSTITUTION, weight_trivial));
    return weight_trivial;
  } // if


  variant.insert(variant.end(), prefix.begin(), prefix.end());

  if (!lcs->reverse_complement)
  {
    variant.push_back(Variant(lcs->reference_index, lcs->reference_index + length, lcs->sample_index, lcs->sample_index + length));
  } // if
  else
  {
    variant.push_back(Variant(lcs->reference_index, lcs->reference_index + length, lcs->sample_index, lcs->sample_index + length, REVERSE_COMPLEMENT, 2 * weight_position + WEIGHT_SEPARATOR + WEIGHT_INVERSION));
  } // else

  variant.insert(variant.end(), suffix.begin(), suffix.end());

  return weight;
} // extractor

size_t extractor_transposition(std::vector<Variant> &variant,
                               char_t const* const   reference,
                               char_t const* const   complement,
                               size_t const          reference_start,
                               size_t const          reference_end,
                               char_t const* const   sample,
                               size_t const          sample_start,
                               size_t const          sample_end,
                               size_t const          weight_trivial)
{
  size_t const sample_length = sample_end - sample_start;

  size_t weight = 0;


#if defined(__debug__)
  fputs("Transposition extraction\n", stderr);
  fprintf(stderr, "  reference %ld--%ld:  ", 0ul, global_reference_length);
  Dprint_truncated(reference, 0, global_reference_length);
  fprintf(stderr, " (%ld)\n", global_reference_length);
  fprintf(stderr, "  complement %ld--%ld: ", 0ul, global_reference_length);
  Dprint_truncated(complement, 0, global_reference_length);
  fprintf(stderr, " (%ld)\n", global_reference_length);
  fprintf(stderr, "  sample %ld--%ld:     ", sample_start, sample_end);
  Dprint_truncated(sample, sample_start, sample_end);
  fprintf(stderr, " (%ld)\n", sample_length);
  fprintf(stderr, "  trivial weight: %ld\n", weight_trivial);
#endif


  if (sample_length <= 2 * weight_position)
  {
    return sample_length * WEIGHT_BASE;
  } // if

  std::vector<Substring> substring;
  size_t const length = LCS(substring, reference, complement, 0, global_reference_length, sample, sample_start, sample_end);

  if (length <= 0 || substring.size() <= 0)
  {
    weight = sample_length * WEIGHT_BASE;
    variant.push_back(Variant(reference_start, reference_end, sample_start, sample_end, SUBSTITUTION, weight));
    return weight;
  } // if

  std::vector<Substring>::iterator const lcs = substring.begin();

  weight += 2 * weight_position + WEIGHT_SEPARATOR;
  if (lcs->reverse_complement)
  {
    weight += WEIGHT_INVERSION;
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
    fputs("\n", stderr);
  } // for
#endif


  std::vector<Variant> prefix;
  weight += extractor_transposition(prefix, reference, complement, reference_start, reference_end, sample, sample_start, lcs->sample_index, weight_trivial - weight);

  if (weight > weight_trivial)
  {
    weight = sample_length * WEIGHT_BASE;
    variant.push_back(Variant(reference_start, reference_end, sample_start, sample_end, SUBSTITUTION, weight));
    return weight;
  } // if

  std::vector<Variant> suffix;
  weight += extractor_transposition(suffix, reference, complement, reference_start, reference_end, sample, lcs->sample_index + length, sample_end, weight_trivial - weight);

  if (weight > weight_trivial)
  {
    weight = sample_length * WEIGHT_BASE;
    variant.push_back(Variant(reference_start, reference_end, sample_start, sample_end, SUBSTITUTION, weight));
    return weight;
  } // if

  variant.insert(variant.end(), prefix.begin(), prefix.end());

  if (!lcs->reverse_complement)
  {
    variant.push_back(Variant(reference_start, reference_end, lcs->sample_index, lcs->sample_index + length, IDENTITY, 2 * weight_position + WEIGHT_SEPARATOR, lcs->reference_index, lcs->reference_index + length));
  } // if
  else
  {
    variant.push_back(Variant(reference_start, reference_end, lcs->sample_index, lcs->sample_index + length, REVERSE_COMPLEMENT, 2 * weight_position + WEIGHT_SEPARATOR + WEIGHT_INVERSION, lcs->reference_index, lcs->reference_index + length));
  } // else

  variant.insert(variant.end(), suffix.begin(), suffix.end());

  return weight + variant.size() - 1;
} // extractor_transposition


size_t LCS(std::vector<Substring> &substring,
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

  static size_t const THRESHOLD = 16000;

  double const a = reference_length >= sample_length ? reference_length : sample_length;
  double const b = reference_length >= sample_length ? sample_length : reference_length;
  size_t const cut_off = (reference_length > THRESHOLD ? ceil((1.0 - b / (a + 0.1 * b)) * b) / 8 : 0) + 1;

  size_t k = reference_length > sample_length ? sample_length / 4 : reference_length / 4;

  while (k > 4 && k > cut_off)
  {


#if defined(__debug__)
  fprintf(stderr, "  k = %ld (cut-off: %ld)\n", k, cut_off);
#endif


    substring.clear();

    size_t const length = LCS_k(substring, reference, complement, reference_start, reference_end, sample, sample_start, sample_end, k);

    // A LCS of sufficient length has been found.
    if (length >= 2 * k && substring.size() > 0)
    {
      return length;
    } // if
    k /= 3;
  } // while

  if (cut_off > 1)
  {


#if defined(__debug__)
  fprintf(stderr, "  cut-off\n");
#endif


    substring.clear();
    return 0;
  } // if


#if defined(__debug__)
  fprintf(stderr, "  k = 1\n");
#endif


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
  bool reverse_complement = false;

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

        if (LCS_line[i % 2][j] >= length)
        {
          if (reverse_complement || LCS_line[i % 2][j] > length)
          {
            length = LCS_line[i % 2][j];
            substring = std::vector<Substring>(1, Substring(j - length + reference_start + 1, i - length + sample_start + 1, length));
          } // if
          else
          {
            substring.push_back(Substring(j - length + reference_start + 1, i - length + sample_start + 1, length));
          } // else
          reverse_complement = false;
        } // if
      } // if
      else
      {
        LCS_line[i % 2][j] = 0;
      } // else

      // If applicable check for a LCS in reverse complement space.
      // The same code is used as before but the complement string is
      // travesed backwards (towards the start).
      if (complement != 0 && complement[reference_end - j - 1] == sample[sample_start + i])
      {
        if (i == 0 || j == 0)
        {
          LCS_line_rc[i % 2][j] = 1;
        } // if
        else
        {
          LCS_line_rc[i % 2][j] = LCS_line_rc[(i + 1) % 2][j - 1] + 1;
        } // else

        if (LCS_line_rc[i % 2][j] > 1 && LCS_line_rc[i % 2][j] > length)
        {
          length = LCS_line_rc[i % 2][j];
          substring = std::vector<Substring>(1, Substring(reference_end - j - 1, i - length + sample_start + 1, length, true));
          reverse_complement = true;
        } // if
      } // if
      else
      {
        LCS_line_rc[i % 2][j] = 0;
      } // else

      if (!reverse_complement && length >= sample_length)
      {
        break;
      } // if

    } // for
  } // for

  delete[] &LCS_line;
  delete[] &LCS_line_rc;

  return length;
} // LCS_1

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
  size_t length = 0;

  // Stop if we cannot partition the strings into k-mers.
  if (k <= 1 || reference_end - reference_start < k || sample_end - sample_start < k)
  {
    return length;
  } // if

  size_t const reference_length = (reference_end - reference_start) / k;
  size_t const sample_length = sample_end - sample_start - k + 1;
  bool reverse_complement = false;

  // Just a fancy way of allocation a continuous (k+1)D array in heap
  // space.
  typedef size_t array[k + 1][reference_length];
  array &LCS_line = *(reinterpret_cast<array*>(new size_t[(k + 1) * reference_length]));
  array &LCS_line_rc = *(reinterpret_cast<array*>(new size_t[(k + 1) * reference_length]));

  // Filling the LCS k-mer matrix (actually only the current and the k
  // previous rows). We count in k-mers.
  for (size_t i = 0; i < sample_length; ++i) // overlapping k-mers
  {
    for (size_t j = 0; j < reference_length; ++j) // non-overlapping
    {
      // A match
      if (string_match(reference + reference_start + j * k, sample + sample_start + i, k))
      {
        if (i < k || j == 0)
        {
          LCS_line[i % (k + 1)][j] = 1;
        } // if
        else
        {
          LCS_line[i % (k + 1)][j] = LCS_line[(i + 1) % (k + 1)][j - 1] + 1;
        } // else

        // Check for a new maximal length.
        if (LCS_line[i % (k + 1)][j] > length)
        {
          length = LCS_line[i % (k + 1)][j];

          // Remove all solutions with a length more than 1 from the
          // new maximal length. And remove also the partial LCS that
          // was extended to the new maximal length because it is
          // guaranteerd to be of maximal length - 1.
          for (std::vector<Substring>::iterator it = substring.begin(); it != substring.end(); ++it)
          {
            if (length - it->length > 1 || (it->reference_index == j - 1 && it->sample_index == i - k && !it->reverse_complement))
            {
              std::vector<Substring>::iterator const temp = it - 1;
              substring.erase(it);
              it = temp;
            } // if
          } // for
          substring.push_back(Substring(j, i, LCS_line[i % (k + 1)][j]));
        } // if
        else if (LCS_line[i % (k + 1)][j] > 0 && length - LCS_line[i % (k + 1)][j] <= 1)
        {
          substring.push_back(Substring(j, i, LCS_line[i % (k + 1)][j]));
        } // if
      } // if
      else
      {
        LCS_line[i % (k + 1)][j] = 0;
      } // else

      if (complement != 0 && string_match_reverse(complement + reference_end - j * k - 1, sample + sample_start + i, k))
      {
        if (i < k || j == 0)
        {
          LCS_line_rc[i % (k + 1)][j] = 1;
        } // if
        else
        {
          LCS_line_rc[i % (k + 1)][j] = LCS_line_rc[(i + 1) % (k + 1)][j - 1] + 1;
        } // else

        // Check for a new maximal length.
        if (LCS_line_rc[i % (k + 1)][j] > length)
        {
          length = LCS_line_rc[i % (k + 1)][j];

          // Remove all solutions with a length more than 1 from the
          // new maximal length. And remove also the partial LCS that
          // was extended to the new maximal length because it is
          // guaranteerd to be of maximal length - 1.
          for (std::vector<Substring>::iterator it = substring.begin(); it != substring.end(); ++it)
          {
            if (length - it->length > 1 || (it->reference_index == j - 1 && it->sample_index == i - k && it->reverse_complement))
            {
              std::vector<Substring>::iterator const temp = it - 1;
              substring.erase(it);
              it = temp;
            } // if
          } // for
          substring.push_back(Substring(j, i, LCS_line_rc[i % (k + 1)][j], true));
        } // if
        else if (LCS_line[i % (k + 1)][j] > 0 && length - LCS_line_rc[i % (k + 1)][j] <= 1)
        {
          substring.push_back(Substring(j, i, LCS_line[i % (k + 1)][j], true));
        } // if
      } // if
      else
      {
        LCS_line_rc[i % (k + 1)][j] = 0;
      } // else

    } // for
  } // for


  // extending LCS
  for (std::vector<Substring>::iterator it = substring.begin(); it != substring.end(); ++it)
  {
    if (!it->reverse_complement)
    {
      it->reference_index = reference_start + (it->reference_index - it->length + 1) * k;
      it->sample_index = sample_start + it->sample_index - (it->length - 1) * k;
      it->length *= k;
      // extending to the right
      {
        size_t i = 0;
        while (i <= k && it->reference_index + it->length + i < reference_end && it->sample_index + it->length + i < sample_end && reference[it->reference_index + it->length + i] == sample[it->sample_index + it->length + i])
        {
          ++i;
        } // while
        it->length += i;
      }
      // extending to the left
      {
        size_t i = 0;
        while (i <= k && it->reference_index - i - 1 >= reference_start && it->sample_index - i - 1 >= sample_start && reference[it->reference_index - i - 1] == sample[it->sample_index - i - 1])
        {
          ++i;
        } // while
        it->reference_index -= i;
        it->sample_index -= i;
        it->length += i;
      }
    } // if
    else
    {
      it->reference_index = reference_end - (it->reference_index + 1) * k;
      it->sample_index = sample_start + it->sample_index - (it->length - 1) * k;
      it->length *= k;
      // extending to the right (sample orientation)
      {
        size_t i = 0;
        while (i <= k && it->reference_index - i - 1 >= reference_start && it->sample_index + it->length + i < sample_end && complement[it->reference_index - i - 1] == sample[it->sample_index + it->length + i])
        {
          ++i;
        } // while
        it->reference_index -= i;
        it->length += i;
      }
      // extending to the left (sample orientation)
      {
        size_t i = 0;
        while (i <= k && it->reference_index + it->length + i < reference_end && it->sample_index - i - 1 >= sample_start && complement[it->reference_index + it->length + i] == sample[it->sample_index - i - 1])
        {
          ++i;
        } // while
        it->sample_index -= i;
        it->length += i;
      }
    } // else

    if (it->length > length)
    {
      length = it->length;
      reverse_complement = it->reverse_complement;
    } // if
    else if (reverse_complement && it->length == length)
    {
      reverse_complement = false;
    } // if
  } // for


  for (std::vector<Substring>::iterator it = substring.begin(); it != substring.end(); ++it)
  {
    if (it->length < length || reverse_complement != it->reverse_complement)
    {
      std::vector<Substring>::iterator const temp = it - 1;
      substring.erase(it);
      it = temp;
    } // if
  } // for

  delete[] &LCS_line;
  delete[] &LCS_line_rc;

  return length;
} // LCS_k


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

bool string_match_reverse(char_t const* const string_1,
                          char_t const* const string_2,
                          size_t const        length)
{
  for (size_t i = 0; i < length; ++i)
  {
    if (string_1[-i] != string_2[i])
    {
      return false;
    } // if
  } // for
  return true;
} // string_match_reverse


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

