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
//   Revision: 1.05a
//   Date:     2014/02/05
// *******************************************************************
// DESCRIPTION:
//   This library can be used to generete HGVS variant descriptions as
//   accepted by the Mutalyzer Name Checker.
// *******************************************************************

#include "extractor.h"

#include <cstdlib>

#if defined(__debug__)
#include <cstdio>
#endif

namespace mutalyzer
{

#if defined(__debug__)
void print(char const* const          reference,
           char const* const          sample,
           std::vector<Variant> const result)
{
  for (size_t i = 0; i < result.size(); ++i)
  {
    if (result[i].type == VARIANT_IDENTITY)
    {
      printf("%ld_%ldidem%ld_%ld\n", result[i].reference_start + 1, result[i].reference_end, result[i].sample_start + 1, result[i].sample_end);
    } // if
    else if (result[i].type == VARIANT_REVERSE_COMPLEMENT)
    {
      printf("%ld_%ldinv\n", result[i].reference_start + 1, result[i].reference_end);
    } // if
    else if (result[i].type == VARIANT_TRANSPOSITION)
    {
      printf("%ld_%ldins%ld_%ld\n", result[i].reference_start + 1, result[i].reference_end, result[i].sample_start, result[i].sample_end);
    } // if
    else if (result[i].type == VARIANT_REVERSE_COMPLEMENT_TRANSPOSITION)
    {
      printf("%ld_%ldinso%ld_%ld\n", result[i].reference_start + 1, result[i].reference_end, result[i].sample_start, result[i].sample_end);
    } // if
    else
    {
      size_t const reference_length = result[i].reference_end - result[i].reference_start;
      size_t const sample_length = result[i].sample_end - result[i].sample_start;
      if (reference_length == 0)
      {
        fprintf(stderr, "%ld_%ldins", result[i].reference_start, result[i].reference_start + 1);
        fwrite(sample + result[i].sample_start, sizeof(char), sample_length, stderr);
        fputc('\n', stderr);
      } // if
      else if (sample_length == 0)
      {
        if (reference_length == 1)
        {
          fprintf(stderr, "%lddel\n", result[i].reference_start + 1);
        } // if
        else
        {
          fprintf(stderr, "%ld_%lddel\n", result[i].reference_start + 1, result[i].reference_end);
        } // else
      } // if
      else if (reference_length == 1 && sample_length == 1)
      {
        fprintf(stderr, "%ld%c>%c\n", result[i].reference_start + 1, reference[result[i].reference_start], sample[result[i].sample_start]);
      } // if
      else
      {
        if (reference_length == 1)
        {
          fprintf(stderr, "%lddelins", result[i].reference_start + 1);
        } // if
        else
        {
          fprintf(stderr, "%ld_%lddelins", result[i].reference_start + 1, result[i].reference_end);
        } // else
        fwrite(sample + result[i].sample_start, sizeof(char), sample_length, stderr);
        fprintf(stderr, "  (%ld--%ld)\n", result[i].sample_start, result[i].sample_end);
      } // else
    } // else
  } // for
  return;
} // print
#endif

// This global variable is a dirty trick used to always have access to
// the complete reference strings even when deep into the recursion.
// This seems necessary to compute transpositions.
size_t global_reference_length = 0;

// This function is more or less equivalent to C's strncmp, but it
// returns true iff both strings are the same.
static inline bool string_match(char const* const string_1,
                                char const* const string_2,
                                size_t const      n)
{
  for (size_t i = 0; i < n; ++i)
  {
    if (string_1[i] != string_2[i])
    {
      return false;
    } // if
  } // for
  return true;
} // string_match

// This function is very similar to C's strncmp, but it traverses
// string_1 from end to start while traversing string_2 from start to
// end (useful for the reverse complement in DNA/RNA), and it returns
// true iff both strings are the same in their respective directions.
static inline bool string_match_reverse(char const* const string_1,
                                        char const* const string_2,
                                        size_t const      n)
{
  for (size_t i = 0; i < n; ++i)
  {
    if (string_1[-i] != string_2[i])
    {
      return false;
    } // if
  } // for
  return true;
} // string_match_reverse

// This function calculates the length (in characters) of the shared
// prefix between two strings. The result of this function is also
// used in the suffix_match function.
static inline size_t prefix_match(char const* const reference,
                                  size_t const      reference_length,
                                  char const* const sample,
                                  size_t const      sample_length)
{
  size_t result = 0;

  // Traverse both strings towards the end as long as their characters
  // are equal. Do not exceed the length of the strings.
  while (result < reference_length && result < sample_length && reference[result] == sample[result])
  {
    ++result;
  } // while
  return result;
} // prefix_match

// This function calculates the length (in characters) of the shared
// suffix between two strings. It needs the calculated shared prefix.
static inline size_t suffix_match(char const* const reference,
                                  size_t const      reference_length,
                                  char const* const sample,
                                  size_t const      sample_length,
                                  size_t const      prefix = 0)
{
  size_t result = 0;

  // Start at the end of both strings and traverse towards the start
  // as long as their characters are equal. Do not exceed the length
  // of the strings.
  while (result < reference_length - prefix && result < sample_length - prefix && reference[reference_length - result - 1] == sample[sample_length - result - 1])
  {
    ++result;
  } // while
  return result;
} // suffix_match

// The main library function. Extract all variants (regions of change)
// from the given strings.
std::vector<Variant> extract(char const* const reference,
                             size_t const      reference_length,
                             char const* const sample,
                             size_t const      sample_length,
                             int const         type)
{
#if defined(__debug__)
  fprintf(stderr, "extractor.cc --- prefix removal\n");
#endif
  // First remove the common prefix and suffix (in that order!).
  size_t const prefix = prefix_match(reference, reference_length, sample, sample_length);
#if defined(__debug__)
  fprintf(stderr, "  prefix: %ld\n", prefix);
#endif

#if defined(__debug__)
  fprintf(stderr, "extractor.cc --- suffix removal\n");
#endif
  size_t const suffix = suffix_match(reference, reference_length, sample, sample_length, prefix);
#if defined(__debug__)
  fprintf(stderr, "  suffix: %ld\n", suffix);
#endif

  // Used to always have access to the complete reference string(s).
  global_reference_length = reference_length;

#if defined(__debug__)
  fprintf(stderr, "extractor.cc --- constructing IUPAC complement\n");
#endif
  // Do NOT construct a complement string for protein strings. All
  // other string types default to DNA/RNA.
  char const* const complement = type != TYPE_PROTEIN ? IUPAC_complement(reference, reference_length) : 0;

  std::vector<Variant> result;
  if (prefix > 0)
  {
    result.push_back(Variant(0, prefix, 0, prefix, VARIANT_IDENTITY));
  } // if

#if defined(__debug__)
  fprintf(stderr, "extractor.cc --- starting variant extraction\n");
#endif
  extractor(reference, complement, prefix, reference_length - suffix, sample, prefix, sample_length - suffix, result);
  if (suffix > 0)
  {
    result.push_back(Variant(reference_length - suffix, reference_length, sample_length - suffix, sample_length, VARIANT_IDENTITY));
  } // if

#if defined(__debug__)
  fprintf(stderr, "extractor.cc --- variant extraction done, cleaning up\n");
#endif
  // do NOT forget to clean up the complement string.
  if (complement != 0)
  {
    delete[] complement;
  } // if
  return result;
} // extract

// This is the recursive extractor function. It works as follows:
// First, determine the ``best fitting'' longest common substring
// (LCS) (possibly as a reverse complement) and discard it from the
// solution. Then apply the same function on the remaining prefixes
// and suffixes. If there is no LCS it is a variant (region of
// change), i.e., an insertion/deletion.
// With regard to the reverse complement: the complement string is, as
// its name suggests, just the complement (DNA/RNA) of the reference
// string but it is NOT reversed.
void extractor(char const* const     reference,
               char const* const     complement,
               size_t const          reference_start,
               size_t const          reference_end,
               char const* const     sample,
               size_t const          sample_start,
               size_t const          sample_end,
               std::vector<Variant> &result)
{
#if defined(__debug__)
  fprintf(stderr, "extractor.cc --- extractor\n  reference: %ld--%ld (%ld)\n  sample:    %ld--%ld (%ld)\n", reference_start, reference_end, reference_end - reference_start, sample_start, sample_end, sample_end - sample_start);
#endif

  // First some base cases to end the recursion.
  // No more reference string.
  if (reference_end - reference_start <= 0)
  {
    // But some of the sample string is remaining: this is an
    // insertion.
    if (sample_end - sample_start > 0)
    {
      // First, we check if we can match the inserted substring
      // somewhere in the complete reference string. This will
      // indicate a possible transposition.
      if (sample_end - sample_start > 64)
      {
        std::vector<Variant> transposition;
#if defined(__debug__)
  fprintf(stderr, "extractor.cc --- starting transposition extraction\n");
#endif
        extractor(reference, complement, 0, global_reference_length, sample, sample_start, sample_end, transposition);
#if defined(__debug__)
  fprintf(stderr, "  transpositions: %ld\n", transposition.size());
  //print(reference, sample, std::vector<Variant>(1, Variant(reference_start, reference_end, sample_start, sample_end)));
  print(reference, sample, transposition);
  fprintf(stderr, "  %ld--%ld\n", reference_start, reference_end);
#endif

        size_t start = reference_start;
        // This variant can be described as a transposition
        for (size_t i = 0; i < transposition.size(); ++i)
        {
          if (transposition[i].type == VARIANT_IDENTITY)
          {
            result.push_back(Variant(start - 1, start + 1, transposition[i].reference_start, transposition[i].reference_end, VARIANT_TRANSPOSITION));
          } // if
          else if (transposition[i].type == VARIANT_REVERSE_COMPLEMENT)
          {
            result.push_back(Variant(start - 1, start + 1, transposition[i].reference_start, transposition[i].reference_end, VARIANT_REVERSE_COMPLEMENT_TRANSPOSITION));
          } // if
          else
          {
            size_t const reference_length = transposition[i].reference_end - transposition[i].reference_start;
            size_t const sample_length = transposition[i].sample_end - transposition[i].sample_start;
            if (reference_length > 0 && sample_length > 0)
            {
              result.push_back(Variant(start, start, transposition[i].sample_start, transposition[i].sample_end, transposition[i].type));
            } // if
            else if (sample_length > 0)
            {
              result.push_back(transposition[i]);
            } // else
          } // else
          ++start;
        } // for
        return;
      } // if

      // Insertion considered to be too short for transposition extraction
      result.push_back(Variant(reference_start, reference_end, sample_start, sample_end));
    } // if
    return;
  } // if

  // Obviously there is a piece of reference string left, but no more
  // sample string: this is a deletion.
  if (sample_end - sample_start <= 0)
  {
    result.push_back(Variant(reference_start, reference_end, sample_start, sample_end));
    return;
  } // if

#if defined(__debug__)
  fprintf(stderr, "extractor.cc --- calculating LCS\n");
#endif
  // Calculate the LCS (possibly in reverse complement) of the two
  // strings.
  std::vector<Substring> LCS_result = LCS(reference, complement, reference_start, reference_end, sample, sample_start, sample_end);

  // No LCS found: this is an insertion/deletion or a substitution.
  if (LCS_result.size() == 0)
  {
    result.push_back(Variant(reference_start, reference_end, sample_start, sample_end));
    return;
  } // if

#if defined(__debug__)
  fprintf(stderr, "  LCS length: %ld\n", LCS_result[0].length);
#endif

  // Pick the ``best fitting'' LCS, i.e., the location of the LCS
  // within their respective strings is close.
  // FIXME: we could extract all non-overlapping LCSs in one go
  size_t difference = (reference_end - reference_start) + (sample_end - sample_start);
  size_t index = 0;
  for (size_t i = 0; i < LCS_result.size(); ++i)
  {
    if (abs(LCS_result[i].reference_index - LCS_result[i].sample_index) < difference)
    {
      difference = abs(LCS_result[i].reference_index - LCS_result[i].sample_index);
      index = i;
    } // if
  } // for

  // Apply this function to the prefixes of the strings.
  extractor(reference, complement, reference_start, LCS_result[index].reference_index, sample, sample_start, LCS_result[index].sample_index, result);

  // We always add a full description including the region that did
  // not change.
  if (!LCS_result[index].reverse_complement)
  {
    result.push_back(Variant(LCS_result[index].reference_index, LCS_result[index].reference_index + LCS_result[index].length, LCS_result[index].sample_index, LCS_result[index].sample_index + LCS_result[index].length, VARIANT_IDENTITY));
  } // if
  else
  {
    result.push_back(Variant(LCS_result[index].reference_index, LCS_result[index].reference_index + LCS_result[index].length, LCS_result[index].sample_index, LCS_result[index].sample_index + LCS_result[index].length, VARIANT_REVERSE_COMPLEMENT));
  } // else
  
  // Apply this function to the suffixes of the strings.
  extractor(reference, complement, LCS_result[index].reference_index + LCS_result[index].length, reference_end, sample, LCS_result[index].sample_index + LCS_result[index].length, sample_end, result);
  return;
} // extractor

// Calculate the LCS in the well-known way using dynamic programming.
// NOT suitable for large strings.
std::vector<Substring> LCS_1(char const* const reference,
                             char const* const complement,
                             size_t const      reference_start,
                             size_t const      reference_end,
                             char const* const sample,
                             size_t const      sample_start,
                             size_t const      sample_end)
{
  size_t const reference_length = reference_end - reference_start;
  size_t const sample_length = sample_end - sample_start;

  // Just a fancy way of allocation a continuous 2D array in heap
  // space.
  typedef size_t array[2][reference_length];
  array &LCS_line = *(reinterpret_cast<array*>(new size_t[2 * reference_length]));
  array &LCS_line_rc = *(reinterpret_cast<array*>(new size_t[2 * reference_length]));

  std::vector<Substring> result;
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
          result = std::vector<Substring>(1, Substring(j - length + reference_start + 1, i - length + sample_start + 1, length));
        } // if
        // Found a LCS of the same maximal length.
        else if (LCS_line[i % 2][j] == length)
        {
          result.push_back(Substring(j - length + reference_start + 1, i - length + sample_start + 1, length));
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
        if (LCS_line_rc[i % 2][j] > length)
        {
          length = LCS_line_rc[i % 2][j];
          result = std::vector<Substring>(1, Substring(reference_end - j - 1, i - length + sample_start + 1, length, true));
        } // if
        else if (LCS_line_rc[i % 2][j] == length)
        {
          result.push_back(Substring(reference_end - j - 1, i - length + sample_start + 1, length, true));
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

  return result;
} // LCS_1

// Calculate the LCS using overlapping and non-overlapping k-mers.
// This function should be suitable for large (similar) strings.
// Be careful: if the resulting LCS is of length <= 2k it might not be
// the actual LCS. Remedy: try again with a reduced k.
std::vector<Substring> LCS_k(char const* const reference,
                             char const* const complement,
                             size_t const      reference_start,
                             size_t const      reference_end,
                             char const* const sample,
                             size_t const      sample_start,
                             size_t const      sample_end,
                             size_t const      k)
{
  size_t const reference_length = (reference_end - reference_start) / k;

  // Stop if we cannot partition the sample string into overlapping
  // k-mers.
  if (sample_end - sample_start < k)
  {
    return std::vector<Substring>();
  } // if
  size_t const sample_length = sample_end - sample_start - k + 1;

  // Just a fancy way of allocation a continuous (k+1)D array in heap
  // space.
  typedef size_t array[k + 1][reference_length];
  array &LCS_line = *(reinterpret_cast<array*>(new size_t[(k + 1) * reference_length]));
  array &LCS_line_rc = *(reinterpret_cast<array*>(new size_t[(k + 1) * reference_length]));

  std::vector<Substring> result;
  size_t length = 0;

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
          for (size_t e = 0; e < result.size(); ++e)
          {
            if (length - result[e].length > 1 || (result[e].reference_index == j - 1 && result[e].sample_index == i - k))
            {
              result.erase(result.begin() + e);
              --e;
            } // if
          } // for

          // The actual indices are calculated afterwards (this is in
          // k-mer space).
          result.push_back(Substring(j, i, length));
        } // if
        // Found a LCS of the same maximal length or maximal
        // length - 1.
        else if (length - LCS_line[i % (k + 1)][j] < 1)
        {
          result.push_back(Substring(j, i, length));
        } // if
      } // if
      else
      {
        LCS_line[i % (k + 1)][j] = 0;
      } // else

      // If applicable check for a LCS in reverse complement space.
      // The same code is used as before but the complement string is
      // travesed backwards (towards the start).
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
        if (LCS_line_rc[i % (k + 1)][j] > length)
        {
          length = LCS_line_rc[i % (k + 1)][j];
          for (size_t e = 0; e < result.size(); ++e)
          {
            if (length - result[e].length > 1 || (result[e].reference_index == j - 1 && result[e].sample_index == i - k))
            {
              result.erase(result.begin() + e);
              --e;
            } // if
          } // for
          result.push_back(Substring(j, i, length, true));
        } // if
        else if (length - LCS_line_rc[i % (k + 1)][j] < 1)
        {
          result.push_back(Substring(j, i, length, true));
        } // if
      } // if
      else
      {
        LCS_line_rc[i % (k + 1)][j] = 0;
      } // else
    } // for
  } // for

  // Now we need to do some extensions of the found LCSs to find their
  // exact lengths. We extend up to k positions to the left (towards
  // the start of the strings) and to the right (towards the end of
  // the string).
  length *= k;
  for (size_t i = 0; i < result.size(); ++i)
  {

    // In reverse complement space matters are complicated: in the
    // complement string we have to extend towards the other
    // direction.
    if (result[i].reverse_complement)
    {
      // The correct positions of the LCS within the strings. Here the
      // indices are updated to their actual value.
      result[i].reference_index = reference_end - (result[i].reference_index + 1) * k;
      result[i].sample_index = result[i].sample_index - (result[i].length - 1) * k + sample_start;
      result[i].length *= k;
      size_t j;

      // Extend to the left.
      for (j = 1; j < k; ++j)
      {
        if (result[i].reference_index + result[i].length + j - 1 >= reference_end || result[i].sample_index - j < sample_start || complement[result[i].reference_index + result[i].length + j - 1] !=  sample[result[i].sample_index - j])
        {
          break;
        } // if
      } // for
      result[i].sample_index -= j - 1;
      result[i].length += j - 1;

      // Extend to the right.
      for (j = 1; j < k; ++j)
      {
        if (result[i].reference_index - j < reference_start || result[i].sample_index + result[i].length + j - 1 >= sample_end || complement[result[i].reference_index - j] != sample[result[i].sample_index + result[i].length + j - 1])
        {
          break;
        } // if
      } // for
      result[i].reference_index -= j - 1;
      result[i].length += j - 1;
      if (result[i].length > length)
      {
        length = result[i].length;
      } // if
    } // if

    // ``Normal'' LCSs (almost the same as reverse complement)
    else
    {
      result[i].reference_index = ((result[i].reference_index + 1) * k + reference_start - 1) - result[i].length * k + 1;
      result[i].sample_index = result[i].sample_index - (result[i].length - 1) * k + sample_start;
      result[i].length *= k;
      size_t j;

      // Extend to the left.
      for (j = 1; j < k; ++j)
      {
        if (result[i].reference_index - j < reference_start || result[i].sample_index - j < sample_start || reference[result[i].reference_index - j] != sample[result[i].sample_index - j])
        {
          break;
        } // if
      } // for
      result[i].reference_index -= j - 1;
      result[i].sample_index -= j - 1;
      result[i].length += j - 1;

      // Extend to the right.
      for (j = 0; j < k - 1; ++j)
      {
        if (result[i].reference_index + result[i].length + j >= reference_end || result[i].sample_index + result[i].length + j >= sample_end || reference[result[i].reference_index + result[i].length + j] != sample[result[i].sample_index + result[i].length + j])
        {
          break;
        } // if
      } // for
      result[i].length += j;
      if (result[i].length > length)
      {
        length = result[i].length;
      } // if
    } // else
  } // for

  // Remove all sub-optimal LCSs from the solution.
  for (size_t i = 0; i < result.size(); ++i)
  {
    if (result[i].length < length)
    {
      result.erase(result.begin() + i);
      --i;
    } // if
  } // for

  delete[] &LCS_line;
  delete[] &LCS_line_rc;

  return result;
} // LCS_k

// This function calculates the LCS using the LCS_k function by
// choosing an initial k and reducing it if necessary until the
// strings represent random strings.
std::vector<Substring> LCS(char const* const reference,
                           char const* const complement,
                           size_t const      reference_start,
                           size_t const      reference_end,
                           char const* const sample,
                           size_t const      sample_start,
                           size_t const      sample_end)
{
  size_t const reference_length = reference_end - reference_start;
  size_t const sample_length = sample_end - sample_start;

  size_t k = reference_length > sample_length ? sample_length / 4 : reference_length / 4;
  size_t const k_initial = k;

  std::vector<Substring> result;

  // FIXME: stop reducing k if the strings appear to be random
  // while (k > log(static_cast<double>(reference_end - reference_start)) / log(static_cast<double>(ALPHABET_SIZE[complement != 0 ? 0 : 1])))
  while (k > 4 && k_initial / k < 16)
  {

#if defined(__debug__)
  fprintf(stderr, "  k = %ld\n", k);
#endif

    result = LCS_k(reference, complement, reference_start, reference_end, sample, sample_start, sample_end, k);

    // A LCS of sufficient length has been found.
    if (result.size() > 0 && result[0].length >= 2 * k)
    {
      return result;
    } // if
    k /= 3;
  } // while

  return std::vector<Substring>();
  // Alternatively, find any LCS using the standard LCS algorithm.
  // Do NOT do this for large strings: instead return an empty set.
  // return LCS_1(reference, complement, reference_start, reference_end, sample, sample_start, sample_end);
} // LCS

} // namespace

