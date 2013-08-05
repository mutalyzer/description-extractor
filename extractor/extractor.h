// *******************************************************************
//   (C) Copyright 2013 Leiden Institute of Advanced Computer Science
//   Universiteit Leiden
//   All Rights Reserved
// *******************************************************************
// Extractor (library)
// *******************************************************************
// FILE INFORMATION:
//   File:     extractor.h (implemented in extractor.cc)
//   Author:   Jonathan K. Vis
//   Revision: 1.03a
//   Date:     2013/07/29
// *******************************************************************
// DESCRIPTION:
//   This library can be used to generete HGVS variant descriptions as
//   accepted by the Mutalyzer Name Checker.
// *******************************************************************

#if !defined(__extractor_h__)
#define __extractor_h__

#include <cstddef>
#include <vector>

namespace mutalyzer
{

// Can be used as cut-off for k-mer reduction. The expected length of
// a longest common substring between to strings of length n is
// log_z(n), were z is the size of the alphabet.
static int const ALPHABET_SIZE[2] =
{
   4, // DNA/RNA
  20  // Protein/other
}; // ALPHABET


// *******************************************************************
// Variant Extraction
//   These functions are used to extract variants (regions of change)
//   between two strings.
// *******************************************************************


// *******************************************************************
// Variant structure
//   This structure describes a variant (region of change).
//
//   @member reference_start: starting position of the variant within
//                            the reference string
//   @member reference_end: ending position of the variant within the
//                          reference string
//   @member sample_start: starting position of the variant within the
//                         sample string
//   @member sample_end: ending position of the variant within the
//                       sample string
//   @member reverse_complement: indicates a reverse complement
//                               variant
// *******************************************************************
struct Variant
{
  size_t reference_start;
  size_t reference_end;
  size_t sample_start;
  size_t sample_end;
  bool   reverse_complement;

  inline Variant(size_t const reference_start,
                 size_t const reference_end,
                 size_t const sample_start,
                 size_t const sample_end,
                 bool   const reverse_complement = false):
         reference_start(reference_start),
         reference_end(reference_end),
         sample_start(sample_start),
         sample_end(sample_end),
         reverse_complement(reverse_complement) { }

  inline Variant(void) { }
}; // Variant


// *******************************************************************
// extract function
//   This function extracts the variants (regions of change) between
//   the reference and the sample string. It automatically constructs
//   the reverse complement string for the reference string if the
//   string type is DNA/RNA.
//
//   @arg reference: reference string
//   @arg reference_length: length of the reference string
//   @arg sample: sample string
//   @arg sample_length: length of the sample string
//   @arg type: type of strings  0 --- DNA/RNA (default)
//                               1 --- Protein
//   @return: list of variants
// *******************************************************************
std::vector<Variant> extract(char const* const reference,
                             size_t const      reference_length,
                             char const* const sample,
                             size_t const      sample_length,
                             int const         type = 0);


// *******************************************************************
// extractor function
//   This function extracts the variants (regions of change) between
//   the reference and the sample string by recursively calling itself
//   on prefixes and suffixes of a longest common substring.
//
//   @arg reference: reference string
//   @arg complement: complement string (can be null for
//                    strings other than DNA/RNA)
//   @arg reference_start: starting position in the reference string
//   @arg reference_end: ending position in the reference string
//   @arg sample: sample string
//   @arg sample_start: starting position in the sample string
//   @arg sample_end: ending position in the sample string
//   @arg result: list of variants
// *******************************************************************
void extractor(char const* const     reference,
               char const* const     complement,
               size_t const          reference_start,
               size_t const          reference_end,
               char const* const     sample,
               size_t const          sample_start,
               size_t const          sample_end,
               std::vector<Variant> &result);


// *******************************************************************
// Longest Common Substring (LCS) calculation
//   These functions are useful for LCS calculation of (large similar)
//   strings.
// *******************************************************************


// *******************************************************************
// Substring structure
//   This structure describes a common substring between two strings.
//
//   @member reference_index: starting position of the substring
//                            within the reference sequence
//   @member sample_index: ending position of the substring within the
//                         sample sequence
//   @member reverse_complement: indicates a reverse complement
//                               substring (only for DNA/RNA)
// *******************************************************************
struct Substring
{
  size_t reference_index;
  size_t sample_index;
  size_t length;
  bool   reverse_complement;

  inline Substring(size_t const reference_index,
                   size_t const sample_index,
                   size_t const length,
                   bool const   reverse_complement = false):
         reference_index(reference_index),
         sample_index(sample_index),
         length(length),
         reverse_complement(reverse_complement) { }

  inline Substring(void) { }
}; // Substring


// *******************************************************************
// LCS_1 function
//   This function calculates the longest common substrings between
//   two (three?) strings. It asumes no similarity between both
//   strings. Not for use for large strings.
//
//   @arg reference: reference string
//   @arg complement: complement string (can be null for strings other
//                    than DNA/RNA)
//   @arg reference_start: starting position in the reference string
//   @arg reference_end: ending position in the reference string
//   @arg sample: sample string
//   @arg sample_start: starting position in the sample string
//   @arg sample_end: ending position in the sample string
//   @return: list of longest common substrings
// *******************************************************************
std::vector<Substring> LCS_1(char const* const reference,
                             char const* const complement,
                             size_t const      reference_start,
                             size_t const      reference_end,
                             char const* const sample,
                             size_t const      sample_start,
                             size_t const      sample_end);


// *******************************************************************
// LCS_k function
//   This function calculates the longest common substrings between
//   two (three?) strings by encoding the reference and complement
//   strings into non-overlapping k-mers and the sample string into
//   overlapping k-mers. This function can be used for large similar
//   strings. If the returned vector is empty or the length of the
//   substrings is less or equal 2k, try again with a smaller k.
//
//   @arg reference: reference string
//   @arg complement: complement string (can be null for strings other
//                    than DNA/RNA)
//   @arg reference_start: starting position in the reference string
//   @arg reference_end: ending position in the reference string
//   @arg sample: sample string
//   @arg sample_start: starting position in the sample string
//   @arg sample_end: ending position in the sample string
//   @arg k: size of the k-mers, must be greater than 1
//   @return list of longest common substrings
// *******************************************************************
std::vector<Substring> LCS_k(char const* const reference,
                             char const* const complement,
                             size_t const      reference_start,
                             size_t const      reference_end,
                             char const* const sample,
                             size_t const      sample_start,
                             size_t const      sample_end,
                             size_t const      k);


// *******************************************************************
// LCS function
//   This function calculates the longest common substrings between
//   two (three?) strings by choosing an initial k and calling the
//   lcs_k function. The k is automatically reduced if necessary until
//   the LCS of the two strings approaches the theoretical length of 
//   the LCS of two random strings. FIXME: cut-off based on alphabet
//   size.
//
//   @arg reference: reference string
//   @arg complement: complement string (can be null for strings other
//                    than DNA/RNA)
//   @arg reference_start: starting position in the reference string
//   @arg reference_end: ending position in the reference string
//   @arg sample: sample string
//   @arg sample_start: starting position in the sample string
//   @arg sample_end: ending position in the sample string
//   @return: list of longest common substrings
// *******************************************************************
std::vector<Substring> LCS(char const* const reference,
                           char const* const complement,
                           size_t const      reference_start,
                           size_t const      reference_end,
                           char const* const sample,
                           size_t const      sample_start,
                           size_t const      sample_end);


// *******************************************************************
// IUPAC Nucleotide Acid Notation functions
//   These functions are useful for calculating the complement of DNA/
//   RNA strings.
// *******************************************************************


// *******************************************************************
// IUPAC_complement function
//   This function converts a IUPAC Nucleotide Acid Notation into its
//   complement.
//
//   @arg base: character from the IUPAC Nucleotide Acid Notation
//              alphabet
//   @return: bit array containing all the possible bases
// *******************************************************************
static inline char IUPAC_complement(char const base)
{
  switch (base)
  {
    case 'A':
      return 'T';
    case 'B':
      return 'V';
    case 'C':
      return 'G';
    case 'D':
      return 'H';
    case 'G':
      return 'C';
    case 'H':
      return 'D';
    case 'K':
      return 'M';
    case 'M':
      return 'K';
    case 'R':
      return 'Y';
    case 'T':
    case 'U':
      return 'A';
    case 'V':
      return 'B';
    case 'Y':
      return 'R';
  } // switch
  return base;
} // IUPAC_complement


// *******************************************************************
// IUPAC_complement function
//   This function converts a string in IUPAC Nucleotide Acid Notation
//   into its complement. A new string is allocated, so deletion is
//   the responsibility of the caller.
//
//   @arg string: string in the IUPAC Nucleotide Acid Notation
//                alphabet
//   @arg n: number of characters in the string to convert (might be
//           less than the actual string length)
//   @return: string containing the complement of the input string
// *******************************************************************
static inline char const* IUPAC_complement(char const* const string,
                                           size_t const      n)
{
  // The caller is responsible for deallocation.
  char* result = new char[n];
  for (size_t i = 0; i < n; ++i)
  {
    result[i] = IUPAC_complement(string[i]);
  } // for
  return result;
} // IUPAC_complement

} // namespace

#endif

