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
//   Revision: 1.02a
//   Date:     2013/07/24
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

static int const ALPHABET_SIZE[2] =
{
   4, // DNA/RNA
  20  // Protein
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
//   the LCS of two random strings. FIXME
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
//   These functions are useful for (inexact) matching of DNA/RNA. All
//   functions are case insensitive.
//   Note: gaps are not part of out alphabet.
// *******************************************************************


// *******************************************************************
// IUPAC_base function
//   This function converts a IUPAC Nucleotide Acid Notation into a
//   bit array for matching.
//
//   @arg base: character from the IUPAC Nucleotide Acid Notation
//              alphabet
//   @return: bit array containing all the possible bases
// *******************************************************************
static inline unsigned int IUPAC_base(char const base);


// *******************************************************************
// IUPAC_complement function
//   This function converts a IUPAC Nucleotide Acid Notation into its
//   complement.
//
//   @arg base: character from the IUPAC Nucleotide Acid Notation
//              alphabet
//   @return: bit array containing all the possible bases
// *******************************************************************
static inline char IUPAC_complement(char const base);


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
                                           size_t const      n);


// *******************************************************************
// IUPAC_match function
//   This function matches two strings in the IUPAC Nucleotide Acid
//   Notation alphabet.
//   Note: some ``weird'' matching is also possible e.g., Y ({C, T})
//   matches S ({C, G}), because they both contain C.
//
//   @arg string_1: string in the IUPAC Nucleotide Acid Notation
//                  alphabet
//   @arg string_2: string in the IUPAC Nucleotide Acid Notation
//                  alphabet
//   @arg n: number of characters in the string to match (might be
//           less than the actual string length)
//   @return: true or false whether the two strings match or not
// *******************************************************************
static inline bool IUPAC_match(char const* const string_1,
                               char const* const string_2,
                               size_t const      n = 1);

// *******************************************************************
// IUPAC_match_reverse function
//   This function matches two strings in the IUPAC Nucleotide Acid
//   Notation alphabet. The first string is matched in reverse order.
//   This function is useful when matching the reverse complement of
//   DNA/RNA.
//   Note: some ``weird'' matching is also possible e.g., Y ({C, T})
//   matches S ({C, G}), because they both contain C.
//
//   @arg string_1: string in the IUPAC Nucleotide Acid Notation
//                  alphabet
//   @arg string_2: string in the IUPAC Nucleotide Acid Notation
//                  alphabet
//   @arg n: number of characters in the string to match (might be
//           less than the actual string length)
//   @return: true or false whether the two strings match or not
// *******************************************************************
static inline bool IUPAC_match_reverse(char const* const string_1,
                                       char const* const string_2,
                                       size_t const      n = 1);


// *******************************************************************
// Implementation Section (not for reference)
// *******************************************************************


static unsigned int const IUPAC_ADENINE  = 1;
static unsigned int const IUPAC_CYTOSINE = 2;
static unsigned int const IUPAC_GUANINE  = 4;
static unsigned int const IUPAC_THYMINE  = 8;

static inline unsigned int IUPAC_base(char const base)
{
  switch (base)
  {
    case 'A':
    case 'a':
      return IUPAC_ADENINE;
    case 'B':
    case 'b':
      return IUPAC_CYTOSINE | IUPAC_GUANINE | IUPAC_THYMINE;
    case 'C':
    case 'c':
      return IUPAC_CYTOSINE;
    case 'D':
    case 'd':
      return IUPAC_ADENINE | IUPAC_GUANINE | IUPAC_THYMINE;
    case 'G':
    case 'g':
      return IUPAC_GUANINE;
    case 'H':
    case 'h':
      return IUPAC_ADENINE | IUPAC_CYTOSINE | IUPAC_THYMINE;
    case 'K':
    case 'k':
      return IUPAC_GUANINE | IUPAC_THYMINE;
    case 'M':
    case 'm':
      return IUPAC_ADENINE | IUPAC_CYTOSINE;
    case 'N':
    case 'n':
      return IUPAC_ADENINE | IUPAC_CYTOSINE |
             IUPAC_GUANINE | IUPAC_THYMINE;
    case 'R':
    case 'r':
      return IUPAC_ADENINE | IUPAC_GUANINE;
    case 'S':
    case 's':
      return IUPAC_CYTOSINE | IUPAC_GUANINE;
    case 'T':
    case 't':
    case 'U':
    case 'u':
      return IUPAC_THYMINE;
    case 'V':
    case 'v':
      return IUPAC_ADENINE | IUPAC_CYTOSINE | IUPAC_GUANINE;
    case 'W':
    case 'w':
      return IUPAC_ADENINE | IUPAC_THYMINE;
    case 'Y':
    case 'y':
      return IUPAC_CYTOSINE | IUPAC_THYMINE;
  } // switch
  return 0;
} // IUPAC_base

static inline char IUPAC_complement(char const base)
{
  switch (base)
  {
    case 'A':
    case 'a':
      return 'T';
    case 'B':
    case 'b':
      return 'V';
    case 'C':
    case 'c':
      return 'G';
    case 'D':
    case 'd':
      return 'H';
    case 'G':
    case 'g':
      return 'C';
    case 'H':
    case 'h':
      return 'D';
    case 'K':
    case 'k':
      return 'M';
    case 'M':
    case 'm':
      return 'K';
    case 'R':
    case 'r':
      return 'Y';
    case 'T':
    case 't':
    case 'U':
    case 'u':
      return 'A';
    case 'V':
    case 'v':
      return 'B';
    case 'Y':
    case 'y':
      return 'R';
  } // switch
  return base;
} // IUPAC_complement

static inline char const* IUPAC_complement(char const* const string,
                                           size_t const      n)
{
  // the caller is responsible for deallocation
  char* result = new char[n];
  for (size_t i = 0; i < n; ++i)
  {
    result[i] = IUPAC_complement(string[i]);
  } // for
  return result;
} // IUPAC_complement

static inline bool IUPAC_match(char const* const string_1,
                               char const* const string_2,
                               size_t const      n)
{
  for (size_t i = 0; i < n; ++i)
  {
    if ((IUPAC_base(string_1[i]) & IUPAC_base(string_2[i])) == 0)
    {
      return false;
    } // if
  } // for
  return true;
} // IUPAC_match

static inline bool IUPAC_match_reverse(char const* const string_1,
                                       char const* const string_2,
                                       size_t const      n)
{
  for (size_t i = 0; i < n; ++i)
  {
    if ((IUPAC_base(string_1[-i]) & IUPAC_base(string_2[i])) == 0)
    {
      return false;
    } // if
  } // for
  return true;
} // IUPAC_match_reverse

} // namespace

#endif

