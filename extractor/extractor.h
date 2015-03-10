// *******************************************************************
//   (C) Copyright 2015 Leiden Institute of Advanced Computer Science
//   Universiteit Leiden
//   All Rights Reserved
// *******************************************************************
// Extractor (library)
// *******************************************************************
// FILE INFORMATION:
//   File:     extractor.h (implemented in extractor.cc)
//   Author:   Jonathan K. Vis
//   Revision: 2.2.0
//   Date:     2015/03/10
// *******************************************************************
// DESCRIPTION:
//   This library can be used to generate HGVS variant descriptions as
//   accepted by the Mutalyzer Name Checker.
// *******************************************************************

#if !defined(__extractor_h__)
#define __extractor_h__

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <vector>


#if defined(__debug__)
#include <cstdio>
#endif


namespace mutalyzer
{

// Version string for run-time identification.
static char const* const VERSION = "2.2.0";


// The character type used for all strings. For now it should just be
// a char.
typedef char char_t;


// *******************************************************************
// Variant Extraction
//   These functions are used to extract variants (regions of change)
//   between two strings.
// *******************************************************************


// These constants can be used to specify the type of string to be
// extracted. The extractor is primarily focussed on DNA/RNA. When
// TYPE_PROTEIN (or another value) is used no complement string is
// constructed and no reverse complement is calculated. For
// TYPE_PROTEIN frame shift detection is applied on
// deletions/insertions.
static int const TYPE_DNA     = 0; // DNA/RNA (default)
static int const TYPE_PROTEIN = 1; // Protein
static int const TYPE_OTHER   = 2; // Other strings


// These constants can be used to deterimine the type of variant.
// Substitution covers most: deletions, insertions, substitutions, and
// insertion/deletions. Indentity is used to describe the unchanged
// (matched) regions. The constants are coded as bitfields and should
// be appropriately combined, e.g., IDENTITY | TRANSPOSITION_OPEN for
// describing a real transposition. Note that some combinations do NOT
// make sense, e.g., SUBSTITUION | REVERSE_COMPLEMENT.
static unsigned int const IDENTITY            = 0x01;
static unsigned int const REVERSE_COMPLEMENT  = 0x02;
static unsigned int const SUBSTITUTION        = 0x04;
static unsigned int const TRANSPOSITION_OPEN  = 0x08;
static unsigned int const TRANSPOSITION_CLOSE = 0x10;
static unsigned int const FRAME_SHIFT         = 0x20;
static unsigned int const FRAME_SHIFT_1       = 0x01;
static unsigned int const FRAME_SHIFT_2       = 0x02;


// These constants are used in calculating the weight of the generated
// description and consequently used to end the description process
// when a certain ``trivial'' weight is exeeded. The weight constants
// are based on their HGVS description lengths, i.e., the amount of
// characters used. The weight_position variable is used to have a
// constant weight for a position description regardless the actual
// position. It is usually set to ceil(log10(|reference| / 4)), and
// its intention is to be constant during an extraction run.
extern size_t       weight_position;

static size_t const WEIGHT_BASE               = 1; // i.e., A, C, G, T
static size_t const WEIGHT_DELETION           = 3; // i.e., del
static size_t const WEIGHT_DELETION_INSERTION = 6; // i.e., delins
static size_t const WEIGHT_INSERTION          = 3; // i.e., ins
static size_t const WEIGHT_INVERSION          = 3; // i.e., inv
static size_t const WEIGHT_SEPARATOR          = 1; // i.e., _, [, ], ;
static size_t const WEIGHT_SUBSTITUTION       = 1; // i.e., >


// Cut-off constants. The threshold is used to specify the maximum
// reference length without any cut-off. Otherwise the weight_position
// is used as a cut-off. For transpositions the cut-off is specified
// as a fraction of the sample_length.
static size_t const THRESHOLD_CUT_OFF     = 16000;
static double const TRANSPOSITION_CUT_OFF =   0.1;


// This global variable is used to have access to the whole reference
// string at any point in the extraction process. Commonly used in
// transposition extraction.
extern size_t global_reference_length;


// Codon to amino acid table assuming the order of A, C, G, T, thus:
// codon_table[0] = AAA, ..., codon_table[63] = TTT.
extern char_t codon_table[64];


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
//   @member type: type of the variant described using the
//                 constants above
//   @member weight: weight of the variant according to the weight
//                   constants above (used internally)
//   @member transposition_start: starting position of a transposition
//                                within the reference string
//   @member transposition_end: ending position of a transposition
//                              within the reference string
// *******************************************************************
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

// *******************************************************************
// Variant_List structure
//   This structure describes a list of variants with associated
//   metadata.
//
//   @member weight_position: weight used for position descriptors
//   @member variants: vector of variants
// *******************************************************************
struct Variant_List
{
  size_t               weight_position;
  std::vector<Variant> variants;
}; // Variant_List

// *******************************************************************
// extract function
//   This function is the interface function for Python. It is just a
//   wrapper for the C++ extract function below.
//
//   @arg reference: reference string
//   @arg reference_length: length of the reference string
//   @arg sample: sample string
//   @arg sample_length: length of the sample string
//   @arg type: type of strings  0 --- DNA/RNA (default)
//                               1 --- Protein
//                               2 --- Other
//   @arg codon_string: serialized codon table: 64 characters
//                      corresponding to the codons AAA, ..., TTT.
//                      Only for protein extraction (frame shifts).
//   @return: variant list with metadata
// *******************************************************************
Variant_List extract(char_t const* const reference,
                     size_t const        reference_length,
                     char_t const* const sample,
                     size_t const        sample_length,
                     int const           type         = TYPE_DNA,
                     char_t const* const codon_string = 0);

// *******************************************************************
// extract function
//   This function extracts the variants (regions of change) between
//   the reference and the sample string. It automatically constructs
//   the reverse complement string for the reference string if the
//   string type is DNA/RNA.
//
//   @arg variant: vector of variants
//   @arg reference: reference string
//   @arg reference_length: length of the reference string
//   @arg sample: sample string
//   @arg sample_length: length of the sample string
//   @arg type: type of strings  0 --- DNA/RNA (default)
//                               1 --- Protein
//                               2 --- Other
//   @arg codon_string: serialized codon table: 64 characters
//                      corresponding to the codons AAA, ..., TTT.
//                      Only for protein extraction (frame shifts).
//   @return: weight of the extracted variants
// *******************************************************************
size_t extract(std::vector<Variant> &variant,
               char_t const* const   reference,
               size_t const          reference_length,
               char_t const* const   sample,
               size_t const          sample_length,
               int const             type         = TYPE_DNA,
               char_t const* const   codon_string = 0);

// *******************************************************************
// extractor function
//   This function extracts the variants (regions of change) between
//   the reference and the sample string by recursively calling itself
//   on prefixes and suffixes of a longest common substring.
//
//   @arg variant: vector of variants
//   @arg reference: reference string
//   @arg complement: complement string (can be null for strings other
//                    than DNA/RNA)
//   @arg reference_start: starting position in the reference string
//   @arg reference_end: ending position in the reference string
//   @arg sample: sample string
//   @arg sample_start: starting position in the sample string
//   @arg sample_end: ending position in the sample string
//   @return: weight of the extracted variants
// *******************************************************************
size_t extractor(std::vector<Variant> &variant,
                 char_t const* const   reference,
                 char_t const* const   complement,
                 size_t const          reference_start,
                 size_t const          reference_end,
                 char_t const* const   sample,
                 size_t const          sample_start,
                 size_t const          sample_end);

// *******************************************************************
// extractor_transposition function
//   This function extracts the variants (regions of change) between
//   a part of the sample string classified as an insertion and the
//   whole reference string.
//
//   @arg variant: vector of variants
//   @arg reference: reference string
//   @arg complement: complement string (can be null for strings other
//                    than DNA/RNA)
//   @arg reference_start: starting position in the reference string
//                         used for the deletion part
//   @arg reference_end: ending position in the reference string used
//                       for the deletion part
//   @arg sample: sample string
//   @arg sample_start: starting position in the sample string
//   @arg sample_end: ending position in the sample string
//   @arg weight_trivial: trivial weight to describe the transposition
//                        as a normal insertion (used for ending the
//                        extraction process)
//   @return: weight of the extracted variants
// *******************************************************************
size_t extractor_transposition(std::vector<Variant> &variant,
                               char_t const* const   reference,
                               char_t const* const   complement,
                               size_t const          reference_start,
                               size_t const          reference_end,
                               char_t const* const   sample,
                               size_t const          sample_start,
                               size_t const          sample_end,
                               size_t const          weight_trivial = 0);

// *******************************************************************
// extractor_protein function
//   This function extracts the variants (regions of change) between
//   the reference and the sample protein string by recursively
//   calling itself on prefixes and suffixes of a longest common
//   substring, calculated by the LCS_1 algorithm (these strings are
//   very short).
//
//   @arg variant: vector of variants
//   @arg reference: reference string
//   @arg reference_start: starting position in the reference string
//   @arg reference_end: ending position in the reference string used
//   @arg sample: sample string
//   @arg sample_start: starting position in the sample string
//   @arg sample_end: ending position in the sample string
//   @return: weight of the extracted variants
// *******************************************************************
size_t extractor_protein(std::vector<Variant> &variant,
                         char_t const* const   reference,
                         size_t const          reference_start,
                         size_t const          reference_end,
                         char_t const* const   sample,
                         size_t const          sample_start,
                         size_t const          sample_end);


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
//   @member length: length of the substring
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
// LCS function
//   This function calculates the longest common substrings between
//   two (three?) strings by choosing an initial k and calling the
//   lcs_k function. The k is automatically reduced if necessary until
//   the LCS of the two strings approaches some cutoff threshold.
//
//   @arg substring: vector of substrings
//   @arg reference: reference string
//   @arg complement: complement string (can be null for strings other
//                    than DNA/RNA)
//   @arg reference_start: starting position in the reference string
//   @arg reference_end: ending position in the reference string
//   @arg sample: sample string
//   @arg sample_start: starting position in the sample string
//   @arg sample_end: ending position in the sample string
//   @arg cut_off: optional cut-off value for the k in LCS_k
//   @return: length of the LCS
// *******************************************************************
size_t LCS(std::vector<Substring> &substring,
           char_t const* const     reference,
           char_t const* const     complement,
           size_t const            reference_start,
           size_t const            reference_end,
           char_t const* const     sample,
           size_t const            sample_start,
           size_t const            sample_end,
           size_t const            cut_off = 1);

// *******************************************************************
// LCS_1 function
//   This function calculates the longest common substrings between
//   two (three?) strings. It asumes no similarity between both
//   strings. Not for use for large strings. This is the classical
//   dynamic programming algorithm.
//
//   @arg substring: vector of substrings
//   @arg reference: reference string
//   @arg complement: complement string (can be null for strings other
//                    than DNA/RNA)
//   @arg reference_start: starting position in the reference string
//   @arg reference_end: ending position in the reference string
//   @arg sample: sample string
//   @arg sample_start: starting position in the sample string
//   @arg sample_end: ending position in the sample string
//   @return: length of the LCS
// *******************************************************************
size_t LCS_1(std::vector<Substring> &substring,
             char_t const* const     reference,
             char_t const* const     complement,
             size_t const            reference_start,
             size_t const            reference_end,
             char_t const* const     sample,
             size_t const            sample_start,
             size_t const            sample_end);

// *******************************************************************
// LCS_k function
//   This function calculates the longest common substrings between
//   two (three?) strings by encoding the reference and complement
//   strings into non-overlapping k-mers and the sample string into
//   overlapping k-mers. This function can be used for large similar
//   strings. If the returned vector is empty or the length of the
//   substrings is less or equal 2k, try again with a smaller k.
//
//   @arg substring: vector of substrings
//   @arg reference: reference string
//   @arg complement: complement string (can be null for strings other
//                    than DNA/RNA)
//   @arg reference_start: starting position in the reference string
//   @arg reference_end: ending position in the reference string
//   @arg sample: sample string
//   @arg sample_start: starting position in the sample string
//   @arg sample_end: ending position in the sample string
//   @arg k: size of the k-mers, must be greater than 1
//   @return: length of the LCS
// *******************************************************************
size_t LCS_k(std::vector<Substring> &substring,
             char_t const* const     reference,
             char_t const* const     complement,
             size_t const            reference_start,
             size_t const            reference_end,
             char_t const* const     sample,
             size_t const            sample_start,
             size_t const            sample_end,
             size_t const            k);


// *******************************************************************
// General string matching functions
//   These functions are useful for string matching.
// *******************************************************************

// *******************************************************************
// string_match function
//   This function is more or less equivalent to C's strncmp.
//
//   @arg string_1: first string to be compared
//   @arg string_2: second string to be compared
//   @arg length: maximum length to be compared

//   @return: true iff string_1 and string_2 match for the given
//            length
// *******************************************************************
bool string_match(char_t const* const string_1,
                  char_t const* const string_2,
                  size_t const        length);

// *******************************************************************
// string_match_reverse function
//   This function is very similar to C's strncmp, but it traverses
//   string_1 from end to start while traversing string_2 from start
//   to end (useful for the reverse complement in DNA/RNA).
//
//   @arg string_1: first string to be compared
//   @arg string_2: second string to be compared
//   @arg length: maximum length to be compared

//   @return: true iff string_1 and string_2 match in their respective
//            directions for the given length
// *******************************************************************
bool string_match_reverse(char_t const* const string_1,
                          char_t const* const string_2,
                          size_t const        length);

// *******************************************************************
// prefix_match function
//   This function calculates the length (in characters) of the common
//   prefix between two strings. The result of this function is also
//   used in the suffix_match function.
//
//   @arg reference: reference string
//   @arg reference_length: reference length
//   @arg sample: sample string
//   @arg sample_length: sample length

//   @return: the length of the common prefix
// *******************************************************************
size_t prefix_match(char_t const* const reference,
                    size_t const        reference_length,
                    char_t const* const sample,
                    size_t const        sample_length);

// *******************************************************************
// suffix_match function
//   This function calculates the length (in characters) of the common
//   suffix between two strings. It needs the calculated common
//   prefix.
//
//   @arg reference: reference string
//   @arg reference_length: reference length
//   @arg sample: sample string
//   @arg sample_length: sample length
//   @arg prefix: length of the common prefix

//   @return: the length of the common suffix
// *******************************************************************
size_t suffix_match(char_t const* const reference,
                    size_t const        reference_length,
                    char_t const* const sample,
                    size_t const        sample_length,
                    size_t const        prefix = 0);


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
//   @return: its corresponding IUPAC complement for single bases only
// *******************************************************************
char_t IUPAC_complement(char_t const base);

// *******************************************************************
// IUPAC_complement function
//   This function converts a string in IUPAC Nucleotide Acid Notation
//   into its complement. A new string is allocated, so deletion is
//   the responsibility of the caller.
//
//   @arg string: string in the IUPAC Nucleotide Acid Notation
//                alphabet
//   @arg length: number of characters in the string to convert (might
//                be less than the actual string length)
//   @return: string containing the complement of the input string
// *******************************************************************
char_t const* IUPAC_complement(char_t const* const string,
                               size_t const        length);


// *******************************************************************
// Amino Acid functions
//   These functions are useful for calculating frame shifts.
// *******************************************************************

// *******************************************************************
// frame_shift function
//   This function calculates the frame shift A reference amino acid
//   is checked against two possible partial overlaps between every
//   combination of two sample (observed) amino acids. Possible
//   results are: 0 (no frame shift), 1, 2, or 3 (both frame shifts
//   are possible).
//
//   @arg reference: reference amino acid
//   @arg sample_1: first sample amino acid
//   @arg sample_2: second sample amino acid
//   @return: frame shift
// *******************************************************************
unsigned int frame_shift(char_t const reference,
                         char_t const sample_1,
                         char_t const sample_2);

// *******************************************************************
// annotate_frame_shift function
//   This function calculates frame shift annotations for
//   deletions/insertions.
//
//   @arg variant: vector of variants
//   @arg reference: reference string
//   @arg reference_start: starting position in the reference string
//   @arg reference_end: ending position in the reference string
//   @arg sample: sample string
//   @arg sample_start: starting position in the sample string
//   @arg sample_end: ending position in the sample string
//   @return: number of frame shift variants
// *******************************************************************
size_t annotate_frame_shift(std::vector<Variant> &variant,
                            char_t const* const   reference,
                            size_t const          reference_start,
                            size_t const          reference_end,
                            char_t const* const   sample,
                            size_t const          sample_start,
                            size_t const          sample_end);


#if defined(__debug__)
// *******************************************************************
// Dprint_truncated function
//   Debug function for printing large strings in truncated form: a
//   prefix of a certain length ... a suffix of the same length.
//
//   @arg string: string to be printed
//   @arg start: starting position in the string
//   @arg end: ending position in the string
//   @arg length: length of the prefix and suffix
//   @arg stream: file stream to print to
//   @return: the length of the printed string
// *******************************************************************
size_t Dprint_truncated(char_t const* const string,
                        size_t const        start,
                        size_t const        end,
                        size_t const        length = 40,
                        FILE*               stream = stderr);
#endif


} // mutalyzer

#endif

