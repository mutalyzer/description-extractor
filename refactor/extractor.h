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

typedef char char_t;


static int const TYPE_DNA     = 0;
static int const TYPE_PROTEIN = 1;


static int const IDENTITY            = 0;
static int const REVERSE_COMPLEMENT  = 1;
static int const SUBSTITUTION        = 2;
static int const TRANSPOSITION_OPEN  = 4;
static int const TRANSPOSITION_CLOSE = 8;


extern size_t       weight_position;
static size_t const WEIGHT_BASE               = 1; // A
static size_t const WEIGHT_DELETION           = 3; // del
static size_t const WEIGHT_DELETION_INSERTION = 6; // delins
static size_t const WEIGHT_INSERTION          = 3; // ins
static size_t const WEIGHT_INVERSE            = 3; // inv
static size_t const WEIGHT_SEPARATOR          = 1; // _[]
static size_t const WEIGHT_SUBSTITUTION       = 1; // >


extern size_t global_reference_length;


struct Variant
{
  size_t reference_start;
  size_t reference_end;
  size_t sample_start;
  size_t sample_end;
  int    type;
  size_t weight;
  size_t transposition_start;
  size_t transposition_end;

  inline Variant(size_t const reference_start,
                 size_t const reference_end,
                 size_t const sample_start,
                 size_t const sample_end,
                 int const    type = IDENTITY,
                 size_t const weight = 0,
                 size_t const transposition_start = 0,
                 size_t const transposition_end = 0):
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


size_t extract(std::vector<Variant> &variant,
               char_t const* const   reference,
               size_t const          reference_length,
               char_t const* const   sample,
               size_t const          sample_start,
               int const             type = TYPE_DNA);

size_t extractor(std::vector<Variant> &variant,
                 char_t const* const   reference,
                 char_t const* const   complement,
                 size_t const          reference_start,
                 size_t const          reference_end,
                 char_t const* const   sample,
                 size_t const          sample_start,
                 size_t const          sample_end,
                 bool const            transposition = false);

size_t extractor_transposition(std::vector<Variant> &variant,
                               char_t const* const   reference,
                               char_t const* const   complement,
                               size_t const          reference_start,
                               size_t const          reference_end,
                               char_t const* const   sample,
                               size_t const          sample_start,
                               size_t const          sample_end);


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


size_t LCS(std::vector<Substring> &substring,
           char_t const* const     reference,
           char_t const* const     complement,
           size_t const            reference_start,
           size_t const            reference_end,
           char_t const* const     sample,
           size_t const            sample_start,
           size_t const            sample_end);

size_t LCS_1(std::vector<Substring> &substring,
             char_t const* const     reference,
             char_t const* const     complement,
             size_t const            reference_start,
             size_t const            reference_end,
             char_t const* const     sample,
             size_t const            sample_start,
             size_t const            sample_end);

/*
size_t LCS_k(std::vector<Substring> &substring,
             char_t const* const     reference,
             char_t const* const     complement,
             size_t const            reference_start,
             size_t const            reference_end,
             char_t const* const     sample,
             size_t const            sample_start,
             size_t const            sample_end,
             size_t const            k);
*/

size_t prefix_match(char_t const* const reference,
                    size_t const        reference_length,
                    char_t const* const sample,
                    size_t const        sample_length);
size_t suffix_match(char_t const* const reference,
                    size_t const        reference_length,
                    char_t const* const sample,
                    size_t const        sample_length,
                    size_t const        prefix = 0);

char_t IUPAC_complement(char_t const base);
char_t const* IUPAC_complement(char_t const* const string,
                               size_t const        length);


#if defined(__debug__)
size_t Dprint_truncated(char_t const* const string,
                        size_t const        start,
                        size_t const        end,
                        size_t const        length = 40,
                        FILE*               stream = stderr);
#endif

} // mutalyzer

#endif

