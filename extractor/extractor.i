// *******************************************************************
// Extractor (library)
// *******************************************************************
// FILE INFORMATION:
//   File:     extractor.i (SWIG interface file)
//   Author:   Jonathan K. Vis
// *******************************************************************
// DESCRIPTION:
//   Defines the SWIG interface for the Extractor library for use in
//   other languages than C/C++.
// *******************************************************************

%include "std_vector.i"

%module extractor
%{
#include "extractor.h"
%}

namespace std
{
%template(VariantVector) vector<mutalyzer::Variant>;
}

namespace mutalyzer
{

static char const* const VERSION;

typedef char char_t;

static int const TYPE_DNA;
static int const TYPE_PROTEIN;
static int const TYPE_OTHER;

static unsigned int const IDENTITY;
static unsigned int const REVERSE_COMPLEMENT;
static unsigned int const SUBSTITUTION;
static unsigned int const TRANSPOSITION_OPEN;
static unsigned int const TRANSPOSITION_CLOSE;
static unsigned int const FRAME_SHIFT;

static unsigned int const FRAME_SHIFT_NONE;
static unsigned int const FRAME_SHIFT_1;
static unsigned int const FRAME_SHIFT_2;
static unsigned int const FRAME_SHIFT_REVERSE;
static unsigned int const FRAME_SHIFT_REVERSE_1;
static unsigned int const FRAME_SHIFT_REVERSE_2;

static size_t const WEIGHT_BASE;
static size_t const WEIGHT_DELETION;
static size_t const WEIGHT_DELETION_INSERTION;
static size_t const WEIGHT_INSERTION;
static size_t const WEIGHT_INVERSION;
static size_t const WEIGHT_SEPARATOR;
static size_t const WEIGHT_SUBSTITUTION;

struct Variant
{
  size_t       reference_start;
  size_t       reference_end;
  size_t       sample_start;
  size_t       sample_end;
  unsigned int type;
  size_t       transposition_start;
  size_t       transposition_end;
};

struct Variant_List
{
  size_t               weight_position;
  std::vector<Variant> variants;
};

Variant_List extract(char_t const* const reference,
                     size_t const        reference_length,
                     char_t const* const sample,
                     size_t const        sample_length,
                     int const           type = TYPE_DNA,
                     char_t const* const codon_string = 0);

}
