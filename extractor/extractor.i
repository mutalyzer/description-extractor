// *******************************************************************
//   (C) Copyright 2014 Leiden Institute of Advanced Computer Science
//   Universiteit Leiden
//   All Rights Reserved
// *******************************************************************
// Extractor (library)
// *******************************************************************
// FILE INFORMATION:
//   File:     extractor.i (SWIG interface file)
//   Author:   Jonathan K. Vis
//   Revision: 2.0.3
//   Date:     2014/07/31
// *******************************************************************
// DESCRIPTION:
//   Defines the SWIG interface for the Extractor library for use in
//   other languages than C/C++.
// *******************************************************************

%include "std_vector.i"

%module extractor
%{
#include "extractor.h"
%} // %module

namespace std
{
%template(VariantVector) vector<mutalyzer::Variant>;
} // namespace

namespace mutalyzer
{

// Version string for run-time identification.
static char const* const VERSION = "2.0.3";

// The character type used for all strings. For now it should just be
// a char.
typedef char char_t;

// These constants can be used to specify the type of string to be
// extracted. The extractor is primarily focussed on DNA/RNA. When
// TYPE_PROTEIN (or another value) is used no complement string is
// constructed and no reverse complement is calculated.
static int const TYPE_DNA     = 0;
static int const TYPE_PROTEIN = 1;

// These constants can be used to deterimine the type of variant.
// Substitution covers most: deletions, insertions, substitutions, and
// insertion/deletions. Indentity is used to describe the unchanged
// (matched) regions. The transposition constants are coded as
// bitfields and should be appropriately combined, e.g.,
// IDENTITY | TRANSPOSITION_OPEN for describing a real transposition.
// Note that some combinations do NOT make sense, e.g.,
// SUBSTITUION | REVERSE_COMPLEMENT.
static int const IDENTITY            = 0x01;
static int const REVERSE_COMPLEMENT  = 0x02;
static int const SUBSTITUTION        = 0x04;
static int const TRANSPOSITION_OPEN  = 0x08;
static int const TRANSPOSITION_CLOSE = 0x10;

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
//                                withing the reference string
//   @member transposition_end: ending position of a transposition
//                              withing the reference string
// *******************************************************************
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
}; // Variant

// *******************************************************************
// extract function
//   This function is the interface function for Python.
//
//   @arg reference: reference string
//   @arg reference_length: length of the reference string
//   @arg sample: sample string
//   @arg sample_length: length of the sample string
//   @arg type: type of strings  0 --- DNA/RNA (default)
//                               1 --- Protein/other
//   @return: vector of variants
// *******************************************************************
std::vector<Variant> extract(char_t const* const reference,
                             size_t const        reference_length,
                             char_t const* const sample,
                             size_t const        sample_length,
                             int const           type = TYPE_DNA);

} // namespace

#include "extractor.h"

