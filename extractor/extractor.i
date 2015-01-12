// *******************************************************************
//   (C) Copyright 2015 Leiden Institute of Advanced Computer Science
//   Universiteit Leiden
//   All Rights Reserved
// *******************************************************************
// Extractor (library)
// *******************************************************************
// FILE INFORMATION:
//   File:     extractor.i (SWIG interface file)
//   Author:   Jonathan K. Vis
//   Revision: 2.1.8
//   Date:     2015/01/12
// *******************************************************************
// DESCRIPTION:
//   Defines the SWIG interface for the Extractor library for use in
//   other languages than C/C++.
// *******************************************************************

%include "std_vector.i"

%module extractor
%{
#include "extractor.h"
%} // extractor

namespace std
{
%template(VariantVector) vector<mutalyzer::Variant>;
} // std

namespace mutalyzer
{

// Version string for run-time identification.
static char const* const VERSION = "2.1.8";

// The character type used for all strings. For now it should just be
// a char.
typedef char char_t;

// These constants can be used to specify the type of string to be
// extracted. The extractor is primarily focussed on DNA/RNA. When
// TYPE_PROTEIN (or another value) is used no complement string is
// constructed and no reverse complement is calculated. For
// TYPE_PROTEIN frame shift detection is applied on
// deletions/insertions.
static int const TYPE_DNA     = 0;
static int const TYPE_PROTEIN = 1;
static int const TYPE_OTHER   = 2;

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
// characters used.
static size_t const WEIGHT_BASE               = 1; // i.e., A, G, T, C
static size_t const WEIGHT_DELETION           = 3; // i.e., del
static size_t const WEIGHT_DELETION_INSERTION = 6; // i.e., delins
static size_t const WEIGHT_INSERTION          = 3; // i.e., ins
static size_t const WEIGHT_INVERSION          = 3; // i.e., inv
static size_t const WEIGHT_SEPARATOR          = 1; // i.e., _, [, ], ;
static size_t const WEIGHT_SUBSTITUTION       = 1; // i.e., >

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
//   @member transposition_start: starting position of a transposition
//                                withing the reference string
//   @member transposition_end: ending position of a transposition
//                              withing the reference string
// *******************************************************************
struct Variant
{
  size_t       reference_start;
  size_t       reference_end;
  size_t       sample_start;
  size_t       sample_end;
  unsigned int type;
  size_t       transposition_start;
  size_t       transposition_end;
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
//   This function is the interface function for Python.
//
//   @arg reference: reference string
//   @arg reference_length: length of the reference string
//   @arg sample: sample string
//   @arg sample_length: length of the sample string
//   @arg type: type of strings  0 --- DNA/RNA (default)
//                               1 --- Protein
//                               2 --- Other
//   @arg codon_string: TODO
//   @return: variant list with metadata
// *******************************************************************
Variant_List extract(char_t const* const reference,
                     size_t const        reference_length,
                     char_t const* const sample,
                     size_t const        sample_length,
                     int const           type = TYPE_DNA,
                     char_t const* const codon_string = 0);

} // mutalyzer

#include "extractor.h"

