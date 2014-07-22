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
//   Revision: 2.01a
//   Date:     2014/07/22
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

typedef char char_t;

static int const TYPE_DNA     = 0;
static int const TYPE_PROTEIN = 1;

static int const IDENTITY            = 0;
static int const REVERSE_COMPLEMENT  = 1;
static int const SUBSTITUTION        = 2;
static int const TRANSPOSITION_OPEN  = 4;
static int const TRANSPOSITION_CLOSE = 8;

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

std::vector<Variant> extract(char_t const* const reference,
                             size_t const        reference_length,
                             char_t const* const sample,
                             size_t const        sample_length,
                             int const           type = TYPE_DNA);

} // namespace

#include "extractor.h"

