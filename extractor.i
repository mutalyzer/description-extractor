// *******************************************************************
//   (C) Copyright 2013 Leiden Institute of Advanced Computer Science
//   Universiteit Leiden
//   All Rights Reserved
// *******************************************************************
// Extractor (library)
// *******************************************************************
// FILE INFORMATION:
//   File:     extractor.i (SWIG interface file)
//   Author:   Jonathan K. Vis
//   Revision: 1.03a
//   Date:     2013/07/25
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
%template(SubstringVector) vector<mutalyzer::Substring>;
} // namespace

namespace mutalyzer
{

struct Variant
{
  size_t reference_start;
  size_t reference_end;
  size_t sample_start;
  size_t sample_end;
  bool   reverse_complement;
}; // Variant

std::vector<Variant> extract(char const* const reference,
                             size_t const      reference_length,
                             char const* const sample,
                             size_t const      sample_length,
                             size_t const      type = 0);

void extractor(char const* const     reference,
               char const* const     complement,
               size_t const          reference_start,
               size_t const          reference_end,
               char const* const     sample,
               size_t const          sample_start,
               size_t const          sample_end,
               std::vector<Variant> &result);

struct Substring
{
  size_t reference_index;
  size_t sample_index;
  size_t length;
  bool   reverse_complement;
}; // Substring

std::vector<Substring> LCS_1(char const* const reference,
                             char const* const complement,
                             size_t const      reference_start,
                             size_t const      reference_end,
                             char const* const sample,
                             size_t const      sample_start,
                             size_t const      sample_end);

std::vector<Substring> LCS_k(char const* const reference,
                             char const* const complement,
                             size_t const      reference_start,
                             size_t const      reference_end,
                             char const* const sample,
                             size_t const      sample_start,
                             size_t const      sample_end,
                             size_t const      k);

std::vector<Substring> LCS(char const* const reference,
                           char const* const complement,
                           size_t const      reference_start,
                           size_t const      reference_end,
                           char const* const sample,
                           size_t const      sample_start,
                           size_t const      sample_end);

char IUPAC_complement(char const base);

char const* IUPAC_complement(char const* const string,
                             size_t const      n);

} // namespace

#include "extractor.h"

