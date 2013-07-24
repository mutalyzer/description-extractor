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

struct Variant
{
  long reference_start;
  long reference_end;
  long sample_start;
  long sample_end;
  bool reverse_complement;
}; // Variant

std::vector<Variant> extract(char const* const     reference,
                             long const            reference_length,
                             char const* const     sample,
                             long const            sample_length,
                             long const            type);

}

#include "extractor.h"

