// *******************************************************************
// Extractor (library)
// *******************************************************************
// FILE INFORMATION:
//   File:     extractor.cc (depends on extractor.h)
//   Author:   Jonathan K. Vis
// *******************************************************************
// DESCRIPTION:
//   This library can be used to generete HGVS variant descriptions as
//   accepted by the Mutalyzer Name Checker.
// *******************************************************************

#include "extractor.h"

namespace mutalyzer
{

// The (average) description length of a position. Depends on the
// reference string length: ceil(log10(|reference| / 4)).
size_t weight_position = 1;

// This global variable is a dirty trick used to always have access to
// the complete reference strings even when deep into the recursion.
// This seems necessary to compute transpositions.
size_t global_reference_length = 0;

static char_t const IUPAC_ALPHA[16] =
{
  'x',  // 0x00
  'A',  // 0x01
  'C',  // 0x02
  'M',  // 0x03  A | C
  'G',  // 0x04
  'R',  // 0x05  A | G
  'S',  // 0x06  C | G
  'V',  // 0x07  A | C | G
  'T',  // 0x08
  'W',  // 0x09  A | T
  'Y',  // 0x0a  C | T
  'H',  // 0x0b  A | C | T
  'K',  // 0x0c  G | T
  'D',  // 0x0d  A | G | T
  'B',  // 0x0e  C | G | T
  'N'   // 0x0f  A | C | G | T
}; // IUPAC_ALPHA

static char_t const IUPAC_BASE[4] =
{
  'A',
  'C',
  'G',
  'T'
}; // IUPAC_BASE

// The actual frame shift map indexed on the lower 127 ASCII
// characters. This map should be precalculated given a codon string
// by the initialize_frame_shift_map function.
static uint8_t frame_shift_map[128][128][128] = {{{FRAME_SHIFT_NONE}}};

// A frequency count of all possible frame shifts (5) for all
// combinations of two amino acids (indexed by the lower 127 ASCII
// characters). Used to calculate the frame shift probability.
static uint8_t frame_shift_count[128][128][5] = {{{0}}};

static uint64_t acid_map[128] = {0x0ull};

static double acid_frequency[128] = {.0f};
static double frame_shift_frequency[128][128][5] = {{{.0f}}};

// This character is always ignored when LCS matching and can be used for
// repeat masking
static char_t const MASK = '$';

// Only used to interface to Python: calls the C++ extract function.
Variant_List extract(char_t const* const reference,
                     size_t const        reference_length,
                     char_t const* const sample,
                     size_t const        sample_length,
                     int const           type,
                     char_t const* const codon_string)
{
  Variant_List variant_list;
  extract(variant_list.variants, reference, reference_length, sample, sample_length, type, codon_string);
  variant_list.weight_position = weight_position;
  return variant_list;
} // extract

// The main library function. Extract all variants (regions of change)
// from the given strings.
size_t extract(std::vector<Variant> &variant,
               char_t const* const   reference,
               size_t const          reference_length,
               char_t const* const   sample,
               size_t const          sample_length,
               int const             type,
               char_t const* const   codon_string)
{
  // The global variables are set for this extraction run.
  global_reference_length = reference_length;
  weight_position = ceil(log10(reference_length / 4));
  if (weight_position <= 0)
  {
    weight_position = 1;
  } // if

  // Common prefix and suffix snooping.
  size_t const prefix = prefix_match(reference, reference_length, sample, sample_length);
  size_t const suffix = suffix_match(reference, reference_length, sample, sample_length, prefix);


#if defined(__debug__)
  fputs("reference: ", stderr);
  Dprint_truncated(reference, 0, reference_length);
  fprintf(stderr, " (%ld)\n", reference_length);
  fputs("sample:    ", stderr);
  Dprint_truncated(sample, 0, sample_length);
  fprintf(stderr, " (%ld)\n", sample_length);
  fputs("prefix: ", stderr);
  Dprint_truncated(reference, 0, prefix);
  fprintf(stderr, " (%ld)\n", prefix);
  fputs("suffix: ", stderr);
  Dprint_truncated(reference, reference_length - suffix, reference_length);
  fprintf(stderr, " (%ld)\n", suffix);
  if (type == TYPE_DNA)
  {
    fputs("Construction IUPAC complement\n", stderr);
  } // if
#endif


  // Do NOT construct a complement string for protein strings. All
  // other string types default to protein strings.
  char_t const* complement = type == TYPE_DNA ? IUPAC_complement(reference, reference_length) : 0;


#if defined(__debug__)
  if (type == TYPE_DNA)
  {
    fputs("complement: ", stderr);
    Dprint_truncated(complement, 0, reference_length);
    fprintf(stderr, " (%ld)\n", reference_length);
  } // if
  fprintf(stderr, "position weight: %ld\n", weight_position);
#endif


  if (prefix > 0)
  {
    variant.push_back(Variant(0, prefix, 0, prefix));
  } // if

  // The actual extraction process starts here.
  size_t weight;
  if (type == TYPE_PROTEIN)
  {
    initialize_frame_shift_map(codon_string);

    weight = extractor_protein(variant, reference, prefix, reference_length - suffix, sample, prefix, sample_length - suffix);
  } // if
  else
  {
    weight = extractor(variant, reference, complement, prefix, reference_length - suffix, sample, prefix, sample_length - suffix);
  } // else

  if (suffix > 0)
  {
    variant.push_back(Variant(reference_length - suffix, reference_length, sample_length - suffix, sample_length));
  } // if


  // Frame shift annotation starts here.
  if (type == TYPE_PROTEIN)
  {
    std::vector<Variant> merged;
    for (std::vector<Variant>::iterator it = variant.begin(); it != variant.end(); ++it)
    {
      merged.push_back(*it);
      if (it->type == SUBSTITUTION)
      {
        std::vector<Variant> annotation;
        extractor_frame_shift(annotation, reference, it->reference_start, it->reference_end, sample, it->sample_start, it->sample_end);
        merged.insert(merged.end(), annotation.begin(), annotation.end());
      } // if
    } // for
    variant = merged;
  } // if


  // Do NOT forget to clean up the complement string.
  delete[] complement;

  return weight;
} // extract

// This is the recursive extractor function. It works as follows:
// First, determine the ``best fitting'' longest common substring
// (LCS) (possibly as a reverse complement) and discard it from the
// solution. Then apply the same function on the remaining prefixes
// and suffixes. If there is no LCS it is a variant (region of
// change), i.e., a deletion/insertion.
// With regard to the reverse complement: the complement string is, as
// its name suggests, just the complement (DNA/RNA) of the reference
// string but it is NOT reversed.
size_t extractor(std::vector<Variant> &variant,
                 char_t const* const   reference,
                 char_t const* const   complement,
                 size_t           reference_start,
                 size_t           reference_end,
                 char_t const* const   sample,
                 size_t           sample_start,
                 size_t           sample_end)
{
  // First do prefix and suffix matching on the MASK character
  size_t i = 0;
  while (reference_start + i < reference_end && reference[reference_start + i] == MASK)
  {
    ++i;
  } // while
  reference_start += i;
  i = 0;
  while (reference_end - i - 1 > reference_start && reference[reference_end - i - 1] == MASK)
  {
    ++i;
  } // while
  reference_end -= i;

  i = 0;
  while (sample_start + i < sample_end && sample[sample_start + i] == MASK)
  {
    ++i;
  } // while
  sample_start += i;
  i = 0;
  while (sample_end - i - 1 > sample_start && sample[sample_end - i - 1] == MASK)
  {
    ++i;
  } // while
  sample_end -= i;


  size_t const reference_length = reference_end - reference_start;
  size_t const sample_length = sample_end - sample_start;

  // Assume this is a deletion/insertion.
  size_t const weight_trivial = weight_position + WEIGHT_DELETION_INSERTION + WEIGHT_BASE * sample_length + (reference_length != 1 ? weight_position + WEIGHT_SEPARATOR : 0);
  size_t weight = 0;


#if defined(__debug__)
  fputs("Extraction\n", stderr);
  fprintf(stderr, "  reference %ld--%ld:  ", reference_start, reference_end);
  Dprint_truncated(reference, reference_start, reference_end);
  fprintf(stderr, " (%ld)\n", reference_length);
  fprintf(stderr, "  complement %ld--%ld: ", reference_start, reference_end);
  Dprint_truncated(complement, reference_start, reference_end);
  fprintf(stderr, " (%ld)\n", reference_length);
  fprintf(stderr, "  sample %ld--%ld:     ", sample_start, sample_end);
  Dprint_truncated(sample, sample_start, sample_end);
  fprintf(stderr, " (%ld)\n", sample_length);
  fprintf(stderr, "  trivial weight: %ld\n", weight_trivial);
#endif


  // First some base cases to end the recursion.
  // No more reference string.
  if (reference_length <= 0)
  {
    // But some of the sample string is remaining: this is an
    // insertion or transposition.
    if (sample_length > 0)
    {
      weight = 2 * weight_position + WEIGHT_SEPARATOR + WEIGHT_INSERTION + WEIGHT_BASE * sample_length;

      // First, we check if we can match the inserted substring
      // somewhere in the complete reference string. This will
      // indicate a possible transposition. Otherwise it is a regular
      // insertion.
      std::vector<Variant> transposition;
      size_t const weight_transposition = extractor_transposition(transposition, reference, complement, reference_start, reference_end, sample, sample_start, sample_end, weight) + 2 * weight_position + 3 * WEIGHT_SEPARATOR + WEIGHT_INSERTION;


#if defined(__debug__)
  fprintf(stderr, "Transpositions: %ld (trivial: %ld)\n", weight_transposition, weight);
  for (std::vector<Variant>::const_iterator it = transposition.begin(); it != transposition.end(); ++it)
  {
    fprintf(stderr, "  %ld--%ld, %ld--%ld, %d, %ld, %ld--%ld\n", it->reference_start, it->reference_end, it->sample_start, it->sample_end, it->type, it->weight, it->transposition_start, it->transposition_end);
  } // for
#endif


      // Add transpositions if any.
      if (weight > weight_transposition && transposition.size() > 0 && !(transposition.size() == 1 && transposition.front().type == SUBSTITUTION))
      {
        transposition.front().type |= TRANSPOSITION_OPEN;
        transposition.back().type |= TRANSPOSITION_CLOSE;
        variant.insert(variant.end(), transposition.begin(), transposition.end());
        return weight_transposition;
      } // if

      // This is an actual insertion.
      variant.push_back(Variant(reference_start, reference_end, sample_start, sample_end, SUBSTITUTION, weight));
    } // if
    return weight;
  } // if

  // Obviously there is a piece of reference string left, but no more
  // sample string: this is a deletion.
  if (sample_length <= 0)
  {
    weight = weight_position + WEIGHT_DELETION + (reference_length > 1 ? weight_position + WEIGHT_SEPARATOR : 0);
    variant.push_back(Variant(reference_start, reference_end, sample_start, sample_end, SUBSTITUTION, weight));
    return weight;
  } // if

  // Single-nucleotide polymorphism (SNP): a special case for HGVS.
  if (reference_length == 1 && sample_length == 1)
  {
    weight = weight_position + 2 * WEIGHT_BASE + WEIGHT_SUBSTITUTION;
    variant.push_back(Variant(reference_start, reference_end, sample_start, sample_end, SUBSTITUTION, weight));
    return weight;
  } // if


  // Calculate the LCS (possibly in reverse complement) of the two
  // strings.
  size_t const cut_off = reference_length < THRESHOLD_CUT_OFF ? 1 : weight_position;
  std::vector<Substring> substring;
  size_t const length = LCS(substring, reference, complement, reference_start, reference_end, sample, sample_start, sample_end, cut_off);


  // No LCS found: this is a transposition or a deletion/insertion.
  if (length <= 0 || substring.size() <= 0)
  {
    weight = weight_trivial;

    // First, we check if we can match the inserted substring
    // somewhere in the complete reference string. This will
    // indicate a possible transposition.
    std::vector<Variant> transposition;
    size_t const weight_transposition = extractor_transposition(transposition, reference, complement, reference_start, reference_end, sample, sample_start, sample_end, weight)  + 2 * weight_position + 3 * WEIGHT_SEPARATOR + WEIGHT_DELETION_INSERTION;


#if defined(__debug__)
  fprintf(stderr, "Transpositions: %ld (trivial: %ld)\n", weight_transposition, weight);
  for (std::vector<Variant>::const_iterator it = transposition.begin(); it != transposition.end(); ++it)
  {
    fprintf(stderr, "  %ld--%ld, %ld--%ld, %d, %ld, %ld--%ld\n", it->reference_start, it->reference_end, it->sample_start, it->sample_end, it->type, it->weight, it->transposition_start, it->transposition_end);
  } // for
#endif


    // Add transpositions if any.
    if (weight > weight_transposition && transposition.size() > 0 && !(transposition.size() == 1 && transposition.front().type == SUBSTITUTION))
    {
      transposition.front().type |= TRANSPOSITION_OPEN;
      transposition.back().type |= TRANSPOSITION_CLOSE;
      variant.insert(variant.end(), transposition.begin(), transposition.end());
      return weight_transposition;
    } // if

    // This is an actual deletion/insertion.
    variant.push_back(Variant(reference_start, reference_end, sample_start, sample_end, SUBSTITUTION, weight_trivial));
    return weight_trivial;
  } // if


  // Pick the ``best fitting'' LCS, i.e., the remaining prefixes and
  // suffixes are of the same length.
  size_t diff = (reference_end - reference_start) + (sample_end - sample_start);
  std::vector<Substring>::const_iterator lcs = substring.begin();
  for (std::vector<Substring>::const_iterator it = substring.begin(); it != substring.end(); ++it)
  {
    size_t const prefix_diff = abs((it->reference_index - reference_start) - (it->sample_index - sample_start));
    size_t const suffix_diff = abs((reference_end - (it->reference_index + it->length)) - (sample_end - (it->sample_index + it->length)));
    if (prefix_diff + suffix_diff < diff)
    {
      // A better fitting LCS.
      diff = prefix_diff + suffix_diff;
      lcs = it;
    } // if
  } // for


  // Add some weight for reverse complement LCS.
  if (lcs->reverse_complement)
  {
    weight = 2 * weight_position + WEIGHT_SEPARATOR + WEIGHT_INVERSION;
  } // if


#if defined(__debug__)
  fprintf(stderr, "  LCS (x%ld)\n", substring.size());
  for (std::vector<Substring>::const_iterator it = substring.begin() ; it != substring.end(); ++it)
  {
    if (!it->reverse_complement)
    {
      fprintf(stderr, "    %ld--%ld: ", it->reference_index, it->reference_index + length);
      Dprint_truncated(reference, it->reference_index, it->reference_index + length);
      fprintf(stderr, " (%ld)\n    %ld--%ld: ", length, it->sample_index, it->sample_index + length);
      Dprint_truncated(sample, it->sample_index, it->sample_index + length);
      fprintf(stderr, " (%ld)", length);
    } // if
    else
    {
      fprintf(stderr, "    %ld--%ld: ", it->reference_index, it->reference_index + length);
      Dprint_truncated(complement, it->reference_index, it->reference_index + length);
      fprintf(stderr, " (%ld), (reverse complement)\n    %ld--%ld: ", length, it->sample_index, it->sample_index + length);
      Dprint_truncated(sample, it->sample_index, it->sample_index + length);
      fprintf(stderr, " (%ld)", length);
    } // else
    fputs("\n", stderr);
  } // for
#endif


  // Recursively apply this function to the prefixes of the strings.
  std::vector<Variant> prefix;
  weight += extractor(prefix, reference, complement, reference_start, lcs->reference_index, sample, sample_start, lcs->sample_index);

  // Stop if the weight of the variant exeeds the trivial weight.
  if (weight > weight_trivial)
  {
    weight = weight_trivial;

    // First, we check if we can match the inserted substring
    // somewhere in the complete reference string. This will
    // indicate a possible transposition.
    std::vector<Variant> transposition;
    size_t const weight_transposition = extractor_transposition(transposition, reference, complement, reference_start, reference_end, sample, sample_start, sample_end, weight)  + 2 * weight_position + 3 * WEIGHT_SEPARATOR + WEIGHT_DELETION_INSERTION;


#if defined(__debug__)
  fprintf(stderr, "Transpositions: %ld (trivial: %ld)\n", weight_transposition, weight);
  for (std::vector<Variant>::const_iterator it = transposition.begin(); it != transposition.end(); ++it)
  {
    fprintf(stderr, "  %ld--%ld, %ld--%ld, %d, %ld, %ld--%ld\n", it->reference_start, it->reference_end, it->sample_start, it->sample_end, it->type, it->weight, it->transposition_start, it->transposition_end);
  } // for
#endif


    // Add transpositions if any.
    if (weight > weight_transposition && transposition.size() > 0 && !(transposition.size() == 1 && transposition.front().type == SUBSTITUTION))
    {
      transposition.front().type |= TRANSPOSITION_OPEN;
      transposition.back().type |= TRANSPOSITION_CLOSE;
      variant.insert(variant.end(), transposition.begin(), transposition.end());
      return weight_transposition;
    } // if

    // This is an actual deletion/insertion.
    variant.push_back(Variant(reference_start, reference_end, sample_start, sample_end, SUBSTITUTION, weight_trivial));
    return weight_trivial;
  } // if

  // Recursively apply this function to the suffixes of the strings.
  std::vector<Variant> suffix;
  weight += extractor(suffix, reference, complement, lcs->reference_index + length, reference_end, sample, lcs->sample_index + length, sample_end);

  // Stop if the weight of the variant exeeds the trivial weight.
  if (weight > weight_trivial)
  {
    weight = weight_trivial;

    // First, we check if we can match the inserted substring
    // somewhere in the complete reference string. This will
    // indicate a possible transposition.
    std::vector<Variant> transposition;
    size_t const weight_transposition = extractor_transposition(transposition, reference, complement, reference_start, reference_end, sample, sample_start, sample_end, weight)  + 2 * weight_position + 3 * WEIGHT_SEPARATOR + WEIGHT_DELETION_INSERTION;


#if defined(__debug__)
  fprintf(stderr, "Transpositions: %ld (trivial: %ld)\n", weight_transposition, weight);
  for (std::vector<Variant>::const_iterator it = transposition.begin(); it != transposition.end(); ++it)
  {
    fprintf(stderr, "  %ld--%ld, %ld--%ld, %d, %ld, %ld--%ld\n", it->reference_start, it->reference_end, it->sample_start, it->sample_end, it->type, it->weight, it->transposition_start, it->transposition_end);
  } // for
#endif


    // Add transpositions if any.
    if (weight > weight_transposition && transposition.size() > 0 && !(transposition.size() == 1 && transposition.front().type == SUBSTITUTION))
    {
      transposition.front().type |= TRANSPOSITION_OPEN;
      transposition.back().type |= TRANSPOSITION_CLOSE;
      variant.insert(variant.end(), transposition.begin(), transposition.end());
      return weight_transposition;
    } // if

    // This is an actual deletion/insertion.
    variant.push_back(Variant(reference_start, reference_end, sample_start, sample_end, SUBSTITUTION, weight_trivial));
    return weight_trivial;
  } // if


  // Add all variants (in order) to the variant vector.
  variant.insert(variant.end(), prefix.begin(), prefix.end());

  if (!lcs->reverse_complement)
  {
    variant.push_back(Variant(lcs->reference_index, lcs->reference_index + length, lcs->sample_index, lcs->sample_index + length));
  } // if
  else
  {
    variant.push_back(Variant(lcs->reference_index, lcs->reference_index + length, lcs->sample_index, lcs->sample_index + length, REVERSE_COMPLEMENT, 2 * weight_position + WEIGHT_SEPARATOR + WEIGHT_INVERSION));
  } // else

  variant.insert(variant.end(), suffix.begin(), suffix.end());

  return weight;
} // extractor

// This function tries to extract transpositions from inserted
// sequences (insertions or deletion/insertions). Again we use a
// recursive method: extract the LCS and apply to the remaining prefix
// and suffix.
size_t extractor_transposition(std::vector<Variant> &variant,
                               char_t const* const   reference,
                               char_t const* const   complement,
                               size_t const          reference_start,
                               size_t const          reference_end,
                               char_t const* const   sample,
                               size_t const          sample_start,
                               size_t const          sample_end,
                               size_t const          weight_trivial)
{
  size_t const sample_length = sample_end - sample_start;

  size_t weight = 0;


#if defined(__debug__)
  fputs("Transposition extraction\n", stderr);
  fprintf(stderr, "  reference %ld--%ld:  ", 0ul, global_reference_length);
  Dprint_truncated(reference, 0, global_reference_length);
  fprintf(stderr, " (%ld)\n", global_reference_length);
  fprintf(stderr, "  complement %ld--%ld: ", 0ul, global_reference_length);
  Dprint_truncated(complement, 0, global_reference_length);
  fprintf(stderr, " (%ld)\n", global_reference_length);
  fprintf(stderr, "  sample %ld--%ld:     ", sample_start, sample_end);
  Dprint_truncated(sample, sample_start, sample_end);
  fprintf(stderr, " (%ld)\n", sample_length);
  fprintf(stderr, "  trivial weight: %ld\n", weight_trivial);
#endif


  // End transposition extraction if no more of the sample string
  // remains.
  if (sample_length <= 0)
  {
    return weight;
  } // if


  // Only consider large enough inserted regions (>> 1), based on
  // (average) description length of a position, otherwise it is just
  // a deletion/insertion.
  if (sample_length <= 2 * weight_position)
  {
    weight = sample_length * WEIGHT_BASE;
    variant.push_back(Variant(reference_start, reference_end, sample_start, sample_end, SUBSTITUTION, weight));
    return weight;
  } // if


  // Extract the LCS (from the whole reference string).
  size_t const cut_off = global_reference_length < THRESHOLD_CUT_OFF ? 1 : TRANSPOSITION_CUT_OFF * sample_length;
  std::vector<Substring> substring;
  size_t const length = LCS(substring, reference, complement, 0, global_reference_length, sample, sample_start, sample_end, cut_off);


  // No LCS found: this is a deletion/insertion.
  if (length <= 0 || substring.size() <= 0)
  {
    weight = sample_length * WEIGHT_BASE;
    variant.push_back(Variant(reference_start, reference_end, sample_start, sample_end, SUBSTITUTION, weight));
    return weight;
  } // if

  // We do NOT have to bother about a ``best fitting'' LCS here.
  std::vector<Substring>::const_iterator const lcs = substring.begin();

  // Update the weight of the transposition.
  weight += 2 * weight_position + WEIGHT_SEPARATOR;
  if (lcs->reverse_complement)
  {
    weight += WEIGHT_INVERSION;
  } // if


#if defined(__debug__)
  fprintf(stderr, "  LCS (x%ld)\n", substring.size());
  for (std::vector<Substring>::const_iterator it = substring.begin() ; it != substring.end(); ++it)
  {
    if (!it->reverse_complement)
    {
      fprintf(stderr, "    %ld--%ld: ", it->reference_index, it->reference_index + length);
      Dprint_truncated(reference, it->reference_index, it->reference_index + length);
      fprintf(stderr, " (%ld)\n    %ld--%ld: ", length, it->sample_index, it->sample_index + length);
      Dprint_truncated(sample, it->sample_index, it->sample_index + length);
      fprintf(stderr, " (%ld)", length);
    } // if
    else
    {
      fprintf(stderr, "    %ld--%ld: ", it->reference_index, it->reference_index + length);
      Dprint_truncated(complement, it->reference_index, it->reference_index + length);
      fprintf(stderr, " (%ld), (reverse complement)\n    %ld--%ld: ", length, it->sample_index, it->sample_index + length);
      Dprint_truncated(sample, it->sample_index, it->sample_index + length);
      fprintf(stderr, " (%ld)", length);
    } // else
    fputs("\n", stderr);
  } // for
#endif


  // Recursively apply this function to the prefixes of the strings
  std::vector<Variant> prefix;
  weight += extractor_transposition(prefix, reference, complement, reference_start, reference_end, sample, sample_start, lcs->sample_index, lcs->sample_index - sample_start) + WEIGHT_SEPARATOR;

  // Stop if the weight of the variant exeeds the trivial weight.
  if (weight > weight_trivial)
  {
    weight = sample_length * WEIGHT_BASE;
    variant.push_back(Variant(reference_start, reference_end, sample_start, sample_end, SUBSTITUTION, weight));
    return weight;
  } // if

  // Recursively apply this function to the suffixes of the strings.
  std::vector<Variant> suffix;
  weight += extractor_transposition(suffix, reference, complement, reference_start, reference_end, sample, lcs->sample_index + length, sample_end, sample_end - (lcs->sample_index + length)) + WEIGHT_SEPARATOR;

  // Stop if the weight of the variant exeeds the trivial weight.
  if (weight > weight_trivial)
  {
    weight = sample_length * WEIGHT_BASE;
    variant.push_back(Variant(reference_start, reference_end, sample_start, sample_end, SUBSTITUTION, weight));
    return weight;
  } // if


  // Add all variants (in order) to the variant vector.
  variant.insert(variant.end(), prefix.begin(), prefix.end());

  if (!lcs->reverse_complement)
  {
    variant.push_back(Variant(reference_start, reference_end, lcs->sample_index, lcs->sample_index + length, IDENTITY, 2 * weight_position + WEIGHT_SEPARATOR, lcs->reference_index, lcs->reference_index + length));
  } // if
  else
  {
    variant.push_back(Variant(reference_start, reference_end, lcs->sample_index, lcs->sample_index + length, REVERSE_COMPLEMENT, 2 * weight_position + WEIGHT_SEPARATOR + WEIGHT_INVERSION, lcs->reference_index, lcs->reference_index + length));
  } // else

  variant.insert(variant.end(), suffix.begin(), suffix.end());

  return weight;
} // extractor_transposition

// This is the recursive protein extractor function. It works as the
// regular extractor function, but no reverse complements nor
// transposion matching is used.
size_t extractor_protein(std::vector<Variant> &variant,
                         char_t const* const   reference,
                         size_t const          reference_start,
                         size_t const          reference_end,
                         char_t const* const   sample,
                         size_t const          sample_start,
                         size_t const          sample_end)
{
  size_t const reference_length = reference_end - reference_start;
  size_t const sample_length = sample_end - sample_start;

  // Assume this is a deletion/insertion.
  size_t const weight_trivial = weight_position + WEIGHT_DELETION_INSERTION + WEIGHT_BASE * sample_length + (reference_length != 1 ? weight_position + WEIGHT_SEPARATOR : 0);
  size_t weight = 0;


#if defined(__debug__)
  fputs("Extraction (protein)\n", stderr);
  fprintf(stderr, "  reference %ld--%ld:  ", reference_start, reference_end);
  Dprint_truncated(reference, reference_start, reference_end);
  fprintf(stderr, " (%ld)\n", reference_length);
  fprintf(stderr, "  sample %ld--%ld:     ", sample_start, sample_end);
  Dprint_truncated(sample, sample_start, sample_end);
  fprintf(stderr, " (%ld)\n", sample_length);
  fprintf(stderr, "  trivial weight: %ld\n", weight_trivial);
#endif


  // First some base cases to end the recursion.
  // No more reference string.
  if (reference_length <= 0)
  {
    // But some of the sample string is remaining: this is an
    // insertion.
    if (sample_length > 0)
    {
      weight = 2 * weight_position + WEIGHT_SEPARATOR + WEIGHT_INSERTION + WEIGHT_BASE * sample_length;
      variant.push_back(Variant(reference_start, reference_end, sample_start, sample_end, SUBSTITUTION, weight));
    } // if
    return weight;
  } // if

  // Obviously there is a piece of reference string left, but no more
  // sample string: this is a deletion.
  if (sample_length <= 0)
  {
    weight = weight_position + WEIGHT_DELETION + (reference_length > 1 ? weight_position + WEIGHT_SEPARATOR : 0);
    variant.push_back(Variant(reference_start, reference_end, sample_start, sample_end, SUBSTITUTION, weight));
    return weight;
  } // if

  // Single substitution: a special case for HGVS.
  if (reference_length == 1 && sample_length == 1)
  {
    weight = weight_position + 2 * WEIGHT_BASE + WEIGHT_SUBSTITUTION;
    variant.push_back(Variant(reference_start, reference_end, sample_start, sample_end, SUBSTITUTION, weight));
    return weight;
  } // if


  // Calculate the LCS of the two strings.
  std::vector<Substring> substring;
  size_t const length = LCS_1(substring, reference, 0, reference_start, reference_end, sample, sample_start, sample_end);


  // No LCS found: this is a deletion/insertion.
  if (length <= 0 || substring.size() <= 0)
  {
    weight = weight_trivial;

    // This is an actual deletion/insertion.
    variant.push_back(Variant(reference_start, reference_end, sample_start, sample_end, SUBSTITUTION, weight_trivial));
    return weight_trivial;
  } // if


  // Pick the ``best fitting'' LCS, i.e., the remaining prefixes and
  // suffixes are of the same length.
  size_t diff = (reference_end - reference_start) + (sample_end - sample_start);
  std::vector<Substring>::const_iterator lcs = substring.begin();
  for (std::vector<Substring>::const_iterator it = substring.begin(); it != substring.end(); ++it)
  {
    size_t const prefix_diff = abs((it->reference_index - reference_start) - (it->sample_index - sample_start));
    size_t const suffix_diff = abs((reference_end - (it->reference_index + it->length)) - (sample_end - (it->sample_index + it->length)));
    if (prefix_diff + suffix_diff < diff)
    {
      // A better fitting LCS.
      diff = prefix_diff + suffix_diff;
      lcs = it;
    } // if
  } // for


#if defined(__debug__)
  fprintf(stderr, "  LCS (x%ld)\n", substring.size());
  for (std::vector<Substring>::const_iterator it = substring.begin() ; it != substring.end(); ++it)
  {
    fprintf(stderr, "    %ld--%ld: ", it->reference_index, it->reference_index + length);
    Dprint_truncated(reference, it->reference_index, it->reference_index + length);
    fprintf(stderr, " (%ld)\n    %ld--%ld: ", length, it->sample_index, it->sample_index + length);
    Dprint_truncated(sample, it->sample_index, it->sample_index + length);
    fprintf(stderr, " (%ld)", length);
    fputs("\n", stderr);
  } // for
#endif


  // Recursively apply this function to the prefixes of the strings.
  std::vector<Variant> prefix;
  weight += extractor_protein(prefix, reference, reference_start, lcs->reference_index, sample, sample_start, lcs->sample_index);

  // Stop if the weight of the variant exeeds the trivial weight.
  if (weight > weight_trivial)
  {
    weight = weight_trivial;

    // This is an actual deletion/insertion.
    variant.push_back(Variant(reference_start, reference_end, sample_start, sample_end, SUBSTITUTION, weight_trivial));
    return weight_trivial;
  } // if

  // Recursively apply this function to the suffixes of the strings.
  std::vector<Variant> suffix;
  weight += extractor_protein(suffix, reference, lcs->reference_index + length, reference_end, sample, lcs->sample_index + length, sample_end);

  // Stop if the weight of the variant exeeds the trivial weight.
  if (weight > weight_trivial)
  {
    weight = weight_trivial;

    // This is an actual deletion/insertion.
    variant.push_back(Variant(reference_start, reference_end, sample_start, sample_end, SUBSTITUTION, weight_trivial));
    return weight_trivial;
  } // if


  // Add all variants (in order) to the variant vector.
  variant.insert(variant.end(), prefix.begin(), prefix.end());
  variant.push_back(Variant(lcs->reference_index, lcs->reference_index + length, lcs->sample_index, lcs->sample_index + length));
  variant.insert(variant.end(), suffix.begin(), suffix.end());

  return weight;
} // extractor_protein

void extractor_frame_shift(std::vector<Variant> &annotation,
                           char_t const* const   reference,
                           size_t const          reference_start,
                           size_t const          reference_end,
                           char_t const* const   sample,
                           size_t const          sample_start,
                           size_t const          sample_end)
{
  size_t const reference_length = reference_end - reference_start;
  size_t const sample_length = sample_end - sample_start;


#if defined(__debug__)
  fputs("Extraction (frame shift)\n", stderr);
  fprintf(stderr, "  reference %ld--%ld:  ", reference_start, reference_end);
  Dprint_truncated(reference, reference_start, reference_end);
  fprintf(stderr, " (%ld)\n", reference_length);
  fprintf(stderr, "  sample %ld--%ld:     ", sample_start, sample_end);
  Dprint_truncated(sample, sample_start, sample_end);
  fprintf(stderr, " (%ld)\n", sample_length);
#endif


  // First the base cases to end the recursion.
  if (reference_length <= 0 || sample_length <= 0)
  {
    return;
  } // if


  // Calculate the frame shift LCS of the two strings.
  std::vector<Substring> substring;
  LCS_frame_shift(substring, reference, reference_start, reference_end, sample, sample_start, sample_end);


  // Pick the ``best fitting'' frame shift LCS, i.e., pushed as far to
  // the start of the reference string as possible. Also update in
  // case of compound frame shifts.
  Substring lcs(0, 0, 0, FRAME_SHIFT_NONE);
  for (size_t i = 0; i < 5; ++i)
  {
    if (substring[i].length > lcs.length ||
        (substring[i].length == lcs.length && substring[i].reference_index < lcs.reference_index))
    {
      lcs = substring[i];
    } // if
    // Update compound frame shifts, e.g.,
    // FRAME_SHIFT_1 | FRAME_SHIFT_2
    else if (substring[i].length == lcs.length &&
             substring[i].reference_index == lcs.reference_index &&
             substring[i].sample_index == lcs.sample_index)
    {
      lcs.type |= substring[i].type;
    } // if
  } // for


  // No LCS found: no frame shift annotation.
  if (lcs.length <= 0)
  {
    return;
  } // if


#if defined(__debug__)
  fprintf(stderr, "  LCS type = %d\n", lcs.type);
  fprintf(stderr, "    %ld--%ld: ", lcs.reference_index, lcs.reference_index + lcs.length);
  Dprint_truncated(reference, lcs.reference_index, lcs.reference_index + lcs.length);
  fprintf(stderr, " (%ld)\n    %ld--%ld: ", lcs.length, lcs.sample_index, lcs.sample_index + lcs.length);
  Dprint_truncated(sample, lcs.sample_index, lcs.sample_index + lcs.length);
  fprintf(stderr, " (%ld)", lcs.length);
  fputs("\n", stderr);
#endif


  // Calculate the frame shift probability.
  double probability = 1.f;
  for (size_t i = 0; i < lcs.length; ++i)
  {
    double probability_compound = .0f;
    if ((lcs.type & FRAME_SHIFT_1) == FRAME_SHIFT_1)
    {
      probability_compound += frame_shift_frequency[reference[lcs.reference_index + i] & 0x7f][reference[lcs.reference_index + i + 1] & 0x7f][0];
    } // if
    if ((lcs.type & FRAME_SHIFT_2) == FRAME_SHIFT_2)
    {
      probability_compound += frame_shift_frequency[reference[lcs.reference_index + i] & 0x7f][reference[lcs.reference_index + i + 1] & 0x7f][1];
    } // if
    if ((lcs.type & FRAME_SHIFT_REVERSE) == FRAME_SHIFT_REVERSE)
    {
      probability_compound += frame_shift_frequency[reference[lcs.reference_index + i] & 0x7f][reference[lcs.reference_index + i] & 0x7f][2];
    } // if
    if ((lcs.type & FRAME_SHIFT_REVERSE_1) == FRAME_SHIFT_REVERSE_1)
    {
      probability_compound += frame_shift_frequency[reference[lcs.reference_index + i] & 0x7f][reference[lcs.reference_index + i + 1] & 0x7f][3];
    } // if
    if ((lcs.type & FRAME_SHIFT_REVERSE_2) == FRAME_SHIFT_REVERSE_2)
    {
      probability_compound += frame_shift_frequency[reference[lcs.reference_index + i] & 0x7f][reference[lcs.reference_index + i + 1] & 0x7f][4];
    } // if
    probability *= probability_compound;
  } // for


  // Recursively apply this function to the prefixes of the strings.
  std::vector<Variant> prefix;
  extractor_frame_shift(prefix, reference, reference_start, lcs.reference_index, sample, sample_start, lcs.sample_index);


  // Recursively apply this function to the suffixes of the strings.
  std::vector<Variant> suffix;
  extractor_frame_shift(suffix, reference, lcs.reference_index + lcs.length, reference_end, sample, lcs.sample_index + lcs.length, sample_end);


  // Add all variants (in order) to the annotation vector.
  annotation.insert(annotation.end(), prefix.begin(), prefix.end());
  Variant variant(lcs.reference_index, lcs.reference_index + lcs.length, lcs.sample_index, lcs.sample_index + lcs.length, FRAME_SHIFT | lcs.type);
  variant.probability = probability;
  annotation.push_back(variant);
  annotation.insert(annotation.end(), suffix.begin(), suffix.end());

  return;
} // extractor_frame_shift

// This function calculates the LCS using the LCS_k function by
// choosing an initial k and reducing it if necessary until the
// strings represent random strings modeled by a threshold value.
size_t LCS(std::vector<Substring> &substring,
           char_t const* const     reference,
           char_t const* const     complement,
           size_t const            reference_start,
           size_t const            reference_end,
           char_t const* const     sample,
           size_t const            sample_start,
           size_t const            sample_end,
           size_t const            cut_off)
{
  size_t const reference_length = reference_end - reference_start;
  size_t const sample_length = sample_end - sample_start;


  // The initial k.
  size_t k = reference_length > sample_length ? sample_length / 8 : reference_length / 8;


  // Reduce k until the cut-off is reached.
  while (k > 8 && k > cut_off)
  {


#if defined(__debug__)
  fprintf(stderr, "  k = %ld (cut-off: %ld)\n", k, cut_off);
#endif


    // Try to find a LCS with k.
    substring.clear();
    size_t const length = LCS_k(substring, reference, complement, reference_start, reference_end, sample, sample_start, sample_end, k);

    // A LCS of sufficient length has been found.
    if (length >= 2 * k && substring.size() > 0)
    {
      return length;
    } // if
    k /= 3;
  } // while

  // Cut-off: no LCS found.
  if (cut_off > 1)
  {


#if defined(__debug__)
  fprintf(stderr, "  cut-off\n");
#endif


    substring.clear();
    return 0;
  } // if


#if defined(__debug__)
  fprintf(stderr, "  k = 1\n");
#endif


  // As a last resort try running the classical algorithm.
  return LCS_1(substring, reference, complement, reference_start, reference_end, sample, sample_start, sample_end);
} // LCS

// Calculate the LCS in the well-known way using dynamic programming.
// NOT suitable for large strings.
size_t LCS_1(std::vector<Substring> &substring,
             char_t const* const     reference,
             char_t const* const     complement,
             size_t const            reference_start,
             size_t const            reference_end,
             char_t const* const     sample,
             size_t const            sample_start,
             size_t const            sample_end)
{
  size_t const reference_length = reference_end - reference_start;
  size_t const sample_length = sample_end - sample_start;
  bool reverse_complement = false;

  // Just a fancy way of allocation a continuous 2D array in heap
  // space.
  typedef size_t array[2][reference_length];
  array &LCS_line = *(reinterpret_cast<array*>(new size_t[2 * reference_length]));
  array &LCS_line_rc = *(reinterpret_cast<array*>(new size_t[2 * reference_length]));

  size_t length = 0;

  // Filling the LCS matrix (actually only the current and the
  // previous row).
  for (size_t i = 0; i < sample_length; ++i)
  {
    for (size_t j = 0; j < reference_length; ++j)
    {
      // A match
      if (reference[reference_start + j] == sample[sample_start + i] && reference[reference_start + j] != MASK)
      {
        if (i == 0 || j == 0)
        {
          LCS_line[i % 2][j] = 1;
        } // if
        else
        {
          LCS_line[i % 2][j] = LCS_line[(i + 1) % 2][j - 1] + 1;
        } // else

        // Check for a new maximal length.
        if (LCS_line[i % 2][j] >= length)
        {
          if (reverse_complement || LCS_line[i % 2][j] > length)
          {
            length = LCS_line[i % 2][j];
            substring = std::vector<Substring>(1, Substring(j - length + reference_start + 1, i - length + sample_start + 1, length));
          } // if
          else
          {
            substring.push_back(Substring(j - length + reference_start + 1, i - length + sample_start + 1, length));
          } // else
          reverse_complement = false;
        } // if
      } // if
      else
      {
        LCS_line[i % 2][j] = 0;
      } // else

      // If applicable check for a LCS in reverse complement space.
      // The same code is used as before but the complement string is
      // travesed backwards (towards the start).
      if (complement != 0 && complement[reference_end - j - 1] == sample[sample_start + i] && complement[reference_end - j - 1] != MASK)
      {
        if (i == 0 || j == 0)
        {
          LCS_line_rc[i % 2][j] = 1;
        } // if
        else
        {
          LCS_line_rc[i % 2][j] = LCS_line_rc[(i + 1) % 2][j - 1] + 1;
        } // else

        if (LCS_line_rc[i % 2][j] > 1 && LCS_line_rc[i % 2][j] > length)
        {
          length = LCS_line_rc[i % 2][j];
          substring = std::vector<Substring>(1, Substring(reference_end - j - 1, i - length + sample_start + 1, length, true));
          reverse_complement = true;
        } // if
      } // if
      else
      {
        LCS_line_rc[i % 2][j] = 0;
      } // else


      // We can stop if the whole sample string is part of the LCS.
      if (!reverse_complement && length >= sample_length)
      {
        break;
      } // if

    } // for
  } // for

  // Cleaning up.
  delete[] &LCS_line;
  delete[] &LCS_line_rc;

  return length;
} // LCS_1

// Calculate the LCS using overlapping and non-overlapping k-mers.
// This function should be suitable for large (similar) strings.
// Be careful: if the resulting LCS is of length <= 2k it might not be
// the actual LCS. Remedy: try again with a reduced k.
size_t LCS_k(std::vector<Substring> &substring,
             char_t const* const     reference,
             char_t const* const     complement,
             size_t const            reference_start,
             size_t const            reference_end,
             char_t const* const     sample,
             size_t const            sample_start,
             size_t const            sample_end,
             size_t const            k)
{
  size_t length = 0;

  // Stop if we cannot partition the strings into k-mers.
  if (k <= 1 || reference_end - reference_start < k || sample_end - sample_start < k)
  {
    return length;
  } // if

  size_t const reference_length = (reference_end - reference_start) / k;
  size_t const sample_length = sample_end - sample_start - k + 1;
  bool reverse_complement = false;

  // Just a fancy way of allocation a continuous (k+1)D array in heap
  // space.
  typedef size_t array[k + 1][reference_length];
  array &LCS_line = *(reinterpret_cast<array*>(new size_t[(k + 1) * reference_length]));
  array &LCS_line_rc = *(reinterpret_cast<array*>(new size_t[(k + 1) * reference_length]));

  // Filling the LCS k-mer matrix (actually only the current and the k
  // previous rows). We count in k-mers.
  for (size_t i = 0; i < sample_length; ++i) // overlapping k-mers
  {
    for (size_t j = 0; j < reference_length; ++j) // non-overlapping
    {
      // A match
      if (string_match(reference + reference_start + j * k, sample + sample_start + i, k))
      {
        if (i < k || j == 0)
        {
          LCS_line[i % (k + 1)][j] = 1;
        } // if
        else
        {
          LCS_line[i % (k + 1)][j] = LCS_line[(i + 1) % (k + 1)][j - 1] + 1;
        } // else

        // Check for a new maximal length.
        if (LCS_line[i % (k + 1)][j] > length)
        {
          length = LCS_line[i % (k + 1)][j];

          // Remove all solutions with a length more than 1 from the
          // new maximal length. And remove also the partial LCS that
          // was extended to the new maximal length because it is
          // guaranteerd to be of maximal length - 1.
          for (std::vector<Substring>::iterator it = substring.begin(); it != substring.end(); ++it)
          {
            if (length - it->length > 1 || (it->reference_index == j - 1 && it->sample_index == i - k && !it->reverse_complement))
            {
              std::vector<Substring>::iterator const temp = it - 1;
              substring.erase(it);
              it = temp;
            } // if
          } // for
          substring.push_back(Substring(j, i, LCS_line[i % (k + 1)][j]));
        } // if
        else if (LCS_line[i % (k + 1)][j] > 0 && length - LCS_line[i % (k + 1)][j] <= 1)
        {
          substring.push_back(Substring(j, i, LCS_line[i % (k + 1)][j]));
        } // if
      } // if
      else
      {
        LCS_line[i % (k + 1)][j] = 0;
      } // else

      // If applicable check for a LCS in reverse complement space.
      // The same code is used as before but the complement string is
      // travesed backwards (towards the start).
      if (complement != 0 && string_match_reverse(complement + reference_end - j * k - 1, sample + sample_start + i, k))
      {
        if (i < k || j == 0)
        {
          LCS_line_rc[i % (k + 1)][j] = 1;
        } // if
        else
        {
          LCS_line_rc[i % (k + 1)][j] = LCS_line_rc[(i + 1) % (k + 1)][j - 1] + 1;
        } // else

        // Check for a new maximal length.
        if (LCS_line_rc[i % (k + 1)][j] > length)
        {
          length = LCS_line_rc[i % (k + 1)][j];

          // Remove all solutions with a length more than 1 from the
          // new maximal length. And remove also the partial LCS that
          // was extended to the new maximal length because it is
          // guaranteerd to be of maximal length - 1.
          for (std::vector<Substring>::iterator it = substring.begin(); it != substring.end(); ++it)
          {
            if (length - it->length > 1 || (it->reference_index == j - 1 && it->sample_index == i - k && it->reverse_complement))
            {
              std::vector<Substring>::iterator const temp = it - 1;
              substring.erase(it);
              it = temp;
            } // if
          } // for
          substring.push_back(Substring(j, i, LCS_line_rc[i % (k + 1)][j], true));
        } // if
        else if (LCS_line[i % (k + 1)][j] > 0 && length - LCS_line_rc[i % (k + 1)][j] <= 1)
        {
          substring.push_back(Substring(j, i, LCS_line[i % (k + 1)][j], true));
        } // if
      } // if
      else
      {
        LCS_line_rc[i % (k + 1)][j] = 0;
      } // else

    } // for
  } // for


  // Now we need to do some extensions of the found LCSs to find their
  // exact lengths. We extend up to k positions to the left (towards
  // the start of the strings) and to the right (towards the end of
  // the string).
  for (std::vector<Substring>::iterator it = substring.begin(); it != substring.end(); ++it)
  {
    if (!it->reverse_complement)
    {
      it->reference_index = reference_start + (it->reference_index - it->length + 1) * k;
      it->sample_index = sample_start + it->sample_index - (it->length - 1) * k;
      it->length *= k;
      // Extending to the right.
      {
        size_t i = 0;
        while (i <= k && it->reference_index + it->length + i < reference_end && it->sample_index + it->length + i < sample_end && reference[it->reference_index + it->length + i] == sample[it->sample_index + it->length + i] && reference[it->reference_index + it->length + i] != MASK)
        {
          ++i;
        } // while
        it->length += i;
      }
      // Extending to the left.
      {
        size_t i = 0;
        while (i <= k && it->reference_index - i - 1 >= reference_start && it->sample_index - i - 1 >= sample_start && reference[it->reference_index - i - 1] == sample[it->sample_index - i - 1] && reference[it->reference_index - i - 1] != MASK)
        {
          ++i;
        } // while
        it->reference_index -= i;
        it->sample_index -= i;
        it->length += i;
      }
    } // if
    else
    {
      it->reference_index = reference_end - (it->reference_index + 1) * k;
      it->sample_index = sample_start + it->sample_index - (it->length - 1) * k;
      it->length *= k;
      // Extending to the right (sample orientation).
      {
        size_t i = 0;
        while (i <= k && it->reference_index - i - 1 >= reference_start && it->sample_index + it->length + i < sample_end && complement[it->reference_index - i - 1] == sample[it->sample_index + it->length + i] && complement[it->reference_index - i - 1] != MASK)
        {
          ++i;
        } // while
        it->reference_index -= i;
        it->length += i;
      }
      // Extending to the left (sample orientation).
      {
        size_t i = 0;
        while (i <= k && it->reference_index + it->length + i < reference_end && it->sample_index - i - 1 >= sample_start && complement[it->reference_index + it->length + i] == sample[it->sample_index - i - 1] && complement[it->reference_index + it->length + i] != MASK)
        {
          ++i;
        } // while
        it->sample_index -= i;
        it->length += i;
      }
    } // else

    if (it->length > length)
    {
      length = it->length;
      reverse_complement = it->reverse_complement;
    } // if
    else if (reverse_complement && it->length == length)
    {
      reverse_complement = false;
    } // if
  } // for


  // Remove all sub-optimal LCSs from the solution.
  for (std::vector<Substring>::iterator it = substring.begin(); it != substring.end(); ++it)
  {
    if (it->length < length || reverse_complement != it->reverse_complement)
    {
      std::vector<Substring>::iterator const temp = it - 1;
      substring.erase(it);
      it = temp;
    } // if
  } // for

  // Cleaning up.
  delete[] &LCS_line;
  delete[] &LCS_line_rc;

  return length;
} // LCS_k

// This function calculates the frame shift LCS. The five possible
// frame shift LCSs are calculated separately. Picking the longest
// ``best'' fitting one is the reponsibility of the caller.
// This function is a version of the LCS_1 function (not suitable for
// very large strings).
void LCS_frame_shift(std::vector<Substring> &substring,
                     char_t const* const     reference,
                     size_t const            reference_start,
                     size_t const            reference_end,
                     char_t const* const     sample,
                     size_t const            sample_start,
                     size_t const            sample_end)
{
  size_t const reference_length = reference_end - reference_start;
  size_t const sample_length = sample_end - sample_start;

  size_t lcs[2][reference_length][5];
  for (size_t i = 0; i < reference_length; ++i)
  {
    lcs[0][i][0] = 0;
    lcs[0][i][1] = 0;
    lcs[0][i][2] = 0;
    lcs[0][i][3] = 0;
    lcs[0][i][4] = 0;
    lcs[1][i][0] = 0;
    lcs[1][i][1] = 0;
    lcs[1][i][2] = 0;
    lcs[1][i][3] = 0;
    lcs[1][i][4] = 0;
  } // for

  Substring fs_substring[5];
  for (size_t i = 0; i < sample_length; ++i)
  {
    uint8_t const shift_reverse = frame_shift(reference[reference_end - 1], reference[reference_end - 1], sample[sample_start + i]);
    if ((shift_reverse & FRAME_SHIFT_REVERSE) == FRAME_SHIFT_REVERSE)
    {
      lcs[i % 2][0][2] = 1;
    } // if
    if (lcs[i % 2][0][2] > fs_substring[2].length)
    {
      fs_substring[2] = Substring(reference_start - lcs[i % 2][0][2] + 1, sample_start + i - lcs[i % 2][0][2] + 1, lcs[i % 2][0][2], FRAME_SHIFT_REVERSE);
    } // if
    for (size_t j = 1; j < reference_length; ++j)
    {
      uint8_t const shift_forward = frame_shift(reference[reference_start + j - 1], reference[reference_start + j], sample[sample_start + i]);
      uint8_t const shift_reverse = frame_shift(reference[reference_end - j - 1], reference[reference_end - j], sample[sample_start + i]);
      if ((shift_forward & FRAME_SHIFT_1) == FRAME_SHIFT_1)
      {
        lcs[i % 2][j][0] = lcs[(i + 1) % 2][j - 1][0] + 1;
      } // if
      else
      {
        lcs[i % 2][j][0] = 0;
      } // else
      if ((shift_forward & FRAME_SHIFT_2) == FRAME_SHIFT_2)
      {
        lcs[i % 2][j][1] = lcs[(i + 1) % 2][j - 1][1] + 1;
      } // if
      else
      {
        lcs[i % 2][j][1] = 0;
      } // else
      if ((shift_reverse & FRAME_SHIFT_REVERSE) == FRAME_SHIFT_REVERSE)
      {
        lcs[i % 2][j][2] = lcs[(i + 1) % 2][j - 1][2] + 1;
      } // if
      else
      {
        lcs[i % 2][j][2] = 0;
      } // else
      if ((shift_reverse & FRAME_SHIFT_REVERSE_1) == FRAME_SHIFT_REVERSE_1)
      {
        lcs[i % 2][j][3] = lcs[(i + 1) % 2][j - 1][3] + 1;
      } // if
      else
      {
        lcs[i % 2][j][3] = 0;
      } // else
      if ((shift_reverse & FRAME_SHIFT_REVERSE_2) == FRAME_SHIFT_REVERSE_2)
      {
        lcs[i % 2][j][4] = lcs[(i + 1) % 2][j - 1][4] + 1;
      } // if
      else
      {
        lcs[i % 2][j][4] = 0;
      } // else

      if (lcs[i % 2][j][0] > fs_substring[0].length)
      {
        fs_substring[0] = Substring(reference_start + j - lcs[i % 2][j][0], sample_start + i - lcs[i % 2][j][0] + 1, lcs[i % 2][j][0], FRAME_SHIFT_1);
      } // if
      if (lcs[i % 2][j][1] > fs_substring[1].length)
      {
        fs_substring[1] = Substring(reference_start + j - lcs[i % 2][j][1], sample_start + i - lcs[i % 2][j][1] + 1, lcs[i % 2][j][1], FRAME_SHIFT_2);
      } // if
      if (lcs[i % 2][j][2] > fs_substring[2].length)
      {
        fs_substring[2] = Substring(reference_end - j - 1, sample_start + i - lcs[i % 2][j][2] + 1, lcs[i % 2][j][2], FRAME_SHIFT_REVERSE);
      } // if
      if (lcs[i % 2][j][3] > fs_substring[3].length)
      {
        fs_substring[3] = Substring(reference_end - j - 1, sample_start + i - lcs[i % 2][j][3] + 1, lcs[i % 2][j][3], FRAME_SHIFT_REVERSE_1);
      } // if
      if (lcs[i % 2][j][4] > fs_substring[4].length)
      {
        fs_substring[4] = Substring(reference_end - j - 1, sample_start + i - lcs[i % 2][j][4] + 1, lcs[i % 2][j][4], FRAME_SHIFT_REVERSE_2);
      } // if
    } // for
  } // for
  substring = std::vector<Substring>(1, fs_substring[0]);
  substring.push_back(fs_substring[1]);
  substring.push_back(fs_substring[2]);
  substring.push_back(fs_substring[3]);
  substring.push_back(fs_substring[4]);
  return;
} // LCS_frame_shift

// This function is more or less equivalent to C's strncmp, but it
// returns true iff both strings are the same.
bool string_match(char_t const* const string_1,
                  char_t const* const string_2,
                  size_t const        length)
{
  for (size_t i = 0; i < length; ++i)
  {
    if (string_1[i] != string_2[i] || string_1[i] == MASK)
    {
      return false;
    } // if
  } // for
  return true;
} // string_match

// This function is very similar to C's strncmp, but it traverses
// string_1 from end to start while traversing string_2 from start to
// end (useful for the reverse complement in DNA/RNA), and it returns
// true iff both strings are the same in their respective directions.
bool string_match_reverse(char_t const* const string_1,
                          char_t const* const string_2,
                          size_t const        length)
{
  for (size_t i = 0; i < length; ++i)
  {
    if (string_1[-i] != string_2[i] || string_1[-i] == MASK)
    {
      return false;
    } // if
  } // for
  return true;
} // string_match_reverse

// This function calculates the length (in characters) of the common
// prefix between two strings. The result of this function is also
// used in the suffix_match function.
size_t prefix_match(char_t const* const reference,
                    size_t const        reference_length,
                    char_t const* const sample,
                    size_t const        sample_length)
{
  size_t i = 0;

  // Traverse both strings towards the end as long as their characters
  // are equal. Do NOT exceed the length of the strings.
  while (i < reference_length && i < sample_length && reference[i] == sample[i] && reference[i] != MASK)
  {
    ++i;
  } // while
  return i;
} // prefix_match

// This function calculates the length (in characters) of the common
// suffix between two strings. It needs the calculated common prefix.
size_t suffix_match(char_t const* const reference,
                    size_t const        reference_length,
                    char_t const* const sample,
                    size_t const        sample_length,
                    size_t const        prefix)

{
  size_t i = 0;

  // Start at the end of both strings and traverse towards the start
  // as long as their characters are equal. Do not exceed the length
  // of the strings.
  while (i < reference_length - prefix && i < sample_length - prefix && reference[reference_length - i - 1] == sample[sample_length - i - 1] && reference[reference_length - i - 1] != MASK)
  {
    ++i;
  } // while
  return i;
} // suffix_match

// This function converts a IUPAC Nucleotide Acid Notation into its
// complement.
char_t IUPAC_base_complement(char_t const base)
{
  switch (base)
  {
    case 'A':
      return 'T';
    case 'C':
      return 'G';
    case 'G':
      return 'C';
    case 'T':
    case 'U':
      return 'A';
  } // switch
  return base;
} // IUPAC_base_complement

// This function converts a string in IUPAC Nucleotide Acid Notation
// into its complement. A new string is allocated, so deletion is
// the responsibility of the caller.
char_t const* IUPAC_complement(char_t const* const string,
                               size_t const        length)
{
  char_t* complement = new char_t[length];
  for (size_t i = 0; i < length; ++i)
  {
    complement[i] = IUPAC_base_complement(string[i]);
  } // for
  return complement;
} // IUPAC_complement

void backtranslation(char_t              ref_DNA[],
                     char_t              alt_DNA[],
                     char_t const* const reference,
                     size_t const        reference_start,
                     char_t const* const sample,
                     size_t const        sample_start,
                     size_t const        length,
                     uint8_t const       type)
{
  size_t reference_DNA[3 * length];
  size_t sample_DNA[3 * length];
  for (size_t i = 0; i < 3 * length; i++)
  {
    reference_DNA[i] = 0;
    sample_DNA[i] = 0;
  } // for

  for (size_t p = 0; p < length; ++p)
  {
    for (size_t i = 0; i < 64; ++i)
    {
      if (((acid_map[reference[reference_start + p] & 0x7f] >> i) & 0x1ull) == 0x1ull)
      {
        size_t const codon_reverse = ((i >> 0x4) | (i & 0xc) | ((i & 0x3) << 0x4)) ^ 0x3f;
        for (size_t k = 0; k < 64; ++k)
        {
          if (((acid_map[sample[sample_start + length - p - 1] & 0x7f] >> k) & 0x1ull) == 0x1ull)
          {
            if ((type & FRAME_SHIFT_REVERSE) == FRAME_SHIFT_REVERSE && codon_reverse == k)
            {
              reference_DNA[p * 3] |= 0x1 << (i >> 4);
              reference_DNA[p * 3 + 1] |= 0x1 << ((i >> 2) & 0x3);
              reference_DNA[p * 3 + 2] |= 0x1 << (i & 0x3);
              sample_DNA[(length - p) * 3 - 3] |= 0x1 << (codon_reverse >> 4);
              sample_DNA[(length - p) * 3 - 2] |= 0x1 << ((codon_reverse >> 2) & 0x3);
              sample_DNA[(length - p) * 3 - 1] |= 0x1 << (codon_reverse & 0x3);
            } // if
          } // if
        } // for

        for (size_t j = 0; j < 64; ++j)
        {
          if (((acid_map[reference[reference_start + p + 1] & 0x7f] >> j) & 0x1ull) == 0x1ull)
          {
            size_t const codon_1 = ((i & 0x3) << 0x4) | ((j & 0x3c) >> 0x2);
            size_t const codon_2 = ((i & 0xf) << 0x2) | (j >> 0x4);
            size_t const codon_reverse_1 = (((i & 0xc) >> 0x2) | ((i & 0x3) << 0x2) | (j & 0x30)) ^ 0x3f;
            size_t const codon_reverse_2 = ((i & 0x3) | ((j & 0x30) >> 0x2) | ((j & 0xc) << 0x2)) ^ 0x3f;

            for (size_t k = 0; k < 64; ++k)
            {
              if (((acid_map[sample[sample_start + p] & 0x7f] >> k) & 0x1ull) == 0x1ull)
              {
                if ((type & FRAME_SHIFT_1) == FRAME_SHIFT_1 && codon_1 == k)
                {
                  reference_DNA[p * 3] |= 0x1 << (i >> 4);
                  reference_DNA[p * 3 + 1] |= 0x1 << ((i >> 2) & 0x3);
                  reference_DNA[p * 3 + 2] |= 0x1 << (i & 0x3);
                  sample_DNA[p * 3] |= 0x1 << (codon_1 >> 4);
                  sample_DNA[p * 3 + 1] |= 0x1 << ((codon_1 >> 2) & 0x3);
                  sample_DNA[p * 3 + 2] |= 0x1 << (codon_1 & 0x3);
                } // if
                if ((type & FRAME_SHIFT_2) == FRAME_SHIFT_2 && codon_2 == k)
                {
                  reference_DNA[p * 3] |= 0x1 << (i >> 4);
                  reference_DNA[p * 3 + 1] |= 0x1 << ((i >> 2) & 0x3);
                  reference_DNA[p * 3 + 2] |= 0x1 << (i & 0x3);
                  sample_DNA[p * 3] |= 0x1 << (codon_2 >> 4);
                  sample_DNA[p * 3 + 1] |= 0x1 << ((codon_2 >> 2) & 0x3);
                  sample_DNA[p * 3 + 2] |= 0x1 << (codon_2 & 0x3);
                } // if
              } // if

              if (((acid_map[sample[sample_start + length - p - 1] & 0x7f] >> k) & 0x1ull) == 0x1ull)
              {
                if ((type & FRAME_SHIFT_REVERSE_1) == FRAME_SHIFT_REVERSE_1 && codon_reverse_1 == k)
                {
                  reference_DNA[p * 3] |= 0x1 << (i >> 4);
                  reference_DNA[p * 3 + 1] |= 0x1 << ((i >> 2) & 0x3);
                  reference_DNA[p * 3 + 2] |= 0x1 << (i & 0x3);
                  sample_DNA[(length - p) * 3 - 3] |= 0x1 << (codon_reverse_1 >> 4);
                  sample_DNA[(length - p) * 3 - 2] |= 0x1 << ((codon_reverse_1 >> 2) & 0x3);
                  sample_DNA[(length - p) * 3 - 1] |= 0x1 << (codon_reverse_1 & 0x3);
                } // if
                if ((type & FRAME_SHIFT_REVERSE_2) == FRAME_SHIFT_REVERSE_2 && codon_reverse_2 == k)
                {
                  reference_DNA[p * 3] |= 0x1 << (i >> 4);
                  reference_DNA[p * 3 + 1] |= 0x1 << ((i >> 2) & 0x3);
                  reference_DNA[p * 3 + 2] |= 0x1 << (i & 0x3);
                  sample_DNA[(length - p) * 3 - 3] |= 0x1 << (codon_reverse_2 >> 4);
                  sample_DNA[(length - p) * 3 - 2] |= 0x1 << ((codon_reverse_2 >> 2) & 0x3);
                  sample_DNA[(length - p) * 3 - 1] |= 0x1 << (codon_reverse_2 & 0x3);
                } // if
              } // if
            } // for
          } // if
        } // for
      } // if
    } // for
  } // for
  for (size_t i = 0; i < 3 * length; i++)
  {
    ref_DNA[i] = IUPAC_ALPHA[reference_DNA[i]];
    alt_DNA[i] = IUPAC_ALPHA[sample_DNA[i]];
  } // for
  return;
} // backtranslation

static void initialize_acid_frequency(void)
{
  acid_frequency['A'] = .09515673f;
  acid_frequency['C'] = .01157279f;
  acid_frequency['D'] = .05151007f;
  acid_frequency['E'] = .05762795f;
  acid_frequency['F'] = .03890338f;
  acid_frequency['G'] = .07374416f;
  acid_frequency['H'] = .02266328f;
  acid_frequency['I'] = .06010209f;
  acid_frequency['K'] = .04406110f;
  acid_frequency['L'] = .10672657f;
  acid_frequency['M'] = .02819341f;
  acid_frequency['N'] = .03945573f;
  acid_frequency['P'] = .04425210f;
  acid_frequency['Q'] = .04439959f;
  acid_frequency['R'] = .05510809f;
  acid_frequency['S'] = .05802322f;
  acid_frequency['T'] = .05398938f;
  acid_frequency['U'] = .00000221f;
  acid_frequency['V'] = .07073316f;
  acid_frequency['W'] = .01531018f;
  acid_frequency['X'] = .00001106f;
  acid_frequency['Y'] = .02845373f;
  return;
} // initialize_acid_frequency

// This function precalculates the frame_shift_map and frequency count
// based on a given codon string.
void initialize_frame_shift_map(char_t const* const codon_string)
{
  for (size_t i = 0; i < 64; ++i)
  {
    acid_map[i] = 0x0ull;
  } // for
  for (size_t i = 0; i < 128; ++i)
  {
    for (size_t j = 0; j < 128; ++j)
    {
      for (size_t k = 0; k < 5; ++k)
      {
        frame_shift_count[i][j][k] = 0;
        frame_shift_frequency[i][j][k] = .05f;
      } // for
    } // for
  } // for
  initialize_acid_frequency();
  for (size_t i = 0; i < 64; ++i)
  {
    acid_map[codon_string[i] & 0x7f] |= (0x1ull << i);
  } // for
  for (size_t i = 0; i < 128; ++i)
  {
    if (acid_map[i] != 0x0ull)
    {
      for (size_t j = 0; j < 128; ++j)
      {
        if (acid_map[j] != 0x0ull)
        {
          for (size_t k = 0; k < 128; ++k)
          {
            if (acid_map[k] != 0x0ull)
            {
              uint8_t const shift = calculate_frame_shift(i, j, k);
              frame_shift_map[i][j][k] = shift;

              if ((shift & FRAME_SHIFT_1) == FRAME_SHIFT_1)
              {
                ++frame_shift_count[i][j][0];
                frame_shift_frequency[i][j][0] += acid_frequency[k];
              } // if
              if ((shift & FRAME_SHIFT_2) == FRAME_SHIFT_2)
              {
                ++frame_shift_count[i][j][1];
                frame_shift_frequency[i][j][1] += acid_frequency[k];
              } // if
              if ((shift & FRAME_SHIFT_REVERSE) == FRAME_SHIFT_REVERSE)
              {
                ++frame_shift_count[i][j][2];
                frame_shift_frequency[i][j][2] += acid_frequency[k];
              } // if
              if ((shift & FRAME_SHIFT_REVERSE_1) == FRAME_SHIFT_REVERSE_1)
              {
                ++frame_shift_count[i][j][3];
                frame_shift_frequency[i][j][3] += acid_frequency[k];
              } // if
              if ((shift & FRAME_SHIFT_REVERSE_2) == FRAME_SHIFT_REVERSE_2)
              {
                ++frame_shift_count[i][j][4];
                frame_shift_frequency[i][j][4] += acid_frequency[k];
              } // if

            } // if
          } // for
        } // if
      } // for
    } // if
  } // for
  return;
} // initialize_frame_shift_map

// Used to precalculate the frame_shift_map. It computes for all
// combinations of two reference amino acids the corresponding DNA
// sequence and the (partial) overlap between all possible DNA
// sequences of the sample amico acid.
uint8_t calculate_frame_shift(size_t const reference_1,
                              size_t const reference_2,
                              size_t const sample)
{
  uint8_t shift = FRAME_SHIFT_NONE;
  for (size_t i = 0; i < 64; ++i)
  {
    if (((acid_map[reference_1] >> i) & 0x1ull) == 0x1ull)
    {
      size_t const codon_reverse = ((i >> 0x4) | (i & 0xc) | ((i & 0x3) << 0x4)) ^ 0x3f;
      for (size_t j = 0; j < 64; ++j)
      {
        if (((acid_map[reference_2] >> j) & 0x1ull) == 0x1ull)
        {
          size_t const codon_1 = ((i & 0x3) << 0x4) | ((j & 0x3c) >> 0x2);
          size_t const codon_2 = ((i & 0xf) << 0x2) | (j >> 0x4);
          size_t const codon_reverse_1 = (((i & 0xc) >> 0x2) | ((i & 0x3) << 0x2) | (j & 0x30)) ^ 0x3f;
          size_t const codon_reverse_2 = ((i & 0x3) | ((j & 0x30) >> 0x2) | ((j & 0xc) << 0x2)) ^ 0x3f;
          for (size_t k = 0; k < 64; ++k)
          {
            if (((acid_map[sample] >> k) & 0x1ull) == 0x1ull)
            {
              if (codon_1 == k)
              {
                shift |= FRAME_SHIFT_1;
              } // if
              if (codon_2 == k)
              {
                shift |= FRAME_SHIFT_2;
              } // if
              if (codon_reverse == k)
              {
                shift |= FRAME_SHIFT_REVERSE;
              } // if
              if (codon_reverse_1 == k)
              {
                shift |= FRAME_SHIFT_REVERSE_1;
              } // if
              if (codon_reverse_2 == k)
              {
                shift |= FRAME_SHIFT_REVERSE_2;
              } // if
            } // if
          } // for
        } // if
      } // for
    } // if
  } // for
  return shift;
} // calculate_frame_shift

// This function calculates the frame shift. A reference amino acid is
// checked against two possible partial overlaps between every
// combination of two sample (observed) amino acids.
uint8_t frame_shift(char_t const reference_1,
                    char_t const reference_2,
                    char_t const sample)
{
  return frame_shift_map[reference_1 & 0x7f][reference_2 & 0x7f][sample & 0x7f];
} // frame_shift


#if defined(__debug__)
// Debug function for printing large strings in truncated form: a
// prefix of a certain length ... a suffix of the same length.
size_t Dprint_truncated(char_t const* const string,
                        size_t const        start,
                        size_t const        end,
                        size_t const        length,
                        FILE*               stream)
{
  if (end - start <= length - 3)
  {
    return fwrite(string + start, sizeof(char_t), end - start, stream);
  } // if
  return fwrite(string + start, sizeof(char_t), length / 2 - 1, stream) +
         fputs("...", stream) +
         fwrite(string + end - (length / 2) + 1, sizeof(char_t), length / 2 - 1, stream);
} // Dprint_truncated

size_t Dprint_codon(size_t const index,
                    FILE*        stream)
{
  return fprintf(stream, "%c%c%c", IUPAC_BASE[index >> 0x4],
                                   IUPAC_BASE[(index >> 0x2) & 0x3],
                                   IUPAC_BASE[index & 0x3]);
} // Dprint_codon
#endif


} // mutalyzer

