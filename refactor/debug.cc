#include "extractor.h"
using namespace mutalyzer;

#include <cstdio>

int main(int argc, char* argv[])
{

  if (argc < 3)
  {
    fprintf(stderr, "usage: %s reference sample\n", argv[0]);
    return 1;
  } // if

  fprintf(stderr, "HGVS description extractor\n");
  FILE* file = fopen(argv[1], "r");
  if (file == 0)
  {
    fprintf(stderr, "ERROR: could not open file `%s'\n", argv[1]);
    return 1;
  } // if
  fseek(file, 0, SEEK_END);
  size_t const reference_length = ftell(file) - 1;
  rewind(file);
  char_t* reference = new char_t[reference_length];
  fread(reference, sizeof(char_t), reference_length, file);
  fclose(file);

  file = fopen(argv[2], "r");
  if (file == 0)
  {
    fprintf(stderr, "ERROR: could not open file `%s'\n", argv[2]);
    delete[] reference;
    return 1;
  } // if
  fseek(file, 0, SEEK_END);
  size_t const sample_length = ftell(file) - 1;
  rewind(file);
  char_t* sample = new char_t[sample_length];
  fread(sample, sizeof(char_t), sample_length, file);
  fclose(file);


  std::vector<Variant> variant;
  size_t const weight = extract(variant, reference, reference_length, sample, sample_length);


  fprintf(stderr, "\nVariants (%ld / %ld):\n", variant.size(), weight);
  for (std::vector<Variant>::iterator it = variant.begin() ; it != variant.end(); ++it)
  {
    //if (it->type != IDENTITY)
    {
      fprintf(stderr, "%ld--%ld, %ld--%ld, %d, %ld, %ld--%ld\n", it->reference_start, it->reference_end, it->sample_start, it->sample_end, it->type, it->weight, it->transposition_start, it->transposition_end);
    } // if
  } // for

  delete[] reference;
  delete[] sample;

  return 0;
} // main

