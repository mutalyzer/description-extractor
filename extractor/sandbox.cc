// *******************************************************************
//   (C) Copyright 2015 Leiden Institute of Advanced Computer Science
//   Universiteit Leiden
//   All Rights Reserved
// *******************************************************************
// Sandbox (testing playground)
// *******************************************************************
// FILE INFORMATION:
//   File:     sandbox.cc
//   Author:   Jonathan K. Vis
//   Revision: 1.1.5
//   Date:     2015/02/02
// *******************************************************************
// DESCRIPTION:
//  General testing playground: automatic tandem repeat detection
//  based on run length encoding with variable run lengths.
// *******************************************************************

#include <cstdlib>
#include <vector>

#include <cstdio>

typedef char char_t;

static size_t const THRESHOLD = 10000;

struct Repeat
{
  inline Repeat(size_t const start,
                size_t const end,
                size_t const count = 0):
         start(start), end(end), count(count) { }

  inline Repeat(void) { }

  size_t start;
  size_t end;
  size_t count;
}; // Repeat

bool string_match(char_t const* const string_1,
                  char_t const* const string_2,
                  size_t const        length)
{
  for (size_t i = 0; i < length; ++i)
  {
    if (string_1[i] != string_2[i])
    {
      return false;
    } // if
  } // for
  return true;
} // string_match

void tandem_repeat_annotation(std::vector<Repeat> &repeat,
                              char_t const* const  string,
                              size_t const         start,
                              size_t const         end)
{
  size_t const length = end - start;
  size_t const k_max = length > THRESHOLD ? THRESHOLD / 2 : length / 2 + 1;

  size_t i = 0;
  size_t last_repeat = i;
  while (i < length)
  {
    size_t max_count = 0;
    size_t max_k = 1;
    for (size_t k = 1; k < k_max; ++k)
    {
      size_t count = 0;
      for (size_t j = i + k; j < length - k + 1; j += k)
      {
        if (!string_match(string + start + i, string + start + j, k))
        {
          break;
        } // if
        ++count;
      } // for
      if (count > 0 && count >= max_count)
      {
        max_count = count;
        max_k = k;
      } // if
    } // for
    if (max_count > 0)
    {
      if (last_repeat < i)
      {
        repeat.push_back(Repeat(start + last_repeat, start + i));
      } // if
      repeat.push_back(Repeat(start + i, start + i + max_k, max_count));
      last_repeat = i + max_k * (max_count + 1);
    } // if
    i += max_k * (max_count + 1);
  } // while
  if (last_repeat < i)
  {
    repeat.push_back(Repeat(start + last_repeat, start + i));
  } // if
  return;
} // tandem_repeat_annotation

int main(int argc, char* argv[])
{
  if (argc < 2)
  {
    fprintf(stderr, "usage: %s string [start] [end]\n", argv[0]);
    return 1;
  } // if
  fprintf(stderr, "Tandem Repeat Annotator\n");

  FILE* file = fopen(argv[1], "r");
  if (file == 0)
  {
    fprintf(stderr, "ERROR: could not open file `%s'\n", argv[1]);
    return 1;
  } // if
  fseek(file, 0, SEEK_END);
  size_t const length = ftell(file);
  rewind(file);
  char_t* string = new char_t[length];
  size_t const ref_length = fread(string, sizeof(char_t), length, file);
  static_cast<void>(ref_length);
  fclose(file);

  size_t start = 0;
  if (argc > 2)
  {
    start = atoi(argv[2]);
    if (start > length)
    {
      start = 0;
    } // if
  } // if
  size_t end = length;
  if (argc > 3)
  {
    end = atoi(argv[3]);
    if (end > length)
    {
      end = length;
    } // if
  } // if

  std::vector<Repeat> repeat;
  tandem_repeat_annotation(repeat, string, start, end);

  for (std::vector<Repeat>::const_iterator it = repeat.begin(); it != repeat.end(); ++it)
  {
    for (size_t i = it->start; i < it->end; ++i)
    {
      printf("%c", string[i]);
    } // for
    if (it->count > 0)
    {
      printf("%ld", it->count + 1);
    } // if
    printf(";");
  } // for
  printf("\n");

  char_t* reverse = new char_t[length];
  for (size_t i = 0; i < length; ++i)
  {
    reverse[i] = string[length - i - 1];
  } // for

  std::vector<Repeat> repeat_reverse;
  tandem_repeat_annotation(repeat_reverse, reverse, length - end, length - start);

  for (std::vector<Repeat>::const_reverse_iterator it = repeat_reverse.rbegin(); it != repeat_reverse.rend(); ++it)
  {
    size_t const repeat_length = it->end - it->start;
    for (size_t i = 0; i < repeat_length; ++i)
    {
      printf("%c", reverse[it->end - 1 - i]);
    } // for
    if (it->count > 0)
    {
      printf("%ld", it->count + 1);
    } // if
    printf(";");
  } // for
  printf("\n");

  delete[] reverse;
  delete[] string;

  return 0;
} // main

