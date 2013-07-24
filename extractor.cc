// *******************************************************************
//   (C) Copyright 2013 Leiden Institute of Advanced Computer Science
//   Universiteit Leiden
//   All Rights Reserved
// *******************************************************************
// Extractor (library)
// *******************************************************************
// FILE INFORMATION:
//   File:     extractor.cc (depends on extractor.h)
//   Author:   Jonathan K. Vis
//   Revision: 1.02a
//   Date:     2013/07/24
// *******************************************************************
// DESCRIPTION:
//   This library can be used to generete HGVS variant descriptions as
//   accepted by the Mutalyzer Name Checker.
// *******************************************************************

#include <cmath>
#include <cstdlib>
#include <cstring>

#include "extractor.h"

namespace mutalyzer
{

std::vector<Variant> extract(char const* const reference,
                             size_t const      reference_length,
                             char const* const sample,
                             size_t const      sample_length,
                             int const         type)
{
  if (type == 1) // Protein
  {
    std::vector<Variant> result;
    extractor(reference, 0, 0, reference_length, sample, 0, sample_length, result);
    return result;
  } // if

  char const* complement = IUPAC_complement(reference, reference_length);
  std::vector<Variant> result;
  extractor(reference, complement, 0, reference_length, sample, 0, sample_length, result);
  delete[] complement;
  return result;
} // extract

void extractor(char const* const     reference,
               char const* const     complement,
               size_t const          reference_start,
               size_t const          reference_end,
               char const* const     sample,
               size_t const          sample_start,
               size_t const          sample_end,
               std::vector<Variant> &result)
{
  if (reference_end - reference_start <= 0)
  {
    if (sample_end - sample_start > 0)
    {
      result.push_back(Variant(reference_start, reference_end, sample_start, sample_end));
    } // if
    return;
  } // if

  if (sample_end - sample_start <= 0)
  {
    result.push_back(Variant(reference_start, reference_end, sample_start, sample_end));
    return;
  } // if

  std::vector<Substring> LCS_result = LCS(reference, complement, reference_start, reference_end, sample, sample_start, sample_end);

  if (LCS_result.size() == 0)
  {
    result.push_back(Variant(reference_start, reference_end, sample_start, sample_end));
    return;
  } // if

  size_t difference = (reference_end - reference_start) + (sample_end - sample_start);
  size_t index = 0;
  for (size_t i = 0; i < LCS_result.size(); ++i)
  {
    if (abs(LCS_result[i].reference_index - LCS_result[i].sample_index) < difference)
    {
      difference = abs(LCS_result[i].reference_index - LCS_result[i].sample_index);
      index = i;
    } // if
  } // for

  extractor(reference, complement, reference_start, LCS_result[index].reference_index, sample, sample_start, LCS_result[index].sample_index, result);
  if (LCS_result[index].reverse_complement)
  {
    result.push_back(Variant(LCS_result[index].reference_index, LCS_result[index].reference_index + LCS_result[index].length, LCS_result[index].sample_index, LCS_result[index].sample_index + LCS_result[index].length, LCS_result[index].length > 1));
  } // if
  extractor(reference, complement, LCS_result[index].reference_index + LCS_result[index].length, reference_end, sample, LCS_result[index].sample_index + LCS_result[index].length, sample_end, result);
  return;
} // extractor

std::vector<Substring> LCS_1(char const* const reference,
                             char const* const complement,
                             size_t const      reference_start,
                             size_t const      reference_end,
                             char const* const sample,
                             size_t const      sample_start,
                             size_t const      sample_end)
{
  size_t const reference_length = reference_end - reference_start;
  size_t const sample_length = sample_end - sample_start;

  size_t** LCS_line = new size_t*[2];
  size_t** LCS_line_rc = new size_t*[2];
  LCS_line[0] = new size_t[reference_length];
  LCS_line[1] = new size_t[reference_length];
  LCS_line_rc[0] = new size_t[reference_length];
  LCS_line_rc[1] = new size_t[reference_length];

  std::vector<Substring> result;

  size_t length = 0;
  for (size_t i = 0; i < sample_length; ++i)
  {
    for (size_t j = 0; j < reference_length; ++j)
    {
      if (complement == 0)
      {
        if (reference[reference_start + j] == sample[sample_start + i])
        {
          if (i == 0 || j == 0)
          {
            LCS_line[i % 2][j] = 1;
          } // if
          else
          {
            LCS_line[i % 2][j] = LCS_line[(i + 1) % 2][j - 1] + 1;
          } // else
          if (LCS_line[i % 2][j] > length)
          {
            length = LCS_line[i % 2][j];
            result = std::vector<Substring>(1, Substring(j - length + reference_start + 1, i - length + sample_start + 1, length));
          } // if
          else if (LCS_line[i % 2][j] == length)
          {
            result.push_back(Substring(j - length + reference_start + 1, i - length + sample_start + 1, length));
          } // if
        } // if
        else
        {
          LCS_line[i % 2][j] = 0;
        } // else
      } // if
      else
      {
        if (IUPAC_match(reference + reference_start + j, sample + sample_start + i))
        {
          if (i == 0 || j == 0)
          {
            LCS_line[i % 2][j] = 1;
          } // if
          else
          {
            LCS_line[i % 2][j] = LCS_line[(i + 1) % 2][j - 1] + 1;
          } // else
          if (LCS_line[i % 2][j] > length)
          {
            length = LCS_line[i % 2][j];
            result = std::vector<Substring>(1, Substring(j - length + reference_start + 1, i - length + sample_start + 1, length));
          } // if
          else if (LCS_line[i % 2][j] == length)
          {
            result.push_back(Substring(j - length + reference_start + 1, i - length + sample_start + 1, length));
          } // if
        } // if
        else
        {
          LCS_line[i % 2][j] = 0;
        } // else
        if (IUPAC_match_reverse(complement + reference_end - j - 1, sample + sample_start + i))
        {
          if (i == 0 || j == 0)
          {
            LCS_line_rc[i % 2][j] = 1;
          } // if
          else
          {
            LCS_line_rc[i % 2][j] = LCS_line_rc[(i + 1) % 2][j - 1] + 1;
          } // else
          if (LCS_line_rc[i % 2][j] > length)
          {
            length = LCS_line_rc[i % 2][j];
            result = std::vector<Substring>(1, Substring(reference_end - j - 1, i - length + sample_start + 1, length, true));
          } // if
          else if (LCS_line_rc[i % 2][j] == length)
          {
            result.push_back(Substring(reference_end - j - 1, i - length + sample_start + 1, length, true));
          } // if
        } // if
        else
        {
          LCS_line_rc[i % 2][j] = 0;
        } // else
      } // else
    } // for
  } // for

  delete[] LCS_line[0];
  delete[] LCS_line[1];
  delete[] LCS_line_rc[0];
  delete[] LCS_line_rc[1];
  delete[] LCS_line;
  delete[] LCS_line_rc;
  return result;
} // LCS_1

std::vector<Substring> LCS_k(char const* const reference,
                             char const* const complement,
                             size_t const      reference_start,
                             size_t const      reference_end,
                             char const* const sample,
                             size_t const      sample_start,
                             size_t const      sample_end,
                             size_t const      k)
{
  size_t const reference_length = (reference_end - reference_start) / k;
  size_t const sample_length = sample_end - sample_start - k + 1;

  size_t** LCS_line = new size_t*[k + 1];
  size_t** LCS_line_rc = new size_t*[k + 1];
  for (size_t i = 0; i < k + 1; ++i)
  {
    LCS_line[i] = new size_t[reference_length];
    LCS_line_rc[i] = new size_t[reference_length];
  } // for

  std::vector<Substring> result;

  size_t length = 0;
  for (size_t i = 0; i < sample_length; ++i)
  {
    for (size_t j = 0; j < reference_length; ++j)
    {
      if (complement == 0)
      {
        if (strncmp(reference + reference_start + j * k, sample + sample_start + i, k) == 0)
        {
          if (i < k || j == 0)
          {
            LCS_line[i % (k + 1)][j] = 1;
          } // if
          else
          {
            LCS_line[i % (k + 1)][j] = LCS_line[(i + 1) % (k + 1)][j - 1] + 1;
          } // else
          if (LCS_line[i % (k + 1)][j] > length)
          {
            length = LCS_line[i % (k + 1)][j];
            for (size_t e = 0; e < result.size(); ++e)
            {
              if (length - result[e].length > 1 || (result[e].reference_index == j - 1 && result[e].sample_index == i - k))
              {
                result.erase(result.begin() + e);
                --e;
              } // if
            } // for
            result.push_back(Substring(j, i, length));
          } // if
          else if (length - LCS_line[i % (k + 1)][j] < 1)
          {
            result.push_back(Substring(j, i, length));
          } // if
        } // if
        else
        {
          LCS_line[i % (k + 1)][j] = 0;
        } // else
      } // if
      else
      {
        if (IUPAC_match(reference + reference_start + j * k, sample + sample_start + i, k))
        {
          if (i < k || j == 0)
          {
            LCS_line[i % (k + 1)][j] = 1;
          } // if
          else
          {
            LCS_line[i % (k + 1)][j] = LCS_line[(i + 1) % (k + 1)][j - 1] + 1;
          } // else
          if (LCS_line[i % (k + 1)][j] > length)
          {
            length = LCS_line[i % (k + 1)][j];
            for (size_t e = 0; e < result.size(); ++e)
            {
              if (length - result[e].length > 1 ||
                  (result[e].reference_index == j - 1 && result[e].sample_index == i - k))
              {
                result.erase(result.begin() + e);
                --e;
              } // if
            } // for
            result.push_back(Substring(j, i, length));
          } // if
          else if (length - LCS_line[i % (k + 1)][j] < 1)
          {
            result.push_back(Substring(j, i, length));
          } // if
        } // if
        else
        {
          LCS_line[i % (k + 1)][j] = 0;
        } // else
        if (IUPAC_match_reverse(complement + reference_end - j * k - 1, sample + sample_start + i, k))
        {
          if (i < k || j == 0)
          {
            LCS_line_rc[i % (k + 1)][j] = 1;
          } // if
          else
          {
            LCS_line_rc[i % (k + 1)][j] = LCS_line_rc[(i + 1) % (k + 1)][j - 1] + 1;
          } // else
          if (LCS_line_rc[i % (k + 1)][j] > length)
          {
            length = LCS_line_rc[i % (k + 1)][j];
            for (size_t e = 0; e < result.size(); ++e)
            {
              if (length - result[e].length > 1 ||
                  (result[e].reference_index == j - 1 && result[e].sample_index == i - k))
              {
                result.erase(result.begin() + e);
                --e;
              } // if
            } // for
            result.push_back(Substring(j, i, length, true));
          } // if
          else if (length - LCS_line_rc[i % (k + 1)][j] < 1)
          {
            result.push_back(Substring(j, i, length, true));
          } // if
        } // if
        else
        {
          LCS_line_rc[i % (k + 1)][j] = 0;
        } // else
      } // else
    } // for
  } // for
  length *= k;
  for (size_t i = 0; i < result.size(); ++i)
  {
    if (complement == 0)
    {
      result[i].reference_index = ((result[i].reference_index + 1) * k + reference_start - 1) - result[i].length * k + 1;
      result[i].sample_index = result[i].sample_index - (result[i].length - 1) * k + sample_start;
      result[i].length *= k;
      size_t j;
      for (j = 1; j < k; ++j)
      {
        if (result[i].reference_index - j < reference_start || result[i].sample_index - j < sample_start || reference[result[i].reference_index - j] != sample[result[i].sample_index - j])
        {
          break;
        } // if
      } // for
      result[i].reference_index -= j - 1;
      result[i].sample_index -= j - 1;
      result[i].length += j - 1;
      for (j = 0; j < k - 1; ++j)
      {
        if (result[i].reference_index + result[i].length + j >= reference_end || result[i].sample_index + result[i].length + j >= sample_end || reference[result[i].reference_index + result[i].length + j] != sample[result[i].sample_index + result[i].length + j])
        {
          break;
        } // if
      } // for
      result[i].length += j;
      if (result[i].length > length)
      {
        length = result[i].length;
      } // if
    } // if
    else if (result[i].reverse_complement)
    {
      result[i].reference_index = reference_end - (result[i].reference_index + 1) * k;
      result[i].sample_index = result[i].sample_index - (result[i].length - 1) * k + sample_start;
      result[i].length *= k;
      size_t j;
      for (j = 1; j < k; ++j)
      {
        if (result[i].reference_index + result[i].length + j - 1 >= reference_end || result[i].sample_index - j < sample_start || !IUPAC_match_reverse(complement + result[i].reference_index + result[i].length + j - 1, sample + result[i].sample_index - j))
        {
          break;
        } // if
      } // for
      result[i].sample_index -= j - 1;
      result[i].length += j - 1;
      for (j = 1; j < k; ++j)
      {
        if (result[i].reference_index - j < reference_start || result[i].sample_index + result[i].length + j - 1 >= sample_end || !IUPAC_match_reverse(complement + result[i].reference_index - j, sample + result[i].sample_index + result[i].length + j - 1))
        {
          break;
        } // if
      } // for
      result[i].reference_index -= j - 1;
      result[i].length += j - 1;
      if (result[i].length > length)
      {
        length = result[i].length;
      } // if
    } // if
    else
    {
      result[i].reference_index = ((result[i].reference_index + 1) * k + reference_start - 1) - result[i].length * k + 1;
      result[i].sample_index = result[i].sample_index - (result[i].length - 1) * k + sample_start;
      result[i].length *= k;
      size_t j;
      for (j = 1; j < k; ++j)
      {
        if (result[i].reference_index - j < reference_start || result[i].sample_index - j < sample_start || !IUPAC_match(reference + result[i].reference_index - j, sample + result[i].sample_index - j))
        {
          break;
        } // if
      } // for
      result[i].reference_index -= j - 1;
      result[i].sample_index -= j - 1;
      result[i].length += j - 1;
      for (j = 0; j < k - 1; ++j)
      {
        if (result[i].reference_index + result[i].length + j >= reference_end || result[i].sample_index + result[i].length + j >= sample_end || !IUPAC_match(reference + result[i].reference_index + result[i].length + j, sample + result[i].sample_index + result[i].length + j))
        {
          break;
        } // if
      } // for
      result[i].length += j;
      if (result[i].length > length)
      {
        length = result[i].length;
      } // if
    } // else
  } // for
  for (size_t i = 0; i < result.size(); ++i)
  {
    if (result[i].length < length)
    {
      result.erase(result.begin() + i);
      --i;
    } // if
  } // for

  for (size_t i = 0; i < k + 1; ++i)
  {
    delete[] LCS_line[i];
    delete[] LCS_line_rc[i];
  } // for
  delete[] LCS_line;
  delete[] LCS_line_rc;

  return result;
} // LCS_k

std::vector<Substring> LCS(char const* const reference,
                           char const* const complement,
                           size_t const      reference_start,
                           size_t const      reference_end,
                           char const* const sample,
                           size_t const      sample_start,
                           size_t const      sample_end)
{
  size_t k = (reference_end - reference_start) / 3;

  std::vector<Substring> result;

// FIXME
//  while (k > log(static_cast<double>(reference_end - reference_start)) / log(static_cast<double>(ALPHABET_SIZE[complement != 0 ? 0 : 1])))
  while (k > 1)
  {
    result = LCS_k(reference, complement, reference_start, reference_end, sample, sample_start, sample_end, k);
    if (result.size() > 0 && result[0].length >= 2 * k)
    {
      return result;
    } // if
    k /= 2;
  } // while
// FIXME
//  return std::vector<Substring>();
  return LCS_1(reference, complement, reference_start, reference_end, sample, sample_start, sample_end);
} // LCS

} // namespace

