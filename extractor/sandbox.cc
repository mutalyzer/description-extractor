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
//   Revision: 1.0.1
//   Date:     2015/01/21
// *******************************************************************
// DESCRIPTION:
//  General testing playground: automatic (small) repeat annotation
//  and protein frame shift construction code (in comments).
// *******************************************************************

#include <cstdio>
#include <cstdlib>

typedef char char_t;

struct Node
{
  inline explicit Node(void): start(0), count(0)
  {
    child[0] = child[1] = child[2] = child[3] = 0;
  } // Node
  Node*  child[4];
  size_t start;
  size_t count;
}; // Node

static char_t const IUPAC_BASE[4] =
{
  'A',
  'C',
  'G',
  'T'
}; // IUPAC_BASE

int idx(char_t const base)
{
  switch (base)
  {
    case 'A':
      return 0;
    case 'C':
      return 1;
    case 'G':
      return 2;
    case 'T':
      return 3;
  } // switch
  return -1;
} // idx

char_t const* random_dna(size_t const length)
{
  char_t* const string = new char_t[length];
  for (size_t i = 0; i < length; ++i)
  {
    string[i] = IUPAC_BASE[rand() % 4];
  } // for
  return string;
} // random_dna

void add_string(Node* const root, char_t const* const string, size_t const start, size_t const length)
{
  Node* node = root;
  for (size_t i = 0; i < length; ++i)
  {
    int const c = idx(string[start + i]);
    if (node->child[c] == 0)
    {
      node->child[c] = new Node;
    } // if
    node = node->child[c];
  } // for
  if (node->count == 0 || start >= node->start + length)
  {
    node->start = start;
    ++node->count;
  } // if
  return;
} // add_string

void traverse(Node const* const root, char_t const* const string, size_t const depth = 0)
{
  if (root != 0)
  {
    traverse(root->child[0], string, depth + 1);
    traverse(root->child[1], string, depth + 1);
    traverse(root->child[2], string, depth + 1);
    traverse(root->child[3], string, depth + 1);
    if (root->count > 1)
    {
      for (size_t i = 0; i < depth; ++i)
      {
        printf("%c", string[root->start + i]);
      } // for
      printf(": %ld\n", root->count);
    } // if
  } // if
  return;
} // traverse

void remove(Node const* const root)
{
  if (root != 0)
  {
    remove(root->child[0]);
    remove(root->child[1]);
    remove(root->child[2]);
    remove(root->child[3]);
    delete root;
  } // if
  return;
} // remove

size_t primary(Node const* const root, char_t const* const string, size_t const depth, size_t &start)
{
  if (root == 0)
  {
    return 0;
  } // if

  if (depth == 0)
  {
    start = root->start;
    return root->count;
  } // if

  size_t max_count = primary(root->child[0], string, depth - 1, start);
  size_t here = 0;
  size_t count = primary(root->child[1], string, depth - 1, here);
  if (count > max_count)
  {
    start = here;
    max_count = count;
  } // if
  count = primary(root->child[2], string, depth - 1, here);
  if (count > max_count)
  {
    start = here;
    max_count = count;
  } // if
  count = primary(root->child[3], string, depth - 1, here);
  if (count > max_count)
  {
    start = here;
    max_count = count;
  } // if
  return max_count;
} // primary

int main(int, char*[])
{
  size_t const length = 16;
  //char_t const* const string = random_dna(length);
  //char_t const* const string = "ATAGATAGATAGATAG";
  // length = 220
  //char_t const* const string = "CATGCTGGCCATATTCACTTGCCCACTTCTGCCCAGGGATCTATTTTTCTGTGGTGTGTATTCCCTGTGCCTTTGGGGGCATCTCTTATACTCATGAAATCAACAGAGGCTTGCATGTATCTATCTGTCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATGAGACAGGGTCTTGCTCTGTCACCCAGATTGGACTGCAGT";
  // length = 229
  //char_t const* const string = "ATATGTGAGTCAATTCCCCAAGTGAATTGCCTTCTATCTATCTATCTATCTATCTGTCTGTCTGTCTGTCTGTCTGTCTATCTATCTATATCTATCTATCATCTATCTATCCATATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATATCTATCGTCTATCTATCCAGTCTATCTACCTCCTATTAGTCTGTCTCTGGAGAACATTGACTAATACA";
  // length = 204
  //char_t const* const string = "GGCGACTGAGCAAGACTCAGTCTCAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAATTGTAAGGAGTTTTCTCAATTAATAACCCAAATAAGAGAATTCTTTCCATGTATCAATCATGATACTAAGCACTTTACACACATGTATGTTATGTAATCATTATATCATGCATGCAAGGTAATGAGT";

  for (size_t i = 0; i < length; ++i)
  {
    printf("%c", string[i]);
  } // for
  printf("\n");

  Node* root = new Node;

  for (size_t sub_length = 2; sub_length < length / 2 + 1; ++sub_length)
  {
    for (size_t i = 0; i < length - sub_length + 1; ++i)
    {
      /*
      for (size_t j = 0; j < sub_length; ++j)
      {
        printf("%c", string[i + j]);
      } // for
      printf("\n");
      */
      add_string(root, string, i, sub_length);
    } // for
  } // for

  traverse(root, string);

  for (size_t sub_length = 2; sub_length < length / 2 + 1; ++sub_length)
  {
    size_t start = 0;
    size_t const count = primary(root, string, sub_length, start);
    if (count > 1)
    {
      printf("primary (%ld) = ", sub_length);
      for (size_t i = 0; i < sub_length; ++i)
      {
        printf("%c", string[start + i]);
      } // for
      printf(" (%ld) --- %ld\n", count, count * sub_length);
    } // if
  } // for

  //delete[] string;
  remove(root);
  return 0;
} // main

/*
#include <map>

static unsigned int const FRAME_SHIFT_1       = 0x01;
static unsigned int const FRAME_SHIFT_2       = 0x02;

static char_t const* const codon_string = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";

char_t codon_table[64] = {'\0'};
std::multimap<char_t const, unsigned int const> codon_map;

unsigned int frame_shift(char_t const reference, char_t const sample_1, char_t const sample_2)
{
  std::pair<std::multimap<char_t const, unsigned int const>::const_iterator, std::multimap<char_t const, unsigned int const>::const_iterator> const range_1 = codon_map.equal_range(sample_1);
  std::pair<std::multimap<char_t const, unsigned int const>::const_iterator, std::multimap<char_t const, unsigned int const>::const_iterator> const range_2 = codon_map.equal_range(sample_2);

  unsigned int shift = 0x0;
  for (std::multimap<char_t const, unsigned int const>::const_iterator it_1 = range_1.first; it_1 != range_1.second; ++it_1)
  {
    for (std::multimap<char_t const, unsigned int const>::const_iterator it_2 = range_2.first; it_2 != range_2.second; ++it_2)
    {
      const unsigned int codon_1 = ((it_1->second << 0x2) | (it_2->second >> 0x4)) & 0x3f;
      const unsigned int codon_2 = ((it_1->second << 0x4) | (it_2->second >> 0x2)) & 0x3f;
      if (reference == codon_table[codon_1])
      {
        shift |= FRAME_SHIFT_1;
        if (shift == (FRAME_SHIFT_1 | FRAME_SHIFT_2))
        {
          return shift;
        } // if
      } // if
      if (reference == codon_table[codon_2])
      {
        shift |= FRAME_SHIFT_2;
        if (shift == (FRAME_SHIFT_1 | FRAME_SHIFT_2))
        {
          return shift;
        } // if
      } // if
      // printf("%d %d -> %d %d\n", it_1->second, it_2->second, codon_1, codon_2);
    } // for
  } // for
  return shift;
} // frame_shift


int main(int, char* [])
{

  for (unsigned int i = 0; i < 64; ++i)
  {
    codon_table[i] = codon_string[i];
    codon_map.insert(std::pair<char_t const, unsigned int const>(codon_table[i], i));
  } // for

  //unsigned int const shift = frame_shift('F', 'L', 'L');
  //printf("%d\n", shift);

  std::map<std::pair<char_t const, char_t const> const, size_t> fs_1;
  std::map<std::pair<char_t const, char_t const> const, size_t> fs_2;
  for (std::multimap<char_t const, unsigned int const>::const_iterator it_1 = codon_map.begin(); it_1 != codon_map.end(); ++it_1)
  {
    for (std::multimap<char_t const, unsigned int const>::const_iterator it_2 = codon_map.begin(); it_2 != codon_map.end(); ++it_2)
    {
      for (std::multimap<char_t const, unsigned int const>::const_iterator it_3 = codon_map.begin(); it_3 != codon_map.end(); ++it_3)
      {
        unsigned int const shift = frame_shift(it_3->first, it_1->first, it_2->first);
        if ((shift & FRAME_SHIFT_1) == FRAME_SHIFT_1)
        {
          ++fs_1[std::pair<char_t const, char_t const>(it_1->first, it_2->first)];
        } // if
        if ((shift & FRAME_SHIFT_2) == FRAME_SHIFT_2)
        {
          ++fs_2[std::pair<char_t const, char_t const>(it_1->first, it_2->first)];
        } // if

        char_t const key = it_3->first;
        while (it_3 != codon_map.end() && key == it_3->first)
        {
          ++it_3;
        } // while
        --it_3;
      } // for
      char_t const key = it_2->first;
      while (it_2 != codon_map.end() && key == it_2->first)
      {
        ++it_2;
      } // while
      --it_2;
    } // for
    char_t const key = it_1->first;
    while (it_1 != codon_map.end() && key == it_1->first)
    {
      ++it_1;
    } // while
    --it_1;
  } // for

  for (std::map<std::pair<char_t const, char_t const> const, size_t>::const_iterator it = fs_1.begin(); it != fs_1.end(); ++it)
  {
    printf("%c%c: %ld\n", it->first.first, it->first.second, it->second);
  } // for
  return 1;

  size_t codon_count = 0;
  for (std::multimap<char_t const, unsigned int const>::const_iterator it = codon_map.begin(); it != codon_map.end(); ++it)
  {
    ++codon_count;
    int count = 0;
    char_t const key = it->first;
    while (it != codon_map.end() && key == it->first)
    {
      ++count;
      ++it;
    } // while
    --it;
    //printf("%c: %d\n", key, count);
  } // for
  //printf("%ld\n", codon_count);

  std::map<std::pair<char_t const, std::pair<char_t const, char_t const> const> const, unsigned int const> fs_map_1;
  std::map<std::pair<char_t const, std::pair<char_t const, char_t const> const> const, unsigned int const> fs_map_2;
  for (std::multimap<char_t const, unsigned int const>::const_iterator it_1 = codon_map.begin(); it_1 != codon_map.end(); ++it_1)
  {
    for (std::multimap<char_t const, unsigned int const>::const_iterator it_2 = codon_map.begin(); it_2 != codon_map.end(); ++it_2)
    {
      for (std::multimap<char_t const, unsigned int const>::const_iterator it_3 = codon_map.begin(); it_3 != codon_map.end(); ++it_3)
      {
        unsigned int const shift = frame_shift(it_3->first, it_1->first, it_2->first);
        if ((shift & FRAME_SHIFT_1) == FRAME_SHIFT_1)
        {
          fs_map_1.insert(std::pair<std::pair<char_t const, std::pair<char_t const, char_t const> const> const, unsigned int const>(std::pair<char_t const, std::pair<char_t const, char_t const> const>(it_1->first, std::pair<char_t const, char_t const>(it_2->first, it_3->first)), FRAME_SHIFT_1));
        } // if
        if ((shift & FRAME_SHIFT_2) == FRAME_SHIFT_2)
        {
          fs_map_2.insert(std::pair<std::pair<char_t const, std::pair<char_t const, char_t const> const> const, unsigned int const>(std::pair<char_t const, std::pair<char_t const, char_t const> const>(it_1->first, std::pair<char_t const, char_t const>(it_2->first, it_3->first)), FRAME_SHIFT_2));
        } // if

        char_t const key = it_3->first;
        while (it_3 != codon_map.end() && key == it_3->first)
        {
          ++it_3;
        } // while
        --it_3;
      } // for
      char_t const key = it_2->first;
      while (it_2 != codon_map.end() && key == it_2->first)
      {
        ++it_2;
      } // while
      --it_2;
    } // for
    char_t const key = it_1->first;
    while (it_1 != codon_map.end() && key == it_1->first)
    {
      ++it_1;
    } // while
    --it_1;
  } // for

  for (std::multimap<char_t const, unsigned int const>::const_iterator it_1 = codon_map.begin(); it_1 != codon_map.end(); ++it_1)
  {
    for (std::multimap<char_t const, unsigned int const>::const_iterator it_2 = codon_map.begin(); it_2 != codon_map.end(); ++it_2)
    {
      printf("%c%c 1: ", it_1->first, it_2->first);
      int count = 0;
      for (std::multimap<char_t const, unsigned int const>::const_iterator it_3 = codon_map.begin(); it_3 != codon_map.end(); ++it_3)
      {
        size_t const is = fs_map_1.count(std::pair<char_t const, std::pair<char_t const, char_t const> const>(it_1->first, std::pair<char_t const, char_t const>(it_2->first, it_3->first)));

        if (is > 0)
        {
          printf("%c ", it_3->first);
          ++count;
        } // if

        char_t const key = it_3->first;
        while (it_3 != codon_map.end() && key == it_3->first)
        {
          ++it_3;
        } // while
        --it_3;
      } // for
      printf("(%d)  2: ", count);
      count = 0;
      for (std::multimap<char_t const, unsigned int const>::const_iterator it_3 = codon_map.begin(); it_3 != codon_map.end(); ++it_3)
      {
        size_t const is = fs_map_2.count(std::pair<char_t const, std::pair<char_t const, char_t const> const>(it_1->first, std::pair<char_t const, char_t const>(it_2->first, it_3->first)));

        if (is > 0)
        {
          printf("%c ", it_3->first);
          ++count;
        } // if

        char_t const key = it_3->first;
        while (it_3 != codon_map.end() && key == it_3->first)
        {
          ++it_3;
        } // while
        --it_3;
      } // for

      printf("(%d)\n", count);
      char_t const key = it_2->first;
      while (it_2 != codon_map.end() && key == it_2->first)
      {
        ++it_2;
      } // while
      --it_2;
    } // for
    char_t const key = it_1->first;
    while (it_1 != codon_map.end() && key == it_1->first)
    {
      ++it_1;
    } // while
    --it_1;
  } // for

  return 0;
} // main
*/
