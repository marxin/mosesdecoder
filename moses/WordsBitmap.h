// $Id$

/***********************************************************************
Moses - factored phrase-based language decoder
Copyright (C) 2006 University of Edinburgh

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
***********************************************************************/

#ifndef moses_WordsBitmap_h
#define moses_WordsBitmap_h

#include <limits>
#include <vector>
#include <iostream>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include "TypeDef.h"
#include "WordsRange.h"
#include <boost/dynamic_bitset.hpp>

#define UINT_SIZE 64

namespace Moses
{
typedef unsigned long WordsBitmapID;

/** vector of boolean used to represent whether a word has been translated or not
*/
class WordsBitmap
{
  friend std::ostream& operator<<(std::ostream& out, const WordsBitmap& wordsBitmap);
protected:
  boost::dynamic_bitset<> m_new_bitmap;

  WordsBitmap(); // not implemented

  //! set all elements to false
  void Initialize() {
  }

  //sets elements by vector
  void Initialize(std::vector<bool> vector) {
    size_t vector_size = vector.size();
    for (size_t pos = 0 ; pos < m_new_bitmap.size(); pos++) {
      if (pos < vector_size && vector[pos] == true) {
	m_new_bitmap[pos] = true;
      }
    }
  }


public:
  //! create WordsBitmap of length size and initialise with vector
  WordsBitmap(size_t size, std::vector<bool> initialize_vector) {
    m_new_bitmap.resize(size);
    Initialize(initialize_vector);
  }
  //! create WordsBitmap of length size and initialise
  WordsBitmap(size_t size) {
    m_new_bitmap.resize(size);
  }
  //! deep copy
  WordsBitmap(const WordsBitmap &copy) {
    m_new_bitmap = copy.m_new_bitmap;
  }
  //! count of words translated
  size_t GetNumWordsCovered() const {
    return m_new_bitmap.count ();
  }

  size_t GetFirstGapPos() const {
    if (m_new_bitmap.size() < UINT_SIZE)
    {
      unsigned long long v = m_new_bitmap.to_ulong();
      unsigned index = (unsigned)__builtin_ctzll(~v);
      return index < m_new_bitmap.size() ? index : 0;
    }

    for (size_t pos = 0 ; pos < m_new_bitmap.size() ; pos++) {
      if (!m_new_bitmap[pos]) {
        return pos;
      }
    }
    // no starting pos
    return NOT_FOUND;
  }


  //! position of last word not yet translated, or NOT_FOUND if everything already translated
  size_t GetLastGapPos() const {
    for (int pos = (int) m_new_bitmap.size() - 1 ; pos >= 0 ; pos--) {
      if (!m_new_bitmap[pos]) {
        return pos;
      }
    }
    // no starting pos
    return NOT_FOUND;
  }


  //! position of last translated word
  size_t GetLastPos() const {
    for (int pos = (int) m_new_bitmap.size() - 1 ; pos >= 0 ; pos--) {
      if (m_new_bitmap[pos]) {
        return pos;
      }
    }
    // no starting pos
    return NOT_FOUND;
  }

  bool IsAdjacent(size_t startPos, size_t endPos) const;

  //! whether a word has been translated at a particular position
  bool GetValue(size_t pos) const {
    return m_new_bitmap[pos];
  }
  //! set value at a particular position
  void SetValue( size_t pos, bool value ) {
    m_new_bitmap[pos] = value;
  }
  //! set value between 2 positions, inclusive
  void SetValue( size_t startPos, size_t endPos, bool value ) {
    for(size_t pos = startPos ; pos <= endPos ; pos++) {
      m_new_bitmap[pos] = value;
    }
  }
  //! whether every word has been translated
  bool IsComplete() const {
    return GetSize() == GetNumWordsCovered();
  }
  //! whether the wordrange overlaps with any translated word in this bitmap
  bool Overlap(const WordsRange &compare) const {
    for (size_t pos = compare.GetStartPos() ; pos <= compare.GetEndPos() ; pos++) {
      if (m_new_bitmap[pos])
        return true;
    }
    return false;
  }
  //! number of elements
  size_t GetSize() const {
    return m_new_bitmap.size();
  }

  inline int Compare (const WordsBitmap &compare) const {
    bool eq = m_new_bitmap == compare.m_new_bitmap;
    if (eq)
      return 0;

    bool lt = m_new_bitmap < compare.m_new_bitmap;
    return lt ? -1 : 1;
  }

  bool operator< (const WordsBitmap &compare) const {
    return m_new_bitmap < compare.m_new_bitmap;
  }

  inline size_t GetEdgeToTheLeftOf(size_t l) const {
    if (l == 0) return l;
    while (l && !m_new_bitmap[l-1]) {
      --l;
    }
    return l;
  }

  inline size_t GetEdgeToTheRightOf(size_t r) const {
    if (r+1 == m_new_bitmap.size()) return r;
    while (r+1 < m_new_bitmap.size() && !m_new_bitmap[r+1]) {
      ++r;
    }
    return r;
  }


  //! TODO - ??? no idea
  int GetFutureCosts(int lastPos) const ;

  //! converts bitmap into an integer ID: it consists of two parts: the first 16 bit are the pattern between the first gap and the last word-1, the second 16 bit are the number of filled positions. enforces a sentence length limit of 65535 and a max distortion of 16
  WordsBitmapID GetID() const {
    assert(m_new_bitmap.size() < (1<<16));

    size_t start = GetFirstGapPos();
    if (start == NOT_FOUND) start = m_new_bitmap.size(); // nothing left

    size_t end = GetLastPos();
    if (end == NOT_FOUND) end = 0; // nothing translated yet

    assert(end < start || end-start <= 16);
    WordsBitmapID id = 0;
    for(size_t pos = end; pos > start; pos--) {
      id = id*2 + (int) GetValue(pos);
    }
    return id + (1<<16) * start;
  }

  //! converts bitmap into an integer ID, with an additional span covered
  WordsBitmapID GetIDPlus( size_t startPos, size_t endPos ) const {
	  assert(m_new_bitmap.size() < (1<<16));

    size_t start = GetFirstGapPos();
    if (start == NOT_FOUND) start = m_new_bitmap.size(); // nothing left

    size_t end = GetLastPos();
    if (end == NOT_FOUND) end = 0; // nothing translated yet

    if (start == startPos) start = endPos+1;
    if (end < endPos) end = endPos;

    assert(end < start || end-start <= 16);
    WordsBitmapID id = 0;
    for(size_t pos = end; pos > start; pos--) {
      id = id*2;
      if (GetValue(pos) || (startPos<=pos && pos<=endPos))
        id++;
    }
    return id + (1<<16) * start;
  }

  TO_STRING();
};

// friend
inline std::ostream& operator<<(std::ostream& out, const WordsBitmap& wordsBitmap)
{
  for (size_t i = 0 ; i < wordsBitmap.m_new_bitmap.size() ; i++) {
    out << (wordsBitmap.GetValue(i) ? 1 : 0);
  }
  return out;
}

}
#endif
