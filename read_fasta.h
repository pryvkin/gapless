/*  gapless - gapless alignment of sequences
    Copyright (C) 2013 Paul Ryvkin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#ifndef _READ_FASTA_H_
#define _READ_FASTA_H_

#include <iostream>
#include <string>

struct FastaSeq {
  std::string hdr;
  std::string seq;

  FastaSeq() : hdr(), seq() { }
};

std::istream& operator>> (std::istream& is, FastaSeq& fa);

#endif
