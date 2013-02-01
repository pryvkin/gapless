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


#include <iostream>
#include <string>

#include "read_fasta.h"

// removes newlines at the end of a string
void chomp(std::string &s) {
  int end = s.size()-1;
  int startchomp = end + 1;

  for(int i=end; i >=0; --i) {
    if (s[i] == '\r' || s[i] == '\n')
      startchomp = i;
    else
      break;
  }

  if (startchomp <= end)
    s.erase(startchomp);
}

std::istream& operator>> (std::istream& is, FastaSeq& fa) {
  std::string line;
  line.clear();

  fa.hdr.clear();
  fa.seq.clear();

  while(!is.eof()) {
    std::getline(is, line);

    if (line[0] == '>') {
      chomp(line);
      fa.hdr = line.substr(1, line.size()-1);
      break;
    }
  }

  if (is.eof() && fa.hdr.empty())
    return(is);

  int prevpos;

  while(!is.eof()) {
    prevpos = is.tellg();

    std::getline(is, line);
    chomp(line);

    if (line.size() > 0) {
      if (line[0] == '>') {
	// found another header. return stream to prev line
	// and stop
	is.seekg(prevpos, std::ios_base::beg);
	is.clear();
	return(is);

      } else {
	// a seq line
	fa.seq += line;
	
      }
    }
  }

  return(is);
}

// #include <fstream>
// int main(int argc, char **argv) {
//   FastaSeq fa;

//   std::ifstream is(argv[1]);

//   while(!is.eof()) {
//     is >> fa;

//     std::cout << "hdr: " << fa.hdr << "\n";
//     std::cout << "seq: " << fa.seq << "\n";
//   }
// }
