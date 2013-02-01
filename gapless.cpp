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

#include <fstream>
#include <sstream>
#include <vector>

#include "read_fasta.h"

using namespace std;

int main(int argc, char **argv) {
  if (argc < 3) {
    cerr << "USAGE: " << argv[0] << " matrix file1 file2\n";
    cerr << "or     " << argv[0] << " matrix file\n";
    cerr << "The second usage is for a symmetric all vs all query.\n";
    return(-1);
  }

  bool one_fasta = (argc == 3);

  ifstream matrix_file(argv[1]);
  ifstream fa_file1(argv[2]);
  ifstream fa_file2;

  if (!one_fasta)
    fa_file2.open(argv[3]);

  ifstream &fa_db = one_fasta ? fa_file1 : fa_file2;

  string line;
  getline(matrix_file, line);
  
  istringstream iss(line);

  int matrix[256][256];
  bool valid_char[256];
  int matrix_min = 65535;

  for(int i=0; i<256; ++i) {
    valid_char[i] = false;

    for(int j=0; j<256; ++j)
      matrix[i][j] = 0;
  }

  string s;
  vector<char> chars;

  // read chars
  while(!iss.eof()) {
    iss >> s;
    chars.push_back(s[0]);
    valid_char[(int)s[0]] = true;
  }
  int nchars = chars.size();

  // read matrix values
  int curr_char_i = 0;
  int curr_char_j = 0;

  while(!matrix_file.eof()) {
    getline(matrix_file, line);
    istringstream iss(line);

    for(curr_char_j=0; curr_char_j < nchars; ++curr_char_j) {
      int v;
      iss >> v;
      matrix[(int)chars[curr_char_i]][(int)chars[curr_char_j]] = v;
      if (v < matrix_min)
	matrix_min = v;
    }

    ++curr_char_i;
  }

  // fill the rest of the matrix (for invalid chars)
  // with the minimal matrix value
  for(int i=0; i<256; ++i)
    if (!valid_char[i])
      for(int j=0; j<256; ++j)
	if (!valid_char[j])
	  matrix[i][j] = matrix_min;

  // display matrix
  // cerr << "\t";
  // for(int i=0; i<nchars; ++i)
  //   cerr << chars[i] << "\t";
  // cerr << "\n";

  // for(int i=0; i<nchars; ++i) {
  //   cerr << chars[i] << "\t";
  //   for(int j=0; j<nchars; ++j)
  //     cerr << matrix[(int)chars[i]][(int)chars[j]] << "\t";
  //   cerr << "\n";
  // }

  // create db of target seqs
  vector<FastaSeq> seqs;
  while(!fa_db.eof()) {
    FastaSeq fa;
    fa_db >> fa;
    
    if (!fa.hdr.empty() && !fa.seq.empty())
      seqs.push_back(fa);
  }

  int nseqs = seqs.size();

  cerr << "Read " << nseqs << " sequences.\n"; 

  // process query seqs
  if (one_fasta) {
    fa_file1.close();
    fa_file1.open(argv[2]);
  }


  int seq_i = 0;

  while(!fa_file1.eof()) {
    FastaSeq fa;
    fa_file1 >> fa;
    
    if (fa.hdr.empty() || fa.seq.empty())
      continue;

    int seq_j_start = one_fasta ? seq_i : 0;

    for(int seq_j = seq_j_start; seq_j < nseqs; ++seq_j) {
      //cout << fa.seq << "\n" << seqs[seq_j].seq << "\n\n";

      string &seq1 = fa.seq;
      string &seq2 = seqs[seq_j].seq;

      string seq2_padded(seq1.size(), '-');

      seq2_padded += seq2;
      seq2_padded += string(seq1.size(), '-');

      int max_offset = seq1.size() + seq2.size();

      vector<int> scores(max_offset+1);
      int best_offset = -1;
      int best_score = -65535;

      for(int offset=0; offset <= max_offset; ++offset) {
	for(unsigned int i=0; i < seq1.size(); ++i) {
	  scores[offset] += matrix[(int)seq1[i]][(int)seq2_padded[offset+i]];
	}
	if (scores[offset] > best_score) {
	  best_offset = offset;
	  best_score = scores[offset];
	}
	//cout << "offset: " << offset << "; score: " << scores[offset] << "\n";
	//cout << seq1 << "\n";
	//cout << seq2_padded.substr(offset, seq1.size()) << "\n\n";
      }

      //// tons of ugly code to make the alignments display correctly
      
      int offset2 = best_offset - seq1.length();
      //cout << "o2: " << offset2 << "\n";
      unsigned int aln_len = best_offset+((seq1.size() < seq2.size()) ? seq2.size() : seq1.size());
      string aln1 = ((offset2 < 0) ? string() : string(offset2, '-')) +  seq1;
      if (seq1.size() < aln_len) {
	aln1 += string(aln_len-seq1.size(), '-');
      }
      string aln2 = seq2_padded.substr(best_offset - ((offset2 < 0) ? 0 : offset2), aln_len);

      unsigned int rtrim;
      string &long_seq = (aln1.size() > aln2.size()) ? aln1 : aln2;
      string &short_seq = (aln1.size() > aln2.size()) ? aln2 : aln1;

      for(rtrim=long_seq.size()-1; rtrim >= 0; --rtrim) {
	if (long_seq[rtrim] != '-' || (rtrim < short_seq.size() && short_seq[rtrim] != '-')) {
	  ++rtrim;
	  break;
	}
      }
      //cout << long_seq.size() << "," << short_seq.size() << "," << rtrim << "\n";
      if (rtrim < long_seq.size())
	long_seq.erase(rtrim);
      if (rtrim < short_seq.size())
	short_seq.erase(rtrim);

      cout << fa.hdr << "\t" << seqs[seq_j].hdr << "\t" << best_score << "\t"
	   << aln1 << "\t"
	   << aln2 << "\n";

    }

    ++seq_i;

  }

  return(0);
}
