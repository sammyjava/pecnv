#include <process_readmappings_common.hpp>
#include <sstream>
#include <iostream>
#include <algorithm>

using namespace std;
using namespace Sequence;

string toSAM(const bamrecord & b,
	     const bamreader & reader)
{
  if(b.refid() > reader.ref_cend()-reader.ref_cbegin())
    {
      cerr << "Error: refID " << b.refid()
	   << " does not exist in BAM header. Line "
	   << __LINE__ << " of " << __FILE__ << '\n';
    }
  if(b.next_refid() > reader.ref_cend()-reader.ref_cbegin())
    {
      cerr << "Error: next_refID " << b.next_refid()
	   << " does not exist in BAM header. Line "
	   << __LINE__ << " of " << __FILE__ << '\n';
    }
  ostringstream buffer;
  samflag f(b.flag());
  bool um = f.query_unmapped,mum=f.mate_unmapped;
  buffer << b.read_name() << '\t'
	 << b.flag() << '\t'
	 << (!um ? (reader.ref_cbegin()+b.refid())->first : std::string("*"))<< '\t'
	 << b.pos()+1 << '\t'
	 << b.mapq() << '\t'
	 << b.cigar() << '\t'
	 << (!mum ? (reader.ref_cbegin()+b.next_refid())->first : std::string("*") )<< '\t'
	 << b.next_pos()+1 << '\t'
	 << ((!um && !mum) ? b.tlen() : 0)<< '\t'
	 << b.seq() << '\t';
  std::for_each(b.qual_cbegin(),b.qual_cend(),[&](const char & ch) {
      buffer << char(ch+33);
    });
  buffer << '\t'
	 << b.allaux() << '\n';
  return buffer.str();
}
