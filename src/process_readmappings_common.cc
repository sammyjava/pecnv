#include <process_readmappings_common.hpp>
#include <common.hpp>
#include <intermediateIO.hpp>
#include <sstream>
#include <iostream>
#include <algorithm>

using namespace std;
using namespace Sequence;

output_files::output_files(const char * structural_base, const char * um_base) :
  structural_fn(structural_base),
  structural_sam_fn(structural_base),
  um_u_fn(um_base),
  um_m_fn(um_base),
  um_sam_fn(um_base),
  structural(nullptr),
  um_u(nullptr),
  um_m(nullptr),
  structural_sam(nullptr),
  um_sam(nullptr)
{
  structural_fn += ".csv.gz";
  structural_sam_fn += ".sam.gz";
  um_u_fn += "_u.csv.gz";
  um_m_fn += "_m.csv.gz";
  um_sam_fn += ".sam.gz";
  
  structural = gzopen(structural_fn.c_str(),"w");
  if ( structural == NULL ) {
    cerr << "Error, could not open " << structural_fn
	 << " for writing\n";
    exit(1);
  }
  
  structural_sam = gzopen(structural_sam_fn.c_str(),"w");
  if ( structural_sam == NULL ) {
    cerr << "Error, could not open " << structural_sam_fn
	 << " for writing\n";
    exit(1);
    }
  
  um_u = gzopen(um_u_fn.c_str(),"w");
  if ( um_u == NULL ) {
    cerr << "Error, could not open " << um_u_fn
	   << " for writing\n";
    exit(1);
  }
  
  um_m = gzopen(um_m_fn.c_str(),"w");
  if ( um_m == NULL ) {
    cerr << "Error, could not open " << um_m_fn
	 << " for writing\n";
    exit(1);
  }
  
  um_sam = gzopen(um_sam_fn.c_str(),"w");
  if ( um_sam == NULL ) {
    cerr << "Error, could not open " << um_sam_fn
	 << " for writing\n";
    exit(1);
  }
}

output_files::~output_files()
{
  gzclose(structural);
  gzclose(structural_sam);
  gzclose(um_u);
  gzclose(um_m);
  gzclose(um_sam);
}

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

void updateBucket( readbucket & rb, bamrecord & b, 
		   gzFile csvfile, gzFile samfile,
		   const char * maptype,
		   const bamreader & reader )
{
  string n = editRname(b.read_name());
  auto i = rb.find(n);
  if(i == rb.end())
    {
      rb.insert( make_pair(n, std::move(b)) );
    }
  else //We've got our pair, so write it out
    {
      //ostringstream o;
      if(i->second.refid() > reader.ref_cend()-reader.ref_cbegin())
	{
	  cerr << "Error: reference ID number : "<< i->second.refid()
	       << " is not present in the BAM file header. "
	       << " Line " << __LINE__ << " of " << __FILE__ << '\n';
	  exit(1);
	}
      if(b.refid() > reader.ref_cend()-reader.ref_cbegin())
	{
	  cerr << "Error: reference ID number : "<< b.refid()
	       << " is not present in the BAM file header. "
	       << " Line " << __LINE__ << " of " << __FILE__ << '\n';
	  exit(1);
	}
      auto REF = reader.ref_cbegin()+i->second.refid();
      if ( gzwriteCstr( csvfile,editRname(i->second.read_name() ) ) <= 0 )
	{
	  cerr << "Error: gzwrite error encountered at line " << __LINE__ 
	       << " of " << __FILE__ << '\n';
	  exit(1);
	}
      if ( gzwriteCstr( csvfile,REF->first ) <= 0 )
	{
	  cerr << "Error: gzwrite error encountered at line " << __LINE__ 
	       << " of " << __FILE__ << '\n';
	  exit(1);
	}
      REF = reader.ref_cbegin()+b.refid();
      if ( gzwriteCstr( csvfile,REF->first ) <= 0 )
	{
	  cerr << "Error: gzwrite error encountered at line " << __LINE__ 
	       << " of " << __FILE__ << '\n';
	  exit(1);
	}
      if( gzputs( csvfile, maptype ) <= 0 )
	{
	  cerr << "Error: gzputs error encountered at line " << __LINE__ 
	       << " of " << __FILE__ << '\n';
	  exit(1);
	}
      alnInfo a1(i->second),a2(b);
      if( a1.write(csvfile) <= 0 )
	{
	  cerr << "Error: gzwrite error encountered at line " << __LINE__ 
	       << " of " << __FILE__ << '\n';
	  exit(1);
	}
      if( a2.write(csvfile) <= 0 )
	{
	  cerr << "Error: gzwrite error encountered at line " << __LINE__ 
	       << " of " << __FILE__ << '\n';
	  exit(1);
	}
      string SAM = toSAM(b,reader);
      if(!gzwrite(samfile,SAM.c_str(),SAM.size()))
	{
	  cerr << "Error: gzwrite error at line "
	       << __LINE__ << " of " << __FILE__ << '\n';
	  exit(1);
	}
      SAM = toSAM(i->second,reader);
      if(!gzwrite(samfile,SAM.c_str(),SAM.size()))
	{
	  cerr << "Error: gzwrite error at line "
	       << __LINE__ << " of " << __FILE__ << '\n';
	  exit(1);
	}
      rb.erase(i);
    }
}
