/*
  pecnv orphan module
  Author: Kevin Thornton
 */

#include <Sequence/bamrecord.hpp>
#include <Sequence/bamreader.hpp>
#include <Sequence/samfunctions.hpp>
#include <boost/program_options.hpp>
#include <file_common.hpp>
#include <common.hpp>
#include <bed.hpp>
#include <string>
#include <sstream>
#include <unordered_map>
#include <map>
#include <zlib.h>

using namespace std;
using namespace Sequence;
using namespace boost::program_options;

int orphan_main(int argc, char ** argv)
{
  string bamfile,outfile,sampleid;

  options_description desc("pecnv orphan: extract orphaned reads from a BAM file.");
  desc.add_options()
    ("help,h","Produce help message")
    ("bamfile,b",value<string>(&bamfile),"BAM file to scan")
    ("outfile,o",value<string>(&outfile),"Output file.  This will be a gzipped bed file, unless stdout is given as the option, in which case the output will be written to the STDOUT stream")
    ("sampleid,s",value<string>(&sampleid)->default_value("sample"),"Identifier for the sample")
    ;

  variables_map vm;
  store(parse_command_line(argc, argv, desc), vm);
  notify(vm);

  if ( argc == 1 ||
       vm.count("help") ||
       !vm.count("bamfile") ||
       !vm.count("outfile") )
    {
      cerr << desc << '\n';
      exit(0);
    }

  if( ! file_exists(bamfile.c_str()) )
    {
      cerr << "Error: input file "
	   << bamfile << " does not exist.\n";
      exit(1);
    }

  bamreader reader(bamfile.c_str());

  if( ! reader )
    {
      cerr << "Error: could not open "
	   << bamfile
	   << " for reading\n";
      exit(1);
    }

  auto pos = reader.tell(); //Store the position after reading the header

  unordered_map<string,bamrecord> records;

  //Pass 1 finds all the candidates
  while(!reader.eof() && !reader.error())
    {
      bamrecord b = reader.next_record();
      if(b.empty()) break; //hit EOF
      samflag flag = b.flag();

      if( !flag.query_unmapped && flag.mate_unmapped )
	{
	  records.insert( make_pair( editRname(b.read_name()), std::move(b) ) );
	}
    }

  reader.seek(pos,SEEK_SET);
  //Pass 2 makes sure that the SAM flags of the mates found in Pass 1 
  //agrees

  //ostringstream out;
  unsigned recordno=0; //for making the name column in the BED
  map<string,vector< BED6 > > b6data; //The map will keep the data sorted lexically by chromosome
  while(!reader.eof() && !reader.error())
    {
      bamrecord b = reader.next_record();
      if(b.empty()) break; //hit EOF
      samflag flag = b.flag();

      if( flag.query_unmapped && !flag.mate_unmapped )
	{
	  auto itr =  records.find(editRname(b.read_name()));
	  if ( itr != records.end() )
	    //Write a BED record for the store
	    {
	      auto REF = reader.ref_cbegin() + itr->second.refid();
	      b6data[REF->first].emplace_back( BED6(REF->first,
						    itr->second.pos(),
						    itr->second.pos() + alignment_length(itr->second),
						    string(sampleid + "_orphan" + to_string(recordno++)),
						    itr->second.mapq(),
						    ( itr->second.flag().qstrand ? '-' : '+' ) ) );
	    }
	}
    }

  //Sort the data w/in chromosome
  for( auto i = b6data.begin() ; i != b6data.end() ; ++i )
    {
      sort(i->second.begin(),i->second.end(),
	   [](const BED6 & lhs, const BED6 & rhs) {
	     return lhs.start < rhs.start;
	   } );
    }

  if( outfile == "stdout" )
    {
      for( auto i = b6data.begin() ; i != b6data.end() ; ++i )
	{
	  copy(i->second.begin(),i->second.end(),
	       ostream_iterator<BED6>(cout,"\n"));
	}
    }
  else
    {
      ostringstream out;
      for( auto i = b6data.begin() ; i != b6data.end() ; ++i )
	{
	  copy(i->second.begin(),i->second.end(),
	       ostream_iterator<BED6>(out,"\n"));
	}
      gzFile gzout = gzopen(outfile.c_str(),"w");
      if(gzout == NULL)
  	{
  	  cerr << "Error: could not open "
  	       << outfile
  	       << " for writing at line " << __LINE__ 
  	       << " of " << __FILE__ << '\n';
  	  exit(1);
  	}
      if( gzwrite( gzout, out.str().c_str(), out.str().size() ) <= 0 )
  	{
  	  cerr << "Error: gzwrite error at line "
  	       << __LINE__ << " of "
  	       << __FILE__ << '\n';
  	  exit(1);
  	}
      gzclose(gzout);
    }
  return 0;
}
