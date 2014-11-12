/*  
  cluster_cnv2.cc

  Faster version that only reads input files once.

  Copyright 2010 Kevin Thornton, University of California Irvine

  This code is released under the terms of the GNU Public License

  This code reads in the *_structural*.csv.gz from a single line,
  and does the following:

  For each chromosome, read pairs are clustered into the same
  CNV call if the following criteria are met:

  1. Two reads have the same lane and read id
  OR
  2. Reads from different read pairs map to the same chromosome
  within mdist base pairs of each other on the same strand
*/

#include <cstdlib>
#include <string>
#include <sstream>
#include <map>
#include <vector>
#include <cassert>
#include <algorithm>
#include <functional>
#include <iostream>
#include <fstream>
#include <limits>
#include <Sequence/IOhelp.hpp>
#include <zlib.h>

using namespace std;

struct linkeddata
{
  mutable unsigned a,aS,b,bS; //positions on strands -- start1,stop1,start2,stop2
  mutable string readname;
  short strand1,strand2;
  linkeddata(const unsigned & __a, 
	     const unsigned & __aS, 
	     const unsigned & __b,
	     const unsigned & __bS,
	     const string & __readname,
	     const short & _strand1,
	     const short & _strand2) : a(__a),
				       aS(__aS),
				       b(__b),
				       bS(__bS),
				       readname( __readname ),
				       strand1(_strand1),strand2(_strand2)
  {
  }
};

auto order_clusters = []( const vector<vector<linkeddata>::const_iterator> & a,
			  const vector<vector<linkeddata>::const_iterator> & b ) 
  {
    unsigned min1 = numeric_limits<unsigned>::max(),
      min2 = numeric_limits<unsigned>::max();
    for(unsigned i=0;i<a.size();++i)
      {
	min1=min(min1,a[i]->a);
      }
    for(unsigned i=0;i<b.size();++i)
      {
	min2=min(min2,b[i]->a);
      }
    return min1 < min2;
  };

using cluster_container = vector< vector<vector<linkeddata>::const_iterator> >;
using lvector = vector<linkeddata>;
using putCNVs = map<string,lvector>;

cluster_container cluster_linked( const lvector & raw,
				  const unsigned & mdist );

void reduce_clusters( cluster_container & clusters,
		      const unsigned & mdist );

unsigned mindist(const unsigned & st1,
		 const unsigned & stp1,
		 const unsigned & st2,
		 const unsigned & stp2);

bool pair_should_cluster( lvector::const_iterator & pair,
			  vector<lvector::const_iterator> & cluster,
			  const unsigned & mdist);

bool unique_positions(const lvector & data,
		      const unsigned & start,
		      const unsigned & start2);

void write_clusters( gzFile o,
		     const string & chrom1,
		     const string & chrom2,
		     const cluster_container & clusters,
		     unsigned * eventid );


void read_data_details(putCNVs & raw_div,
		       putCNVs & raw_par,
		       map<string,putCNVs > & raw_ul,
		       gzFile lin,
		       const unsigned & min_mqual,
		       const unsigned & max_mm,
		       const unsigned & max_gap)
{
  string chrom,chrom2,pairname,pairname2;
  unsigned mqual,strand,mm,gap,
    mqual2,strand2,mm2,gap2;

  int start,stop,start2,stop2;
  std::string type,type2;

  do
    {
      //Very lazy input method...
      auto data = Sequence::IOhelp::gzreadline(lin);
      if(!data.second) break;
      istringstream pdata(data.first);
      pdata >> pairname
	    >> mqual >> chrom >> start >> stop >> strand >> mm >> gap >> type
	    >> mqual2 >> chrom2 >> start2 >> stop2 >> strand2 >> mm2 >> gap2 >> type2 >> ws;

      if( mqual >= min_mqual && mqual2 >= min_mqual &&
	  mm <= max_mm && gap <= max_gap &&
	  mm2 <= max_mm && gap2 <= max_gap )
	{
	  if(type == "DIV")
	    {
	      if ( unique_positions(raw_div[chrom],
				    (strand==0) ? start2 : start,
				    (strand==0) ? start : start2) )
		{
		  assert( (strand==0) ? (strand2 == 1) : (strand == 1) );
		  raw_div[chrom].push_back( linkeddata( (strand==0) ? start2 : start,
							(strand==0) ? stop2 : stop,
							(strand==0) ? start : start2,
							(strand==0) ? stop : stop2,
							pairname,1,0 ) );
		}
	    }
	  else if (type == "PAR")
	    {
	      if ( unique_positions(raw_par[chrom],
				    (start<start2) ? start : start2,
				    (start<start2) ? start2 : start) )
		{
		  raw_par[chrom].push_back( linkeddata( (start<start2) ? start : start2,
							(start<start2) ? stop : stop2,
							(start<start2) ? start2 : start,
							(start<start2) ? stop2 : stop,
							pairname,
							(start<start2) ? strand : strand2,
							(start<start2) ? strand2 : strand ));
		}
	    }
	  else if (type == "UL")
	    {
	      if(chrom==chrom2)
		{
		  cerr << chrom << ' ' <<chrom2 << '\n';
		  cerr << data.first << '\n';
		}
	      assert(chrom != chrom2);
	      if( chrom > chrom2 )
		{
		  swap(chrom,chrom2);
		  swap(start,start2);
		  swap(stop,stop2);
		  swap(strand,strand2);
		}
	      if ( unique_positions(raw_ul[chrom][chrom2],start,start2) )
		{
		  raw_ul[chrom][chrom2].push_back( linkeddata(start,stop,start2,stop2,
							      pairname,strand,strand2) );
		}
	    }
#ifndef NDEBUG
	  else
	    {
	      abort();
	    }
#endif
	}
    } while(!gzeof(lin));
}

void read_data( putCNVs & raw_div,
		putCNVs & raw_par,
		map<string,putCNVs > & raw_ul,
		const char * left,
		const unsigned & min_mqual,
		const unsigned & max_mm,
		const unsigned & max_gap)
{
  gzFile input = gzopen(left,"r");
  if(input == NULL) {
    cerr << "Error: could not open "
	 << left
	 << " for reading\n";
    exit(1);
  }
  read_data_details( raw_div,raw_par,raw_ul , input, min_mqual,max_mm,max_gap );
  gzclose(input);
}
		


int main(int argc, char ** argv)
{
  int argn=1;
  if( argc < 8 )
    {
      cerr << "usage: "
	   << argv[0]
	   << " min_quality max_mm max_gap insert_size outfile_div outfile_par outfile_ul "
	   << "structural_file1a structural_file1b ... structural_fileNa structural_fileNb\n";
      exit(0);
    }
  const unsigned min_mqual = atoi(argv[argn++]);
  const unsigned max_mm = atoi(argv[argn++]);
  const unsigned max_gap = atoi(argv[argn++]);
  const unsigned mdist = atoi(argv[argn++]);
  const char * divfile = argv[argn++];
  const char * parfile = argv[argn++];
  const char * ulfile = argv[argn++];

  //make sure output files are writable
  const string header = "id\tchrom1\tcoverage\tstrand1\tstart1\tstop1\tchrom2\tstrand2\tstart2\tstop2\treads";

  gzFile divstream = gzopen(divfile,"wb");
  if(divstream==NULL) {
    cerr << "Error: could not open "
	 << divfile
	 << " for writing\n";
  }
 if( gzprintf(divstream,"%s\n",header.c_str()) <= 0 )
    {
      cerr << "Error: gzprintf error encountered at line " << __LINE__ 
	   << " of " << __FILE__ << '\n';
      exit(1);
    }

  gzFile parstream = gzopen(parfile,"wb");
  if(parstream == NULL) {
    cerr << "Error: could not open "
	 << parfile << " for writing\n";
    exit(1);
  }
  if (gzprintf(parstream,"%s\n",header.c_str()) <= 0 )
    {
      cerr << "Error: gzprintf error encountered at line " << __LINE__ 
	   << " of " << __FILE__ << '\n';
      exit(1);
    }

  gzFile ulstream = gzopen(ulfile,"wb");
  if(ulstream == NULL)
    {
      cerr << "Error: could not open "
	   << ulfile
	   << " for writing\n";
      exit(1);
    }
  if (gzprintf(ulstream,"%s\n",header.c_str()) <= 0 )
    {
      cerr << "Error: gzprintf error encountered at line " << __LINE__ 
	   << " of " << __FILE__ << '\n';
      exit(1);
    }

  map<string, lvector > raw_div;
  map<string, lvector > raw_par;
  map<string, putCNVs > raw_ul;
  for(int i = argn;i<argc;++i)//i+=2)
    {
      cerr << "processing " << argv[i] << '\n';
      read_data(raw_div,raw_par,raw_ul,
		argv[i],min_mqual,max_mm,max_gap);
    }

  unsigned eventid=0;
  cerr << "clustering div\n";
  for(putCNVs::iterator itr = raw_div.begin();
      itr != raw_div.end();++itr)
    {
      sort(itr->second.begin(),
	   itr->second.end(),
	   [](const linkeddata & lhs, const linkeddata & rhs){
	     return lhs.a < rhs.a && lhs.b < rhs.b;
	   });
      cluster_container clusters = cluster_linked(itr->second,mdist);
      cluster_container clusters2(clusters);
      sort(clusters.begin(),clusters.end(),order_clusters);
      write_clusters( divstream, 
		      itr->first,itr->first,
		      clusters,&eventid );
    }

  cerr << "clustering par\n";
  eventid=0;
  for(putCNVs::iterator itr = raw_par.begin();
      itr != raw_par.end();++itr)
    {
      sort(itr->second.begin(),
	   itr->second.end(),
	   [](const linkeddata & lhs, const linkeddata & rhs){
	     return lhs.a < rhs.a && lhs.b < rhs.b;
	   });
      cluster_container clusters = cluster_linked(itr->second,mdist);
      sort(clusters.begin(),clusters.end(),order_clusters);
      write_clusters( parstream, 
		      itr->first,
		      itr->first,
		      clusters,&eventid );
    }

  cerr << "clustering ul\n";
  eventid=0;
  for( map<string, putCNVs >::iterator itr = raw_ul.begin() ;
       itr != raw_ul.end() ; ++itr )
    {
      for( putCNVs::iterator itr2 = itr->second.begin() ; 
	   itr2 != itr->second.end() ; ++itr2 )
	{
	  assert(itr->first < itr2->first);
	  sort(itr2->second.begin(),itr2->second.end(),
	       [](const linkeddata & lhs, const linkeddata & rhs){
		 return lhs.a < rhs.a && lhs.b < rhs.b;
	       });
	  cluster_container clusters = cluster_linked(itr2->second,mdist);
	  sort(clusters.begin(),clusters.end(),order_clusters);
	  write_clusters( ulstream, 
			  itr->first,
			  itr2->first,
			  clusters,&eventid );
	}
    }
}

bool unique_positions(const lvector & data,
		      const unsigned & start,
		      const unsigned & start2)
{
  for(unsigned i=0;i<data.size();++i)
    {
      if( start == data[i].a && 
	  start2 == data[i].b )
	{
	  return false;
	}
    }
  return true;
}

/*
  OLD VERSION, USED IN DGRP
bool pair_should_cluster( lvector::const_iterator & pair,
			  vector<lvector::const_iterator> & cluster,
			  const unsigned & mdist)
{
  for( unsigned i = 0 ; i < cluster.size() ; ++i )
    {
      if( pair->strand1 == cluster[i]->strand1
	  && pair->strand2 == cluster[i]->strand2 )
	{
	  if ( (max(pair->a,cluster[i]->a)-min(pair->a,cluster[i]->a) <= mdist || //a1 vs a2
		max(pair->a,cluster[i]->aS)-min(pair->a,cluster[i]->aS) <= mdist ||// a1 vs aS2
		max(pair->aS,cluster[i]->a)-min(pair->aS,cluster[i]->a) <= mdist ||// aS1 vs a2
		max(pair->aS,cluster[i]->aS)-min(pair->aS,cluster[i]->aS) <= mdist )//as1 vs aS2
	       &&
	       (max(pair->b,cluster[i]->b)-min(pair->b,cluster[i]->b) <= mdist ||//b1 vs b2
		max(pair->b,cluster[i]->bS)-min(pair->b,cluster[i]->bS) <= mdist ||//b1 vs bS2
		max(pair->bS,cluster[i]->b)-min(pair->bS,cluster[i]->b) <= mdist ||//bS vs b2
		max(pair->bS,cluster[i]->bS)-min(pair->bS,cluster[i]->bS) <= mdist ) )//bS vs bS2
	    {
	      return true;
	    }
	}
      else
	{
	  return false;
	}
    }
  return false;
}
*/

unsigned mindist(const unsigned & st1,
		 const unsigned & stp1,
		 const unsigned & st2,
		 const unsigned & stp2)
{
  unsigned a = max(st1,st2) - min(st1,st2);
  unsigned b = max(st1,stp2) - min(st1,stp2);
  unsigned c = max(stp1,stp2) - min(stp1,stp2);
  unsigned d = max(stp1,st2) - min(stp1,st2);

  return( min(a,min(min(b,c),d)) );
}



bool pair_should_cluster( lvector::const_iterator & pair,
			  vector<lvector::const_iterator> & cluster,
			  const unsigned & mdist)
{
  for( unsigned i = 0 ; i < cluster.size() ; ++i )
    {
      if( pair->strand1 == cluster[i]->strand1
	  && pair->strand2 == cluster[i]->strand2 )
	{
	  if( mindist(pair->a,pair->aS,cluster[i]->a,cluster[i]->aS) <= mdist &&
	      mindist(pair->b,pair->bS,cluster[i]->b,cluster[i]->bS) <= mdist )
	    {
	      return true;
	    }
	}
      else if(pair->strand1 == cluster[i]->strand2 &&
	      pair->strand2 == cluster[i]->strand1)
	{
	  if( mindist(pair->a,pair->aS,cluster[i]->b,cluster[i]->bS) <= mdist &&
	      mindist(pair->b,pair->bS,cluster[i]->a,cluster[i]->aS) <= mdist )
	    {
	      return true;
	    }
	}
    }
  return false;
}

cluster_container cluster_linked( const lvector & raw,
				  const unsigned & mdist )
{
  using citr = lvector::const_iterator;
  cluster_container clusters;

  clusters.push_back( vector<citr>(1,raw.begin()) );
  for( citr i = raw.begin()+1; i < raw.end();++i )
    {
      bool clustered=false;
      for(unsigned j=0;j<clusters.size();++j)
	{
	  if( pair_should_cluster(i,clusters[j],mdist) )
	    {
	      clusters[j].push_back(i);
	      clustered=true;
	      j=clusters.size();
	    } 
	}
      if(!clustered)
	{
	  clusters.push_back( vector<citr>(1,i) );
	}
    }
  reduce_clusters(clusters,mdist);
  return clusters;
}

void reduce_clusters( cluster_container & clusters,
		      const unsigned & mdist )
{
  typedef cluster_container::iterator citr;
  citr i=clusters.end()-1,j,beg=clusters.begin();
  while( i>beg )
    {
      bool merged = false;
      for( j = i - 1 ; !merged&&j >= beg ; --j )
	{
	  //do any reads in i cluster with any data in j?
	  for(unsigned k=0 ; !merged && k < i->size() ; ++k )
	    {
	      if( pair_should_cluster(*(i->begin()+k),*j,mdist ) )
		{
		  //answer is yes, so we merge i into j, 
		  //delete i, and take care of iterator invalidation
		  copy( i->begin(),i->end(),std::back_inserter(*j) );
		  clusters.erase(i);
		  i=clusters.end()-1;
		  beg=clusters.begin();
		  merged=true;
		}
	    }
	}
      if(!merged)
	{
	  --i;
	}
    }
}

//void write_clusters( filtering_ostream & o,
void write_clusters( gzFile gzout,
		     const string & chrom1,
		     const string & chrom2,
		     const cluster_container & clusters,
		     unsigned * eventid )
{
  for(unsigned i=0;i<clusters.size();++i)
    {
      //get the boundaries of each event
      unsigned min1=numeric_limits<unsigned>::max(),max1=0,min2=numeric_limits<unsigned>::max(),max2=0;
      string readnames;
      for(unsigned j=0;j<clusters[i].size();++j)
	{
	  min1 = min(min1,clusters[i][j]->a);
	  max1 = max(max1,clusters[i][j]->aS);
	  min2 = min(min2,clusters[i][j]->b);
	  max2 = max(max2,clusters[i][j]->bS);
	  ostringstream t;
	  t << ';' 
	    << clusters[i][j]->a << ',' 
	    << clusters[i][j]->aS << ','
	    << clusters[i][j]->strand1 << ','
	    << clusters[i][j]->b << ',' 
	    << clusters[i][j]->bS << ','
	    << clusters[i][j]->strand2;
	  if ( readnames.empty() )
	    {
	      readnames += clusters[i][j]->readname;
	    }
	  else
	    {
	      readnames += "|";
	      readnames += clusters[i][j]->readname;
	    }
	  readnames += t.str();
	}
      ostringstream o;
      o << *eventid << '\t'
	<< chrom1 << '\t'
	<< clusters[i].size() << '\t'
	<< clusters[i][0]->strand1 << '\t'
	<< min1 << '\t' << max1 << '\t'
	<< chrom2 << '\t'
	<< clusters[i][0]->strand2 << '\t'
	//<< min2 << '\t' << max2 << '\t';
	<< min2 << '\t' << max2 << '\t'
	<< readnames << '\n';
      if(!gzwrite(gzout,o.str().c_str(),o.str().size()))
	{
	  cerr << "Error: gzwrite error encountered at line " << __LINE__ 
	       << " of " << __FILE__ << '\n';
	  exit(1);
	}
      ++(*eventid);
    } 
}
