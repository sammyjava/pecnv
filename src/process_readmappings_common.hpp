#ifndef __PECNV_PROECSS_READMAPPINGS_COMMON_HPP__
#define __PECNV_PROECSS_READMAPPINGS_COMMON_HPP__

#include <Sequence/bamrecord.hpp>
#include <Sequence/bamreader.hpp>
#include <unordered_map>
#include <string>
#include <zlib.h>

using readbucket = std::unordered_map<std::string, Sequence::bamrecord>; //name, alignment

std::string toSAM(const Sequence::bamrecord & b,
		  const Sequence::bamreader & reader);

void updateBucket( readbucket & rb, Sequence::bamrecord & b, 
		   gzFile csvfile, gzFile samfile,
		   const char * maptype,
		   const Sequence::bamreader & reader );
#endif
