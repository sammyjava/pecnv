#ifndef __PECNV_PROECSS_READMAPPINGS_COMMON_HPP__
#define __PECNV_PROECSS_READMAPPINGS_COMMON_HPP__

#include <Sequence/bamrecord.hpp>
#include <Sequence/bamreader.hpp>
#include <unordered_map>
#include <string>
#include <zlib.h>

struct output_files
{
  //enum MAPTYPE {DIV,PAR,UL,UMU,UMM};
  std::string structural_fn,structural_sam_fn,
    um_u_fn,um_m_fn,um_sam_fn;

  gzFile structural,um_u,um_m,
    structural_sam,um_sam;

  output_files(const char * structural_base, const char * um_base);
  //closes all the gzFiles
  ~output_files();
};

using readbucket = std::unordered_map<std::string, Sequence::bamrecord>; //name, alignment

std::string toSAM(const Sequence::bamrecord & b,
		  const Sequence::bamreader & reader);

void updateBucket( readbucket & rb, Sequence::bamrecord & b, 
		   gzFile csvfile, gzFile samfile,
		   const char * maptype,
		   const Sequence::bamreader & reader );
#endif
