#ifndef __PECNV_PROECSS_READMAPPINGS_COMMON_HPP__
#define __PECNV_PROECSS_READMAPPINGS_COMMON_HPP__

#include <string>
#include <Sequence/bamrecord.hpp>
#include <Sequence/bamreader.hpp>

std::string toSAM(const Sequence::bamrecord & b,
		  const Sequence::bamreader & reader);

#endif
