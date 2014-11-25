#ifndef __PECNV_BED_HPP__
#define __PECNV_BED_HPP__

//Very simple structures for handling BED-format records
#include <string>
#include <cstdint>

struct BED6
{
  std::string chrom,name;
  std::int32_t start,stop,score;
  char strand;
};


#endif
