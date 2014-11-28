#ifndef __PECNV_BED_HPP__
#define __PECNV_BED_HPP__

//Very simple structures for handling BED-format records
#include <string>
#include <cstdint>
#include <iosfwd>

struct BED6
{
  std::string chrom,name;
  std::int32_t start,stop,score;
  char strand;
  //In order, these are what you exped in a BED6
  BED6( const std::string &,
	std::int32_t,
	std::int32_t,
	const std::string &, 
	std::int32_t,
	char );
  std::ostream & print(std::ostream & o) const;
};

std::ostream & operator<< (std::ostream & s, const BED6 & b6);

#endif
