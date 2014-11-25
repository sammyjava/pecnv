#include <bed.hpp>
#include<iostream>
using namespace std;

BED6::BED6( const std::string __chrom,
	    const std::int32_t __start,
	    const std::int32_t __stop,
	    const std::string __name, 
	    const std::int32_t __score,
	    const char __strand ) : chrom( move(__chrom) ),
  name( move(__name) ),
  start( move(__start) ),
  stop( move(__stop) ),
  score( move(__score) ),
  strand( move(__strand) )
{
}

ostream & BED6::print( std::ostream & s ) const
{
  s << chrom << '\t'
    << start << '\t'
    << stop << '\t'
    << name << '\t'
    << score << '\t'
    << strand;
  return s;
}

ostream & operator<< (ostream & s, const BED6 & b6)
{
  return b6.print(s);
}
