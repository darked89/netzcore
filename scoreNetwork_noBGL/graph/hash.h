#ifndef HASH_H_
#define HASH_H_

#include "globals.h"
#include "hash_map.h"

struct eqstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return strcmp(s1, s2) == 0;
  }
};

struct isEqualString
{
  bool operator()(string s1, string s2) const
  {
    return s1.compare(s2) == 0;
  }
};

/*
class stringhasher : public stdext::hash_compare <std::string>
{
public:
    // Inspired by the java.lang.String.hashCode() algorithm somewhat processor cache-friendly), @param The string to be hashed, @return The hash value of s
    size_t operator() (const std::string& s) const {
        size_t h = 0;
        std::string::const_iterator p, p_end;
        for(p = s.begin(), p_end = s.end(); p != p_end; ++p) {
            h = 31 * h + (*p);
        }
        return h;
    }
    // @param s1 The first string, @param s2 The second string @return true if the first string comes before the second in lexicographical order
    bool operator() (const std::string& s1, const std::string& s2) const {
        return s1 < s2;
    }
};

typedef stdext::hash_map<std::string, std::string, stringhasher> HASH_S_S;

    hm[std::string("novembro")] = std::string("November");
    it = hm.find(std::string("março"));
    if(it != hm.end())
        std::cout << "The value corresponding to the key 'março' is " << it->second << std::endl;
    else
        std::cout << "The value corresponding to the key 'março' was not found" << std::endl;
*/

#endif /*HASH_H_*/

