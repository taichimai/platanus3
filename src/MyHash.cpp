#ifndef MYHASH_CPP
#define MYHASH_CPP
#include <cstdint>
#include <functional>

#define BIG_CONSTANT(x) (x##LLU)

inline uint64_t rotl64 ( uint64_t x, int8_t r ){
  return (x << r) | (x >> (64 - r));
}

inline uint64_t fmix64 ( uint64_t k ){
  k ^= k >> 33;
  k *= BIG_CONSTANT(0xff51afd7ed558ccd);
  k ^= k >> 33;
  k *= BIG_CONSTANT(0xc4ceb9fe1a85ec53);
  k ^= k >> 33;
  return k;
}

template<typename Key>
void GetDoubleHash_64bit(const Key *hashed_item,void* out){
    const uint64_t c1 = BIG_CONSTANT(0x87c37b91114253d5);
    const uint64_t c2 = BIG_CONSTANT(0x4cf5ad432745937f);
    uint64_t h1=std::hash<Key>()(*hashed_item); 
    uint64_t h2=h1^c2; 
    h2=rotl64(h2,31); h1^=c1;  h1=rotl64(h1,33);
    h1+=h2; h2+=h1; h1^=c2;h2^=c1;
    h1=fmix64(h1); h2=fmix64(h2);
    h1+=h2; h2+=h1; 

    ((uint64_t*)out)[0] = h1;
    ((uint64_t*)out)[1] = h2;
    return;
}
#endif