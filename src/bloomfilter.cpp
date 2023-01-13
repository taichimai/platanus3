#ifndef BLOOMFILTER_CPP
#define BLOOMFILTER_CPP
#include"common.h"
#include"MyHash.cpp"

template<typename Key>
class BF {
    public:
      BF(uint64_t size, uint8_t numHashes);

      void add(const Key *data, std::size_t len);

      bool possiblyContains(const Key *data, std::size_t len);

      std::array<uint64_t, 2> hash(const Key *data,std::size_t len);

      inline uint64_t  nthHash(uint8_t n,
                        uint64_t hashA,
                        uint64_t hashB,
                        uint64_t filterSize);

    private:
      uint8_t m_numHashes;
      std::vector<bool> m_bits; 
};

template<typename Key>
BF<Key>::BF(uint64_t size, uint8_t numHashes)
      : m_bits(size),
        m_numHashes(numHashes) {}

template<typename Key>
std::array<uint64_t, 2> BF<Key>::hash(const Key *data,std::size_t len){
  std::array<uint64_t, 2> hashValue;
  GetDoubleHash_64bit<Key>(data,hashValue.data());
  return hashValue;
}

template<typename Key>
inline uint64_t BF<Key>::nthHash(uint8_t n,
                        uint64_t hashA,
                        uint64_t hashB,
                        uint64_t filterSize) {
    

    return (hashA + n * hashB) % filterSize;
}

template<typename Key>
void BF<Key>::add(const Key *data, std::size_t len) {
  auto hashValues = hash(data, len);
  for (int n = 0; n < m_numHashes; n++) {
      m_bits[nthHash(n, hashValues[0], hashValues[1], m_bits.size())] = true;
  }
}

template<typename Key>
bool  BF<Key>::possiblyContains(const Key *data, std::size_t len) {
  auto hashValues = hash(data, len);
  for (int n = 0; n < m_numHashes; n++) {
      if (!m_bits[nthHash(n, hashValues[0], hashValues[1], m_bits.size())]) {
          return false;
      }
  }

  return true;
}
#endif 