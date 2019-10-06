#pragma once

#include <cassert>
#include <fstream>
#include <vector>

namespace succinctrits {

template <uint8_t>
class rank_support;

class trit_vector {
  public:
    class builder {
      public:
        builder() = default;

        void reserve(uint64_t capa) {
            m_trytes.reserve(capa / TRITS_PER_BYTE + 1);
        }
        void push_back(uint8_t t) {
            assert(t < 3);

            ++m_count;
            ++m_num_trits;

            if (m_count == 1) {
                m_trytes.emplace_back(t);
            } else if (m_count == 2) {
                m_trytes.back() += t * 3;
            } else if (m_count == 3) {
                m_trytes.back() += t * 9;
            } else if (m_count == 4) {
                m_trytes.back() += t * 27;
            } else if (m_count == 5) {
                m_trytes.back() += t * 81;
                m_count = 0;
            }
        }

      private:
        std::vector<uint8_t> m_trytes;
        uint64_t m_num_trits = 0;
        uint8_t m_count = 0;

        friend class trit_vector;
    };

  public:
    trit_vector() = default;

    template <class Iterator>
    trit_vector(Iterator it, uint64_t num_trits) {
        build(it, num_trits);
    }
    explicit trit_vector(builder* b) {
        build(b);
    }

    template <class Iterator>
    void build(Iterator it, uint64_t num_trits) {
        builder b;
        b.reserve(num_trits);
        for (uint64_t i = 0; i < num_trits; ++i) {
            b.push_back(*it);
            ++it;
        }
        build(&b);
    }
    void build(builder* b) {
        m_trytes = std::move(b->m_trytes);
        m_num_trits = b->m_num_trits;
    }

    uint8_t get(uint64_t i) const {
        assert(i < m_num_trits);

        const uint64_t pos = i / TRITS_PER_BYTE;
        const uint64_t mod = i % TRITS_PER_BYTE;

        switch (mod) {
            case 0:
                return m_trytes[pos] % 3;
            case 1:
                return m_trytes[pos] / 3 % 3;
            case 2:
                return m_trytes[pos] / 9 % 3;
            case 3:
                return m_trytes[pos] / 27 % 3;
            case 4:
                return m_trytes[pos] / 81 % 3;
        }

        // should not come
        assert(false);
        return uint8_t(-1);
    }

    uint8_t operator[](uint64_t i) const {
        return get(i);
    }

    uint64_t get_num_trits() const {
        return m_num_trits;
    }
    uint64_t size_in_bytes() const {
        return m_trytes.size() * sizeof(uint8_t) + sizeof(m_num_trits);
    }

    void save(std::ostream& os) const {
        size_t n = m_trytes.size();
        os.write(reinterpret_cast<const char*>(&n), sizeof(size_t));
        os.write(reinterpret_cast<const char*>(m_trytes.data()), sizeof(uint8_t) * n);
        os.write(reinterpret_cast<const char*>(&m_num_trits), sizeof(uint64_t));
    }

    void load(std::istream& is) {
        size_t n = 0;
        is.read(reinterpret_cast<char*>(&n), sizeof(size_t));
        m_trytes.resize(n);
        is.read(reinterpret_cast<char*>(m_trytes.data()), sizeof(uint8_t) * n);
        is.read(reinterpret_cast<char*>(&m_num_trits), sizeof(uint64_t));
    }

  private:
    static constexpr uint64_t TRITS_PER_BYTE = 5;

    std::vector<uint8_t> m_trytes;  // each of 5 trits
    uint64_t m_num_trits = 0;

    friend class rank_support<0>;
    friend class rank_support<1>;
    friend class rank_support<2>;
};

}  // namespace succinctrits