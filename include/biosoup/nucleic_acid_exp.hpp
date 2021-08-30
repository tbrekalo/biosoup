// Copyright (c) 2020 Robert Vaser

#ifndef BIOSOUP_NUCLEIC_ACID_EXP_HPP_
#define BIOSOUP_NUCLEIC_ACID_EXP_HPP_

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace biosoup {
namespace exp {

/* clang-format off */
constexpr static std::uint8_t kNucleotideCoder[] = {
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    0,   255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 0,   1,   1,   0,   255, 255, 2,   3,   255, 255,
    2,   255, 1,   0,   255, 255, 255, 0,   1,   3,   3,   2,   0,   255, 3,
    255, 255, 255, 255, 255, 255, 255, 0,   1,   1,   0,   255, 255, 2,   3,
    255, 255, 2,   255, 1,   0,   255, 255, 255, 0,   1,   3,   3,   2,   0,
    255, 3,   255, 255, 255, 255, 255, 255};

constexpr static char kNucleotideDecoder[] = {'A', 'C', 'G', 'T'};
/* clang-format on */

class NucleicAcid {
 public:
  NucleicAcid() = default;

  NucleicAcid(const std::string& name, const std::string& data)
      : NucleicAcid(name.c_str(), name.size(), data.c_str(), data.size()) {}

  NucleicAcid(const char* name, std::uint32_t name_len, const char* data,
              std::uint32_t data_len)
      : id(num_objects++),
        name(name, name_len),
        deflated_data(),
        inflated_len(data_len),
        is_reverse_complement(0) {
    deflated_data.reserve(std::ceil(data_len / 32.));
    std::uint64_t block = 0;
    for (std::uint32_t i = 0; i < data_len; ++i) {
      std::uint64_t c = kNucleotideCoder[static_cast<std::uint8_t>(data[i])];
      if (c == 255ULL) {
        throw std::invalid_argument(
            "[biosoup::NucleicAcid::NucleicAcid] error: not a nucleotide");
      }
      block |= c << ((i << 1) & 63);
      if (((i + 1) & 31) == 0 || i == data_len - 1) {
        deflated_data.emplace_back(block);
        block = 0;
      }
    }
  }

  NucleicAcid(const std::string& name, const std::string& data,
              const std::string& quality)
      : NucleicAcid(name.c_str(), name.size(), data.c_str(), data.size(),
                    quality.c_str(), quality.size()) {}

  NucleicAcid(const char* name, std::uint32_t name_len, const char* data,
              std::uint32_t data_len, const char* quality,
              std::uint32_t quality_len)
      : NucleicAcid(name, name_len, data, data_len) {
    deflated_quality.reserve(std::ceil(quality_len / 32.));
    qlvl.reserve(std::ceil(quality_len / 128.));

    std::vector<std::uint8_t> levels;
    std::vector<std::uint8_t> frequencies(128, 0);
    for (std::uint32_t i = 0; i < quality_len; i += 128U) {
      std::uint32_t j = std::min(i + 128U, quality_len);
      std::fill(frequencies.begin(), frequencies.end(),
                0);  // TODO: consider memset

      double mean = 0;
      std::uint8_t min = -1;
      std::uint8_t max = 0;

      std::uint8_t mode = 0;
      std::uint32_t mode_freq = 0;
      for (std::uint32_t k = i; k < j; ++k) {
        std::uint8_t v = quality[k] - '!';
        if (v < min) {
          min = v;
        }

        if (v > max) {
          max = v;
        }

        mean += v;
        if (++frequencies[v] >= mode_freq) {
          mode = v;
        }
      }

      mean /= j - i;

      if (mean < mode) {
        double step = static_cast<double>(mode - min) / 3.0;
        levels.emplace_back(std::round(mode - 2 * step));
        levels.emplace_back(std::round(mode - step));
        levels.emplace_back(mode);
        step = static_cast<double>(max - mode) / 2.0;
        levels.emplace_back(std::round(mode + step));
      } else {
        double step = static_cast<double>(mode - min) / 2;
        levels.emplace_back(std::round(mode - step));
        levels.emplace_back(mode);
        step = static_cast<double>(max - mode) / 3.0;
        levels.emplace_back(std::round(mode + step));
        levels.emplace_back(std::round(mode + 2 * step));
      }

      std::uint64_t block = 0;
      for (std::uint32_t k = i; k < j; ++k) {
        std::uint64_t m = 0;
        for (unsigned l = 1; l < levels.size(); ++l) {
          if (std::abs((quality[k] - '!') - levels[l]) <
              std::abs((quality[k] - '!') - levels[m])) {
            m = l;
          }
        }
        block |= (3ULL - m) << ((k << 1) & 63);
        if (((k + 1) & 31) == 0 || k == j - 1) {
          deflated_quality.emplace_back(block);
          block = 0;
        }
      }

      std::uint32_t level = 0;
      for (const auto& it : levels) {
        level <<= 8;
        level |= it;
      }
      qlvl.emplace_back(level);
      levels.clear();
    }
  }

  NucleicAcid(const NucleicAcid&) = default;
  NucleicAcid& operator=(const NucleicAcid&) = default;

  NucleicAcid(NucleicAcid&&) = default;
  NucleicAcid& operator=(NucleicAcid&&) = default;

  ~NucleicAcid() = default;

  std::uint64_t Code(std::uint32_t i) const {
    std::uint64_t x = 0;
    if (is_reverse_complement) {
      i = inflated_len - i - 1;
      x = 3;
    }
    return ((deflated_data[i >> 5] >> ((i << 1) & 63)) & 3) ^ x;
  }

  std::uint8_t Score(std::uint32_t i) const {
    if (is_reverse_complement) {
      i = inflated_len - i - 1;
    }
    return (qlvl[i >> 7] >>
            ((deflated_quality[i >> 5] >> ((i << 1) & 63)) << 3)) &
           127;
  }

  std::string InflateData(std::uint32_t i = 0, std::uint32_t len = -1) const {
    if (i >= inflated_len) {
      return std::string{};
    }
    len = std::min(len, inflated_len - i);

    std::string dst{};
    dst.reserve(len);
    for (; len; ++i, --len) {
      dst += kNucleotideDecoder[Code(i)];
    }
    return dst;
  }

  std::string InflateQuality(std::uint32_t i = 0,
                             std::uint32_t len = -1) const {  // NOLINT
    if (deflated_quality.empty() || i >= inflated_len) {
      return std::string{};
    }
    len = std::min(len, inflated_len - i);

    std::string dst{};
    dst.reserve(len);
    for (; len; ++i, --len) {
      dst += Score(i) + '!';
    }
    return dst;
  }

  void ReverseAndComplement() {  // Watson-Crick base pairing
    is_reverse_complement ^= 1;
  }

  static std::atomic<std::uint32_t> num_objects;

  std::uint32_t id;  // (optional) initialize num_objects to 0
  std::string name;
  std::vector<std::uint64_t> deflated_data;
  std::vector<std::uint64_t>
      deflated_quality;  // (optional) Phred quality scores
  std::vector<std::uint32_t> qlvl;
  std::uint32_t inflated_len;
  bool is_reverse_complement;
};

}  // namespace exp
}  // namespace biosoup

#endif  // BIOSOUP_NUCLEIC_ACID_EXP_HPP_
