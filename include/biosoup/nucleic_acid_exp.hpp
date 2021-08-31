// Copyright (c) 2020 Robert Vaser

#ifndef BIOSOUP_NUCLEIC_ACID_EXP_HPP_
#define BIOSOUP_NUCLEIC_ACID_EXP_HPP_

#include <algorithm>
#include <array>
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
    static auto freqs = FreqsBufType{};
    static auto levels = LvlsBuftype{};

    deflated_quality.reserve(std::ceil(quality_len / 32.f));
    qlvl.reserve(std::ceil(quality_len / static_cast<float>(kBlockSize)));

    for (std::uint32_t i = 0; i < quality_len; i += kBlockSize) {
      std::uint32_t j = std::min(i + kBlockSize, quality_len);
      std::fill(freqs.begin(), freqs.end(), 0);

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
        if (++freqs[v] > mode_freq) {
          mode = v;
        }
      }

      mean /= j - i;

      if (mean < mode) {
        float step = static_cast<float>(mode - min) / 3.f;
        levels[0] = std::round(mode - 2 * step);
        levels[1] = std::round(mode - step);
        levels[2] = mode;
        step = static_cast<float>(max - mode) / 2.f;
        levels[3] = std::round(mode + step);
      } else {
        float step = static_cast<float>(mode - min) / 2;
        levels[0] = std::round(mode - step);
        levels[1] = mode;
        step = static_cast<float>(max - mode) / 3.f;
        levels[2] = std::round(mode + step);
        levels[3] = std::round(mode + 2 * step);
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
    }
  }

  NucleicAcid(const NucleicAcid&) = default;
  NucleicAcid& operator=(const NucleicAcid&) = default;

  NucleicAcid(NucleicAcid&&) = default;
  NucleicAcid& operator=(NucleicAcid&&) = default;

  ~NucleicAcid() = default;

  std::uint64_t Code(std::uint32_t i) const noexcept {
    std::uint64_t x = 0;
    if (is_reverse_complement) {
      i = inflated_len - i - 1;
      x = 3;
    }
    return ((deflated_data[i >> 5] >> ((i << 1) & 63)) & 3) ^ x;
  }

  std::uint8_t Score(std::uint32_t i) const noexcept {
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

  void ReverseAndComplement() noexcept {  // Watson-Crick base pairing
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

 private:
  static constexpr std::uint32_t kBlockSize = 128;
  static constexpr std::uint32_t kLvlsCap = 4;

  using FreqsBufType = std::array<std::uint8_t, kBlockSize>;
  using LvlsBuftype = std::array<std::uint8_t, kLvlsCap>;
};
}  // namespace exp
}  // namespace biosoup

#endif  // BIOSOUP_NUCLEIC_ACID_EXP_HPP_
