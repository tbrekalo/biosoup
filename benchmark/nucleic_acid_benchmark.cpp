#include <cstdint>
#include <string>

#include "benchmark/benchmark.h"
#include "biosoup/nucleic_acid.hpp"
#include "biosoup/nucleic_acid_exp.hpp"

namespace biosoup {
namespace benchmark {
namespace detail {

/* clang-format off */
constexpr auto kMockDataUnit =    "GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT";
constexpr auto kMockQualityUnit = "!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65";
/* clang-format on */

constexpr std::size_t kMockMinLen = 1e4;

struct MockData {
  MockData() : id("@kMockSeq") {
    while (data.size() < kMockMinLen) {
      data += kMockDataUnit;
      quality += kMockQualityUnit;
    }
  }

  std::string id;
  std::string data;
  std::string quality;
};

}  // namespace detail
}  // namespace benchmark
}  // namespace biosoup

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};
std::atomic<std::uint32_t> biosoup::exp::NucleicAcid::num_objects{0};

static void BM_constructor(::benchmark::State& state) {
  using namespace biosoup;
  namespace bmd = biosoup::benchmark::detail;
  const auto kMockData = bmd::MockData();
  for (auto _ : state) {
    const auto obj =
        NucleicAcid(kMockData.id, kMockData.data, kMockData.quality);
  }
}

BENCHMARK(BM_constructor);

static void BM_constructor_exp(::benchmark::State& state) {
  using namespace biosoup;
  namespace bmd = biosoup::benchmark::detail;
  const auto kMockData = bmd::MockData();
  for (auto _ : state) {
    const auto boj =
        exp::NucleicAcid(kMockData.id, kMockData.data, kMockData.quality);
  }
}

BENCHMARK(BM_constructor_exp);

BENCHMARK_MAIN();
