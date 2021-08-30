#include <cstdint>
#include <memory>
#include <string>
#include <vector>

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
  MockData(std::string name, std::string data, std::string quality)
      : name(std::move(name)),
        data(std::move(data)),
        quality(std::move(quality)) {}

  std::string name;
  std::string data;
  std::string quality;
};

const MockData& GetFatMock() {
  static std::unique_ptr<MockData> obj_ptr;
  if (obj_ptr) {
    return *obj_ptr;
  } else {
    std::string data;
    std::string quality;
    while (data.size() < kMockMinLen) {
      data += kMockDataUnit;
      quality += kMockQualityUnit;
    }

    obj_ptr = std::unique_ptr<MockData>(
        new MockData("kMockFatSeq", std::move(data), std::move(quality)));
    return *obj_ptr;
  }
}

}  // namespace detail
}  // namespace benchmark
}  // namespace biosoup

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};
std::atomic<std::uint32_t> biosoup::exp::NucleicAcid::num_objects{0};

/* CONSTRUCTORS */
static void BM_Constructor(benchmark::State& state) {
  using namespace biosoup;
  namespace bmd = biosoup::benchmark::detail;
  const auto kMockData = bmd::GetFatMock();
  for (auto _ : state) {
    // const auto sequence =
    ::benchmark::DoNotOptimize(
        NucleicAcid(kMockData.name, kMockData.data, kMockData.quality));
    ::benchmark::ClobberMemory();
  }
}

BENCHMARK(BM_Constructor);

static void BM_Constructor_exp(benchmark::State& state) {
  using namespace biosoup;
  namespace bmd = biosoup::benchmark::detail;
  const auto kMockData = bmd::GetFatMock();
  for (auto _ : state) {
    // const auto sequence =
    ::benchmark::DoNotOptimize(
        exp::NucleicAcid(kMockData.name, kMockData.data, kMockData.quality));
    ::benchmark::ClobberMemory();
  }
}

BENCHMARK(BM_Constructor_exp);

/* QUALITY INFLATION */
static void BM_Quality(benchmark::State& state) {
  using namespace biosoup;
  namespace bmd = biosoup::benchmark::detail;

  const auto kMockSeqData = bmd::GetFatMock();
  const auto kSequence =
      NucleicAcid(kMockSeqData.name, kMockSeqData.data, kMockSeqData.quality);

  for (auto _ : state) {
    auto const quality = kSequence.InflateQuality(0);
    ::benchmark::DoNotOptimize(quality.data());
    ::benchmark::ClobberMemory();
  }
}

BENCHMARK(BM_Quality);

static void BM_Quality_exp(benchmark::State& state) {
  using namespace biosoup;
  namespace bmd = biosoup::benchmark::detail;

  const auto kMockSeqData = bmd::GetFatMock();
  const auto kSequence =
      NucleicAcid(kMockSeqData.name, kMockSeqData.data, kMockSeqData.quality);

  for (auto _ : state) {
    auto const quality = kSequence.InflateQuality(0);
    ::benchmark::DoNotOptimize(quality.data());
    ::benchmark::ClobberMemory();
  }
}

BENCHMARK(BM_Quality_exp);

BENCHMARK_MAIN();
