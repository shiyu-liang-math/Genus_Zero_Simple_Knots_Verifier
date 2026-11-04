#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

struct Row {
  std::int32_t q;
  std::int32_t m;
  std::int32_t c;
  std::int32_t min_g;
  std::int32_t max_g;
  std::int32_t argmin;
  std::int32_t argmax;
  std::int32_t width;
  std::int32_t genus;
};
static_assert(sizeof(Row) == 9 * sizeof(std::int32_t));

std::int64_t extended_gcd(std::int64_t a, std::int64_t b, std::int64_t& x,
                          std::int64_t& y) {
  if (b == 0) {
    x = 1;
    y = 0;
    return a;
  }
  std::int64_t x1 = 0;
  std::int64_t y1 = 0;
  const auto divisor = extended_gcd(b, a % b, x1, y1);
  x = y1;
  y = x1 - (a / b) * y1;
  return divisor;
}

int inverse_mod(int value, int modulus) {
  std::int64_t x = 0;
  std::int64_t y = 0;
  if (extended_gcd(value, modulus, x, y) != 1) {
    throw std::runtime_error("m is not invertible modulo q^2");
  }
  x %= modulus;
  if (x < 0) x += modulus;
  return static_cast<int>(x);
}

struct BlockSummary {
  int total = 0;
  int min_prefix = 0;
  int argmin = 0;
  int max_prefix = 0;
  int argmax = 0;
};

class Workspace {
 public:
  explicit Workspace(int q_value) : q(q_value) {
    size = 1;
    while (size < q) size *= 2;
    const int tree_count = 2 * size;
    total.resize(tree_count);
    min_prefix.resize(tree_count);
    max_prefix.resize(tree_count);
    argmin.resize(tree_count);
    argmax.resize(tree_count);
    length.resize(tree_count);
    states.resize(q);
    values.resize(q);
    heads.resize(q);
    event_column.resize(2 * q);
    event_state.resize(2 * q);
    event_next.resize(2 * q);
    blocks.resize(q);
  }

  Row compute(int m, int c) {
    const int p = q * q;
    const int m_star = inverse_mod(m, p);
    const int ell_residue = static_cast<int>(
        (static_cast<std::int64_t>(m_star) * c * q) % p);
    const int ell = ell_residue == 0 ? p : ell_residue;
    if (ell % q != 0) throw std::runtime_error("ell is not divisible by q");
    const int t = ell / q;
    const int s0 = m_star % q;

    std::fill(heads.begin(), heads.end(), -1);
    event_count = 0;
    for (int b = 1; b < q; ++b) {
      const int column = b - 1;
      const int shift = static_cast<int>(
          ((static_cast<std::int64_t>(m_star) * b) / q) % q);
      states[column] = shift < t;
      values[column] = states[column] ? t - q : t;
      add_event((q - shift) % q, column, true);
      add_event((t - shift + q) % q, column, false);
    }
    const int endpoint_column = q - 1;
    states[endpoint_column] = 1 <= s0 && s0 <= t;
    values[endpoint_column] = states[endpoint_column] ? t - q : t;
    add_event((1 - s0 + q) % q, endpoint_column, true);
    add_event((t + 1 - s0 + q) % q, endpoint_column, false);

    build_tree();
    const int m0 = m % q;
    for (int h = 0; h < q; ++h) {
      if (h > 0) {
        for (int event = heads[h]; event >= 0; event = event_next[event]) {
          const int column = event_column[event];
          states[column] = event_state[event] != 0;
          update(column, states[column] ? t - q : t);
        }
      }
      const int a = static_cast<int>((static_cast<std::int64_t>(h) * m0) % q);
      blocks[a] = root_summary();
    }

    int offset = 0;
    int min_g = 0;
    int max_g = 0;
    int first_min = 0;
    int first_max = 0;
    for (int a = 0; a < q; ++a) {
      const auto& block = blocks[a];
      const int min_value = offset + block.min_prefix;
      const int min_position = a * q + block.argmin;
      if (min_position < p &&
          (min_value < min_g ||
           (min_value == min_g && min_position < first_min))) {
        min_g = min_value;
        first_min = min_position;
      }
      const int max_value = offset + block.max_prefix;
      const int max_position = a * q + block.argmax;
      if (max_position < p &&
          (max_value > max_g ||
           (max_value == max_g && max_position < first_max))) {
        max_g = max_value;
        first_max = max_position;
      }
      offset += block.total;
    }
    if (offset != 0) throw std::runtime_error("nonzero total increment");
    const int width = max_g - min_g;
    const int genus_numerator = width - 2 * q + 2;
    if (genus_numerator < 0 || genus_numerator % 2 != 0) {
      throw std::runtime_error("width does not give a nonnegative integral genus");
    }
    return Row{q, m, c, min_g, max_g, first_min, first_max, width,
               genus_numerator / 2};
  }

 private:
  int q;
  int size = 0;
  int event_count = 0;
  std::vector<int> total;
  std::vector<int> min_prefix;
  std::vector<int> max_prefix;
  std::vector<int> argmin;
  std::vector<int> argmax;
  std::vector<int> length;
  std::vector<std::uint8_t> states;
  std::vector<int> values;
  std::vector<int> heads;
  std::vector<int> event_column;
  std::vector<std::uint8_t> event_state;
  std::vector<int> event_next;
  std::vector<BlockSummary> blocks;

  void add_event(int h, int column, bool state) {
    if (h == 0) return;
    const int index = event_count++;
    event_column[index] = column;
    event_state[index] = state ? 1 : 0;
    event_next[index] = heads[h];
    heads[h] = index;
  }

  void set_leaf(int node, int value) {
    length[node] = 1;
    total[node] = value;
    if (value < 0) {
      min_prefix[node] = value;
      argmin[node] = 1;
    } else {
      min_prefix[node] = 0;
      argmin[node] = 0;
    }
    if (value > 0) {
      max_prefix[node] = value;
      argmax[node] = 1;
    } else {
      max_prefix[node] = 0;
      argmax[node] = 0;
    }
  }

  void pull(int node) {
    const int left = 2 * node;
    const int right = left + 1;
    const int left_length = length[left];
    length[node] = left_length + length[right];
    total[node] = total[left] + total[right];

    const int alternate_min = total[left] + min_prefix[right];
    const int alternate_min_position = left_length + argmin[right];
    min_prefix[node] = min_prefix[left];
    argmin[node] = argmin[left];
    if (alternate_min < min_prefix[node] ||
        (alternate_min == min_prefix[node] &&
         alternate_min_position < argmin[node])) {
      min_prefix[node] = alternate_min;
      argmin[node] = alternate_min_position;
    }

    const int alternate_max = total[left] + max_prefix[right];
    const int alternate_max_position = left_length + argmax[right];
    max_prefix[node] = max_prefix[left];
    argmax[node] = argmax[left];
    if (alternate_max > max_prefix[node] ||
        (alternate_max == max_prefix[node] &&
         alternate_max_position < argmax[node])) {
      max_prefix[node] = alternate_max;
      argmax[node] = alternate_max_position;
    }
  }

  void build_tree() {
    std::fill(total.begin(), total.end(), 0);
    std::fill(min_prefix.begin(), min_prefix.end(), 0);
    std::fill(max_prefix.begin(), max_prefix.end(), 0);
    std::fill(argmin.begin(), argmin.end(), 0);
    std::fill(argmax.begin(), argmax.end(), 0);
    std::fill(length.begin(), length.end(), 0);
    for (int index = 0; index < size; ++index) {
      const int node = size + index;
      if (index < q) set_leaf(node, values[index]);
    }
    for (int node = size - 1; node > 0; --node) pull(node);
  }

  void update(int index, int value) {
    int node = size + index;
    set_leaf(node, value);
    for (node /= 2; node > 0; node /= 2) pull(node);
  }

  BlockSummary root_summary() const {
    return BlockSummary{total[1], min_prefix[1], argmin[1], max_prefix[1],
                        argmax[1]};
  }
};

std::vector<int> read_presentations(const std::string& path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("could not open presentation input");
  std::vector<int> presentations;
  int m = 0;
  while (input >> m) presentations.push_back(m);
  if (presentations.empty()) throw std::runtime_error("presentation input is empty");
  return presentations;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc != 6) {
      std::cerr << "usage: g0sv-native Q PRESENTATIONS.txt OUTPUT.bin THREADS REPEAT\n";
      return 2;
    }
    const int q = std::stoi(argv[1]);
    const std::string input_path = argv[2];
    const std::string output_path = argv[3];
    const int requested_threads = std::max(1, std::stoi(argv[4]));
    const int repeat = std::max(1, std::stoi(argv[5]));
    if (q < 2) throw std::runtime_error("q must be at least 2");
    const auto presentations = read_presentations(input_path);
    std::vector<int> units;
    for (int c = 1; c < q; ++c) {
      if (std::gcd(c, q) == 1) units.push_back(c);
    }
    std::vector<Row> rows(presentations.size() * units.size());
    const int thread_count = std::min<int>(requested_threads, presentations.size());
    const auto started = std::chrono::steady_clock::now();
    for (int iteration = 0; iteration < repeat; ++iteration) {
      std::atomic<std::size_t> next_presentation{0};
      std::atomic<bool> failed{false};
      std::string failure_message;
      std::vector<std::thread> workers;
      workers.reserve(thread_count);
      for (int thread_index = 0; thread_index < thread_count; ++thread_index) {
        workers.emplace_back([&]() {
          try {
            Workspace workspace(q);
            while (!failed.load()) {
              const auto index = next_presentation.fetch_add(1);
              if (index >= presentations.size()) break;
              for (std::size_t unit_index = 0; unit_index < units.size(); ++unit_index) {
                rows[index * units.size() + unit_index] =
                    workspace.compute(presentations[index], units[unit_index]);
              }
            }
          } catch (const std::exception& error) {
            failed.store(true);
            failure_message = error.what();
          }
        });
      }
      for (auto& worker : workers) worker.join();
      if (failed.load()) throw std::runtime_error(failure_message);
    }
    const auto finished = std::chrono::steady_clock::now();
    const double seconds =
        std::chrono::duration_cast<std::chrono::duration<double>>(finished - started)
            .count();
    std::ofstream output(output_path, std::ios::binary | std::ios::trunc);
    if (!output) throw std::runtime_error("could not open binary output");
    output.write(reinterpret_cast<const char*>(rows.data()),
                 static_cast<std::streamsize>(rows.size() * sizeof(Row)));
    if (!output) throw std::runtime_error("failed while writing binary output");
    std::cerr << "rows=" << rows.size() << " repeat=" << repeat
              << " seconds=" << seconds
              << " rows_per_second=" << (rows.size() * repeat / seconds) << "\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "g0sv-native: " << error.what() << "\n";
    return 1;
  }
}
