# IGoR Math Module

A high-performance C++23 N-dimensional tensor library designed for zero-copy data manipulation, probabilistic modeling, and log-space stability. Built on top of standard `std::mdspan`, this module provides the mathematical foundation for IGoR's recombination inference engine.

## 🚀 Key Features

- **Standardized**: Fully built on C++23 `std::mdspan` (with a compatibility layer for C++20/Linux/Windows).
- **Zero-Copy Architecture**: Perform complex slicing and broadcasting operations without memory allocations or data replication.
- **Multidimensional Tensor**: A flexible `Tensor<T>` wrapper that supports runtime rank determination and dynamic shapes while maintaining `mdspan` performance.
- **Linalg Suite**:
    - **Element-wise Ops**: Add, subtract, multiply, divide, scale.
    - **Reductions**: Sum, Min/Max, ArgMax.
    - **Probabilistic**: Matrix multiplication, dot products, and normalization.
    - **Log-Space Stability**: High-precision `log_add_exp`, `log_sum_exp`, and `log_normalize` for numerically stable inference.
- **Zero-Copy Broadcasting**: NumPy-style broadcasting support using `layout_stride` to repeat dimensions without memory allocation.

---

## ⚡ Quick Start

```cpp
// Single include for all math functionality
#include <igor/Math.h>

using namespace igor::math;

// Create probability distribution
Tensor<double> probs({100});
for (size_t i = 0; i < 100; ++i)
    probs(i) = 1.0 + i;

// Normalize to sum = 1
linalg::normalize(probs, probs);

// Compute in log-space for stability
Tensor<double> log_probs({100});
for (size_t i = 0; i < 100; ++i)
    log_probs(i) = std::log(probs(i));

double stable_lse = linalg::log_sum_exp(log_probs);
// Returns ~0.0 (since sum(probs) = 1.0)
```

---

## 📦 Core Components

### 1. `Tensor<T>`
The primary container for N-dimensional data. It manages ownership of a contiguous memory buffer (`HybridBuffer`) while providing runtime-determined rank access.

```cpp
#include <igor/Math/Tensor.h>

// Create a 3D Tensor [64, 25, 12]
igor::math::Tensor<double> prob({64, 25, 12});

// Access elements
prob(10, 5, 2) = 0.5;

// Get a compile-time rank-3 view (std::mdspan)
auto view = prob.view<3>();

// Query shape
std::cout << "Shape: " << prob.shape()[0] << "×"
          << prob.shape()[1] << "×" << prob.shape()[2] << std::endl;
```

### 2. `Linalg`
A suite of free functions for operating on `mdspan` views or `Tensor` objects.

```cpp
#include <igor/Math/Linalg.h>
namespace linalg = igor::math::linalg;

// Normalize a probability distribution
linalg::normalize(prob, prob);

// Stable reduction in log-space
double lse = linalg::log_sum_exp(prob);
```

### 3. Broadcasting
Zero-copy dimension expansion for efficient tensor math using stride manipulation.

**Example: Shape [12] broadcasts to [64, 25, 12]**
```
Strides:  [0,  0,  1]  ← Zero strides repeat data
Memory:   [a, b, c, d, e, f, g, h, i, j, k, l]  (12 elements)
Result:   All 64×25 slices see the same 12 elements
```

```cpp
// Multiply a 3D tensor by a 1D marginal along the last dimension
// P(V,D,J) * P(J)
linalg::broadcast_multiply(p_vdj.view<3>(), p_j.view<1>(), out.view<3>());
```

### 4. Slicing
Efficiently extract sub-tensors using the `slice` method, which returns a zero-copy `std::mdspan` view.

```cpp
// Create a 3x4 Tensor
Tensor<double> matrix({3, 4});

// Get the 2nd row (index 1) as a rank-1 view
auto row_view = matrix.slice<2>(0, 1); 

// row_view is now a 1D mdspan of length 4
row_view[0] = 5.0; // Modifies original matrix
```

---

## 🛠 Compatibility Layer

The `MdspanCompat.h` header ensures consistent behavior across compilers and operating systems:

| Platform | Compiler | Implementation | Min Version |
| :--- | :--- | :--- | :--- |
| **macOS** | Clang | Native `std::mdspan` | 19.0+ |
| **Windows** | MSVC | Native `std::mdspan` | 19.38+ |
| **Linux** | GCC | Native `std::mdspan` | 14.0+ |
| **Legacy** | GCC/Clang | Kokkos `mdspan` | via Pixi |

---

## 🧩 Module Structure

- **`Math.h`**: Umbrella header involving all math functionality.
- **`Tensor.h/tpp`**: High-level N-dim wrapper.
- **`Linalg.h/tpp`**: Mathematical operations (Core & Tensor overloads).
- **`HybridBuffer.h/tpp`**: Efficient storage with small-buffer optimization.
- **`MdspanCompat.h`**: cross-platform `mdspan` polyfill.
- **`CMakeLists.txt`**: Interface library definition.

---

## 📋 Requirements

- **C++23** compiler (Clang 19+, GCC 14+, MSVC 19.38+)
- **CMake** 3.28+
- **Pixi** package manager (for dependencies)

### Building

```bash
pixi install
pixi run build
```

---

## 🎯 Design Philosophy

1. **Algorithm-Data Separation**: Operations work on views (`mdspan`), not containers
2. **Compile-Time Guarantees**: Rank checking at compile-time via templates
3. **Standard Library Only**: Zero external dependencies (except Kokkos fallback)
4. **IGoR-Specific**: Optimized for probabilistic inference workflows

---

## 🧪 Testing

The module is covered by a robust suite of unit tests verifying rank dispatch, numerical stability, and broadcasting logic.

```bash
# Build and run math tests
pixi run build
./build/bin/math_tests

# Expected output: All tests passed
# - Element-wise operations
# - Probabilistic operations
# - Log-space stability checks
# - Broadcasting logic verification
# - Tensor Slicing
```

---

*Developed by the IGoR Team for High-Performance Immunogenomics.*
