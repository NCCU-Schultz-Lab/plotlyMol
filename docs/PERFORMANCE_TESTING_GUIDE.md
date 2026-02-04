# Performance Testing Guide

**Date:** 2026-02-03
**Purpose:** Quantitative performance testing and optimization for plotlyMol

---

## Overview

This guide provides methods for **quantitatively measuring** plotlyMol's performance, particularly for GUI/Streamlit applications. Instead of subjectively saying "the GUI is laggy," you can now measure exact rendering times, memory usage, and identify specific bottlenecks.

---

## Available Tools

### 1. Performance Testing Script

**Location:** [tests/test_performance.py](../tests/test_performance.py)

A standalone Python script that runs comprehensive benchmarks and saves results.

**Usage:**
```bash
# Run all benchmarks
python tests/test_performance.py

# Results saved to: benchmark_results/
```

**What it measures:**
- ✅ Rendering time vs molecule size
- ✅ Resolution impact on performance
- ✅ Vibration file parsing speed
- ✅ Vibration visualization modes (arrows, heatmap, animation)
- ✅ Animation frame count impact
- ✅ Memory usage for all operations

**Output:**
- CSV files with detailed results
- Console summary with recommendations
- Performance ratios and speedup calculations

---

### 2. Performance Benchmarking Notebook

**Location:** [examples/performance_benchmarking.ipynb](../examples/performance_benchmarking.ipynb)

Interactive Jupyter notebook for in-depth performance analysis.

**Features:**
- Interactive visualizations (matplotlib charts)
- Customizable benchmarks
- Statistical analysis (mean, std, min, max)
- Memory profiling with `tracemalloc`
- Real-time results display

**Usage:**
```bash
jupyter notebook examples/performance_benchmarking.ipynb
```

---

## Key Metrics

### 1. Rendering Time

**What it measures:** Time to generate 3D molecular visualization

**Factors affecting performance:**
- **Molecule size** (number of atoms)
- **Rendering mode** (`ball+stick`, `stick`, `vdw`)
- **Resolution** (sphere tessellation quality)
- **Bond count** (especially for complex molecules)

**Typical values (from benchmarks):**
- Small molecules (<20 atoms): 50-200 ms
- Medium molecules (20-50 atoms): 200-500 ms
- Large molecules (>50 atoms): 500-2000 ms

**Optimization tips:**
- Use `resolution=16` for preview (2-3x faster than default 32)
- Use `stick` mode for molecules >50 atoms
- Use `vdw` mode for very large systems (>100 atoms)

---

### 2. Memory Usage

**What it measures:** Peak memory consumption during rendering

**Typical values:**
- Small molecules: 5-20 MB
- Medium molecules: 20-50 MB
- Large molecules: 50-200 MB
- Animations: 2-5x base memory (depends on frame count)

**Memory profiling:**
```python
import tracemalloc

tracemalloc.start()
fig = draw_3D_rep(smiles="CCO", mode="ball+stick")
current, peak = tracemalloc.get_traced_memory()
tracemalloc.stop()

print(f"Peak memory: {peak / 1024 / 1024:.1f} MB")
```

---

### 3. Vibration Parsing Speed

**What it measures:** Time to parse vibrational data from quantum chemistry files

**Typical values (water molecule):**
- Gaussian (.log): 10-30 ms
- ORCA (.out): 10-30 ms
- Molden (.molden): 5-15 ms

**Scaling:** Roughly linear with file size and number of modes

**Optimization:**
- Cache parsed results with `@st.cache_resource` in Streamlit
- Parse once, visualize many times

---

### 4. Vibration Visualization Performance

**What it measures:** Time to add vibration overlays to molecular figure

**Typical values (water molecule):**
- Static arrows: +20-50 ms
- Heatmap coloring: +30-60 ms
- Animation (20 frames): 500-1500 ms

**Optimization tips:**
- Use static arrows for fastest visualization
- Limit animations to 20-30 frames during development
- Use lower resolution (16) for animation preview

---

## Quantifying GUI Lag

### Problem: "The Streamlit GUI feels laggy"

### Solution: Measure Specific Operations

#### Step 1: Identify the Slow Operation

Run the performance script to get baseline measurements:

```bash
python tests/test_performance.py
```

Look at the console output:

```
BENCHMARK: Rendering Performance vs Molecule Size
==================================================================

Water (3 atoms):
  ball+stick  :   85.3 ± 12.1 ms
  stick       :   45.2 ±  8.4 ms
  vdw         :   92.1 ± 15.3 ms

Cholesterol (74 atoms):
  ball+stick  : 1823.5 ± 145.2 ms  <-- THIS IS SLOW!
  stick       :  892.3 ±  67.8 ms
  vdw         : 2145.7 ± 178.9 ms
```

**Conclusion:** Large molecules in `ball+stick` mode are the bottleneck.

---

#### Step 2: Test Optimization Strategies

Compare different settings:

```python
# Benchmark different resolutions
python tests/test_performance.py

# Look at resolution results:
# Resolution 16: 456.2 ms (2.0x faster)
# Resolution 32: 912.3 ms (baseline)
# Resolution 64: 1834.5 ms (2.0x slower)
```

**Conclusion:** Using `resolution=16` gives 2x speedup with minimal visual quality loss.

---

#### Step 3: Measure in Your Actual GUI

Add timing to your Streamlit app:

```python
import time
import streamlit as st

# Before rendering
start_time = time.perf_counter()

# Your rendering code
fig = draw_3D_rep(smiles=smiles, mode=mode, resolution=resolution)

# After rendering
end_time = time.perf_counter()
render_time_ms = (end_time - start_time) * 1000

# Display to user
st.sidebar.metric("Render Time", f"{render_time_ms:.0f} ms")
```

Now you have **quantitative data** showing exactly how long operations take!

---

## Performance Testing Workflow

### For Development

1. **Run baseline benchmarks**
   ```bash
   python tests/test_performance.py
   ```

2. **Identify bottlenecks** from the output

3. **Make optimizations** (change resolution, mode, etc.)

4. **Re-run benchmarks** to measure improvement

5. **Compare results** using saved CSV files

---

### For Production

1. **Profile real-world molecules** that your users will visualize

2. **Add custom benchmark** in `performance_benchmarking.ipynb`:
   ```python
   # Your specific molecule
   my_smiles = "CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C"

   stats = benchmark_multiple_runs(
       draw_3D_rep,
       n_runs=5,
       smiles=my_smiles,
       mode="ball+stick",
       resolution=32
   )

   print(f"Your molecule: {stats['mean_time_ms']:.1f} ms")
   ```

3. **Set performance budgets**:
   - Interactive preview: <200 ms target
   - Final render: <2000 ms acceptable
   - Animation: <5000 ms total

4. **Monitor performance** over time as codebase evolves

---

## Streamlit-Specific Optimizations

### 1. Caching

**Problem:** Re-parsing files on every widget interaction

**Solution:** Use `@st.cache_resource`

```python
@st.cache_resource(show_spinner=False)
def cached_parse_vibrations(vib_bytes: bytes, filename: str):
    """Cache vibration parsing results."""
    # ... parsing logic ...
    return parse_vibrations(temp_path)
```

**Speedup:** 10-100x for repeated access

---

### 2. Performance Mode

**Problem:** High-resolution rendering is slow during exploration

**Solution:** Implement a performance toggle

```python
# In sidebar
perf_mode = st.sidebar.selectbox(
    "Mode",
    ["Balanced", "Performance"],
    help="Performance uses lower resolution for faster rendering"
)

# Use lower resolution in performance mode
resolution_used = 16 if perf_mode == "Performance" else 32
```

**Speedup:** ~2x faster rendering

---

### 3. Lazy Loading

**Problem:** Generating expensive visualizations upfront

**Solution:** Generate only when needed

```python
# Don't generate animation until user requests it
if vib_display_type == "Animation":
    with st.spinner("Creating animation..."):
        fig = create_vibration_animation(...)
else:
    # Regular figure (faster)
    fig = draw_3D_rep(...)
```

---

### 4. Profiling in Streamlit

Add a performance profiler to your app:

```python
import time
import streamlit as st

class PerformanceProfiler:
    def __init__(self):
        self.timings = {}

    def measure(self, name):
        """Context manager for timing operations."""
        class Timer:
            def __enter__(self_):
                self_.start = time.perf_counter()
                return self_

            def __exit__(self_, *args):
                end = time.perf_counter()
                self.timings[name] = (end - self_.start) * 1000

        return Timer()

    def display_metrics(self):
        """Show performance metrics in sidebar."""
        st.sidebar.markdown("### ⚡ Performance")
        for name, time_ms in self.timings.items():
            color = "🟢" if time_ms < 200 else "🟡" if time_ms < 1000 else "🔴"
            st.sidebar.metric(name, f"{time_ms:.0f} ms", delta=None)

# Usage in your app
profiler = PerformanceProfiler()

with profiler.measure("Molecule Rendering"):
    fig = draw_3D_rep(smiles=smiles, mode=mode)

with profiler.measure("Vibration Overlay"):
    fig = add_vibrations_to_figure(fig, vib_data, mode_number=1)

profiler.display_metrics()
```

Now you have **real-time performance monitoring** in your GUI!

---

## Interpreting Results

### Good Performance

```
Water (3 atoms):
  ball+stick  :   85.3 ms  ✅ Fast

Benzene (12 atoms):
  ball+stick  :  234.5 ms  ✅ Acceptable

Glucose (24 atoms):
  ball+stick  :  456.8 ms  ✅ Acceptable
```

**Response:** No optimization needed, feels snappy.

---

### Borderline Performance

```
Cholesterol (74 atoms):
  ball+stick  : 1823.5 ms  ⚠️ Noticeable lag
```

**Response:** Consider performance mode or optimize:
- Switch to `resolution=16`: ~900 ms ✅ Better
- Switch to `stick` mode: ~900 ms ✅ Better

---

### Poor Performance

```
Large Protein (200 atoms):
  ball+stick  : 8234.5 ms  ❌ Very slow (8+ seconds)
```

**Response:** Aggressive optimization needed:
- Use `stick` mode: ~3500 ms ⚠️ Still slow
- Use `vdw` mode: ~4200 ms ⚠️ Still slow
- Use `resolution=16`: ~4000 ms ⚠️ Still slow
- **Consider:** Downsampling, progressive loading, or WebGL optimization

---

## Advanced Profiling

### Line-by-Line Profiling

For deep analysis, use `line_profiler`:

```bash
pip install line-profiler
```

```python
# Add @profile decorator
@profile
def draw_3D_rep(smiles, mode, resolution):
    # ... function code ...

# Run with line profiler
kernprof -l -v your_script.py
```

---

### Memory Profiling

Use `memory_profiler` for detailed memory analysis:

```bash
pip install memory-profiler
```

```python
from memory_profiler import profile

@profile
def render_large_molecule():
    fig = draw_3D_rep(smiles="CC(C)CCCC(C)...", mode="ball+stick")
    return fig

render_large_molecule()
```

Run with:
```bash
python -m memory_profiler your_script.py
```

---

## Performance Budgets

Set targets for your application:

| Operation | Target | Acceptable | Needs Optimization |
|-----------|--------|------------|-------------------|
| Small molecule render | <100 ms | <300 ms | >500 ms |
| Medium molecule render | <300 ms | <800 ms | >1500 ms |
| Large molecule render | <800 ms | <2000 ms | >5000 ms |
| Vibration parsing | <50 ms | <200 ms | >500 ms |
| Animation generation (30 frames) | <2000 ms | <5000 ms | >10000 ms |
| GUI interaction | <100 ms | <300 ms | >500 ms |

**Use these benchmarks to guide optimization priorities!**

---

## Common Bottlenecks

### 1. Fibonacci Sphere Generation

**Issue:** Generating sphere vertices for atoms

**Impact:** Scales with resolution² and number of atoms

**Optimization:**
- Pre-compute sphere meshes for common resolutions
- Cache generated meshes
- Use lower resolution

---

### 2. Cylinder Mesh Generation

**Issue:** Creating bond cylinders

**Impact:** Scales with number of bonds and resolution

**Optimization:**
- Use `stick` mode (smaller cylinders)
- Pre-compute cylinder templates
- Reduce resolution

---

### 3. RDKit Coordinate Generation

**Issue:** 3D embedding from SMILES

**Impact:** Varies greatly with molecule size

**Optimization:**
- Use pre-computed coordinates when available
- Cache RDKit molecules
- Skip UFF optimization for preview

---

### 4. Plotly Figure Serialization

**Issue:** Converting figure to JSON for display

**Impact:** Scales with number of traces

**Optimization:**
- Reduce number of separate traces
- Combine meshes when possible
- Use Plotly's WebGL mode for large datasets

---

## Recommended Workflow

### Day-to-Day Development

1. Run benchmarks weekly to catch regressions
2. Add timing metrics to GUI for user visibility
3. Use performance mode by default during development
4. Profile slow operations with the notebooks

### Before Release

1. Run full benchmark suite
2. Test with real-world molecules from target users
3. Verify all operations meet performance budgets
4. Document expected performance characteristics

### After User Reports Lag

1. Ask user for specific molecule (SMILES or file)
2. Benchmark that exact molecule
3. Identify bottleneck (size, mode, resolution)
4. Provide specific recommendation
5. Consider adding automatic optimization suggestions to GUI

---

## Example: Diagnosing Lag

**User Report:** "The GUI is laggy when I visualize my molecules"

**Your Response:**

1. **Gather information:**
   ```
   Can you share:
   - The molecule (SMILES or file)
   - What rendering mode you're using
   - What operations feel slow
   ```

2. **Run benchmarks with their molecule:**
   ```python
   # In performance_benchmarking.ipynb
   user_smiles = "CC(C)CCCC(C)..."  # Their molecule

   for mode in ["ball+stick", "stick"]:
       for resolution in [16, 32]:
           stats = benchmark_multiple_runs(
               draw_3D_rep,
               smiles=user_smiles,
               mode=mode,
               resolution=resolution
           )
           print(f"{mode}, res={resolution}: {stats['mean_time_ms']:.1f} ms")
   ```

3. **Provide data-driven recommendation:**
   ```
   Results for your molecule (74 atoms):
   - ball+stick, resolution=32: 1823 ms (slow)
   - ball+stick, resolution=16: 912 ms (2x faster)
   - stick, resolution=16: 456 ms (4x faster)

   Recommendation: Use "stick" mode with Performance setting
   for this size molecule. This will reduce lag from 1.8s to 0.5s.
   ```

**Now you have objective, quantitative data to support your recommendation!**

---

## Summary

### ✅ Tools Available

1. **Performance test script** (`tests/test_performance.py`)
2. **Benchmark notebook** (`examples/performance_benchmarking.ipynb`)
3. **Profiling utilities** in notebooks
4. **This guide** for interpretation

### ✅ Metrics to Track

- Rendering time (ms)
- Memory usage (MB)
- Parsing speed (ms)
- Frame generation rate (ms/frame)

### ✅ Optimization Strategies

- Resolution adjustment (8-64)
- Mode selection (ball+stick, stick, vdw)
- Caching (@st.cache_resource)
- Performance mode toggle
- Lazy loading

### ✅ Next Steps

1. Run baseline benchmarks on your system
2. Test with your specific molecules
3. Add performance monitoring to your GUI
4. Set performance budgets for your use case
5. Profile and optimize bottlenecks

---

**Last Updated:** 2026-02-03
**Status:** Production Ready ✅
