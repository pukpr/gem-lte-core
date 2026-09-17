import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Load the data
data_path = 'amo/lte_results.csv'
data = np.loadtxt(data_path, delimiter=',')

# Column 1: decimal year (index 0)
# Column 3: amplitude (index 2)
t = data[:, 0]
y = data[:, 2]

# Filter region after 1885.0
mask = t > 1885.0
t_filtered = t[mask]
y_filtered = y[mask]

# Constants
t_min = t_filtered.min()
t_max = t_filtered.max()
T = t_max - t_min
t_prime = t_filtered - t_min

# Find minimum N for correlation >= 0.85
best_N = None
best_corr = 0.0

for N in range(1, 100):
    # Construct design matrix
    X = [np.ones_like(t_prime)]
    for n in range(1, N + 1):
        X.append(np.cos(2 * np.pi * n * t_prime / T))
        X.append(np.sin(2 * np.pi * n * t_prime / T))
    X = np.column_stack(X)
    
    # Solve least squares
    coeffs, residuals, rank, s = np.linalg.lstsq(X, y_filtered, rcond=None)
    
    # Calculate fit and correlation
    y_fit = X @ coeffs
    corr = np.corrcoef(y_filtered, y_fit)[0, 1]
    
    if corr >= 0.85:
        best_N = N
        best_corr = corr
        break

print(f"Optimal N (harmonics): {best_N}")
print(f"Number of terms: {2 * best_N + 1}")
print(f"Correlation coefficient: {best_corr:.6f}")

# Re-run for the best N to get coefficients
X = [np.ones_like(t_prime)]
for n in range(1, best_N + 1):
    X.append(np.cos(2 * np.pi * n * t_prime / T))
    X.append(np.sin(2 * np.pi * n * t_prime / T))
X = np.column_stack(X)
coeffs, residuals, rank, s = np.linalg.lstsq(X, y_filtered, rcond=None)
y_fit = X @ coeffs

# Extract coefficients
a0 = coeffs[0]
a_cos = []
b_sin = []
for n in range(1, best_N + 1):
    a_cos.append(coeffs[2 * n - 1])
    b_sin.append(coeffs[2 * n])

# Save coefficients to a CSV file for easy access/reusability
coeff_df = pd.DataFrame({
    'Harmonic_n': range(1, best_N + 1),
    'Frequency_Hz_per_year': [n / T for n in range(1, best_N + 1)],
    'Period_years': [T / n for n in range(1, best_N + 1)],
    'a_n_cos': a_cos,
    'b_n_sin': b_sin
})
coeff_df.to_csv('fourier_coefficients.csv', index=False)

# Write a summary text file with the explicit formula and details
summary_content = f"""Fourier Series Fit Summary
===========================
Target: Fit column 3 (amplitude) vs column 1 (decimal year) for 'amo/lte_results.csv' for t > 1885.0.
Condition: Correlation Coefficient >= 0.85.

Results:
--------
- Number of harmonics (N): {best_N}
- Total number of terms: {2 * best_N + 1} (1 constant, {best_N} cosine terms, {best_N} sine terms)
- Fit Time Range: {t_min:.6f} to {t_max:.6f}
- Fundamental Period (T): {T:.6f} years
- Constant term (a0): {a0:.6e}
- Pearson Correlation Coefficient: {best_corr:.6f}

Mathematical Formula:
---------------------
f(t) = a0 + sum_{{n=1}}^{{{best_N}}} [ a_n * cos(2 * pi * n * (t - {t_min:.6f}) / {T:.6f}) + b_n * sin(2 * pi * n * (t - {t_min:.6f}) / {T:.6f}) ]

The top 5 most significant harmonics by magnitude (sqrt(a_n^2 + b_n^2)):
"""

magnitudes = np.sqrt(np.array(a_cos)**2 + np.array(b_sin)**2)
top_indices = np.argsort(magnitudes)[::-1][:5]
for idx in top_indices:
    n = idx + 1
    mag = magnitudes[idx]
    period = T / n
    summary_content += f"  - Harmonic n={n:2d} (Period = {period:5.2f} years): cos_coeff = {a_cos[idx]:.4e}, sin_coeff = {b_sin[idx]:.4e}, magnitude = {mag:.4e}\n"

with open('fourier_fit_summary.txt', 'w') as f:
    f.write(summary_content)

print("Saved fourier_fit_summary.txt and fourier_coefficients.csv")

# Plotting the results
plt.figure(figsize=(12, 6))
plt.plot(t_filtered, y_filtered, label='Original Data (Amplitude)', color='blue', alpha=0.6, linewidth=1.5)
plt.plot(t_filtered, y_fit, label=f'Fourier Fit (N={best_N}, r={best_corr:.3f})', color='red', linestyle='--', linewidth=1.5)
plt.title(f'Fourier Series Fit of AMO Amplitude (t > 1885.0)\nCorrelation Coefficient: {best_corr:.4f}', fontsize=14)
plt.xlabel('Decimal Year', fontsize=12)
plt.ylabel('Amplitude', fontsize=12)
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend(fontsize=11)
plt.tight_layout()
plt.savefig('amo_fourier_fit.png', dpi=150)
plt.close()

print("Generated and saved amo_fourier_fit.png")
