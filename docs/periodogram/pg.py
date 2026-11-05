import argparse
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import signal

# Define region metadata (column mapping for NOAA CPC ERSSTv5 ASCII format)
REGION_MAP = {
    "nino12": {"col": "ANOM",   "name": "Niño 1+2 (0-10S, 90W-80W)"},
    "nino3":  {"col": "ANOM.1", "name": "Niño 3 (5N-5S, 150W-90W)"},
    "nino34": {"col": "ANOM.2", "name": "Niño 3.4 (5N-5S, 170W-120W)"},
    "nino4":  {"col": "ANOM.3", "name": "Niño 4 (5N-5S, 160E-150W)"}
}

def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute and plot the Power Spectral Density (PSD) periodogram for El Niño SST indices."
    )
    
    # Target region selection argument
    parser.add_argument(
        "-r", "--region",
        choices=["nino4", "nino34", "nino3", "nino12"],
        default="nino4",
        help="Select the El Niño SST index region (default: nino4)"
    )
    
    # Spectral estimation options
    parser.add_argument(
        "--nperseg",
        type=int,
        default=256,
        help="Length of each segment for Welch's PSD calculation (default: 256)"
    )
    
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Save plot to a file path (e.g., plot.png) instead of displaying"
    )

    return parser.parse_args()

def main():
    args = parse_args()
    region_info = REGION_MAP[args.region]
    
    # Download dataset from NOAA CPC
    url = "https://www.cpc.ncep.noaa.gov/data/indices/ersst5.nino.mth.91-20.ascii"
    df = pd.read_csv(url, sep=r'\s+')
    
    # Extract SST anomaly vector
    anom_data = df[region_info["col"]].values
    
    # Compute Welch PSD (fs = 12 samples/year for monthly data)
    fs = 12.0
    freqs, psd = signal.welch(anom_data, fs=fs, nperseg=args.nperseg, scaling='density')
    
    # Convert frequencies (cycles/year) to period (years)
    periods = 1.0 / freqs[1:]
    psd_vals = psd[1:]
    
    # Plot Periodogram
    plt.figure(figsize=(10, 5))
    plt.plot(periods, psd_vals, color='#004d40', lw=2, label=f'{region_info["name"]} PSD')
#    plt.xscale('log')
#    plt.yscale('log')
    plt.xlim(0.5, 20)
    plt.xlabel('Period (Years)', fontsize=12)
    plt.ylabel('Power Spectral Density (deg C^2 / Year)', fontsize=12)
    plt.title(f'Power Spectrum: {region_info["name"]}', fontsize=14)
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend()
    plt.tight_layout()
    
    if args.output:
        plt.savefig(args.output, dpi=300)
        print(f"Saved plot to {args.output}")
    else:
        plt.show()

if __name__ == "__main__":
    main()
