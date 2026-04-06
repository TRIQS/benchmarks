#!/usr/bin/env python
"""Download heavy benchmark data (chi3, chi4, > 50 MB) from external host.

This script is a skeleton for future use. The external data host
(Zenodo or institutional) will be decided when the first heavy dataset
is generated.

Usage:
  python scripts/download_heavy_data.py                    # download all
  python scripts/download_heavy_data.py Hubbard_Atom       # one model
  python scripts/download_heavy_data.py --list             # list available
"""

import argparse
import os
import sys


# Placeholder: will be filled with actual URLs when external hosting is set up
DATA_MANIFEST = {
    # 'Hubbard_Atom': {
    #     'pyed_chi4.h5': 'https://zenodo.org/record/XXXXX/files/pyed_chi4.h5',
    # },
}


def download_file(url, dest):
    """Download a file from url to dest."""
    import urllib.request
    print(f"  Downloading {os.path.basename(dest)} ...")
    urllib.request.urlretrieve(url, dest)
    print(f"  -> {dest}")


def main():
    parser = argparse.ArgumentParser(description="Download heavy benchmark data")
    parser.add_argument('model', nargs='?', help="Model to download (default: all)")
    parser.add_argument('--list', action='store_true', help="List available datasets")
    args = parser.parse_args()

    if not DATA_MANIFEST:
        print("No heavy data available for download yet.")
        print("This script will be populated when external hosting is configured.")
        return

    if args.list:
        for model, files in DATA_MANIFEST.items():
            print(f"\n{model}:")
            for fname, url in files.items():
                print(f"  {fname}")
        return

    models = {args.model: DATA_MANIFEST[args.model]} if args.model else DATA_MANIFEST

    for model, files in models.items():
        print(f"\n=== {model} ===")
        results_dir = os.path.join(model, 'results')
        os.makedirs(results_dir, exist_ok=True)
        for fname, url in files.items():
            dest = os.path.join(results_dir, fname)
            if os.path.exists(dest):
                print(f"  {fname} already exists, skipping")
                continue
            download_file(url, dest)


if __name__ == '__main__':
    main()
