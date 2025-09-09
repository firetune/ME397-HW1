#!/usr/bin/env python3
"""Compute atomic masses for any element using a CSV of stable isotopes.

Usage:
    python atomic_weight.py --symbol Sn
    python atomic_weight.py --symbol Cu --csv isotopes.csv

If --csv is provided and exists, it will be used; otherwise the built-in seed
(Tin only) in isotopes_db.py is used.
"""
import argparse
from pathlib import Path
from isotopes_db import get_database, atomic_weight

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--symbol", required=True, help="Element symbol (e.g., Sn, Cu, O)")
    ap.add_argument("--csv", default="isotopes.csv", help="Path to isotopes CSV")
    args = ap.parse_args()

    db = get_database(args.csv)
    aw = atomic_weight(args.symbol, db)
    print(f"Atomic weight (natural) for {args.symbol}: {aw:.6f} u")

if __name__ == "__main__":
    main()
