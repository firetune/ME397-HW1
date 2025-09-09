"""
Isotope database loader for stable isotopes and natural abundances.
- Prefers loading from a CSV file (default: 'isotopes.csv' in the same folder)
- Falls back to a small built-in seed dataset (Sn) if CSV not found.

CSV schema (header required):
    element,symbol,A,mass_u,abundance_percent,stable

Notes:
- 'abundance_percent' should be atom % (number fraction × 100) for natural composition.
- Include only stable (or observationally stable) isotopes if you want the natural atomic weight.
- For elements without stable isotopes (e.g., Tc, Pm), leave them out (the calculator will warn).
"""
from __future__ import annotations
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

@dataclass(frozen=True)
class Isotope:
    element: str     # e.g., 'Tin'
    symbol: str      # e.g., 'Sn'
    A: int           # mass number
    mass_u: float    # isotopic mass in u
    abundance_percent: float  # natural abundance in % (by number)
    stable: bool = True

# Built-in seed (Tin) so the module works out of the box
_SEED: Dict[str, list[Isotope]] = {
    "Sn": [
        Isotope("Tin", "Sn", 112, 111.90482387, 0.97, True),
        Isotope("Tin", "Sn", 114, 113.9027827, 0.66, True),
        Isotope("Tin", "Sn", 115, 114.903344699, 0.34, True),
        Isotope("Tin", "Sn", 116, 115.90174280, 14.54, True),
        Isotope("Tin", "Sn", 117, 116.90295398, 7.68, True),
        Isotope("Tin", "Sn", 118, 117.90160657, 24.22, True),
        Isotope("Tin", "Sn", 119, 118.90331117, 8.59, True),
        Isotope("Tin", "Sn", 120, 119.90220163, 32.58, True),
        Isotope("Tin", "Sn", 122, 121.9034438,  4.63, True),
        Isotope("Tin", "Sn", 124, 123.9052766, 5.79, True),
    ]
}

def load_csv(path: str | Path = "isotopes.csv") -> Dict[str, list[Isotope]]:
    p = Path(path)
    db: Dict[str, list[Isotope]] = {}
    if not p.exists():
        return {}
    with p.open(newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Skip rows that aren't marked stable if you only want stable isotopes
            if str(row.get("stable", "True")).strip().lower() not in {"true","1","yes","y"}:
                continue
            try:
                iso = Isotope(
                    element = row["element"].strip(),
                    symbol  = row["symbol"].strip(),
                    A       = int(row["A"]),
                    mass_u  = float(row["mass_u"]),
                    abundance_percent = float(row["abundance_percent"]),
                    stable  = True
                )
            except Exception as e:
                raise ValueError(f"Bad row in CSV: {row!r}\n{e}")
            db.setdefault(iso.symbol, []).append(iso)
    # Keep isotopes sorted by mass number
    for sym in db:
        db[sym] = sorted(db[sym], key=lambda x: x.A)
    return db

def get_database(csv_path: str | Path = "isotopes.csv") -> Dict[str, list[Isotope]]:
    db = load_csv(csv_path)
    if db:
        return db
    # Fallback to seed
    return _SEED.copy()

def atomic_weight(symbol: str, db: Optional[Dict[str, list[Isotope]]] = None) -> float:
    symbol = symbol.strip().capitalize() if len(symbol)<=2 else symbol.strip()
    # Capitalize properly for one-/two-letter symbols
    if len(symbol) == 1:
        symbol = symbol.upper()
    elif len(symbol) == 2:
        symbol = symbol[0].upper() + symbol[1].lower()
    db = db or get_database()
    if symbol not in db:
        raise KeyError(f"No stable isotope data for element symbol '{symbol}'. "
                       "Provide 'isotopes.csv' with stable isotopes and abundances.")
    isotopes = db[symbol]
    total_abund = sum(i.abundance_percent for i in isotopes)
    if abs(total_abund - 100.0) > 0.5:
        # Warn via exception so users fix their data when building a complete table
        raise ValueError(f"Abundances for {symbol} sum to {total_abund:.3f}%, not ~100%. "
                         "Check your input CSV for that element.")
    return sum((i.abundance_percent/100.0) * i.mass_u for i in isotopes)
