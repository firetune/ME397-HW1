#!/usr/bin/env python3
"""Build a full 'isotopes.csv' from NIST's
'Atomic Weights and Isotopic Compositions for All Elements' table.

Source:
  NIST PML — Atomic Weights and Isotopic Compositions for All Elements
  https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?all=all&ascii=ascii2&ele=&isotype=all

This script parses that ASCII table and writes a CSV with columns:
  element,symbol,A,mass_u,abundance_percent,stable

Notes:
  - We include naturally occurring isotopes with listed 'Isotopic Composition',
    which for a few elements can include very long-lived radioisotopes (e.g., 40K).
    These belong in the natural composition used for atomic-weight calculations.
  - 'stable' is set to True for these entries (meaning "include in natural weight").
  - Abundance is written in atom percent (composition × 100).
  - Isotopic masses are the 'Relative Atomic Mass' per-isotope values (u).
"""
from __future__ import annotations
import csv
import re
import sys
from pathlib import Path

try:
    import requests
except Exception as e:
    sys.exit("This builder requires the 'requests' package. Install with:\n  pip install requests")

NIST_URL = "https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?all=all&ascii=ascii2&ele=&isotype=all"

ELEMENTS = {
    1: ("Hydrogen","H"), 2: ("Helium","He"), 3: ("Lithium","Li"), 4: ("Beryllium","Be"), 5: ("Boron","B"),
    6: ("Carbon","C"), 7: ("Nitrogen","N"), 8: ("Oxygen","O"), 9: ("Fluorine","F"), 10: ("Neon","Ne"),
    11: ("Sodium","Na"), 12: ("Magnesium","Mg"), 13: ("Aluminium","Al"), 14: ("Silicon","Si"), 15: ("Phosphorus","P"),
    16: ("Sulfur","S"), 17: ("Chlorine","Cl"), 18: ("Argon","Ar"), 19: ("Potassium","K"), 20: ("Calcium","Ca"),
    21: ("Scandium","Sc"), 22: ("Titanium","Ti"), 23: ("Vanadium","V"), 24: ("Chromium","Cr"), 25: ("Manganese","Mn"),
    26: ("Iron","Fe"), 27: ("Cobalt","Co"), 28: ("Nickel","Ni"), 29: ("Copper","Cu"), 30: ("Zinc","Zn"),
    31: ("Gallium","Ga"), 32: ("Germanium","Ge"), 33: ("Arsenic","As"), 34: ("Selenium","Se"), 35: ("Bromine","Br"),
    36: ("Krypton","Kr"), 37: ("Rubidium","Rb"), 38: ("Strontium","Sr"), 39: ("Yttrium","Y"), 40: ("Zirconium","Zr"),
    41: ("Niobium","Nb"), 42: ("Molybdenum","Mo"), 43: ("Technetium","Tc"), 44: ("Ruthenium","Ru"), 45: ("Rhodium","Rh"),
    46: ("Palladium","Pd"), 47: ("Silver","Ag"), 48: ("Cadmium","Cd"), 49: ("Indium","In"), 50: ("Tin","Sn"),
    51: ("Antimony","Sb"), 52: ("Tellurium","Te"), 53: ("Iodine","I"), 54: ("Xenon","Xe"), 55: ("Cesium","Cs"),
    56: ("Barium","Ba"), 57: ("Lanthanum","La"), 58: ("Cerium","Ce"), 59: ("Praseodymium","Pr"), 60: ("Neodymium","Nd"),
    61: ("Promethium","Pm"), 62: ("Samarium","Sm"), 63: ("Europium","Eu"), 64: ("Gadolinium","Gd"), 65: ("Terbium","Tb"),
    66: ("Dysprosium","Dy"), 67: ("Holmium","Ho"), 68: ("Erbium","Er"), 69: ("Thulium","Tm"), 70: ("Ytterbium","Yb"),
    71: ("Lutetium","Lu"), 72: ("Hafnium","Hf"), 73: ("Tantalum","Ta"), 74: ("Tungsten","W"), 75: ("Rhenium","Re"),
    76: ("Osmium","Os"), 77: ("Iridium","Ir"), 78: ("Platinum","Pt"), 79: ("Gold","Au"), 80: ("Mercury","Hg"),
    81: ("Thallium","Tl"), 82: ("Lead","Pb"), 83: ("Bismuth","Bi"), 84: ("Polonium","Po"), 85: ("Astatine","At"),
    86: ("Radon","Rn"), 87: ("Francium","Fr"), 88: ("Radium","Ra"), 89: ("Actinium","Ac"), 90: ("Thorium","Th"),
    91: ("Protactinium","Pa"), 92: ("Uranium","U"), 93: ("Neptunium","Np"), 94: ("Plutonium","Pu"), 95: ("Americium","Am"),
    96: ("Curium","Cm"), 97: ("Berkelium","Bk"), 98: ("Californium","Cf"), 99: ("Einsteinium","Es"), 100: ("Fermium","Fm"),
    101: ("Mendelevium","Md"), 102: ("Nobelium","No"), 103: ("Lawrencium","Lr"), 104: ("Rutherfordium","Rf"), 105: ("Dubnium","Db"),
    106: ("Seaborgium","Sg"), 107: ("Bohrium","Bh"), 108: ("Hassium","Hs"), 109: ("Meitnerium","Mt"), 110: ("Darmstadtium","Ds"),
    111: ("Roentgenium","Rg"), 112: ("Copernicium","Cn"), 113: ("Nihonium","Nh"), 114: ("Flerovium","Fl"), 115: ("Moscovium","Mc"),
    116: ("Livermorium","Lv"), 117: ("Tennessine","Ts"), 118: ("Oganesson","Og"),
}

def fetch_text(url: str) -> str:
    import requests
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    return r.text

def parse_all(text: str):
    lines = text.splitlines()
    current_Z = None
    current_symbol = None
    current_A = None
    current_mass = None
    current_comp = None

    def flush():
        if (current_Z is not None and current_symbol and current_A is not None
            and current_mass is not None and current_comp is not None and current_comp != ""):
            try:
                mass = float(re.sub(r"[()#].*", "", str(current_mass)).strip())
            except Exception:
                return None
            try:
                comp = float(re.sub(r"[()#].*", "", str(current_comp)).strip()) * 100.0
            except Exception:
                return None
            name, sym = ELEMENTS.get(current_Z, (None, current_symbol))
            return {
                "element": name or current_symbol,
                "symbol": current_symbol,
                "A": int(current_A),
                "mass_u": mass,
                "abundance_percent": comp,
                "stable": True
            }
        return None

    rows = []
    for ln in lines:
        s = ln.strip()
        if s.startswith("Atomic Number ="):
            row = flush()
            if row:
                rows.append(row)
            current_Z = int(s.split("=")[1].strip())
            current_symbol = None
            current_A = None
            current_mass = None
            current_comp = None
        elif s.startswith("Atomic Symbol ="):
            current_symbol = s.split("=")[1].strip()
        elif s.startswith("Mass Number ="):
            row = flush()
            if row:
                rows.append(row)
            current_A = int(s.split("=")[1].strip())
            current_mass = None
            current_comp = None
        elif s.startswith("Relative Atomic Mass ="):
            current_mass = s.split("=")[1].strip()
        elif s.startswith("Isotopic Composition ="):
            current_comp = s.split("=")[1].strip()
    row = flush()
    if row:
        rows.append(row)
    return rows

def write_csv(rows, path: str = "isotopes.csv") -> None:
    import csv
    from pathlib import Path
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["element","symbol","A","mass_u","abundance_percent","stable"])
        w.writeheader()
        for r in rows:
            w.writerow(r)

def main(out="isotopes.csv"):
    txt = fetch_text(NIST_URL)
    rows = parse_all(txt)
    write_csv(rows, out)
    print(f"Wrote {len(rows)} isotopes to {out}")
    from collections import defaultdict
    sums = defaultdict(float)
    for r in rows:
        sums[r["symbol"]] += r["abundance_percent"]
    examples = ["H","O","Sn","Pb","W","Xe","K","Cl"]
    for sym in examples:
        if sym in sums:
            print(f"  {sym}: total abundance ~ {sums[sym]:.3f}% (should be ~100%)")

if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "isotopes.csv"
    main(out)
