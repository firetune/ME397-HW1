# Isotope Atomic Weight Toolkit

Files:
- build_isotopes_csv.py — downloads & parses NIST's 'All Elements' table and writes a full isotopes.csv
- isotopes_db.py       — loader + API (`atomic_weight(symbol)`) with a Tin seed fallback
- atomic_weight.py     — CLI: `python atomic_weight.py --symbol Sn --csv isotopes.csv`

Steps:
1) Ensure you have Python 3.8+ and `pip install requests`.
2) Run `python build_isotopes_csv.py` to generate `isotopes.csv` (or `python build_isotopes_csv.py all_isotopes.csv`).
3) Compute any element's atomic weight:
   `python atomic_weight.py --symbol Xe --csv isotopes.csv`

Data notes:
- Abundances are atom percent for naturally occurring isotopes (long-lived included where part of natural mix).
- For elements with no stable isotopes (e.g., Tc, Pm), natural atomic weights aren't defined the same way; the CSV will not include them.
