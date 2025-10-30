# Progress

---

## Oct-29-2025

1. **Loaded Light Curve Data:**
   Opened a FITS file and read in the columns `TIME`, `PDCSAP_FLUX`, and `QUALITY`.

2. **Filtered Clean Data:**
   Kept only rows where `QUALITY == 0` and both `TIME` and `FLUX` are finite.
   Sorted the resulting data by time to ensure monotonic order.

3. **Normalized the Flux:**
Computed the median flux and divided all flux values by it, subtracting 1 to express variations as fractional deviations (`f/median - 1`).

---
