# Progress

---

## Oct-29-2025

1. **Loaded Light Curve Data:**
   Opened a FITS file and read in the columns `TIME`, `PDCSAP_FLUX`, and `QUALITY`.

2. **Filtered Clean Data:**
   Kept only rows where `QUALITY == 0` and both `TIME` and `FLUX` are finite.
   Sorted the resulting data by time to ensure monotonic order.

---

## Oct-31-2025

1. **Checked Clean Data**:

   Added some `assert!` to test out the `t` (time) and `f` (flux) arrays.


2. **Segmented Gaps**:

   Found and segmented gaps in data, through comparing differences within `t`, by seeing if,
   $t_{i} - t_{i-1} > 5 \Delta t$. Then, we can be somewhat certain of a gap in which the satellite was looking
   at that star. $\Delta t$ is the median cadence, wherein it is the median of all difference pairs within `t`.  

---
