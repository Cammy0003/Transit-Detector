# A simple transit-detection recipe 
This recipe is meant to be less in depth with the process of detrending, so as to not convolute the project. 

## 0) Inputs & hygiene (one screen of code)

* Use columns: **TIME**, **PDCSAP_FLUX**, **PDCSAP_FLUX_ERR**, filter **QUALITY == 0**.
* Sort by time, drop NaNs/Infs.
* Let ($\Delta t$) = median cadence (robust).
* **Segment on gaps**: start a new segment whenever ($t_{i}-t_{i-1} > 5,\Delta t$). This prevents gap edges from screwing with the search.

## 1) Minimal normalization (exactly one line per segment)

For each segment (s):
$$
F^{(s)}_{\text{norm}}(t) = \frac{F^{(s)}_{\text{raw}}(t)}{\operatorname{median}\big(F^{(s)}_{\text{raw}}\big)} - 1
$$
That’s it. Now your baseline is ~0 and “depth” is a **negative fraction** (e.g., ($-500\ \text{ppm} \approx -5\times10^{-4})$).

Optional but helpful (still simple): **clip wild outliers** with a one-sided rule to protect dips:

* Estimate noise via MAD: ($\sigma \approx 1.4826 \cdot \operatorname{MAD}(F^{(s)}_{\text{norm}})$).
* Replace points where ($F^{(s)}_{\text{norm}} > +6\sigma$) with the segment median (ignore negative side to keep transits intact).

> That’s your whole “detrend”: **none**. PDCSAP already removes common trends; segment + normalize + clip is enough for this project.

## 2) Box Least Squares (BLS) on a fixed grid (the core)

You’re fitting a simple **box model** to a **folded** light curve for each trial period. Very small math footprint, easy to implement.

**Grids**

* Period: ($P \in [0.5, 20]$) days (or whatever you want)
* Duration: ($d \in {0.5, 1, 2, 3, 4}$) hours (convert to days)
* Phase: continuous in ($[0, P)$) but we’ll find phase efficiently (below)

**Fold**
For a trial ($P$), compute phases ($\phi_i = \big(t_i \bmod P\big)/P \in [0,1)$).

**Fast in-transit search with prefix sums**

1. Bin phases into ($B$) bins (e.g., $(B=200)$), keep:

    * mean flux per bin ($y_b$) and weights ($w_b = 1/\sigma_b^2$) (or equal weights if you want it simpler).
2. For a trial duration ($d$) (in phase units ($q = d/P)$), you want the **contiguous arc** ($[\phi_0, \phi_0+q)$) whose average flux is most negative.

    * On the circular array ($y_b$), use a **sliding window** over ($L=\lfloor qB \rfloor$) bins with prefix sums to find the minimum weighted average in ($O(B)$).
3. The “depth” estimate for that arc is simply the mean inside window minus the mean outside:
   $$
   \hat{\delta} ;=; \frac{\sum_{\text{in}} w_b y_b}{\sum_{\text{in}} w_b} ;-; \frac{\sum_{\text{out}} w_b y_b}{\sum_{\text{out}} w_b}
   $$
   It should be **negative** for transits.

**Power / SNR**
Estimate per-point noise (\sigma) from the **out-of-window** bins (MAD → ($\sigma$)). Then a simple score:
$$
\text{SNR}(P,d) = \frac{-\hat{\delta}}{\sigma}\sqrt{N_{\text{in}}}
$$
where ($N_{\text{in}}$) is the number of points (or total weight) inside the best window. Keep the **best SNR** over phases (found by the sliding window) for each ($(P,d)$), and then keep the best over ($d$) for each ($P$). Your **periodogram** is ($P \mapsto \max_d \text{SNR}(P,d)$).

This is the essence of BLS, implemented with basic arrays and prefix sums. No external libs required, no spline anything.

## 3) Pick the candidate & refine

* Select the **top peak** in the BLS periodogram: ($(P^*, d^*, \phi^*)$) with max SNR.
* **Fine-tune** around ($P^*$) with a smaller step (e.g., ±2% range, 5–10× finer resolution). Re-estimate ($(P,d,\phi)$).
* Report:

    * ($P$) (days), duration ($d$) (hours), depth ($\delta$) (ppm), SNR
    * Number of observed transits (how many windows land on data)
    * A quick **phase-folded plot** (bin and show in-transit vs out-of-transit)

## 4) Sanity checks 
**To be done for `v2.0` of the project!

* **Odd/even**: Fold on ($P$) and ($2P$). If the two dips at (2P) differ in depth by ($\gtrsim 3\sigma$), it smells like an eclipsing binary → flag.
* **Secondary search** at phase ~0.5 on ($P$). If there’s a shallower “secondary” that’s still strong → flag as EB.
* **Depth stability** across segments: compute depth per segment; if wildly inconsistent → flag (likely systematics).





