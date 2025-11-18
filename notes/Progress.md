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

## Nov-03-2025

1. **Normalized and Centred Segments**:

   Within each segment of flux data, the fluxes were normalized and centred, wherein a `NormSegments`
   data structure was created for clarity.


2. **Cleaned Normalized Segments**:

   Cleaned the data, by computing the median absolute deviation (MAD) of each segment of flux
   data, wherein a $\sigma = 1.4286 \times \text{MAD}$ was computed, to create a threshold
   for all data within each segment, wherein data that exceeded $6\sigma$, is replaced by
   `NaN`, to flag it as bad data. $k = 6$ for $k \sigma$, was chosen to ensure a conservative
   threshold, that really only flags data from possible outliers beyond noise, such as cosmic
   ray hits, detector glitch, etc. 

---

## Nov-05-2025

1. **Binned the flux (Not finished)**
   Created a function, that given a segment of fluxes, folds the data and bins it into so far,
   200 points. This doesn't technically work right now, but I've played around with it enough to fix
   it for next time.

---

## Nov-16-2025

1. **Deleted Binning functions**:

   Deleted the binning functions, because I had realized that I was starting that process too soon. I first 
   needed to consider the folding data.


2. **Folding Data**:

   Folded the time data into phases, for various trial periods. So far, for testing purposes, I am just going to 
   consider, a random single trial period, until I establish all the functions and processes needed to pull out
   the transit data.


3. **Created Plots**:

   Through interoping with python, in order to use `matplotlib`, I've created a function that can plot
   `Vec<f64>` plots, which I will routinely use, just to get a good look at my data.


---

## Nov-17-2025

1. **Binned Flux Data**

   Fluxes were binned into `N_BINS = 200` number of bins. So far, only did this for a particular trial
   period.

2. **Searched Bins to find Deepest Dip**
   
   Created a function that iterated over the 200 binned fluxes, and through using a sliding window that
   had a size corresponding to a given trial duration, I pulled out the deepest dip in the data, to be 
   used for SNR.


