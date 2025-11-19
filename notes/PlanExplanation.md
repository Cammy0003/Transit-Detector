### The Big Picture: What Are We Doing?

Our goal is to find a "transit"—a tiny, repeating dip in a star's brightness caused by a planet passing in front of it.

Your input is a file from the TESS space telescope, which is basically a long list of two-column data: **Time** and **Brightness (Flux)**.

The "recipe" is a simple, 5-step program to sift through this data and shout "Found one!" when it sees a pattern that looks like a planet.

---

### 0. Data Cleaning (Hygiene)

Before you can do any science, you have to clean the data. This is just like data cleaning in any other software field.

* **Use Good Data:** The file tells you which data points are "good" (`QUALITY == 0`). We ignore everything else (e.g., points where the telescope was jittering). We also ignore any `NaN` (Not a Number) or `Inf` (Infinity) values, as they'll break our math.
* **Sort by Time:** Obvious, but essential. We need the data in chronological order.
* **Segment on Gaps:** TESS doesn't watch a star 24/7. It takes data for ~27 days, then stops, moves, and might look at it again later. This creates big "gaps" in the data.
    * This step finds the *typical* time between measurements (the "median cadence," maybe 2 minutes).
    * If it finds a gap that's *way* bigger (e.g., 5 times the median), it splits the data into separate **segments**.
    * **Why?** So you don't try to do math that connects data from Tuesday with data from two Fridays from now. You'll analyze `Segment 1` and `Segment 2` as independent chunks.



---

### 1. Normalization (Setting a "Zero" Point)

A star's brightness might slowly change for its *own* reasons (e.g., star spots). We only care about the *planet dips*. This step removes the slow, boring changes.

* **The Formula:** For each segment, you find the *median* (middle) brightness value. Then, for every data point, you calculate:
  `Normalized_Flux = (Raw_Flux / Median_Flux) - 1`
* **Why this works:**
    * A normal data point (at the median) becomes: `(Median / Median) - 1 = 0`.
    * A dip (say, 0.1% darker) becomes: `(0.999 * Median / Median) - 1 = -0.001`.
* **The Result:** Your data is now "normalized." The star's normal brightness is `0`, and a planet dip is a small *negative number*. This is much easier to work with.
* **Optional Clipping (Highly Recommended):**
    * Sometimes a cosmic ray hits the detector, creating a *huge* spike in brightness. This is junk data and can mess up your search.
    * This step finds the typical "noise" level (`sigma`).
    * It then finds any points that are way too *bright* (e.g., `> 6 * sigma`) and just resets them to `0`.
    * **Crucially:** It *only* clips the positive (bright) spikes. We **leave the negative dips alone** because those might be our planets!

---

### 2. The "Box" Search (The Core Algorithm)

This is the main event: **Box Least Squares (BLS)**.

We're looking for a repeating, *box-shaped* dip. To do this, we have to "brute force" a bunch of guesses. We need to guess:
1.  **The Period:** How often does the dip repeat? (e.g., every 3.5 days?)
2.  **The Duration:** How long does the dip last? (e.g., for 2 hours?)

Here's the process for *each* `Period` we guess (e.g., we might test 4,000 periods between 0.5 and 20 days):

* **A. Fold the Data (Stacking)**
    * Let's test a `Period` of **3 days**.
    * You take your data and "fold" it. Imagine wrapping your data (a long string of numbers) around a drum that is 3 days in circumference.
    * `Day 1` stacks on top of `Day 0`.
    * `Day 2` stacks on top of `Day 1`.
    * `Day 3` stacks on top of `Day 0`.
    * `Day 4` stacks on top of `Day 1`.
    * ...and so on.
    * **The Magic:** If a planet *also* has a 3-day period, all its dips will stack up in the *exact same spot*. The "planet" signal gets stronger, and all the random noise averages out. This "stacked" view is called a "phase-folded" light curve.



* **B. Bin & Search (Finding the Dip)**
    * To make this *fast*, we don't look at all 50,000 data points at once. We "bin" the folded data into, say, 200 time-slices. We just find the *average* brightness in each bin.
    * Now we test our `Duration` guesses (e.g., 2 hours). A 2-hour duration might be 5 bins wide.
    * We use a "sliding window" to check every 5-bin-wide "box" (Bins 1-5, Bins 2-6, Bins 3-7...) to find the one that is the *darkest* (most negative).
    * **Prefix Sums:** This is a classic computer science trick to make the sliding window search *extremely* fast (O(N) instead of O(N\*W)). You can find the sum of *any* window in constant time (O(1)) after one initial pass. This is a great thing to mention in an interview.

* **C. Score the Guess (SNR)**
    * Once you find the best "box" for this `Period`, you give it a score: the **Signal-to-Noise Ratio (SNR)**.
    * In simple terms, the formula is:
      `SNR = (How Deep is the Dip?) / (How Noisy is the Data?) * sqrt(Number of Dips Seen)`
    * This is the *most important* concept. A high SNR means:
        1.  The dip is **deep** (easy to see).
        2.  The data *outside* the dip is **quiet** (low noise).
        3.  The dip was seen **many times** (e.g., 8 transits, not just 1). This `sqrt(N_in)` part is key—it rewards repeating signals.

You repeat A, B, and C for all 4,000 of your `Period` guesses. Your final output from this step is a list of `(Period, SNR)` scores. This list, when plotted, is your **Periodogram**.

---

### 3. Pick the Candidate (Find the Winner)

This part is easy. You look at your Periodogram (the "SNR vs. Period" chart) and find the **highest spike**.

* The `Period` at that highest spike is your **best guess** for the planet's period.
* The `Duration` and `Phase` (where the dip is) associated with that spike are your other best guesses.
* You can then "fine-tune" the search by re-running it, but this time only testing periods *very* close to your best guess.

---

### 4. Sanity Checks (Is it an Impostor?)
**This is to be done for `V2.0` of the project...**

Not all dips are planets. The most common "impostor" is an **Eclipsing Binary (EB)**—two stars orbiting each other. This is how you check for them:

* **Odd/Even Check:** You fold the data at *twice* the period you found.
    * **Planet:** A planet dip is the same depth every time. The dip at `Period 1` will look *exactly* the same as the dip at `Period 2`.
    * **EB:** When two stars orbit, the dip when Star A blocks Star B is *different* from the dip when Star B blocks Star A. One dip will be deeper than the other.
    * **The Test:** If you see a "deep-shallow-deep-shallow" pattern, you flag it as a probable EB.
* **Secondary Eclipse Check:** You look at the "phase-folded" data exactly *opposite* your main dip (at phase 0.5).
    * If you see *another*, shallower dip there, it's the "secondary eclipse" (the smaller star going *behind* the bigger one).
    * Planets are tiny and dark; they don't cause a secondary eclipse you can see. This is a dead giveaway for an EB.
