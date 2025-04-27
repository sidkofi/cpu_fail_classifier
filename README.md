# Finding CPU Fails: A Wavelet Approach

*By Sidney Taylor*

---

## What's the Problem? Why Care About CPU Failures?

Modern CPUs are incredibly complex pieces of engineering. Before they end up in our computers, they undergo rigorous testing under intense conditions (think complex simulations, high-end gaming, video rendering) to ensure they perform reliably. Manufacturers collect vast amounts of data, including operating frequency, temperature, power usage, and more, during these tests.

Most CPUs pass. But sometimes, they exhibit faulty behavior – maybe the frequency gets stuck, drops unexpectedly, or becomes unstable. Catching these failures *before* the CPU ships is crucial. It saves manufacturers money by preventing faulty products from reaching consumers, and it saves users the headache of dealing with an unreliable computer.

[Image: data_privacy.png]
*(Caption: Real-world CPU test data is often vast and proprietary, making analysis challenging.)*

However, analyzing the raw data from these tests presents challenges:

1.  **Scale:** Test runs generate enormous amounts of data (gigabytes or more), making manual inspection impossible.
2.  **Complexity:** Multiple parameters interact in non-obvious ways.
3.  **Subtlety:** Some failures might be transient or hidden within noisy data.
4.  **Privacy:** Detailed performance data is often proprietary and not publicly shared, hindering broad research.

Current methods can be slow or expensive, and sometimes failures slip through. Can we use Digital Signal Processing (DSP) to automatically and efficiently classify different types of CPU failures based on their frequency behavior?

---

## The Approach: From Raw Signals to Failure Types

Our goal is to take the raw CPU frequency data over time (a time-domain signal) and classify it into one of several categories: "None" (healthy), "Thermal Throttling", "Stuck Frequency", or "Frequency Oscillation".

We developed a workflow based on DSP principles:

[Image: workflow_diagram.png] *(Use the diagram from Slide 4)*
*(Caption: Our DSP-based workflow for CPU failure classification.)*

1.  **Generate Data:** Due to the difficulty in obtaining real-world failed-CPU data, we created a realistic simulation to generate frequency profiles mimicking normal operation and the different failure modes.
2.  **Filter Data:** Raw frequency data is inherently noisy. We apply preprocessing filters to remove irrelevant high-frequency noise while preserving the underlying trends and failure signatures.
3.  **Wavelet Analysis:** We use the Discrete Wavelet Transform (DWT) to decompose the filtered signal. Wavelets are excellent tools for analyzing signals with transient events (like failures) because they provide information about frequency content *localized in time*.
4.  **Create Features:** Instead of classifying the entire time series, we extract a small number of numerical "features" from the wavelet decomposition that specifically capture the characteristics distinguishing different failure types.
5.  **Classify Failures:** We build a simple, rule-based classifier (a decision tree) using thresholds on the extracted features to make the final classification.

---

## Understanding the Signals: What Do Failures Look Like?

Let's look at examples of the simulated frequency data (sampled at 10 Hz over 60 seconds) for each category:

* **None (No Failure):** The frequency dynamically adjusts to simulated load changes, staying within expected bounds (typically 3.7-4.6 GHz for our simulated CPU). It looks "choppy" due to normal, rapid adjustments.
    [Image: fail_mode_none.png] *(Use plot from Slide 7)*
* **Thermal Throttling:** A sharp, deep drop in frequency occurs when a simulated thermal limit is hit, followed by a recovery. Performance severely degrades during the drop.
    [Image: fail_mode_thermal.png] *(Use plot from Slide 8)*
* **Stuck Frequency:** The frequency gets locked at a constant level for a period, unable to adjust to load changes. Performance suffers if stuck low; power is wasted if stuck high.
    [Image: fail_mode_stuck.png] *(Use plot from Slide 9)*
* **Frequency Oscillation:** The frequency becomes unstable and fluctuates rapidly around a target level, indicating power delivery or control issues.
    [Image: fail_mode_oscillation.png] *(Use plot from Slide 10)*

---

## Step 1: Filtering the Noise

As mentioned, raw CPU frequency data is noisy due to rapid dynamic adjustments and measurement effects. To focus on the failure signatures, we apply **cascaded preprocessing filters**:

1.  **Moving Average Filter (Window=3):** A simple filter that slightly smooths the data by averaging each point with its immediate neighbors. It reduces minor jitter.
2.  **Butterworth Low-Pass Filter (Order=2, Cutoff=3.0 Hz):** A more advanced filter that removes frequencies above 3 Hz while preserving the slower changes associated with load variations and failure events. It has a smooth response, avoiding distortion of the frequencies we want to keep.

[Image: filter_response_ma.png] *(Use plot from Slide 13, Left)*
*(Caption: Frequency response of the simple Moving Average filter.)*

[Image: filter_response_butter.png] *(Use plot from Slide 13, Right)*
*(Caption: Frequency response of the 2nd-order Butterworth low-pass filter.)*

Filtering helps reveal the underlying structure and makes subsequent analysis more reliable.

---

## Step 2: Wavelet Analysis - Zooming In on Changes

Why use wavelets instead of something like the Fourier Transform? While Fourier analysis tells you *what* frequencies are present overall, it struggles with signals where the frequency content changes *over time* – exactly like our CPU failures, which are transient events!

Wavelets are different. They are small waves that are localized in both time and frequency. This allows them to act like a mathematical microscope, zooming in on specific time intervals to see what frequencies are active *right then*.

We chose the **Daubechies 4 ('db4') wavelet** for several reasons:
* **Compact Support:** It's localized in time, good for analyzing transient failure events.
* **Vanishing Moments (2):** Helps the analysis ignore smooth sections of the signal and focus on abrupt changes or oscillations, making failure signatures stand out in the detail coefficients.
* **Balance:** It offers a good trade-off between time and frequency localization for this type of signal.

[Image: wavelet_db4.png] *(Use plot from Slide 15)*
*(Caption: The 'db4' scaling function (phi) and wavelet function (psi).)*

We decompose the filtered signal into multiple levels using `wavedec`. The "Detail" coefficients at each level capture information within specific frequency bands. We primarily focus on **Detail Level 1 (D1)**, which contains the highest frequency information (approx. 2.5-5 Hz in our case, after filtering) where the sharp transients associated with failures are most prominent.

Here's how the D1 and D2 reconstructions look for our examples:

[Image: wavelet_details_all.png] *(Use combined plot from Slide 16)*
*(Caption: Detail Level 1 & 2 reconstructions for each failure type. Note the different amplitude scales and patterns.)*

Observe how different failure types create distinct patterns in the D1 signal: 'None' has low amplitude noise, 'Thermal' and 'Stuck' have isolated large spikes, and 'Oscillation' has sustained high-frequency activity.

---

## Step 3: Feature Engineering - Distilling Information

Analyzing the entire D1 time series (600 points) for every run is still computationally intensive. We need to distill the essential information into just a few numbers – **features** – that capture the defining characteristics of each failure type. This process, called feature engineering, is key to building an effective and efficient classifier.

[Image: chaos_to_order_visualization.png] *(Use the generated abstract graphic)*
*(Caption: Feature engineering helps bring order to complex data, moving from raw signals to concise, informative metrics through iterative refinement.)*

Based on observing the D1 plots and iterating through different ideas, we developed four key features:

1.  **Significant Peak Count (`pk_count`):** Counts how many high-amplitude spikes occur in the D1 signal. *Helps identify Frequency Oscillation (many peaks).*
    [Image: feature_pk_count_vis.png] *(Use plot from `matlab_vis_pk_count`)*
2.  **Magnitude Ratio (`mag_ratio`):** Measures how much the largest D1 peaks dominate the rest of the signal. *Helps identify Thermal Throttling (very high ratio).*
    [Image: feature_mag_ratio_vis.png] *(Use plot from `matlab_vis_mag_ratio` - the one using separated peaks)*
3.  **Separated Peak Distance (`pk_dist_sep`):** Measures the time between the two largest D1 peaks that are at least 3 seconds apart. *Intended to help identify Stuck Frequency (large separation), but wasn't used in the final logic.*
    [Image: feature_sep_dist_vis.png] *(Use plot from `matlab_vis_sep_dist`)*
4.  **Mean Absolute Deviation (`mad_val`):** Measures the average variability of the D1 signal's *magnitude*. *Excellent for identifying the 'None' state (very low MAD).*
    [Image: feature_mad_vis.png] *(Use plot from `matlab_vis_mad`)*

---

## Step 4: Classification - Making the Decision

With these features calculated for each run, we can build a simple rule-based classifier. This works like a flowchart or **decision tree**, checking the feature values against tuned thresholds in a specific order.

[Image: decision_tree.png] *(Use the generated decision tree PNG)*
*(Caption: The final decision tree logic used for classification.)*

Here's the logic flow:

1.  **Check for 'None':** Is the `mad_val` extremely low (< 0.0006) AND is the `pk_count` high (> 40)? If yes -> Classify as **None**.
2.  **Check for 'Frequency Oscillation':** If not 'None', is `pk_count` high (> 6) AND `mag_ratio` low (< 20.0) AND `mad_val` slightly elevated (> 0.00065)? If yes -> Classify as **Frequency Oscillation**.
3.  **Check for 'Thermal Throttling':** If not 'None' or 'Freq Osc', is `mag_ratio` high (> 39) AND `mad_val` also high (> 0.002)? If yes -> Classify as **Thermal Throttling**.
4.  **Default to 'Stuck Frequency':** If none of the above conditions are met, classify as **Stuck Frequency**.

This sequence leverages the most reliable features first (like MAD for 'None') and uses combinations of features to separate the trickier cases.

---

## Results and Conclusion

We tested this classification logic on 1000 simulated runs (250 of each type). The results were very promising:

[Image: confusion_matrix.png] *(Use the final confusion matrix plot, e.g., from Slide 28)*
*(Caption: Confusion matrix showing classification performance on 1000 simulated runs. Diagonal elements represent correct classifications.)*

* **Overall Accuracy:** Achieved **96.6%** accuracy across all failure types.
* **Key Successes:** The MAD feature perfectly identified all 'None' runs. 'Thermal Throttling' was also classified with 100% accuracy using the Magnitude Ratio and MAD combination.
* **Minor Confusion:** Most errors involved misclassifying some 'Frequency Oscillation' as 'Stuck Frequency' (7 runs) or vice-versa (14 runs), indicating some overlap in their feature characteristics based on the current thresholds.

**Takeaways:**

* Applying DSP techniques (filtering, wavelets) effectively isolates failure signatures in noisy CPU frequency data.
* Carefully engineered features (MAD, Peak Count, Magnitude Ratio) derived from wavelet details can capture the defining characteristics of different failure modes.
* A simple, interpretable rule-based classifier using these features can achieve high classification accuracy on simulated data.
* The MAD feature was surprisingly powerful for identifying normal operation.

**Limitations & Future Work:**

* This analysis used **simulated data**. Validation on diverse real-world CPU telemetry is needed.
* The thresholds were **tuned manually** based on observed distributions; automated optimization could be explored.
* Only **frequency data** was used. Incorporating other metrics (temperature, voltage, utilization) could improve accuracy and potentially identify more failure types.
* Future work could involve **predicting the *time*** of failure onset or using these features as input to more complex **machine learning models**.

---

## Supplemental Material (Code & Data)

All MATLAB code used for data generation, filtering, wavelet analysis, feature extraction, classification, and visualization, along with the generated benchmark data file, can be found in this GitHub repository:

* **Code:** [`code/`](./code/) directory
* **Data:** [`data/`](./data/) directory (`benchmark_data.csv`)

*(Make sure these links point correctly within your repository structure)*

---
