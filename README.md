# EEGDojo🥋  
**Your Personal EEG Processing Toolkit**

![EEGDojo Cover](cover.png)

Welcome to **EEGDojo**, a powerful and versatile toolbox for EEG signal preprocessing and artifact detection. Whether you're a beginner or an expert in EEG analysis, EEGDojo provides clean, structured, and efficient workflows to detect bad channels, run ICA, and clean EEG data like a pro.

---

## Features 🧠
- **Modular and Customizable:** Run modular bad channel detection methods (e.g., Kurtosis, Spectrum, ICLabel).
- **ICA Automation:** Integrates ICLabel, FASTER, and spatial filtering methods.
- **User-Friendly:** Output reports and logs for transparency and reproducibility.
- **Visualization-Ready:** Automatically saves bad channel/component properties for review.
- **Flexible Integration:** Suitable for EEGLAB pipelines and EEG research workflows.

---

## Table of Contents  
1. [Installation](#installation)  
2. [Usage](#usage)  
3. [Example Workflow](#example-workflow)  
4. [Contributing](#contributing)  
5. [License](#license)

## Installation ⚙️  
Clone the repository:  
```bash
git clone https://github.com/yourusername/EEGDojo.git
```

## Usage 🚀  
Example pipeline for bad channel detection and ICA artifact removal:  

```matlab
% Load EEG data
EEG = pop_loadset('filename', 'example_data.set');

% Set parameters
params = load('params_example.mat');

% Run bad channel detection
[EEG, BadChan] = detect_badchannels(EEG, params);

% Run ICA and remove bad components
[EEG, BadICs] = EEGdojo_remove_badICs(EEG, params, 1);

% Save cleaned data
pop_saveset(EEG, 'filename', 'example_data_cleaned.set');
```

## Example Workflow 📊
Below is a sample workflow for EEG data preprocessing:

Detect bad channels using modular methods (e.g., Spectrum, Kurtosis, CleanRawData).
Run ICA to identify and remove bad components.
Save results and export logs/visualizations for further analysis.

## Contributing 🤝
We welcome contributions to EEGDojo! Whether it's fixing a bug, adding new features, or improving documentation, feel free to fork this repository and submit a pull request.


## License 📄
This project is licensed under the MIT License.

Enjoy EEG processing like a ninja in the EEG dojo! 🥋

