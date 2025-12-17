# Adult-Spine-Models
Which single neuron models best predict the somatic tuning and response distribution of a cell from the spine responses


## Installation instructions

1. Adult-Spine-Models requires a recent Python distribution and conda.

2. Open a terminal and type:

```bash
cd path/to/your/code/directory #or similar as relevant for your machine
git clone git@github.com:GreggHeller1/adult-spine-models.git
cd adult-spine-models
conda env create -f environment.yml
conda activate adult-spine-models
```

3. The project should now be installed and the adult-spine-models environment should still be activated.

4. To run the algorithm, copy the default configuration file and set the paths in `src/config.py` to point to your data directory containing soma and spine .mat files:

```bash
cd path/to/adult-spine-models
cp src/default_config.py src/config.py
# Edit src/config.py to set data_path and other paths as desired
python scripts/main.py
```



## Usage instructions
### Input Files:
The code expects to find MATLAB (.mat) files containing experimental data:

1. **Soma data files**: .mat files whose filenames contain the string "ASC" or "BM" that contain somatic calcium imaging traces and metadata. These files should be located in a directory specified by `data_path` in `src/config.py`.

2. **Spine data files**: Corresponding .mat files containing spine calcium imaging traces. The code automatically finds spine data files by replacing "Soma" with "Spine" in the soma data path (e.g., if soma data is in "2023-08 Soma Data/ASC26.mat", it will look for "2023-08 Spine Data/ASC26.mat").

3. **Model specifications**: A CSV file (`model_specifications.csv`) that defines which models to run. Each row specifies a model configuration including:
   - `trials_included`: Which trials to use (e.g., 'all_trials', 'no_bap_trials', 'baps_trials_only')
   - `weighted_by`: How to weight spine contributions (e.g., 'democratic_weights', 'spine_size', 'distance_from_soma')
   - `spines_included`: Which spines to include (e.g., 'all_spines', 'only', 'top_50_percentile')
   - `inclusion_based_on`: Parameter used for spine inclusion filtering
   - `integration`: Integration function (e.g., 'summation', 'two_layer_by_fov_1')
   - `somatic_function`: Somatic transformation function (e.g., 'identity', 'sigmoid', 'unrectified_linear')
   - `run_model`: Boolean flag indicating whether to run this model

The model specifications file path is set in `src/config.py` via the `csv_config_path` variable.

### Usage steps
1. Copy `src/default_config.py` to `src/config.py` (First time only)
2. Edit `src/config.py` to set the following paths:
   - `data_path`: Path to directory containing soma .mat files (or path to a single .mat file)
   - `collect_summary_at_path`: Path where output files will be saved
   - `csv_config_path`: Path to `model_specifications.csv` file
3. Optionally adjust other parameters in `src/config.py` such as:
   - `simulated_trials_per_stim`: Number of simulated trials per stimulus (default: 1000)
   - `walk_dirs`: Whether to walk through subdirectories (default: True)
   - `re_run`: Whether to re-run on already processed files (default: True)
   - `save_traces`: Whether to save trace images (default: True)
4. In a terminal with the conda environment active, run:
```bash
python scripts/main.py
```

Alternatively, you can run the notebook demo in Google Colab by opening `scripts/notebook.ipynb`.

### Output Files:
The code generates several CSV files containing the simulation results and model comparisons:

**Aggregate summary file** (saved at `collect_summary_at_path`):
1. `single_neuron_simulation_scores_{simulated_trials_per_stim}_{datestring}.csv` - CSV containing model performance scores and metadata for all processed cells. Each row represents one model run on one cell. Columns include:
   - `experiment_id`: Identifier for the cell/experiment
   - `full_model_name`: Complete model specification name
   - `model_nickname`: Short model name
   - `soma_responsive_status`: Whether the soma was responsive
   - `trials_included`: Which trials were used
   - `weighted_by`: Weighting scheme used
   - `spines_included`: Spine inclusion criteria
   - `integration_function`: Integration method
   - `somatic_function`: Somatic transformation
   - Various similarity scores comparing model output to actual soma traces (e.g., `model_soma_similarity_score`, `model_correlation_to_soma_r`)

**Per-experiment files** (saved in subdirectories named by `experiment_id` at `collect_summary_at_path`):
1. `{experiment_id}_simulation_mean_stim_response.csv` - CSV containing normalized mean response tuning curves. Each row represents one model (or 'soma' for actual measurements). Columns are stimulus directions, with values representing the normalized mean response amplitude for each direction.

2. `{experiment_id}_temporal_loglikelihood.csv` - CSV containing temporal log-likelihood scores for each model. Each row represents one model. Columns correspond to time samples, with values representing the log-likelihood of the model prediction at each time point.

**Error tracking files** (saved at `collect_summary_at_path`):
1. `failed_dirs_list_{datestring}.json` - JSON file containing a list of directories/files that failed to process
2. `failed_dirs_errors_{datestring}.json` - JSON file containing error messages for failed processing attempts

Please see `scripts/demo_data` for example input files.

## Running the Demo/End-to-End Test

To run a demo with the provided example data:

1. Ensure the conda environment is activated:
```bash
conda activate adult-spine-models
```

2. Set up the configuration to point to the demo data:
   - Edit `src/config.py` to set `data_path` to point to `scripts/demo_data/ASC26_cell_3_soma.mat` (or the full path to this file)
   - Ensure `collect_summary_at_path` points to a directory where you want the output saved
   - Ensure `csv_config_path` points to `model_specifications.csv` in the project root

3. Run the main script:
```bash
python scripts/main.py
```

Alternatively, you can run the interactive notebook demo:
- Open `scripts/notebook.ipynb` in Jupyter or Google Colab
- The notebook can be run directly in Colab using the provided Colab badge link

For unit tests of individual functions, run:
```bash
pytest src/tests.py
```

Code was written by Gregg Heller with advice from Kyle Jenks
