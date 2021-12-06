from util_func import params_search, params_search_summary, sofie2_pipe, submit
from pathlib import Path

# Set the task to be performed
PARAMS_SEARCH = True
PARAMS_SEARCH_SUMMARY = True
SOFIA2_PIPE = True
SUBMIT = True

# Parameters for task `parameters search` and `parameters search summary`.
search_data_file = '../sdc2_data_v2/development/sky_dev_v2.fits'
truth_cat = '../sdc2_data_v2/development/sky_dev_truthcat_v2.txt'
params_file = './sofia_ska.par'

params_search_output_prefix = 'sofia_test_output'
params_search_output_dir = 'params_exploration_1'

search_sub_region = ['0, 300, 0, 300, 0, 3000'] # format: ['xs1, xe1, ys1, ye1, zs1, ze1', 'xs2, xe2, ys2, ye2, zs2, ze2']
# search_sub_region = [' ']

# Parameters for task `SoFiA2-pipe` and `submit`
pipe_data_file = '../sdc2_data_v2/development/sky_dev_v2.fits'
pipe_output_dir = './sofia2_pipe'
pipe_output_prefix = 'sofia_test_output'
pipe_input_params = params_file
pipe_nsub = (2, 2, 2)
pipe_overlap = (10, 10, 10)

# The number of threads
n_workers = 8

#######################################################
#######################################################
if PARAMS_SEARCH:
        
        # The values of each key should be iterable.
        params = {
        "input.data": [search_data_file],
        "input.region": search_sub_region,

        "scaleNoise.enable": ["true"],
        "scaleNoise.mode": ["local"],
        "scaleNoise.fluxRange": ["negative"],
        "scaleNoise.windowXY": [35],
        "scaleNoise.windowZ": [35],
        "scfind.kernelsXY": [[0, 3, 7], [0, 5, 11], [0, 3, 9, 15]],
        "scfind.kernelsZ": [[0, 3, 7, 15, 31], [0, 5, 11, 23, 41]],
        "scfind.threshold": [3.0, 4.0],
        "scfind.replacement": [2.0],
        "scfind.statistic": ["mad"],
        "scfind.fluxRange": ["negative"],
        "linker.minSizeXY": [5],
        "linker.minSizeZ": [25],
        # "linker.positivity": [],
        "reliability.threshold": [0.7, 0.9],
        "reliability.scaleKernel": [0.4],
        "reliability.minSNR": [3.0, 5.0],
        "dilation.enable": ["false"],
        "dilation.iterationsXY": [10],
        "dilation.iterationsZ": [5],
        "dilation.threshold": [0.001],
        "output.directory": [params_search_output_dir],
        "output.filename": [params_search_output_prefix],
        }
        
        params_search(params_list=params,
                      data_file=search_data_file,
                      truth_cat=truth_cat,
                      params_file=params_file,
                      output_dir=params_search_output_dir,
                      output_prefix=params_search_output_prefix,
                      n_workers=n_workers)
    

if PARAMS_SEARCH_SUMMARY:
    
    # ['matched', 'found', 'matched_rate', 'score', 'score_corrected', 'score/source']
    sorted_cond = {'found': 1000,
                   'matched_rate': 500, 
                   'score_corrected': 50}
    
    if Path(params_search_output_dir).exists():
        optimal_params_path = params_search_summary(param_dir=params_search_output_dir,
                                                    data_file=search_data_file,
                                                    sorted_by=sorted_cond,
                                                    return_params_path=True)
        print(optimal_params_path)
        
        
if SOFIA2_PIPE:
    try:
        pipe_input_params = optimal_params_path
    except:
        print("=== Be carefull...you might set optimal parameters ===")
    
    print("\n")
    print(f'=== The parameter file to process the data is {pipe_input_params} ===')
    sofie2_pipe(data_file=pipe_data_file,
                sourcename=pipe_input_params,
                output_dir=pipe_output_dir,
                output_prefix=pipe_output_prefix,
                nsub=pipe_nsub,
                overlap=pipe_overlap,
                n_workers=n_workers)


if SUBMIT:
    submit(res_dir=pipe_output_dir,
           data_file=pipe_data_file,
           output_prefix=pipe_output_prefix)
