import os
import sys
import numpy as np
import pandas as pd
from itertools import product
from pathos.multiprocessing import Pool
from astropy.io import fits
from astropy import units as u
import subprocess
import time
from pathlib import Path


SOFIA_COLUMNS = "name id x y z x_min x_max y_min y_max z_min z_max n_pix f_min f_max f_sum rel flag rms w20 w50 wm50 ell_maj ell_min ell_pa ell3s_maj ell3s_min ell3s_pa kin_pa err_x err_y err_z err_f_sum ra dec freq x_peak y_peak z_peak ra_peak dec_peak freq_peak".split()


def read_true_cat(filename):

    columns = "id ra dec hi_size line_flux_integral central_freq pa i w20"
    col_name = columns.split(" ")

    data = np.loadtxt(filename, skiprows=1)
    data_df = pd.DataFrame(data, columns=col_name)
    data_df = data_df.astype({'id': int})

    return data_df


def read_sofia_cat(filename):

    sofia_df = pd.read_csv(filename, skiprows=13,
                           delim_whitespace=True, header=None)
    sofia_df[0] = sofia_df[0] + ' ' + sofia_df[1]
    sofia_df.drop(1, axis=1, inplace=True)
    sofia_df.columns = SOFIA_COLUMNS
    sofia_df.astype({'id': int})

    return sofia_df


def writeParfile(filename, sourcefile, pars=None):

    with open(sourcefile, 'r') as f:
        lines = f.readlines()

    # Remove leading/trailing whitespace and empty lines
    lines = [line.strip() for line in lines]

    # Remove comments
    lines = [line for line in lines if len(line) > 0 and line[0] != "#"]
    pardict = {}

    for line in lines:
        parameter, value = line.split("=", 1)
        parameter = parameter.strip()
        value = value.split("#")[0].strip()
        module, parname = tuple(parameter.split(".", 1))

        if module not in pardict.keys():
            pardict[module] = {}

        pardict[module][parname] = value

    if pars is not None:
        if isinstance(pars, dict):
            for (par, value) in pars.items():
                module, parname = tuple(par.split(".", 1))
                pardict[module][parname] = value
        else:
            for par in pars:
                p, v = par.split("=", 1)
                p = p.strip()
                v = v.split("#")[0].strip()
                module, parname = tuple(p.split(".", 1))
                pardict[module][parname] = v

    with open(filename, 'w') as f:
        for task in pardict.keys():
            f.writelines("# " + task + "\n")
            for (p, v) in pardict[task].items():
                f.writelines(task + '.' + p + " = " + str(v) + "\n")
            f.writelines("\n")


def sofia_paras_trans(source_sofia):

    pixel_size = 2.8  # arcsec / pixel for sdc2 dataset
    hi_size = source_sofia.ell_maj * pixel_size

    rest_freq = 1420.40575177e6  # Hz
    c = 299792458.0  # m/s
    # sofia - parameter.physical = true. w20 is in the unit of Hz.
    w20 = source_sofia.w20 * c / rest_freq * 1e-3

    i = np.arccos(source_sofia.ell_min / source_sofia.ell_maj) * 180.0 / np.pi

    source_sofia_2 = pd.Series({'ra': source_sofia.ra,
                                'dec': source_sofia.dec,
                                'hi_size': hi_size,
                                'line_flux_integral': source_sofia.f_sum,
                                'central_freq': source_sofia.freq,
                                'pa': source_sofia.kin_pa,
                                'i': i,
                                'w20': w20})

    return source_sofia_2


def weight_score(source_true, source_sub):

    threshold_pos = 0.3
    threshold_hi_size = 0.3
    threshold_flux = 0.1
    threshold_freq = 0.3
    threshold_pa = 10.0
    threshold_inc_angle = 10.0
    threshold_line_width = 0.3

    if source_true.hi_size <= 0.0:
        err_pos = 1e5
        err_hi_size = 1e5
    else:
        err_pos = ((source_sub.ra - source_true.ra)**2 + (source_sub.dec -
                                                          source_true.dec)**2)**0.5 / (source_true.hi_size / 3600)
        err_hi_size = np.abs(source_sub.hi_size -
                             source_true.hi_size) / source_true.hi_size
    if source_true.line_flux_integral <= 0.0:
        err_flux = 1e5
    else:
        err_flux = np.abs(source_sub.line_flux_integral -
                          source_true.line_flux_integral) / source_true.line_flux_integral
    err_freq = np.abs(source_sub.central_freq -
                      source_true.central_freq) / source_true.central_freq
    if source_sub.pa <= 0.0:
        err_pa = 1e5
    else:
        err_pa = abs(np.arctan2(np.sin(abs(source_sub.pa - source_true.pa)),
                                np.cos(abs(source_sub.pa - source_true.pa))))
    err_inc_angle = np.abs(source_sub.i - source_true.i)
    if source_true.w20 <= 0.0:
        err_line_width = 1e5
    else:
        err_line_width = np.abs(
            source_sub.w20 - source_true.w20) / source_true.w20

    w_pos = min(1.0, threshold_pos / err_pos)
    w_hi_size = min(1.0, threshold_hi_size / err_hi_size)
    w_flux = min(1.0, threshold_flux / err_flux)
    w_freq = min(1.0, threshold_freq / err_freq)
    w_pa = min(1.0, threshold_pa / err_pa)
    w_inc_angle = min(1.0, threshold_inc_angle / err_inc_angle)
    w_line_width = min(1.0, threshold_line_width / err_line_width)

    w_total = (w_pos + w_hi_size + w_flux + w_freq +
               w_pa + w_inc_angle + w_line_width) / 7

    return w_total


def scoring(filename, true_cat, outdir, idx):

    true_cat_df = read_true_cat(true_cat)
    sofia_cat_df = read_sofia_cat(filename)

    pixel_size = 0.000777777777778  # 2.8 arcsec
    beam_maj = 0.00194444449153  # deg...7.0 arcsec
    beam_min = 0.00194444449153  # deg...7.0 arcsec
    beam_area = np.pi * beam_maj * beam_min / (4.0 * np.log(2.0) * pixel_size * pixel_size)

    matched_sources_idx = set()
    matched_sources_score = {}

    for idx_true, source_true in true_cat_df.iterrows():

        for idx_sofia, source_sofia in sofia_cat_df.iterrows():
            # spatial_range = source_sofia.ell_maj * pixel_size * np.sqrt(beam_area) # beam convolution? Not sure...
            # spectral_range = source_sofia.w20 * c / rest_freq * 1e-3 # whether this also should be beam-convoluted?
            spatial_range = source_sofia.ell_maj * pixel_size

            if ((source_sofia.ra-spatial_range <= source_true.ra <= source_sofia.ra+spatial_range) and
                (source_sofia.dec-spatial_range <= source_true.dec <= source_sofia.dec+spatial_range) and
                    (source_sofia.freq-source_sofia.w20 <= source_true.central_freq <= source_sofia.freq+source_sofia.w20)):

                # matched_sources.append(source_sofia)
                matched_sources_idx.add(idx_sofia)
                source_sofia_2 = sofia_paras_trans(source_sofia)
                score_i = weight_score(source_true, source_sofia_2)

                if idx_sofia not in matched_sources_score.keys():
                    matched_sources_score[idx_sofia] = [score_i]
                else:
                    matched_sources_score[idx_sofia].append(score_i)

                print(f"idx_true = {idx_true}, idx_sofia = {idx_sofia}, score = {score_i}")

    total_score = sum([max(i) for i in matched_sources_score.values()])

    # print(f"{matched_sources_idx = }")
    matched_sources = pd.DataFrame(
        sofia_cat_df.iloc[list(matched_sources_idx)])
    matched_sources = matched_sources.astype({'id': int})

    matched_sources.to_csv(outdir + "matched_source.txt", sep=' ', index=False)

    print(f"My index is {idx}, the matched_sources_idx = {matched_sources_idx}...There are {len(matched_sources)} / {len(sofia_cat_df)} sources matched.")
    print(f"The total_score = {total_score}...")
    with open(outdir + str(len(matched_sources)) + "_" + str(len(sofia_cat_df)) + '_' + str(total_score) + "_sources", "w") as f:
        pass


def parameters_range(kw_params):
    keys = list(kw_params.keys())
    params_sets = list(product(*kw_params.values()))
    return (keys, params_sets)


def cube_split(x, y, z, nsubx=1, nsuby=1, nsubz=1, overlapx=0, overlapy=0, overlapz=0):

    def split_m(n, m):
        base = n // m
        rem = n % m

        part = base * np.ones(m, dtype=int) + (np.arange(m) < rem).astype(int)
        bound = np.cumsum(np.insert(part, 0, 0))

        return np.array([part, bound[:m], bound[1: (m + 1)]])

    n_subcubes = nsubx * nsuby * nsubz

    n_x, s_x, e_x = split_m(x, nsubx)
    n_y, s_y, e_y = split_m(y, nsuby)
    n_z, s_z, e_z = split_m(z, nsubz)

    if ((overlapx >= n_x[-1]) or (overlapy >= n_y[-1]) or (overlapz >= n_z[-1])):
        raise ValueError(f"The value of overlapx = {overlapx}, overlapy = {overlapy}, overlapz = {overlapz} should not be larger than the size of subcube [{n_x[-1]}, {n_y[-1]}, {n_z[-1]}].")

    s_x_new = s_x
    e_x_new = e_x
    s_y_new = s_y
    e_y_new = e_y
    s_z_new = s_z
    e_z_new = e_z

    s_x_new[1:] = s_x_new[1:] - overlapx // 2
    e_x_new[:-1] = e_x_new[:-1] + overlapx // 2
    s_y_new[1:] = s_y_new[1:] - overlapy // 2
    e_y_new[:-1] = e_y_new[:-1] + overlapy // 2
    s_z_new[1:] = s_z_new[1:] - overlapz // 2
    e_z_new[:-1] = e_z_new[:-1] + overlapz // 2

    return (np.array([s_x_new, e_x_new]), np.array([s_y_new, e_y_new]), np.array([s_z_new, e_z_new]))


def params_search(params_list,
                  data_file,
                  truth_cat,
                  params_file,
                  output_dir="./sofia_params_explorasion/",
                  output_prefix = 'sofia_test_output',
                  n_workers=None):
    """Set a list of values for control parameters of SoFiA, and evaluate the result of each combination with a less precise scoring scheme.

    Parameters
    ----------
    params_list : dict
        Dictionary with parameters names as keys and lists of parameters settings as values.
    data_file : str
        The path of development datacube whose true catalog is known.
    truth_cat : str
        The path of the true catalog of *data_file*.
    params_file : str
        The templete SoFiA-2 parameter file.
    output_dir : str, optional
        The path to save results.
    n_workers : int, optional
        The number of threads to be used.
    """

    params = params_list
    keys, params_set = parameters_range(params)

    dataname = os.path.split(data_file)[1].split('.')[0]

    true_cat = truth_cat
    sourcename = params_file

    def run(idx):
        # print(f"My index is {idx} / {len(params_set)}")
        params = dict(zip(keys, params_set[idx]))

        filename_idx = "/sofia_" + dataname + "_params_" + str(idx)

        outdir = output_dir + filename_idx + '/'
        par_filename = outdir + filename_idx + ".par"

        sofia_cat_out = outdir + output_prefix +"_cat.txt"

        if os.path.exists(outdir):
            pass
        else:
            os.makedirs(outdir)

        params["output.directory"] = outdir
        writeParfile(par_filename, sourcename, pars=params)

        log_file = open(outdir + "sofia_" + dataname +
                        "_params_" + str(idx) + ".out", "w")
        subprocess.run(["sofia", par_filename], stdout=log_file)

        if os.path.exists(sofia_cat_out):
            scoring(sofia_cat_out, true_cat, outdir, idx)

    if n_workers is None:
        n_workers = os.cpu_count()

    p = Pool(n_workers)
    p.map_async(run, range(len(params_set)))
    p.close()
    p.join()

    print("Finish parameters search...")


def params_search_summary(param_dir,
                          data_file,
                          sorted_by='matched_rate',
                          return_params_path=False):
    """Summarize and rank the result of `params_search`.

    Parameters
    ----------
    param_dir : str
        The path where function `params_search` save its output. Might be same with *output_dir* in `params_search`.
    sorted_by : str or dict, optional
        The column to sort by, it can be a string or a ordered dictionary. The dictionary can contain a integer as value of each key, and select only part of sorted result for next sorting. The columns *matched*, *found*, *score* are got from existed files, and columns *matched_rate*, *score_corrected* and *score/source* will be calculated additionally in this function. By default 'matched_rate'.
    return_params_path : bool, optional
        Whether return the path of 'optimal' parameter according to your sorting.
    """

    p = Path(param_dir)
    dataname = os.path.split(data_file)[1].split('.')[0]

    result = pd.DataFrame(
        columns=["param_idx", "matched", "found", "matched_rate", "score"])

    for i in p.glob("_".join(["sofia", dataname, "params", "*"])):
        param_idx = int(i.name.split('_')[-1])
        res = list(i.glob('*_sources'))
        if len(res) != 0:
            res_i = res[0].name.split('_')
            num_matched, num_found, score = int(
                res_i[0]), int(res_i[1]), float(res_i[2])
            if num_found == 0:
                matched_rate = 0
            else:
                matched_rate = num_matched / num_found
            result = result.append({"param_idx": param_idx, "matched": num_matched, "found": num_found,
                                    "matched_rate": matched_rate, "score": score}, ignore_index=True)
        else:
            result = result.append({"param_idx": param_idx, "matched": None,
                                    "found": None, "matched_rate": None, 'score': None}, ignore_index=True)

    result_dropna = result.dropna().copy()

    result_dropna['score_corrected'] = result_dropna.loc[:, 'score'] - (result_dropna.loc[:, 'found'] - result_dropna.loc[:, 'matched'])
    result_dropna['score/source'] = result_dropna.loc[:, 'score_corrected'] / result_dropna.loc[:, 'matched']

    # print(result_dropna)
    result_dropna.replace([np.inf, -np.inf], np.nan, inplace=True)
    result_dropna.dropna(subset=['score/source'], inplace=True)

    if isinstance(sorted_by, dict):
        for k, v in sorted_by.items():
            if v is None:
                v = len(result_dropna)
            # print(result_dropna)
            result_dropna.sort_values(k, inplace=True)
            result_dropna = result_dropna.tail(v)

        print(result_dropna)

    elif isinstance(sorted_by, str):
        print(result_dropna.sort_values(sorted_by, inplace=True).tail(50))

    if return_params_path:
        return_idx = str(int(result_dropna.iloc[-1].param_idx))
        filename_idx = "sofia_" + dataname + "_params_" + return_idx
        return p.joinpath(filename_idx, filename_idx + '.par')


def sofie2_pipe(data_file,
                sourcename,
                output_dir="./sofia2_pipe",
                output_prefix="sofia_test_output",
                nsub=(2, 2, 2),
                overlap=(10, 10, 10),
                ignore_contaminated_spectral=False,
                n_workers=None):
    """Split the large datacube and process each subcube with the *optimal* parameters.

    Parameters
    ----------
    data_file : str
        The path of data cube.
    sourcename : str
        The path of parameter file of SoFiA-2.
    output_dir : str, optional
        The path to save results.
    nsub : tuple or list, optional
        The number of parts along each axis (x, y, z) that will be split into.
    overlap : int, tuple or list, optional
        The pixel number that each subcube shares with its neighbour subcubes.
    ignore_contaminated_spectral : bool, optional
        Whether discard the first few frequencies that are noise-high seen in the SDC2 data.
    n_workers : int, optional
        The number of threads to be used.
    """

    dataname = os.path.split(data_file)[1].split('.')[0]
    print(dataname)

    cube = fits.open(data_file)[0].data
    z_size, y_size, x_size = cube.shape

    nsubx, nsuby, nsubz = nsub
    overlap_x, overlap_y, overlap_z = overlap

    (x_start, x_end), (y_start, y_end), (z_start, z_end) = cube_split(x_size, y_size, z_size,
                                                                      nsubx=nsubx, nsuby=nsuby, nsubz=nsubz,
                                                                      overlapx=overlap_x, overlapy=overlap_y, overlapz=overlap_z)

    cube_idx = list(
        enumerate(product(range(nsubz), range(nsuby), range(nsubx))))

    def run(idx):

        i, (zi, yi, xi) = idx

        if ignore_contaminated_spectral and z_start[zi] == 0:
            sub_region = [x_start[xi], x_end[xi],
                          y_start[yi], y_end[yi], 150, z_end[zi]]
        else:
            sub_region = [x_start[xi], x_end[xi], y_start[yi],
                          y_end[yi], z_start[zi], z_end[zi]]

        sub_region = ', '.join(str(e) for e in sub_region)
        print(f"My region is {sub_region}.")

        filename_idx = "/sofia_" + dataname + "_sub_" + str(i)
        outdir = output_dir + filename_idx + '/'

        par_filename = outdir + filename_idx + '.par'

        if os.path.exists(outdir + output_prefix + "_cat.txt"):
            print(f"The file ||{outdir}/{output_prefix}_cat.txt|| exists...pass...")
        else:
            if not os.path.exists(outdir):
                os.makedirs(outdir)

            pars = {
                "pipeline.threads": 1,
                "input.data": data_file,
                "input.region": sub_region,
                "output.directory": outdir,
                "output.filename": output_prefix,
            }

            writeParfile(par_filename, sourcename, pars=pars)

            log_file = open(outdir + "/sofia_" + dataname +
                            "_sub_" + str(i) + ".out", "w")

            print("running..." + str(i))
            subprocess.run(["sofia", par_filename], stdout=log_file)

    st = time.time()

    if n_workers is None:
        n_workers = os.cpu_count()

    p = Pool(n_workers)
    p.map_async(run, cube_idx)
    p.close()
    p.join()

    print(f"sofie2_pipe - Spend {time.time() - st:.2f}s")


def submit(res_dir,
           data_file,
           output_prefix="sofia_test_output"):
    """Format the output of `sofie2_pipe`.
    """

    p = Path(res_dir)
    dataname = os.path.split(data_file)[1].split('.')[0]
    result_dirs = p.glob("sofia_" + dataname + "_sub_*")

    df = pd.DataFrame()
    for result_dir_i in result_dirs:
        result_catelog_i = result_dir_i.joinpath(output_prefix + "_cat.txt")
        if result_catelog_i.exists():
            sofia_df_i = read_sofia_cat(result_catelog_i)
        else:
            continue
        df = pd.concat([df, sofia_df_i], ignore_index=True)

    # ell_maj(pix -> arcsec)
    pixel_size = 2.8  # arcsec
    df['ell_maj_arcsec'] = df['ell_maj'] * pixel_size

    # Line width (Hz -> km/s) unit conversion.
    rest_freq = 1420.40575177E6  # Hz
    c = 299792458.0  # m/s
    df['w20_velocity'] = df['w20'] * c / rest_freq * 1e-3

    # Inclination angle (degrees) calculation.
    df['i'] = np.arccos(df['ell_min'] / df['ell_maj']) * 180.0 / np.pi

    # kinematic position angle.
    # There might be some unsolved cases that `kin_pa` = -1, and drop these.
    df2 = df.loc[df['kin_pa'] >= -1.0].copy()

    df2['id_new'] = range(len(df2))

    sub_columns = ['id', 'ra', 'dec', 'hi_size',
                   'line_flux_integral', 'central_freq', 'pa', 'i', 'w20']
    
    df_submit = df2[['id_new', 'ra', 'dec', 'ell_maj_arcsec',
                     'f_sum', 'freq', 'kin_pa', 'i', 'w20_velocity']]

    df_submit.columns = sub_columns

    df_submit.to_csv('./result_' + dataname + '.txt', sep=' ', index=False)
