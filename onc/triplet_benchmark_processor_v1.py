"""
This script if for ingesting ECO Triplet classic benchmark files taken on
the BC Ferries. It is intended to be run at the command line and requires
user inputs with appropriate metadata in order to create a processed file. The output
filetype is netCDf4, which contains raw counts, mean values, and median values, along
with key information in the root and variable level attributes.

The following are required:
1 - The factory serial number.
2 - The ONC Device ID.
3 - The date of calibration.
4 - The calibration type.
5 - The calibration solution.
6 - The attempt number of the calibration.
7 - The selection of the filepath to the benchmark file.
8 - The selection of the filepath to the device file.
"""


# pip install netcdf4
# pip install pandas
# pip install xarray

from datetime import datetime
import os
from os import PathLike
import pandas as pd
import re
from typing import NamedTuple
import xarray as xr
import tkinter as tk
from tkinter import filedialog
tk.Tk().withdraw()

SCRIPT_VER = "0.0.1"

NUM_PAT = r"[+-]?[0-9]*[.]?[0-9]+" # REGEX for any number.
SPACE_PAT = r"\s+" # Regex for any number of spaces.

def main():
    print('\n')
    print(f'UWS Triplet Benchmark File Processor {SCRIPT_VER}')
    print('--------------------------------------------------')
    device_sn = input('Enter The Factory Serial Number: ')
    device_sn = device_sn.replace('\n', '')
    if 'BB2FL-' not in device_sn.upper():
        device_sn = 'BB2FL-' + device_sn.upper()

    device_id = input('Enter The ONC Device ID: ')
    device_id = device_id.replace('\n', '')
    if 'DI' not in device_id.upper():
        device_id = 'DI' + device_id.upper()

    cal_date = input('Enter the Calibration Date (YYYY-MM-DD): ')

    cal_type = input('Enter the Calibration Type [PRE/POST]: ')
    cal_type = cal_type.replace('\n', '')

    cal_sol = input('Enter the Calibration Solution [DC,SZ]: ')
    cal_sol = cal_sol.replace('\n', '')

    attempt_number = input('Enter the Cal Attempt Number: ')
    attempt_number = int(attempt_number.replace('\n', ''))

    save_fn = f"{device_id}_{device_sn}_{cal_date}_{cal_type}_{cal_sol}_{attempt_number}"
    save_fn = save_fn.replace('\n', '')
    save_fn = save_fn.upper()
    save_fn = '.'.join((save_fn,'nc'))

    print('Use the pop-up file explorer select a benchmark file.')
    bfp = filedialog.askopenfilename(title='Select Triplet Benchmark File',
                                     initialdir=os.getcwd(),
                                     filetypes=(('Raw Files', '*.raw'),
                                                ('Text Files', '*.txt'),
                                                ('All Files', '*.*')))

    print('Use the pop-up file explorer select a device file.')
    dfp = filedialog.askopenfilename(title='Select Triplet Device File',
                                     initialdir=os.getcwd(),
                                     filetypes=(('Device Files', '*.dev'),
                                                ('All Files', '*.*')))

    bm = import_benchmark(bfp)
    dev = parse_dev(dfp)

    # Assign data based on dev file.
    if 'NTU' in dev.keys() and 700 in bm.wavelength:
        bm['raw_turbidity'] = bm.raw_signal.sel(wavelength=700)
        bm['turbidity'] = raw2val(bm.raw_turbidity, dev['NTU']['scale_factor'],
                                  dev['NTU']['offset'])
        for var in ['raw_turbidity', 'turbidity']:
            bm[var].attrs['wavelength'] = 700

    if 'CHL' in dev.keys() and 695 in bm.wavelength:
        bm['raw_chlorophyll'] = bm.raw_signal.sel(wavelength=695)
        bm['chlorophyll'] = raw2val(bm.raw_chlorophyll, dev['CHL']['scale_factor'],
                                    dev['CHL']['offset'])
        for var in ['raw_chlorophyll', 'chlorophyll']:
            bm[var].attrs['wavelength'] = 695
    if 'CDOM' in dev.keys() and 460 in bm.wavelength:
        bm['raw_cdom'] = bm.raw_signal.sel(wavelength=460)
        bm['cdom'] = raw2val(bm.raw_cdom, dev['CDOM']['scale_factor'],
                             dev['CDOM']['offset'])
        for var in ['raw_cdom', 'cdom']:
            bm[var].attrs['wavelength'] = 460

    # Subset data to derive stats.
    num_times = len(bm.instrument_time.values)
    tp = int(num_times * 0.1)
    edt = bm.instrument_time.values[-10]
    bdt = bm.instrument_time.values[-10 - (tp * 5)]
    subbm = bm.sel(instrument_time=slice(bdt, edt))

    bm['turbidity_mean'] = subbm.turbidity.mean()
    bm['chlorophyll_mean'] = subbm.chlorophyll.mean()
    bm['cdom_mean'] = subbm.cdom.mean()

    bm['turbidity_median'] = subbm.turbidity.median()
    bm['chlorophyll_median'] = subbm.chlorophyll.median()
    bm['cdom_median'] = subbm.cdom.median()

    # Assign factory device file info to root.
    bm.attrs['factory_cal_date'] = dev['cal_date']
    bm.attrs['serial_number'] = device_sn
    bm.attrs['onc_device_id'] = device_id
    bm.attrs['field_calibration_date'] = cal_date
    bm.attrs['calibration_type'] = cal_type
    bm.attrs['calibration_solution'] = cal_sol
    bm.attrs['attempt_number'] = attempt_number
    bm.attrs['device_file'] = dfp
    bm.attrs['raw_benchmark_file'] = bfp
    for k, v in dev.items():
        if k == 'NTU':
            bm['turbidity'].attrs = v
        elif k == 'CHL':
            bm['chlorophyll'].attrs = v
        elif k == 'CDOM':
            bm['cdom'].attrs = v


    NC_ENCODING = {'instrument_time': {'units': 'seconds since 1970-01-01'}}
    fn_verification = input(f"Is this an appropriate filename? {save_fn} [Y/N]: ")
    if fn_verification.upper() == 'Y':
        bm.to_netcdf(save_fn,
                     encoding=NC_ENCODING)
        print(f"Benchmark data saved to {save_fn}.")
    else:
        new_fn = input("Enter a new filename: ")
        bm.to_netcdf(new_fn,encoding = NC_ENCODING)


# Structure for an ECO Triplet Line
class ECOTripletLine(NamedTuple):
    instrument_time: datetime
    wvl1: int
    wvl1_counts: int
    wvl2: int
    wvl2_counts: int
    wvl3: int
    wvl3_counts: int
    raw_temp: int


def parse_data_line(line: str) -> ECOTripletLine:
    """
    Parse a single data line from and ECO Triplet .raw file.

    :param line: A string containing data from the ECO Triplet, typically tab-delimited.

    :return: An ECOTripletLine class containing line data.
    """
    (d, t,
     wvl1, wvl1_counts,
     wvl2, wvl2_counts,
     wvl3, wvl3_counts, raw_temp) = line.split('\t')
    dt = pd.to_datetime(' '.join((d, t))).strftime('%Y-%m-%dT%H:%M:%SZ')
    etl = ECOTripletLine(instrument_time=dt,
                         wvl1=int(wvl1),
                         wvl1_counts=int(wvl1_counts),
                         wvl2=int(wvl2),
                         wvl2_counts=int(wvl2_counts),
                         wvl3=int(wvl3),
                         wvl3_counts=int(wvl3_counts),
                         raw_temp=int(raw_temp))
    return etl


def etl2ds(parsed_line: ECOTripletLine) -> xr.Dataset:
    """
    Convert an ECOTripletLine object to an xarray dataset.

    :param parsed_line: An ECOTripletLine object.
    :return: An xarray Dataset containing parsed data.
    """
    _dict = parsed_line._asdict()
    _dt = pd.to_datetime(_dict['instrument_time'], format='%Y-%m-%dT%H:%M:%SZ')
    _ds = xr.Dataset()
    _ds = _ds.assign_coords({'instrument_time': [_dt],
                             'wavelength': [_dict['wvl1'], _dict['wvl2'],
                                            _dict['wvl3']]})
    _ds['raw_signal'] = (['instrument_time', 'wavelength'],
                         [[_dict['wvl1_counts'], _dict['wvl2_counts'],
                           _dict['wvl3_counts']]])
    _ds['raw_temp'] = (['instrument_time'], [_dict['raw_temp']])
    return _ds


def parse_header_lines(header_lines: list) -> dict:
    """
    Parse a list of header lines, if they exist in the raw file.

    :param header_lines: A list of lines containing header data from the ECO Triplet.
    :return: A dictionary containing parsed data.
    """

    header = {}
    for line in header_lines:
        if 'Ser' in line:
            sn = re.findall(f'(BBFL2{NUM_PAT})', line)[0]
            header['ser'] = sn
        elif 'Ver' in line:
            version = re.findall(f'(Triplet{SPACE_PAT}{NUM_PAT})', line)[0]
            header['ver'] = version
        elif 'Ave' in line:
            average = re.findall(f'Ave{SPACE_PAT}({NUM_PAT})', line)[0]
            header['ave'] = int(average)
        elif 'Pkt' in line:
            packet = re.findall(f'Pkt{SPACE_PAT}({NUM_PAT})', line)[0]
            header['pkt'] = int(packet)
        elif 'Set' in line:
            set = re.findall(f'Set{SPACE_PAT}({NUM_PAT})', line)[0]
            header['set'] = int(set)
        elif 'Rec' in line:
            rec = re.findall(f'Rec{SPACE_PAT}({NUM_PAT})', line)[0]
            header['rec'] = int(rec)
        elif 'Int' in line:
            interval = \
            re.findall(f'Int{SPACE_PAT}({NUM_PAT}:{NUM_PAT}:{NUM_PAT})', line)[0]
            header['int'] = interval
        elif 'Dat' in line:
            _date = re.findall(f'Dat{SPACE_PAT}({NUM_PAT}/{NUM_PAT}/{NUM_PAT})', line)[
                0]
            header['dat'] = _date
        elif 'Clk' in line:
            _time = re.findall(f'Clk{SPACE_PAT}({NUM_PAT}:{NUM_PAT}:{NUM_PAT})', line)[
                0]
            header['clk'] = _time
        elif 'Mem' in line:
            memory = re.findall(f'Mem{SPACE_PAT}({NUM_PAT})', line)[0]
            header['mem'] = int(memory)
    return header


def import_benchmark(filepath: PathLike) -> xr.Dataset:
    """
    Import an ECO Triplet Benchmark file and export it as an Xarray Dataset.

    :param filepath: A content or absolute filepat to the benchmark file.
    :return: An xarray Dataset.
    """
    with open(filepath) as _file:
        lines = _file.readlines()

    # Remove errant pad bytes from each line.
    lines = [line.replace('\x00', '') for line in lines]

    # Remove new line characters from each line.
    lines = [line.replace('\n', '') for line in lines]

    # Remove empty lines.
    lines = [line for line in lines if line]

    # Split lines into header and data lines based on header keywords.
    hkws = ['Ser', 'Ver', 'Ave', 'Pkt', 'Set', 'Rec', 'Int', 'Dat', 'Clk', 'Mem']
    header_lines = [line for line in lines if any(hkw in line for hkw in hkws)]
    data_lines = [line for line in lines if not any(hkw in line for hkw in hkws)]

    header = parse_header_lines(header_lines)

    etls = [parse_data_line(line) for line in data_lines]
    dss = [etl2ds(etl) for etl in etls]
    ds = xr.combine_by_coords(dss, combine_attrs='drop_conflicts')

    ds['instrument_time'] = ds['instrument_time'].astype('datetime64[s]')
    ds['wavelength'] = ds['wavelength'].astype('int32')
    ds['raw_signal'] = ds['raw_signal'].astype('int32')
    ds['raw_temp'] = ds['raw_temp'].astype('int32')

    # Assign root level attributes.
    ds.attrs = header  # Assign header info found in the file, if it exists.

    # Uncomment the following to attempt to automatically identify calibration conditions based on the filename.
    # fp = os.path.basename(filepath.lower())

    # # Assign a calibration solution.
    # if 'sz' in fp or 'spritezero' in fp or 'sprite_zero' in fp or 'sprite' in fp or 'zero' in fp:
    #     ds.attrs['cal_solution'] = 'SZ'
    # elif 'dc' in fp or 'dietcoke' in fp or 'diet_coke' in fp or 'coke' in fp or 'diet' in fp:
    #     ds.attrs['cal_solution'] = 'DC'
    # else:
    #     raise UserWarning(f'Filename does not contain calibration solution info: {fp}')
    #
    # # Assign calibration type.
    # if 'post' in fp or 'after' in fp or 'dirty' in fp:
    #     ds.attrs['cal_type'] = 'post_cal'
    # elif 'pre' in fp or 'before' in fp or 'clean' in fp:
    #     ds.attrs['cal_type'] = 'pre_cal'
    # else:
    #     raise UserWarning(f'Filename does not contain calibration type info: {fp}')

    return ds


def raw2val(output, scale_factor, dark):
    """
    Convert raw counts to something more understandable.

    :param output: The sensor output during the benchmark capture.
    :param scale_factor: The factory supplied scale factor found in the *.dev file or characterization sheet.
    :param dark: The factory supplied dark offset value found in the *.dev file or characterization sheet.
    :return:
    """

    val = scale_factor * (output - dark)
    return val


def parse_dev(filepath: PathLike) -> dict:
    """
    Parse a factory device file.

    :param filepath: The content or absolute file path to the device file.
    :return: A dictionary containing parsed data.
    """
    with open(filepath, 'r') as _file:
        lines = _file.readlines()

    dev = {}
    for line in lines:
        if 'BBFL2' in line:
            sn = re.findall(f'(BBFL2{NUM_PAT})', line)[0]
            dev['ser'] = sn
        elif 'Created' in line:
            cal_date = \
            re.findall(f":{SPACE_PAT}({NUM_PAT})/({NUM_PAT})/({NUM_PAT})", line)[0]
            cal_date = [s.zfill(2) for s in cal_date]
            dev['cal_date'] = datetime.strptime('/'.join(cal_date),
                                                '%m/%d/%y').strftime('%Y-%m-%d')
        elif 'Columns' in line:
            num_cols = int(re.findall(NUM_PAT, line)[0])
            dev['num_cols'] = num_cols
            cidx = lines.index(line)

    header_desc = {}
    col_lines = lines[cidx + 1:]
    for line in col_lines:
        field = re.findall(f"(.*?)=", line)[0]
        nums = re.findall(NUM_PAT, line)
        if len(nums) == 1:
            pos = int(nums[0])
            sf, off = (None, None)
        elif isinstance(nums, list):
            pos, sf1, sf2, off = nums
            sf = 'e'.join((sf1, sf2))
        header_desc[field] = {'column': int(pos),
                              'scale_factor': float(sf) if sf else None,
                              'offset': int(off) if off else None}
        header_desc = {k: v for k, v in header_desc.items() if k != 'N/U'}
    dev_desc = dev | header_desc
    return dev_desc

if __name__ == '__main__':
    main()