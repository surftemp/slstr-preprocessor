#!/usr/bin/env python
"""
Script to regrid SLSTR Vis/NIR channels to the 1 km IR "i stripe"
Also removes fire channels and met_tx data to reduce file size.

Usage:
    s3regrid.py [-h] [--cut-dirs CUT_DIRS] [--prefix PREFIX] file [file ...]

Arguments:
  file                 Input SLSTR zip file(s)

Optional arguments:
  -h, --help           show this help message and exit
  --cut-dirs CUT_DIRS  Cut this number of directory components
  --prefix PREFIX      Output directory prefix (defult current directory)

Example:

  s3regrid.py /neodc/sentinel3a/data/SLSTR/L1_RBT/2020/01/01/S3A_*.zip --cut-dirs=5

  Will save output files into the local directory 2020/01/01 as it is
  removing the first 5 path components: neodc...L1_RBT


"""

import argparse
import atexit
import datetime
import os
from pathlib import Path
import signal
import shutil
import subprocess
import sys
import tempfile
import time
import zipfile
import netCDF4

_workdir = tempfile.mkdtemp()   # Working directory for temporary files

_citation = "Embury, O., & McCarroll, N. (2023). surftemp/slstr-preprocessor: s3regrid (v1.1.0). Zenodo. https://doi.org/10.5281/zenodo.8017576"
_history = "s3regrid.py"

@atexit.register
def _final_cleanup():
    """
    atexit cleanup code to ensure we delete our tempdir
    """
    if _workdir:
        shutil.rmtree(_workdir)


def _catch_signals():
    """
    Set the standard signal handlers so we will always call the atexit code and
    perform appropriate cleanup of the runtime environment.
    """
    def _sigtermhandler(signum, frame):
        sys.exit("Caught signal: {0:d}".format(signum))

    for sig in ['SIGINT', 'SIGTERM', 'SIGUSR1', 'SIGUSR2', 'SIGXCPU']:
        try:
            signal.signal(getattr(signal, sig), _sigtermhandler)
        except AttributeError:
            pass


_catch_signals()


def findexec(name):
    """
    Look for the s3regrid executable. Default location is same as
    this script. Otherwise search the PATH environment variable
    """
    exec = os.path.join(sys.path[0], name)
    if os.path.isfile(exec):
        return exec
    for path in os.environ['PATH'].split(os.pathsep):
        exec = os.path.join(path, name)
        if os.path.isfile(exec):
            return exec
    return name


def replacenc(filename, vars=None, drop=[]):
    """
    Replace the specified netCDF file with a copy which only contains the
    listed variables

    Parameters
    ----------
    filename :  str
        netCDF file to replace

    vars : list
        list of variables to copy from input file. Setting vars to empty list `[]`
        will create an empty file with just dimensions and global attributes of
        the original. Setting vars to None will copy all variables to the output
        (unless explicitly excluded with `drop`)

    drop : list
        list of variables to exclude from input file. Will only be used if
        `vars=None`
    """
    src = netCDF4.Dataset(filename)
    dst = netCDF4.Dataset(filename+'_copy', 'w')
    dst.setncatts(src.__dict__)
    for dim in src.dimensions:
        if src.dimensions[dim].isunlimited():
            dst.createDimension(dim, None)
        else:
            dst.createDimension(dim, len(src.dimensions[dim]))
    if vars is None:
        vars = list(src.variables)
        vars = [v for v in vars if v not in drop]
    for v in vars:
        var = src.variables[v]
        opts = var.filters() or {}
        # netCDF4-python API has changed, so need to remove
        # compression filters from dict.
        comp = ['zlib', 'szip', 'zstd', 'bzip2', 'blosc']
        for c in comp:
            if opts.get(c):
                opts['compression'] = c
            if c in opts:
                del opts[c]
        if var.chunking() and var.chunking() != 'contiguous':
            opts['chunksizes'] = var.chunking()
        if hasattr(var, '_FillValue'):
            opts['fill_value'] = var._FillValue
        out = dst.createVariable(v, var.dtype, var.dimensions, **opts)
        out.setncatts(var.__dict__)
        out[:] = var[:]
    src.close()
    dst.close()
    os.replace(filename+'_copy', filename)


def renamenc(src):
    """
    Rename the specified netCDF file from a stripe to i. This should only be
    used on the quality files as it does not regrid data.
    """
    if src.name.endswith('_an.nc'):
        subs = '_an', '_in'
    elif src.name.endswith('_ao.nc'):
        subs = '_ao', '_io'
    else:
        raise Exception(f'Unrecognised file {src}')
    with netCDF4.Dataset(src.path, 'a') as nc:
        for v in list(nc.variables):
            nc.renameVariable(v, v.replace(*subs))
    os.rename(src, src.path[:-6] + subs[1] + '.nc')


def processzip(filename, useb=True, spt_chans=None):
    """
    Process an input SLSTR zip file to generate an SLSTR-Lite product.
    Parameters
    ----------
    filename :  str
        SLSTR L1b zip files to process

    useb : bool
        Include b-stripe data for channels S4-S6

    spt_chans : str
        Sub pixel texture argument to pass to s3regrid
    """
    zipname = os.path.basename(filename)
    safename = zipname.replace('.zip', '.SEN3')
    with zipfile.ZipFile(filename) as zip:
        # All members should be inside the safedir. If there are any
        # other files in the zip archive then abort.
        if [f for f in zip.namelist() if not f.startswith(safename)]:
            raise Exception('[slstr] Invalid SAFE - files outside safe dir')
        zip.extractall(_workdir)

    safedir = os.path.join(_workdir, safename)
    s3cmd = [findexec('s3regrid')]
    if useb:
        s3cmd.append('-b')
    if spt_chans:
        s3cmd.append('-s')
        s3cmd.append(spt_chans)
    s3cmd.append(safedir)
    subprocess.run(s3cmd, check=True)

    # Copy global attributes including resolution from S9_BT
    for v in 'no':
        with netCDF4.Dataset(os.path.join(safedir, f'S9_BT_i{v}.nc')) as nc:
            attrs = nc.__dict__
        # Append s3regrid info to relevant attributes
        attrs['comment'] = "SLSTR-Lite product"
        if 'references' in attrs:
            attrs['references'] += ", " + _citation
        else:
            attrs['references'] = _citation
        if 'history' in attrs:
            if isinstance(attrs['history'], str):
                attrs['history'] = [attrs['history'], _history]
            else:
                attrs['history'].append(_history)
        else:
            attrs['history'] = _history
        for c in '123456':
            with netCDF4.Dataset(os.path.join(safedir, f'S{c}_radiance_i{v}.nc'), 'a') as nc:
                nc.setncatts(attrs)

    # Remove any files not used for SST-CCI processing
    for f in os.scandir(safedir):
        # Need to keep the Vis/NIR quality files
        if 'quality_a' in f.name:
            # If needed we can rename the quality file/variables to match the
            # "i" grid used for data
            # renamenc(f)
            pass
        # All other a/b stripe files can be removed
        elif f.name[-5:] in ['an.nc', 'ao.nc', 'bn.nc', 'bo.nc', 'cn.nc', 'co.nc', 'fn.nc', 'fo.nc']:
            os.remove(f)
        # Also remove fire channel and met_tx
        elif f.name in ['met_tx.nc', 'F2_BT_in.nc', 'F2_BT_io.nc',
                        'F2_quality_in.nc', 'F2_quality_io.nc']:
            os.remove(f)

    # Cartesian and geodetic are the largest remaining files so remove any
    # variables that are considered unnecessary for 1-km products
    replacenc(os.path.join(safedir, 'cartesian_in.nc'), [])
    replacenc(os.path.join(safedir, 'cartesian_io.nc'), [])
    replacenc(os.path.join(safedir, 'geodetic_in.nc'), ['latitude_in', 'longitude_in', 'elevation_in'])
    replacenc(os.path.join(safedir, 'geodetic_io.nc'), ['latitude_io', 'longitude_io', 'elevation_io'])

    # Drop the orphan pixels from S7-S9
    for b in [7,8,9]:
        for v in ['n', 'o']:
            replacenc(os.path.join(safedir, f'S{b}_BT_i{v}.nc'), [f'S{b}_BT_i{v}', f'S{b}_exception_i{v}'])

    # Make SLSTR-Lite zip file
    subprocess.run(['zip', '-qr', zipname, safename], cwd=_workdir, check=True)
    shutil.rmtree(safedir)
    s1 = os.path.getsize(filename) /1024/1024
    s2 = os.path.getsize(os.path.join(_workdir, zipname)) /1024/1024
    print(f'Compressed {s1:.1f}M to {s2:.1f}M ({100*s2/s1:.1f}%)')
    return os.path.join(_workdir, zipname)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file', nargs='+', help='Input SLSTR zip file(s)')
    parser.add_argument('--cut-dirs', type=int, default=0, help='Cut this number of directory components')
    parser.add_argument('--prefix', default='.', help='Output directory prefix')
    parser.add_argument('--no-bstripe', dest='bstripe', action='store_false', help="Do not include b-stripe data in S4-S6")
    parser.add_argument('-s', '--standard-deviation', dest='chans', help='Channels for sub-pixel texture')
    args = parser.parse_args()

    # Do not include filename / path arguments in history as they can
    # be very long and are not useful to end users.
    _history = datetime.datetime.now(datetime.UTC).strftime('  %Y-%m-%dT%H:%M:%SZ') + ": s3regrid.py"
    if not args.bstripe:
        _history += ' --no-bstripe'
    if args.chans:
        _history += f' -s {args.chans}'
    time0 = time.monotonic()
    nfile = 0
    for src in args.file:
        print(f'Processing {src}')
        time1 = time.monotonic()
        spath = Path(src)
        n = args.cut_dirs + 1 if spath.is_absolute() else args.cut_dirs
        dst = Path(args.prefix).joinpath(*spath.parts[n:])
        try:
            f = processzip(src, args.bstripe, args.chans)

        except zipfile.BadZipFile as e:
            print(f'ERROR {e}: {src}')

        else:
            os.makedirs(dst.parent, exist_ok=True)
            shutil.move(f, dst)
            runtime = time.monotonic() - time1
            print(f'Output {dst} time: {datetime.timedelta(seconds=runtime)}')
            nfile += 1

    runtime = time.monotonic() - time0
    print(f'Processed {nfile}/{len(args.file)} files in {datetime.timedelta(seconds=runtime)}')

