# sectorclass

### Determine whether the TESS telescope was pointed at a known supernova on or around its discovery date. 

## Python Scripts
- `sector.py`

#### Arguments
- `-c`, --configfile`: Specify the file name of the .ini file with settings for this script.
    - Type: str
    - Default: `config.ini`
    - Usage: `-c config.ini` or `--configfile config.ini`
- `-i`, `--inputfile`: A .csv file containing 3 colums: SNid, RA, and Dec.
    - Type: str
    - Default: `None`
    - Usage: `-i input.csv` or `--inputfile input.csv`
- `-o`, `--outputfile`: The .csv or .txt file to write output to. 
    - Type: str
    - Default: `None` 
    - Usage: `-o results.csv` or `--outputfile results.csv`
- `-v`, `--verbose`: Print information about the code as it runs.
    - Type: int
    - Default: `0` (i.e., nothing/very little printed as code runs)
    - Usage: `-v` or `-vv` or `-vvv`
- `--sectoroutput`: If specified, a .csv file to store raw TESS data about when it is pointing where; needed if updating the sectorlist.
    - Type: str
    - Default: `None` (i.e., use unupdated sectorlist.csv from the local directory)
    - Usage: `--sectoroutput newsecs.csv`
- `-u`, `--update`: If specified, retrieve most recent sector pointing data from MIT's TESS website and store in file passed with --sectoroutput.
    - Type: bool
    - Default: `False`
    - Usage: `-u` or `--update`

#### Example commands
- Check if TESS observed given supernovae near its discovery date: `./sector.py -i snlist.csv -o`output.csv
- Update the sectorlist.csv file: `./sector.py --sectoroutput newsectors.csv -u

###  File formats 
- snlist.csv must have 3 columns named SNid, RA, and Dec respectively.

  SNid | RA  |  Dec
--- | --- |  ---
2020ghq | 221.335125 | 38.738419
2020oat | 349.604917 | 58.444131
2017fdl | 325.99489958 | -11.3445704648
2019vxm | 299.618917 | 62.137731
2019yvr | 191.283897036 | -0.459120025663
2017gjn | 40.951792 | 32.526069
2020jfo | 185.460333 | 4.481681
2017glq | 221.335125 | 38.738419

- an output file will contain columns SNid, RA, Dec, Sector, Orbit, Orbit_Start, Orbit_End, Camera, Cam_RA, Cam_Dec, Cam_Roll, ID, Covers, MJD_Difference, where Covers is true if TESS observed the supernova on or around its discovery date.

## Dependencies
- Python 3.11.6 or higher
- math
- pandas
- numpy
- yaml
- pyparsing
- astropy
- requests
- tess_stars2px
- json
- tabulate # sectorclass
