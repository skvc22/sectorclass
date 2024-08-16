#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 13:48:45 2022

@author: arest, svoncoelln
"""

import argparse,glob,re,sys,os,time,copy,math
from xxlimited import new

from pyparsing import col
from pdastro import pdastroclass,pdastrostatsclass,makepath4file,unique,AnotB,AorB,AandB,rmfile
import yaml
import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
import requests
from astropy.time import Time
import tess_stars2px
import json
from tabulate import tabulate

coords = pd.DataFrame()

class sectorclass(pdastroclass):
    def __init__(self):
        pdastroclass.__init__(self)

        self.verbose=0
        
        # self.params will be populated with the arguments from the config file and options.
        # set default values here!
        self.params = {}

        self.params['verbose'] = 0

        self.params['update'] = False

        self.params['inputfile'] = None
        self.params['outputfile'] = None 
        self.secs4SNe = pdastroclass()

        self.t = pd.read_csv('sectorlist.csv')

    def define_options(self,parser=None,usage=None,conflict_handler='resolve'):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage,conflict_handler=conflict_handler)

        # default for config file, if available
        if 'TESSLC_CFGFILE' in os.environ and os.environ['TESSLC_CFGFILE'] != '':
            cfgfilename = os.environ['TESSLC_CFGFILE']
        else:
            cfgfilename = None
        
        parser.add_argument('-c','--configfile', type=str, default=cfgfilename, help='optional config file. default is set to $JWST_QUERY_CFGFILE. Use -vvv to see the full list of all set parameters.')

        parser.add_argument('-i', '--inputfile', default=None, help='input file. assumed to contain a list of SNe')

        parser.add_argument('-o', '--outputfile', type=str, help='name of output file')

        parser.add_argument('-v','--verbose', default=0, action='count')

        parser.add_argument('--sectoroutput', default=None, help='output file to contain sector information. default is None.')

        parser.add_argument('-u','--update', default=False, action='store_true', help='update the input list of sectors and cameras. pass output file name with --outputfile.')

        return(parser)

    def get_arguments(self, args, configfile = None):
        '''

        Parameters
        ----------
        args : list
            pass the command line arguments to self.params.
        configfile : string, optional
            Config filename. The default is None. If None, then
            $TESSLC_CFGFILE is used if exists.

        Returns
        -------
        None.

        '''

        def subenvvarplaceholder(paramsdict):
            """ Loop through all string parameters and substitute environment variables. environment variables have the form $XYZ """
            envvarpattern=re.compile('\$(\w+)')

            for param in paramsdict:
                if isinstance(paramsdict[param], str):
                    envvarnames=envvarpattern.findall(paramsdict[param])
                    if envvarnames:
                        for name in envvarnames:
                            if not (name in os.environ):
                                raise RuntimeError("environment variable %s used in config file, but not set!" % name)
                            envval=os.environ[name]
                            subpattern='\$%s' % (name)
                            paramsdict[param] = re.sub(subpattern,envval,paramsdict[param])
                elif isinstance(paramsdict[param], dict):
                #elif type(dict[param]) is types.DictType:
                    # recursive: sub environment variables down the dictiionary
                    subenvvarplaceholder(paramsdict[param])
            return(0)

        # get the parameters from the config file
        if args.configfile is not None:
            #cfgparams = yaml.load_file(args.configfile)
            if not os.path.isfile(args.configfile):
                raise RuntimeError('config file %s does not exist!' % (args.configfile))
            print(f'Loading config file {args.configfile}')
            cfgparams = yaml.load(open(args.configfile,'r'), Loader=yaml.FullLoader)
            self.params.update(cfgparams)

            # substitute environment variables
            subenvvarplaceholder(self.params)

            if args.verbose>2:
                print('\n### CONFIG FILE PARAMETERS:')
                for p in cfgparams:
                    print('config file args: setting %s to' % (p),cfgparams[p])

        # Go through optional parameters.
        # 'None' does not overwrite previously set parameters (either default or config file)
        if args.verbose>2:
            print('\n### OPTIONAL COMMAND LINE PARAMETERS:')
        argsdict = vars(args)
        for arg in argsdict:

            # skip config file
            if arg=='configfile': continue

            if argsdict[arg] is not None:
                if args.verbose>2:
                    print('optional args: setting %s to %s' % (arg,argsdict[arg]))
                self.params[arg]=argsdict[arg]
            
            else:
                if arg not in  self.params:
                    self.params[arg]=None
        
        # self.params['input'] = args.inputfile

        return(0)

    def readInTables(self):
        """
        Set self.t equal to a csv of sectors, orbits, and orbit start and end times. Read 
        the tables on each year's observation page on the Tess website into the global DataFrame 
        coords, looping through until the pages no longer exist to account for updates. Drop 
        the Dates and Sector columns for redundancy, as well as the Spacecraft column, because 
        the relevant coordinates are where the cameras targeted. 
        """
        global coords
        exists = True
        i = 1
        self.t = pd.read_csv('https://tess.mit.edu/public/files/TESS_orbit_times.csv', comment='#')

        while exists == True:
            link = "http://tess.mit.edu/tess-year-"+str(i)+"-observations/"

            raw = pd.read_html(link)
            temp = raw[0].drop(['Dates','Spacecraft','Sector'], axis=1)

            coords = pd.concat([coords, temp], axis=0)

            i += 1 
            response = requests.get("http://tess.mit.edu/tess-year-"+str(i)+"-observations/")
            if response.status_code != 200:
                exists = False

        coords = coords.loc[:, ~coords.columns.str.contains('^Unnamed')]
        coords = coords.drop_duplicates()
        coords = coords.reset_index(drop=True)

    def toMJD(self):
        """
        Rename the orbit start and end columns to avoid spaces in column names, then convert
        the contents of the columns to MJD.
        """
        self.t.rename(columns = {'Start of Orbit':'Orbit_Start', 'End of Orbit':'Orbit_End'}, inplace = True)
        startlist = self.t['Orbit_Start'].astype(str).tolist()
        endlist = self.t['Orbit_End'].astype(str).tolist()

        starttime = Time(startlist, format='iso', scale='utc')
        endtime = Time(endlist, format='iso', scale='utc')

        self.t['Orbit_Start'] = starttime.mjd
        self.t['Orbit_End'] = endtime.mjd

    def commaSplit(self):
        """
        Split the columns for cameras 1-4 by comma into columns for the RA, Dec, and Roll. 
        Save those columns in a temporary DataFrame, and concatenate it vertically with a 
        second temporary DataFrame. Reorganize the rows to be grouped by sector, with each 
        sector containing four cameras, and save this to coords. 
        """
        print('Splitting comma-separated columns...')
        global coords
        temp = pd.DataFrame()
        newDF = pd.DataFrame()

        for i in range(1,5):
            camName = "Camera " + str(i)
            temp[['Cam_RA', 'Cam_Dec', 'Cam_Roll']] = coords[camName].str.split(', ', expand=True)
            temp['Camera'] = i
            temp = temp[['Camera','Cam_RA','Cam_Dec','Cam_Roll']]
            newDF = pd.concat([newDF, temp], axis=0)

        rows = len(coords)

        coords = newDF.iloc[::rows, :]
        coords = pd.concat([coords, coords])

        for i in range(1,rows):
            coords = pd.concat([coords, newDF.iloc[i::rows, :]])
            coords = pd.concat([coords, newDF.iloc[i::rows, :]])
        coords.index = range(len(coords.index))

    def finalize_table(self):
        """
        Replicate rows as necessary to allow for a successful merge of the two dataframes:
        self.t originally has one row per orbit and thus two per sector, and coords has 
        one row per camera and thus four per sector. Merge rather than concatenate to ensure 
        both lists are the same length and prevent NaNs in the final result.
        """
        global coords

        self.t = pd.DataFrame(np.repeat(self.t.values, 4, axis=0), columns=self.t.columns)
        self.t = pd.merge(self.t, coords, left_index=True, right_index=True)

        orbitList = self.t['Orbit'].tolist()
        cameraList = self.t['Camera'].tolist()
        ID = []

        for i in range(len(orbitList)):
            temp = str(orbitList[i]) + '_' + str(cameraList[i])
            ID.append(temp)

        self.t['ID'] = ID

        # display the dataframe in the file passed as the output file
        self.t.to_csv(self.params['sectoroutput'], index=False)

    def update(self):
        """
        Update self.t with data from the website instead of reading in a stored csv.
        """
        self.readInTables()
        self.toMJD()
        self.commaSplit()
        self.finalize_table()

    def get_separation(self, inputra, inputdec, maxSep=17, sepcol='Separation'):
        """
        Find separation between input coordinates and central coordinates of each orbit+camera combination.

        @param inputra: right ascension to be used for the first skycoord object
        @param inputdec: declination to be used for the first skycoord object
        @param maxSep: maximum separation to allow in the returned DataFrame
        @param sepcol: name of the column to save the separation in
        @return: int array of sectors where the separation is within the given range
        """
        c1 = SkyCoord(ra=inputra, dec=inputdec, unit=u.deg)
        c2 = SkyCoord(ra=self.t['Cam_RA'], dec=self.t['Cam_Dec'], unit=u.deg)

        separation = c1.separation(c2)
        
        self.t[sepcol] = separation.degree
        # print(self.t.to_string())

        ixs = self.ix_inrange(sepcol,None,maxSep)
        #return self.write(indices=ixs,columns=[sepcol,'RA','Dec'])
        # return ixs.tolist()
        result = np.array(self.t.loc[ixs,:]['Sector'].drop_duplicates())
        return result.astype(int)

    def check_sectors_TESS(self, starIDs, starRas, starDecs):
        """
        Based on the tess_point.py tess_stars2px_function_entry function, modified
        to search through a given list of sectors rather than all of them to determine
        in which sectors a given SN was observed.

        @param starIDs: holdover from the original function
        @param starRas: right ascension of the target
        @param starDecs: declination of the target
        @return: list of sectors that contain the target
        """
        inSec = self.get_separation(starRas, starDecs)
        # Instantiate Spacecraft position info
        scinfo = tess_stars2px.TESS_Spacecraft_Pointing_Data()
        # Now make list of the star objects
        starList = tess_stars2px.make_target_objects(np.atleast_1d(starIDs), \
                                    np.atleast_1d(starRas), np.atleast_1d(starDecs))
        # Make rough determination as to which pointing camera combos are worth
        # Checking in detail and then do detailed checking
        findAny=False
        outSec = np.array([-1], dtype=int)
        for i, curTarg in enumerate(starList):
            for curSec in inSec:
                starRas = np.array([curTarg.ra])
                starDecs =  np.array([curTarg.dec])
                idxSec = curSec-1

                starInCam, starCcdNum, starFitsXs, starFitsYs, starCcdXs, starCcdYs = scinfo.fpgObjs[idxSec].radec2pix(\
                        starRas, starDecs)
                for jj, cam in enumerate(starInCam):
                    # SPOC calibrated FFIs have 44 collateral pixels in x and are 1 based  
                    xUse = starCcdXs[jj] + 45.0
                    yUse = starCcdYs[jj] + 1.0
                    xMin = 44.0
                    ymaxCoord = 2049
                    xmaxCoord = 2093
                    if xUse>xMin and yUse>0 and xUse<xmaxCoord and yUse<ymaxCoord:
                        if findAny==False:
                            outSec[0] = curSec
                            findAny=True
                        else:
                            outSec = np.append(outSec, curSec)
        return outSec

    def observation_info(self, name, time='disc', lower=30, upper=100):
        """
        Adapted from tessreduce sn_lookup function. Call check_sectors_TESS to 
        determine which sectors a given SN was observed in, and then determine
        whether it was observed in the given time range.

        @param name: name of the SN to be searched for
        @param time: set whether operations should be performed relative to SN discovery or max MJD
        @param lower: lower bound of the time range
        @param upper: upper bound of the time range
        @return outSecs: list of sectors in which TESS has observed the SN
        @return tab: list of lists containing the an SN was observed in, whether it was observed in 
        the given time range, and margin by which it missed the relevant time range
        """
        url = 'https://api.astrocats.space/{}'.format(name)
        response = requests.get(url)
        json_acceptable_string = response.content.decode("utf-8").replace("'", "").split('\n')[0]
        d = json.loads(json_acceptable_string)
        if list(d.keys())[0] == 'message':
            print(d['message'])
            return None
        else:
            disc_t = d[name]['discoverdate'][0]['value']
            disc_t = Time(disc_t.replace('/','-'))
            
        max_t = d[name]['maxdate'][0]['value']
        max_t = Time(max_t.replace('/','-'))

        ra = d[name]['ra'][-1]['value']
        dec = d[name]['dec'][-1]['value']
        c = SkyCoord(ra,dec, unit=(u.hourangle, u.deg))
        ra = c.ra.deg
        dec = c.dec.deg
                
        outSecs = self.check_sectors_TESS(0, ra, dec)

        if len(outSecs) > 0:
            secs = pd.DataFrame()
            for i in outSecs:
                secs = pd.concat([secs, self.t.loc[self.t['Sector'] == i]])
            if (time.lower() == 'disc') | (time.lower() == 'discovery'):
                disc_start = secs['Orbit_Start'].values - disc_t.mjd
                disc_end = secs['Orbit_End'].values - disc_t.mjd
            elif (time.lower() == 'max') | (time.lower() == 'peak'):
                disc_start = secs['Orbit_Start'].values - max_t.mjd
                disc_end = secs['Orbit_End'].values - max_t.mjd

            covers = []
            differences = []
            tab = []
            for i in range(len(disc_start)):
                ds = disc_start[i]
                de = disc_end[i]
                if (ds-lower < 0) & (de + upper> 0):
                    cover = True
                    dif = 0
                elif (de+upper < 0):
                    cover = False
                    dif = de
                elif (ds-lower > 0):
                    cover = False
                    dif = ds
                covers += [cover]
                differences += [dif]
                tab += [[secs.Sector.values[i], cover, dif]]
            return outSecs, tab
        else:
            print('No TESS coverage')
            return None
            
    def get_secs4SN(self, SNid):
        """
        @param SNid: name of the SN to be searched for
        @return: indices of the table containing the SN
        """
        ixs_sn = self.secs4SNe.ix_equal('SNid', SNid)
        return(ixs_sn)

    def load_secs4SNe(self, filename=None):
        """
        Load the table containing the SNe.

        @param filename: name of the file to be loaded
        """
        if filename is None:
            filename = self.params['inputfile']
        self.secs4SNe.t = pd.read_csv(filename)

        parts = filename.split('.')
        if parts[1] == 'csv':
            self.secs4SNe.t = pd.read_csv(filename)
        else: 
            self.secs4SNe.load(filename)
            
    def save_secs4SNe(self, filename=None):
        """
        Save the table containing the SNe.

        @param filename: name of the file to be saved to
        """
        if filename is None:
            filename = self.params['outputfile']

        parts = filename.split('.')
        if parts[1] == 'csv':
            self.secs4SNe.t.to_csv(filename, index=False)
        else: 
            self.secs4SNe.write(filename)        

    def populate_secs4SNe(self, low=None, high=None):
        """
        Populate the secs4SNe table with the information from the TESS sector files. 
        Replicate rows in the dataframe containing SN names and coordinates as necessary
        to allow for successful concatenation with relevant rows from self.t containing
        sector information, as well as the information on whether the SN was observed in
        the given time range. 

        @param low: lower bound of the time range
        @param high: upper bound of the time range
        """
        self.load_secs4SNe(self.params['inputfile'])

        temp = pd.DataFrame(columns=self.secs4SNe.t.columns)
        templst = []
        secslst = []
        secdf = pd.DataFrame()

        for i in range(len(self.secs4SNe.t)):
            if self.observation_info(self.secs4SNe.t['SNid'][i]) is not None:
                if low is not None and high is not None:
                    outsecs, lst = self.observation_info(self.secs4SNe.t['SNid'][i], lower=low, upper=high)
                else:
                    outsecs, lst = self.observation_info(self.secs4SNe.t['SNid'][i])
                templst += lst
                secslst += outsecs.astype(float).tolist()
                for j in range(len(lst)):
                    temp.loc[len(temp.index)] = self.secs4SNe.t.values[i]

        self.t = self.t.drop(['Separation'], axis=1)
        
        for i in secslst:
            secdf = pd.concat([secdf, self.t.loc[self.t['Sector'] == i]])

        cov = pd.Series(dtype=bool)
        dif = pd.Series(dtype=float)

        for i in range(len(templst)):
            cov = pd.concat([cov, pd.Series(templst[i][1])])
            dif = pd.concat([dif, pd.Series(templst[i][2])])
         
        secdf.index = range(len(secdf.index))
        cov.index = range(len(cov.index))
        dif.index = range(len(dif.index))

        secs = pd.concat([temp, secdf, cov.rename('Covers'), dif.rename("MJD_Difference")], axis=1)
        self.secs4SNe.t = secs

        self.save_secs4SNe(self.params['outputfile'])

if __name__ == '__main__':
    test = sectorclass()
    parser = test.define_options(usage='dummy.py is a dummy module')
    args = parser.parse_args()

    # This reads the config file into self.params, and then
    # adds the command line paramters to it.
    test.get_arguments(args)

    if test.params['verbose']>1:
        print('params:', test.params)

    if test.params['update']:
        test.update()
    else: 
        test.populate_secs4SNe()