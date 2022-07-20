#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 13:48:45 2022

@author: arest, svoncoelln
"""

import argparse,glob,re,sys,os,time,copy,math
from xxlimited import new
from pdastro import pdastroclass,pdastrostatsclass,makepath4file,unique,AnotB,AorB,AandB,rmfile
import yaml
import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from tessreduce import sn_lookup
import requests

coords = pd.DataFrame()

class sectorclass(pdastroclass):
    def __init__(self):
        pdastroclass.__init__(self)

        self.verbose=0
        
        # self.params will be populated with the arguments from the config file and options.
        # set default values here!
        self.params = {}

        self.params['verbose'] = 0

        self.params['outrootdir'] = '.'
        self.params['outsubdir'] = None

        self.params['overwrite'] = False

        self.t = pd.read_csv("https://tess.mit.edu/wp-content/uploads/orbit_times_20220505_1410.csv")


    def define_options(self,parser=None,usage=None,conflict_handler='resolve'):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage,conflict_handler=conflict_handler)

        # default for config file, if available
        if 'TESSLC_CFGFILE' in os.environ and os.environ['TESSLC_CFGFILE'] != '':
            cfgfilename = os.environ['TESSLC_CFGFILE']
        else:
            cfgfilename = None

        parser.add_argument('myinputfile',  help='input file. assumed to contain a list of SNe')
        
        parser.add_argument('-c','--configfile', type=str, default=cfgfilename, help='optional config file. default is set to $JWST_QUERY_CFGFILE. Use -vvv to see the full list of all set parameters.')

        parser.add_argument('--outrootdir', default=None, help='output root directory. The output directory the output root directory + the outsubdir if not None (default=%(default)s)')
        parser.add_argument('--outsubdir', default=None, help='outsubdir added to output root directory (default=%(default)s)')

        # boolean option
        parser.add_argument('--overwrite', default=False, action='store_true', help='overwrite files if they exist.')

        parser.add_argument('-v','--verbose', default=0, action='count')

        parser.add_argument('--outputfile', default=None, help='output file. default is None.')

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

        return(0)

    def set_outdir(self,outrootdir=None,outsubdir=None):
        # if no outrootdir is passed: use the one from the config file and options
        if outrootdir is None:
            outrootdir = self.params['outrootdir']

        # if no outrootdir is passed: use the one from the config file and options
        if outsubdir is None:
            outsubdir = self.params['outsubdir']

        if outrootdir is not None:
            outdir=outrootdir
        else:
            outdir='.'
        # strip superfluous '/' at the end
        outdir = outdir.rstrip("/")

        if outsubdir is not None:
            outdir+=f'/{outsubdir}'
        # strip superfluous '/' at the end
        outdir = outdir.rstrip("/")
             
        self.outdir = outdir
        return(outdir)

    def readInTables(self):
        """
        Read the tables on each year's observation page on the Tess website into the
        global DataFrame coords, looping through until the pages no longer exist to account
        for updates. Drop the Dates and Sector columns for redundancy, as well as the S
        pacecraft column, because the relevant coordinates are where the cameras targeted. 
        """
        global coords
        exists = True
        i = 1
        
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

    def commaSplit(self):
        """
        Split the columns for cameras 1-4 by comma into columns for the RA, Dec, and Roll. 
        Save those columns in a temporary DataFrame, and concatenate it vertically with a 
        second temporary DataFrame. Reorganize the rows to be grouped by sector, with each 
        sector containing four cameras, and save this to coords. 
        
        Replicate rows as necessary to allow for a successful merge of the two dataframes:
        self.t originally has one row per orbit and thus two per sector, and coords has 
        one row per camera and thus four per sector. Merge rather than concatenate to ensure 
        both lists are the same length and prevent NaNs in the final result.
        """
        self.readInTables()
        global coords
        temp = pd.DataFrame()
        newDF = pd.DataFrame()

        for i in range(1,5):
            camName = "Camera " + str(i)
            temp[['Central RA', 'Central Dec', 'Roll']] = coords[camName].str.split(', ', expand=True)
            temp['Camera'] = i
            temp = temp[['Camera','Central RA','Central Dec','Roll']]
            newDF = pd.concat([newDF, temp], axis=0)

        rows = len(coords)

        coords = newDF.iloc[::rows, :]
        coords = pd.concat([coords, coords])

        for i in range(1,rows):
            coords = pd.concat([coords, newDF.iloc[i::rows, :]])
            coords = pd.concat([coords, newDF.iloc[i::rows, :]])
        coords.index = range(len(coords.index))

        self.t = pd.DataFrame(np.repeat(self.t.values, 4, axis=0), columns=self.t.columns)
        self.t = pd.merge(self.t, coords, left_index=True, right_index=True)

    def getSeparation(self, inputra, inputdec, maxSep=17, sepcol='Separation'):
        """
        Create two skycoord objects, one for the input coordinates and one for each row
        in the dataframe. Calculate the separation between the two skycoord objects, and 
        save the result in a new column in the dataframe. Return a list of indices where the
        separation is within a certain range.

        @param inputra: right ascension to be used for the first skycoord object
        @param inputdec: declination to be used for the first skycoord object
        @param maxSep: maximum separation to allow in the returned DataFrame
        @param sepcol: name of the column to save the separation in
        @return: DataFrame made up of rows where the separation is within the given range
        """
        c1 = SkyCoord(ra=inputra, dec=inputdec, unit=u.deg)
        c2 = SkyCoord(ra=self.t['Central RA'], dec=self.t['Central Dec'], unit=u.deg)

        separation = c1.separation(c2)
        
        self.t[sepcol] = separation.degree

        observedlst = (np.where(self.t[sepcol]<maxSep)[0]).tolist()
        observedlst = observedlst[::2]
        
        result = self.t.loc[observedlst]
        return result

    
    def checkSNe(self, snName):
        """
        Call tessreduce on the inputted SN. Return a list of the sectors in which
        TESS has observed the SN.

        @param snName: name of the SN to be searched for
        @return: list of sectors in which TESS has observed the SN
        """
        sectors = []
        tessOutput = sn_lookup(snName, print_table=False)
        for i in range(len(tessOutput)):
            sectors.append(tessOutput[i][2])

        return sectors


if __name__ == '__main__':
    test = sectorclass()
    parser = test.define_options(usage='dummy.py is a dummy module')
    args = parser.parse_args()

    snListDF = pd.read_csv(args.myinputfile)

    # This reads the config file into self.params, and then
    # adds the command line paramters to it.
    test.get_arguments(args)
    if test.params['verbose']>1:
        print('params:', test.params)

    test.set_outdir()

    test.commaSplit()

    for i in range(len(snListDF)):
        print(snListDF['SN'][i] + ':')
        print(test.getSeparation(snListDF.iloc[i]['RA'], snListDF.iloc[i]['Dec']))

    # display the dataframe in the file passed as the output file
    if test.params['outputfile'] is not None:
        f = open(test.params['outputfile'], "w")
        f.write(test.t.to_string())
        f.close()
    
    print(f'This is the output root directory: {test.outdir}')