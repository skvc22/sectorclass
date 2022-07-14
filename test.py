import pandas as pd

maxRA = []
minRA = []

maxDec =[]
minDec = []

class testclass():
    def readInTables(self, link):
        """
        Read in a list of all tables from the link. Isolate the relevant table from the list.
        @param link: link to the TESS Year X Observations webpage.
        @return: dataframe of the data drawn from the webpage.
        """
        raw = pd.read_html(link)
                
        return raw[0]

    def listRA(self, link, camera):
        col = "Camera " + str(camera)
        result = [x.split(",")[0] for x in self.readInTables(link)[col]]
        result = [float(x) for x in result]
        return result
    
    def listDec(self, link, camera):
        col = "Camera " + str(camera)
        result = [x.split(",")[1] for x in self.readInTables(link)[col]]
        result = [float(x) for x in result]
        return result

    def testPandas(self):
        result = pd.DataFrame()
        for i in range(0,5):
            result['num'] = i
            result['other'] = "number" + str(i)
        print(result)
                    
if __name__ == "__main__":
    test = testclass()

    # test.maxAndMinRA()
    # test.maxAndMinDec()

    # print(maxRA)

    # testing = ["90.00, -66.56, 224.19", "90.00, -66.57, 196.63", "90.00, -66.57, 169.06","90.00, -66.56, 141.52","90.00, -66.56, 114.14"]
    # testing = ["a,b,c,d,e,f,g"]
    # df = pd.Series(testing)
    # df = df.str.split(",", expand=True)
    # print(df.to_string())

    # create a dataframe
    # df = pd.DataFrame({
    #     'Address': ['4860 Sunset Boulevard,San Francisco,California',
    #                 '3055 Paradise Lane,Salt Lake City,Utah',
    #                 '682 Main Street,Detroit,Michigan',
    #                 '9001 Cascade Road,Kansas City,Missouri']
    # })
    # # display the dataframe
    # print(df)

    # # split column into multiple columns by delimiter 
    # # split column and add new columns to df
    # df[['Street', 'City', 'State']] = df['Address'].str.split(',', expand=True)
    # # display the dataframe
    # print(df)

    mydataset = {
      'SN': ["2020ghq", "2020oat", "2017fdl", "2019vxm", "2019yvr","2017gjn", "2020jfo","2017glq"],
      'RA': [221.335125, 349.604917, 325.99489958, 299.618917, 191.283897036, 40.951792, 185.460333,221.335125],
      'Dec': [38.738419, 58.444131, -11.3445704648, 62.137731, -0.459120025663, 32.526069, 4.481681,38.738419]
    }

    myvar = pd.DataFrame(mydataset)

    myvar.to_csv("snlist.csv", index=False)
    

    ################################################################################



# frame = pd.DataFrame({"Max RA":maxRA, "Min RA":minRA, "Max Dec":maxDec, "Min Dec":minDec})
# print(frame.to_string())


# do w for loop

# if coord.lower() == "ra":
#     for i in range(1, 6): 
#         for j in range(1, 5):
#             print(test.listRA("http://tess.mit.edu/tess-year-"+str(i)+"-observations/", j))
# elif coord.lower() == "dec":
#     for i in range(1, 6): 
#         for j in range(1, 5):
#             print(test.listDec("http://tess.mit.edu/tess-year-"+str(i)+"-observations/", j))

# for i in range(len(sndf)):
#     coord1 = SkyCoord(ra=sndf.loc[i, 'RA'], dec=sndf.loc[i, 'Dec'], unit=u.deg)
            
#     for i in range(len(coords)):
#         coord2 = SkyCoord(ra=coords.loc[i, 'Central RA'], dec=coords.loc[i, 'Central Dec'], unit=u.deg)
#         separation = coord1.separation(coord2)
#         sepList.append(separation)

# radius = math.sqrt(2*(camWidth ** 2))

# coords['Minimum Dec'] = coords['Central Dec'] - (radius/2)
# coords['Maximum Dec'] = coords['Central Dec'] + (radius/2)

# coords['Minimum RA'] = (coords['Central RA'] - (radius/2)) 
# coords['Maximum RA'] = (coords['Central RA'] + (radius/2)) 

# for i in range(len(coords)):
#     coords.loc[i, 'Minimum RA'] = (coords.loc[i, 'Minimum RA'] / (math.cos(coords['Central Dec'][i]))) % 360
#     coords.loc[i, 'Maximum RA'] = (coords.loc[i, 'Maximum RA'] / (math.cos(coords['Central Dec'][i]))) % 360

# def minAndMax(self):
# """
# Call commaSplit to populate global DataFrame coords and update it with separate
# columns for each coordinate. Split coords into two dataframes, one for each 
# coordinate, and convert the contents from strings to floats. Clear coords and
# repopulate it with the maximum and minimum RA and Dec for each sector. Merge 
# coords with self.t, which contains the time each sector was observed.
# """
# self.commaSplit()

# global coords 

# RAdf = coords.loc[:, ['Camera 1 RA', 'Camera 2 RA', 'Camera 3 RA', 'Camera 4 RA']]
# decdf = coords.loc[:, ['Camera 1 Dec', 'Camera 2 Dec', 'Camera 3 Dec', 'Camera 4 Dec']]

# RAdf = RAdf.apply(pd.to_numeric)
# decdf = decdf.apply(pd.to_numeric)

# coords = coords.drop(['Camera 1 RA', 'Camera 1 Dec', 'Camera 2 RA', 'Camera 2 Dec', 'Camera 3 RA', 'Camera 3 Dec', 'Camera 4 RA', 'Camera 4 Dec'], axis=1)
# coords['Minimum RA'] = RAdf.min(axis=1)
# coords['Maximum RA'] = RAdf.max(axis=1)
# coords['Minimum Dec'] = decdf.min(axis=1)
# coords['Maximum Dec'] = decdf.max(axis=1)