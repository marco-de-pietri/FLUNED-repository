import h5py
import re
import sys
import os

def get_dataset_keys(f):
    ks = []
    f.visit(lambda key : ks.append(key) if isinstance(f[key], h5py.Dataset) 
        else None)
    return ks

def lookup_h5_regex(file, regPat):
    """
    this function just returns the matching path
    """

    allMatches = []
    with h5py.File(file, "r") as fi:
        
        paths = get_dataset_keys(fi)

        for path in paths:
            matches = regPat.findall(path)
            if len(matches) == 1:
                allMatches.append(path)

        if len(allMatches) != 1:
            print ("ERROR with regex search:")
            print (regPat)
            print ("matches found: ", len(allMatches))
            for val in allMatches:
                print (val)

            sys.exit()

        else:

            matchPath = allMatches[0]

    return  matchPath

def get_h5_path_dataset(file, regPat):
    """
    this function open the h5 file, it searches for
    an array and returns the path and its array
    """
    
    with h5py.File(file, "r") as fi:

        allMatches = []
            
        paths = get_dataset_keys(fi)

        for path in paths:
            matches = regPat.findall(path)
            if len(matches) == 1:
                allMatches.append(path)

        if len(allMatches) != 1:
            print ("ERROR with regex search:")
            print (regPat)
            print ("matches found: ", len(allMatches))
            for val in allMatches:
                print (val)

            sys.exit()

        else:

            matchPath = allMatches[0]
            matchArray = fi[matchPath][()]

    return  matchPath, matchArray



def get_h5_dataset(file, regPat):
    """
    this function open the h5 file, it searches for
    an array and returns the path and its array
    """
    
    with h5py.File(file, "r") as fi:

        allMatches = []
            
        paths = get_dataset_keys(fi)

        for path in paths:
            matches = regPat.findall(path)
            if len(matches) == 1:
                allMatches.append(path)

        if len(allMatches) != 1:
            print ("ERROR with regex search:")
            print (regPat)
            print ("matches found: ", len(allMatches))
            for val in allMatches:
                print (val)

            sys.exit()

        else:

            matchPath = allMatches[0]
            matchArray = fi[matchPath][()]

    return   matchArray

def get_h5_dataset_multi(file, regPat):
    """
    this function open the h5 file, it searches for
    an array and returns arrays in a vector of dictionary
    """
    
    resultArray = []

    with h5py.File(file, "r") as fi:

        allMatches = []
            
        paths = get_dataset_keys(fi)

        for path in paths:
            matches = regPat.findall(path)
            if len(matches) == 1:
                allMatches.append(path)

        for patMatch in allMatches:
            newDict = {}
            newDict['minID'] = fi[patMatch].attrs['minId']
            newDict['maxID'] = fi[patMatch].attrs['maxId']
            newDict['values'] = fi[patMatch][()]
            
            resultArray.append(newDict)


    return   resultArray

def extract_multiblock(minID,maxID,dataBlock):

    found = False

    
    for block in dataBlock:

        if (minID >= block['minID']) and (maxID <= block['maxID']):

            
            found = True

            blockIDmin = int(minID - block['minID'])
            blockIDmax = int(maxID - block['minID']+1)

            array = block['values'][blockIDmin:blockIDmax]



    if not found:
        print ("ERROR range not found in the datablock")
        sys.exit()



    return array
