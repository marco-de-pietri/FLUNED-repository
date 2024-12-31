
import re
import struct
import binascii
import numpy as np


def get_fluent_parse_regions(filename):
    """
    this function looks for the string where the hierarchy of the region/faces
    is defined
    """

    headerPtrn = re.compile('^\(cfd-post-mesh-info .*\n')

    



    with open(filename,'rb') as fi:
        fileLines = fi.readlines()
        for line in fileLines:
            matchH = headerPtrn.findall(line.decode("utf-8", "ignore"))
            if len(matchH) != 0:
                break




    return matchH[0]


def get_fluent_parse_headers(filename):
    """
    this function returns a list of dictionaries that contains the headers
    of the ascii/bin files
    """

    headerPtrn = re.compile('^\(\d+\s+.*\(.*\).*\n')

    headerDefPat = re.compile('\(.*?\)')

    indexPat = re.compile('\d+')

    dictVect = []

    with open(filename,'rb') as fi:
        fileLines = fi.readlines()
        for line in fileLines:
            matchH = headerPtrn.findall(line.decode("utf-8", "ignore"))
            if len(matchH) != 0:
                newD = {}
                rawString = matchH[0]
                headerDef = headerDefPat.findall(rawString[1:].strip(' (\n'))
                
                index = int(indexPat.findall(rawString)[0])
                newD['raw'] = rawString
                newD['index'] = index
                newD['header'] = headerDef[0].strip('( )')
                dictVect.append(newD)


    for v in dictVect:
        print(v)


    return dictVect

def get_fluent_binarray_mid(filename,checkStr,index):
    """ this function read a binary/ascii file of the cas/dat format and look
    for the line that start with index and that has the checkStr in the middle 
    and returns the header string"""

    headerPat = re.compile(' (.*).*?\Z')

    strStart = '(' + str(index) + ' ' 
    strStartHex = binascii.hexlify(strStart.encode())
    strStartPat = re.compile(strStartHex)

    strMidHex = binascii.hexlify(checkStr.encode())
    strMidPat = re.compile(strMidHex)

    strEnd = '\n' 
    strEndHex = binascii.hexlify(strEnd.encode())
    strEndPat = re.compile(strEndHex)

    with open(filename,'rb') as fi:
        data = fi.read()
        txtHex = binascii.hexlify(data)

        startMatches = strStartPat.finditer(txtHex)
        startMatches =  list(startMatches)

        endMatches = strEndPat.finditer(txtHex)
        endMatches = list(endMatches)


        midMatches = strMidPat.finditer(txtHex)
        midMatches = list(midMatches)
        midMatch = midMatches[-1] 



        for startMatch in startMatches:


            if startMatch.start() > midMatch.start():

                stringStart = temp
                break
            else:
                temp = startMatch.start()

        for endMatch in endMatches:
            if endMatch.start() > midMatch.end():
                stringEnd = endMatch.end()
                break

        matchString = txtHex[stringStart:stringEnd]

    pString=binascii.unhexlify(matchString).decode('ascii').strip(' \n')[1:-1]

    pString = headerPat.findall(pString)[0][1:-1]

    return pString

def get_fluent_binarray_header(filename,index):
    """ this function read a binary/ascii file of the cas/dat format and look
    for the line that starts with the index and returns the header string"""

    strStart = '(' + str(index) + ' '
    strStartHex = binascii.hexlify(strStart.encode())
    strStartPat = re.compile(strStartHex)

    strEnd = ')\n' 
    strEndHex = binascii.hexlify(strEnd.encode())
    strEndPat = re.compile(strEndHex)

    headerPat = re.compile(' (.*).*?\Z')

    with open(filename,'rb') as fi:
        data = fi.read()
        txtHex = binascii.hexlify(data)

        startMatch = strStartPat.search(txtHex)

        endMatches = strEndPat.finditer(txtHex)


        for endMatch in endMatches:
            if endMatch.start() > startMatch.end():
                stringEnd = endMatch.end()
                break


        matchString = txtHex[startMatch.start():stringEnd]

    pString=binascii.unhexlify(matchString).decode('ascii').strip(' \n')[1:-1]

    pString = headerPat.findall(pString)[0][1:-1]


    return pString


def get_fluent_binarray_double(filename,index,nEntry):
    """ this function read a binary/ascii file of the cas/dat format and look
    for the data that has the passed header"""

    strStart = '(' + str(index) + ' '
    strStartHex = binascii.hexlify(strStart.encode())
    strStartPat = re.compile(strStartHex)
    strEndBinary = ')\nEnd '
    strEndBinaryHex = binascii.hexlify(strEndBinary.encode())
    strEnd = str(index) +  ')' 
    strEndHex = binascii.hexlify(strEnd.encode())
    strEndPat = re.compile(strEndHex)
    strEndBinPat = re.compile(strEndBinaryHex)
    #strEndHex = hex(strEnd)

    with open(filename,'rb') as fi:
        data = fi.read()
        txtHex = binascii.hexlify(data)

        startMatch = strStartPat.search(txtHex)
        endMatch = strEndPat.search(txtHex)

        matchString = txtHex[startMatch.start():endMatch.end()]


        temp = b''
        for i in range(0, len(matchString),2):
            check = matchString[i:i+2]

            if check == b'28':
                temp = b''
            elif check == b'29':
                startData = i + 6 # move after: )\n(
                break
            else:
                temp += check
        


        headerStr = binascii.unhexlify(temp).decode('ascii')

        headerValues = headerStr.split()

        idMin = headerValues[1]

        nVals = int(headerValues[2],16)

        endBinMatch = strEndBinPat.search(matchString)

        match2String = matchString[endBinMatch.start():]
        
        rawData=binascii.unhexlify(matchString[startData:endBinMatch.start()])

        

        fString = '<' + str(nVals*nEntry) + 'd'

        valuesNP = np.frombuffer(rawData)



    return valuesNP


        


