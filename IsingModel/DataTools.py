#! /usr/bin/env python

import numpy as np

def writeDataToFile(filename, data,
            fieldNames=[],
            appendFile=False,
            addTimeField=False,
            dataFormat='%10.5f'):
    commentsStr = '#! '
    delimiterStr =' '
    if(addTimeField):
        time = np.arange(1, len(data[0])+1, dtype=np.float64)
        data.insert(0,time)
    if len(fieldNames)==len(data)-1:
        fieldNames.insert(0,'time')
    fieldsStr = ''
    if len(fieldNames)==len(data):
        fieldsStr += 'FIELDS '
        for f in fieldNames: fieldsStr += f + ' '
    data2 = np.column_stack(data)
    if appendFile:
        file = open(filename,'a')
    else:
        file = open(filename,'w')
    np.savetxt(file, data2 , header=fieldsStr, delimiter=delimiterStr, fmt=dataFormat, comments=commentsStr)
    file.close()
#-------------------------------------------------------------------------------

def getCommentsFromFile(filename,findString=''):
    # assume that the comments are in the first 100 lines
    MaxLines = 100
    commentsPrefix='#'
    comments = []
    file = open(filename,'r')
    for i in range(MaxLines):
        line = file.readline()
        if line[0:1]==commentsPrefix and line.find(findString) != -1:
            comments.append(line)
    file.close()
    return comments
#-------------------------------------------------------------------------------

def getFieldNames(filename):
    fieldsLine = getCommentsFromFile(filename,findString='FIELDS')[0]
    fields = fieldsLine.split()[2:]
    return fields
#-------------------------------------------------------------------------------

def getDataFromFile(filename,ignoreFieldNames=False):
    if ignoreFieldNames:
        fields = []
    else:
        fields = getFieldNames(filename)
    data = np.loadtxt(filename)
    numColumn = data.shape[1]
    data = np.hsplit(data,numColumn)
    for i in range(numColumn):
        data[i] = data[i].reshape( (data[i].size,) )
    return (data, fields)
#-------------------------------------------------------------------------------

def calculateAutocorrelation(data):
    # assert len(data.shape) == 1
    mean = np.mean(data)
    NumSamples = data.size
    autocorr = np.zeros(NumSamples)
    data = data-mean
    for i in range(NumSamples):
        sum = 0.0
        for k in range(NumSamples-i):
            sum += data[k]*data[k+i]
        autocorr[i] = sum/np.float(NumSamples-i)
    autocorr = autocorr/autocorr[0]
    return autocorr
#-------------------------------------------------------------------------------

def estimated_autocorrelation(x):
    """
    http://stackoverflow.com/q/14297012/190597
    http://en.wikipedia.org/wiki/Autocorrelation#Estimation
    """
    n = len(x)
    variance = x.var()
    x = x-x.mean()
    r = np.correlate(x, x, mode = 'full')[-n:]
    assert np.allclose(r, np.array([(x[:n-k]*x[-(n-k):]).sum() for k in range(n)]))
    result = r/(variance*(np.arange(n, 0, -1)))
    return result
#-------------------------------------------------------------------------------
