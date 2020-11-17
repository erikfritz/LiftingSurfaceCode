# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 11:43:18 2020

@author: fritzek
"""

import pandas as pd

def readInput(sfolder,sfile):
#function that reads the input as specified in standard ECN AeroModule input files
#the information is stored inthe dictionary dcinput
    finput = open(sfolder+'/'+sfile,'r')
    l1input = finput.readlines()
    dcinput = {}
    for i in range (0,len(l1input)):
        if l1input[i][0] == '!':
            pass
        else:
            l1line = l1input[i].split()
            if l1line[0] == 'AEROPROPS':
                dcinput['AEROPROPS'] = pd.DataFrame(columns=['zB','chord','t/c','twist','C14','xB','yB'])
            else:
                try: 
                    float(l1line[0])
                    dcinput['AEROPROPS'] = dcinput['AEROPROPS'].append(pd.Series(l1line,index=dcinput['AEROPROPS'].columns),ignore_index=True)
                except:
                    skey = l1line[0]
                    dcinput[skey] = l1line[1]
                    
    if 'INCLUDE' in dcinput:
        sspecialist = dcinput['INCLUDE']
        dcinput = readSpecialist(sfolder,sspecialist,dcinput)
        
    return dcinput

def readSpecialist(sfolder,sspecialist,dcinput):
#function that reads the specialist input as specified in standard ECN AeroModule specialist input files
#the information is added to the dictionary dcinput
    fspecialist = open(sfolder+'/'+sspecialist,'r')
    l1input = fspecialist.readlines()
    for i in range (0,len(l1input)):
        if l1input[i][0] == '!':
            pass
        else:
            l1line = l1input[i].split()
            skey = l1line[0]
            dcinput[skey] = l1line[1]
            
    return dcinput