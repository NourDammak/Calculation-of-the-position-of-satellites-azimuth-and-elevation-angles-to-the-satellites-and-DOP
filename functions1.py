# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 15:51:07 2021

@author: Nour Dammak
"""

MU = 3.986004415e14
omega = 7.2921151467e-5
#a = 6378137
#e2=0.00669438002290


from datetime import date
import numpy as np
import math

def read_yuma(almanac_file):
  
    alm=open(almanac_file)
    alm_lines=alm.readlines()
    all_sats=[]
    for idx, value in enumerate(alm_lines):
        print(idx,value)
        
        if value[0:3]=="ID:":
            one_sat_block = alm_lines[idx:idx+13]
            one_sat=[]
            for line in one_sat_block:
                data=line.split(':')
                one_sat.append(float(data[1].strip()))
            all_sats.append(one_sat)
    alm.close()
    all_sats=np.array(all_sats)
    return(all_sats) 
    
    
    
    
    
def date2gpstime(date1):
    
    dif_of_days = date.toordinal(date(date1[0], date1[1],date1[2])) - date.toordinal(date(2019,4,7))
    
    weeks = dif_of_days//7
    
    day_of_week = dif_of_days-weeks*7
    
    #tow = time of the week
    
    tow = day_of_week*86400 + date1[3]*3600 + date1[4]*60 + date1[5]
     
    return [weeks, tow]





def satpos_alm(week,t,nav1sat):
    PRN = nav1sat[0]
    e = nav1sat[2]
    toa = nav1sat[3]
    i = nav1sat[4]
    OMGdot = nav1sat[5]
    a = nav1sat[6]**2
    OMG = nav1sat[7]
    omg = nav1sat[8]
    M = nav1sat[9]
    gps_week = nav1sat[-1]
    
    
    # The time since reference epoch
    #tk = t - toa
    Tk = (t + week*7*86400) - (toa + gps_week*7*86400)
    
    # Mean motion
    n=np.sqrt(MU/a**3)
    
    # Mean anomaly
    Mk = M + n*Tk
    
    # Eccentric anomaly
    Ek = Mk
    while 1:
        Eprev = Ek
        Ek = Mk+e*math.sin(Eprev)
        if abs(Eprev-Ek)<10**(-12): break
    
    
    # True anomaly 
    Vk=math.atan2(math.sqrt(1-e**2)*np.sin(Ek),np.cos(Ek)-e)
    
    # Arguent of latitude
    phik = Vk + omg
    
    # Orbit radius
    Rk = a*(1-e*math.cos(Ek))
    
    # Orbital coordinates
    xk = Rk*math.cos(phik)
    yk = Rk*math.sin(phik)
    
    # Longitude of ascending noce
    OMGk = OMG + (OMGdot - omega) * Tk - omega*toa
    
    # ECEF coordinates
    Xk = xk*math.cos(OMGk) - yk*math.cos(i)*math.sin(OMGk)
    Yk = xk*math.sin(OMGk) + yk*math.cos(i)*math.cos(OMGk)
    Zk = yk*math.sin(i)
    
    XYZ = np.array([Xk, Yk, Zk])
    return XYZ




