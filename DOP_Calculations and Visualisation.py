# -*- coding: utf-8 -*-
"""
Created on Tue Nov  23 16:55:30 2021

@author: Nour Dammak
"""

import numpy as np
import math 
import statistics
from functions import *



def satpos(t,week,nav1sat):
    MU=3.986004415e14
    omge=7.2921151467e-5
   
# t=172800
    prn =nav1sat[0]
    e=nav1sat[2]
    toa=nav1sat[3]
    i=nav1sat[4]
    OMG_dot=nav1sat[5]
    A=nav1sat[6]**2
    OMG=nav1sat[7]
    omega=nav1sat[8]
    M=nav1sat[9]
    gps_week=nav1sat[-1]

    tk= (t+week*7*86400) - (toa + gps_week*7*86400)
    
#mean motion
    n=math.sqrt(MU/A**3)

#mean anomaly
    Mk =M+n*tk

#eccentric anomaly
    Ek=Mk

    while 1:
        Eprev=Ek
        Ek=Mk+e*math.sin(Eprev)
        if abs(Eprev-Ek)<1e-12: break

#true anomaly
    vk = math.atan2(math.sqrt(1-e**2)*math.sin(Ek), math.cos(Ek)-e)

#argument of Latitude
    phik = vk + omega

#radius of the orbit
    rk = A*(1-e*math.cos(Ek))

#orbital coordinates
    xk = rk*math.cos(phik)
    yk = rk*math.sin(phik)

# Longitude of ascending node
    OMGk= OMG + (OMG_dot - omge)*tk - omge*toa

#Calculation of the ECEF coordinates of the satellite:
    Xk = xk * math.cos(OMGk) - yk*math.cos(i)*math.sin(OMGk)
    Yk = xk * math.sin(OMGk) + yk*math.cos(i)*math.cos(OMGk)
    Zk = yk *math.sin(i)
    Xs=np.array([Xk,Yk,Zk])
    return Xs

almanac_file='almanac.yuma.week0131.319488.txt'
nav = read_yuma(almanac_file)
# nav1sat = nav[1,:]

date1 = [2021,  10, 12, 0, 0, 0]
week, t1 = date2gpstime(date1)

date2= [2021,  10, 12, 23, 59, 59]
week2, t2 = date2gpstime(date2)

#Calculations of 
psi=np.deg2rad(52)
lambda1=np.deg2rad(21)
h=100 #height above sea level
e2=0.00669438002290

a=6378137
N= a/ math.sqrt(1-e2*(math.sin(psi))**2*lambda1)

X= (N+h)*math.cos(psi)*math.cos(lambda1)
Y= (N+h)*math.cos(psi)*math.sin(lambda1)
Z= (N*(1-e2)+h)*math.sin(psi)
Xr=np.array([X,Y,Z])
mask=10

n_sat_all=[]
a=[]
GDOP_0=[]
PDOP_0=[]
TDOP_0=[]
PDOPneu_0=[]
HDOP_0=[]
VDOP_0=[]
dt=15*60

for t in range(t1,t2,dt):
    A=np.zeros((0,4))
    n_sat=0
    for sat in nav: 
        Xs=satpos(t,week,sat)
        
        r=np.linalg.norm(Xs-Xr)
        uv=np.array([[math.cos(psi)*math.cos(lambda1),math.cos(psi)*math.sin(lambda1),math.sin(psi)]])

        nv=np.array([[-math.sin(psi)*math.cos(lambda1),-math.sin(psi)*math.sin(lambda1), math.cos(psi)]])

        ev=np.array([[-math.sin(lambda1), math.cos(lambda1),0]])

        R=np.hstack([nv.T, ev.T, uv.T])

        Xrs=Xs-Xr

        Xrneu=np.dot (R.T,Xrs.T)

        #Calculation of the azimuth(Az) and elevation(el) angles
        n=Xrneu[0]
        e=Xrneu[1]
        u=Xrneu[2]

        Az=np.rad2deg(math.atan2(e,n))
        el=np.rad2deg(math.atan2(u,math.sqrt(e**2+n**2)))
        if el > mask:
            A1row=np.array([[-(Xs[0]-Xr[0])/r, -(Xs[1]-Xr[1])/r, -(Xs[2]-Xr[2])/r,1 ]])
            A=np.append(A,A1row, axis=0)
            n_sat+=1
    Q = np.linalg.inv(np.dot(A.T,A))
    Qxyz=np.array([[Q[0,0],Q[0,1],Q[0,2]],[Q[1,0],Q[1,1],Q[1,2]],[Q[2,0],Q[2,1],Q[3,2]]])
    Qneu=np.dot(R.T,Qxyz,R)
    GDOP=np.sqrt(Q[0,0] + Q[1,1] + Q[2,2] + Q[3,3])
    PDOP=np.sqrt(Q[0,0] + Q[1,1] + Q[2,2])
    TDOP=np.sqrt(Q[3,3])
    PDOPneu=np.sqrt(Qxyz[0,0] + Qxyz[1,1] + Qxyz[2,2])
    HDOP=np.sqrt(Qxyz[0,0] + Qxyz[1,1])
    VDOP=np.sqrt(Qxyz[2,2])
    
    #to verify that PDOPneu equal to sqrt((HDOP**2+VDOP**2)))
    #a=PDOPneu-(np.sqrt((HDOP**2+VDOP**2)))
    #print("a=",a)
    
    GDOP_0.append(GDOP)
    PDOP_0.append(PDOP)
    TDOP_0.append(TDOP)
    PDOPneu_0.append(PDOP)
    HDOP_0.append(HDOP)
    VDOP_0.append(VDOP)
    
    n_sat_all.append(n_sat)
   
print(GDOP_0)
print(PDOP_0)
print(TDOP_0)
print(PDOPneu_0)
print(HDOP_0)
print(VDOP_0)

print("maximum of GDOP values equals to :", max(GDOP_0))
print("maximum of PDOP values equals to :", max(PDOP_0))
print("maximum of TDOP values equals to :", max(TDOP_0))
print("maximum of PDOPneu values equals to :", max(PDOPneu_0))
print("maximum of HDOP values equals to :", max(HDOP_0))
print("maximum of VDOP values equals to :", max(VDOP_0))

print("the mean value of GDOP values equals to :",statistics.mean(GDOP_0))
print("the mean value of PDOP values equals to :",statistics.mean(PDOP_0))
print("the mean value of TDOP values equals to :",statistics.mean(TDOP_0))
print("the mean value of PDOPneu values equals to :",statistics.mean(PDOP_0))
print("the mean value of HDOP values equals to :",statistics.mean(HDOP_0))
print("the mean value of VDOP values equals to :",statistics.mean(VDOP_0))




print(n_sat_all)




#Visualisation of the number of visible satellites during the day by 15 minutes
time=(np.arange(t1,t2,15*60)-2*86400)/3600
import matplotlib.pyplot as plt 
fig = plt.figure(figsize=(14,5))
ax = fig.add_subplot()
ax.plot(time,n_sat_all)
ax.set_xticks([0,4,8,12,16,20,24])
ax.set_xlim([0,24])
ax.set_yticks(np.arange(0,16,1))
ax.grid(True)
ax.set_title("Number of visible satellites above the elevation mask",fontdict={"fontsize" :18,"fontweight" : "bold"})
ax.set_xlabel("time[h]")
ax.set_ylabel("number of sattelite")





