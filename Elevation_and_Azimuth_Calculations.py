# -*- coding: utf-8 -*-
"""
Created on Sun Nov  7 16:03:30 2021

@author: Nour Dammak
"""

import math 
import numpy as np



#phi, lamda, h : Geodetic coordinates of the receiver
phi = 52
lamda = 21
h = 100
# s : distance between the sattelite and the receiver
#s = 
# z : zenith distance of thesatellite (z = el − 90o)
#z = 

e2= 0.00669438002290
a=6378137

#Xs : coordinates of the satellite
X = -20626210.63751
Y= -14177649.82168
Z = -9649997.6934
Xs = np.array([X, Y,Z]) 



N= a/(math.sqrt(1-(e2*math.sin(np.deg2rad(phi)**2))))

Xk = (N+h)*(math.cos(np.deg2rad(phi)))*(math.cos(np.deg2rad(lamda)))
Yk = (N+h)*math.cos(np.deg2rad(phi))*(math.sin(np.deg2rad(lamda)))
Zk=(N*(1-e2)+h)*math.sin(np.deg2rad(phi))


#Xr : coordinates of the receiver
Xr = np.array([Xk,Yk,Zk])



Xrs= Xs - Xr


u = np.array([math.cos(np.deg2rad(phi))*math.cos(np.deg2rad(lamda)),math.cos(np.deg2rad(phi))*math.sin(np.deg2rad(lamda)),math.sin(np.deg2rad(phi))])

n = np.array([-math.sin(np.deg2rad(phi))*math.cos(np.deg2rad(lamda)),-math.sin(np.deg2rad(phi))*math.sin(np.deg2rad(lamda)),math.cos(np.deg2rad(phi))])

e = (1/math.cos(np.deg2rad(phi)))* np.array([-math.sin(np.deg2rad(lamda)),math.cos(np.deg2rad(lamda)),0])

Rneu=np.array([[-math.sin(np.deg2rad(phi))*math.cos(np.deg2rad(lamda)),-math.sin(np.deg2rad(lamda)),(math.cos(np.deg2rad(lamda))*math.cos(np.deg2rad(phi)))],\
                [-math.sin(np.deg2rad(phi))*math.sin(np.deg2rad(lamda)),math.cos(np.deg2rad(lamda)),(math.cos(np.deg2rad(phi))*math.sin(np.deg2rad(lamda)))],\
                [math.cos(np.deg2rad(phi)),0,math.sin(np.deg2rad(phi))]])

#Rneu = np.array([n,e,u])
Xr_neu=np.dot(Rneu.T,Xrs)
print(Xr_neu)

u1=Xr_neu[2]
n1=Xr_neu[0]
e1=Xr_neu[1]
#Xr_neu = np.array([s*math.sin(z)*math.cos(Az),s*math.sin(z)*math.sin(Az),s*math.cos(z)])

Az = math.atan2(e1,n1)*180/math.pi
if Az<0:
    Az = Az+360


el=math.atan2(u1,(math.sqrt((e1**2)+(n1**2))))

el_in_degrees=el*180/math.pi

print(el_in_degrees)








