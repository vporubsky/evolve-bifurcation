# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 13:48:54 2019

@author: vporu

Title: Counting Number of Oscillators
"""
import os
numOscillator = 0
numModels = 0
numFiles = 0
for file in os.listdir(os.getcwd()):
    filename = os.fsdecode(file)
    numFiles += 1
    if filename.endswith("1.ant"): 
        numOscillator += 1
        continue
    elif filename.endswith("0.ant"):
        numModels += 1
        continue
print('number of oscillators: ' + str(numOscillator))
print('number of models: ' + str(numModels))
print('number of files: ' + str(numFiles))

# number of oscillators: 328
# number of total random networks: 9877
# frequency of oscillators: 328/9877 = 0.03320846410853498
# Lightly damped and sustained oscillators, some severely damped oscillators