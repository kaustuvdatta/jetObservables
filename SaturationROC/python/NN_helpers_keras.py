import os, sys
import array
import glob
import numpy as np
import h5py

from tensorflow import keras
from tensorflow.keras import layers,callbacks,models,optimizers,layers,initializers,regularizers,constraints
from tensorflow.keras.models import model_from_json

import random
import math
import sklearn
from sklearn.metrics import roc_curve, roc_auc_score, auc
from sklearn.model_selection import train_test_split

import copy
import matplotlib.cm as cmx
import matplotlib.colors as colors

np.random.RandomState(seed=42)

class BinaryClassifier:
    def __init__(self, input_shape, hidden_size1, hidden_size2, hidden_size3, hidden_size4, name):
        print("input shape = ", input_shape)
        self.input_shape = input_shape
        self.hidden_size1 = hidden_size1
        self.hidden_size2 = hidden_size2
        self.hidden_size3 = hidden_size3
        self.hidden_size4 = hidden_size4
        self.modelname = name
    def build_model(self):
        model = keras.Sequential([
            layers.Dense(self.hidden_size1, activation='relu', input_shape=(self.input_shape,)),
            layers.Dropout(0.2),
            layers.Dense(self.hidden_size2, activation='relu'),
            layers.Dropout(0.2),
            layers.Dense(self.hidden_size3, activation='relu'),
            layers.Dropout(0.1),
            layers.Dense(self.hidden_size4, activation='relu'),
            layers.Dropout(0.1),
            layers.Dense(1, activation='sigmoid'),
        ],  name=self.modelname)

        model.compile(optimizer='adam',
                      loss='binary_crossentropy',
                      metrics=['accuracy'])

        return model

def split_data(inputs, target, train_size=0.8, test_size=0.2, shuffle=True):
    # Split the data into training, validation, and testing sets, by default 60:20:20 splitting (change to make more generic later)
    if train_size+test_size>1.:
        print('Error: Cannot have train and test splits add up to greater than 1. Aborting.')
        return 1
    
    if shuffle:     
        inds_shuffle = np.random.permutation(inputs.shape[0])
        inputs = inputs[inds_shuffle]
        target = target[inds_shuffle]
    
    inp, inputs_test, target, target_test = train_test_split( inputs, target, 
                                                              test_size=test_size,
                                                              train_size=train_size,
                                                              random_state=42)
    
    inputs_train, inputs_validate, target_train, target_validate = train_test_split(inp,target,test_size = 0.2,train_size =0.8)
    
    data_dict = dict()
    data_dict['inputs'] = {'train':inputs_train,
                           'test':inputs_test,
                           'validate':inputs_validate
                          }
    data_dict['targets'] = {'train':target_train,
                            'test':target_test,
                            'validate':target_validate
                           } 
    
    return data_dict


def get_cmap(N):
    
    '''Returns a function that maps each index in 0, 1, ... N-1 to a distinct 
    RGB color.'''
    color_norm  = colors.Normalize(vmin=0, vmax=N-1)
    scalar_map = cmx.ScalarMappable(norm=color_norm, cmap='hsv') 
    def map_index_to_rgb_color(index):
        return scalar_map.to_rgba(index)
    return map_index_to_rgb_color

def loadmodel(name, inDir='', weights = True):
    json_file = open(f'{inDir}/{name}_m.json', 'r')
    loaded_model_json = json_file.read()
    json_file.close()
    model = model_from_json(loaded_model_json)
    #load weights into new model
    if weights==True:
        model.load_weights(f'{inDir}/{name}_w.h5')
    #print (model.summary())
    print("Loaded model %s from disk"%name)
    return model

def savemodel(model,outDir='',name="neural network"):
    #print "Saving model:"
    model_name = name
    #model.summary()
    model.save_weights(f"{outDir}/{model_name}_w.h5", overwrite=True)
    model_json = model.to_json()
    with open(f"{outDir}/{model_name}_m.json", "w") as json_file:
        json_file.write(model_json)
        
def savelosses(hist,outDir='',name="neural network"):  
    print ("Saving losses:")
    loss = np.array(hist.history['loss'])
    valoss = np.array(hist.history['val_loss'])
    f = h5py.File(f"{outDir}/{name}_h.h5","w")
    f.create_dataset('loss',data=loss)
    f.create_dataset('val_loss',data=valoss)
    f.close()
    
    
def makeAUC(sig, bkg):
    
    sig = np.sort(sig.flatten())
    bkg = np.sort(bkg.flatten())
    
    #print sig, bkg
    
    min_sb = np.minimum(sig[0], bkg[0])
    max_sb = np.maximum(sig[sig.shape[0]-1], bkg[bkg.shape[0]-1])
    max_sb = max_sb+max_sb/10001.
    
    interval = (max_sb-min_sb)/10001.
    
    #print max_sb, min_sb
    
    hs, se = np.histogram(a=sig, bins=np.arange(min_sb,max_sb,interval))
    hb, be = np.histogram(a=bkg, bins=np.arange(min_sb,max_sb,interval))
      
    #print hs.shape[0], hb.shape[0]
    
    #plt.plot(se[0:100],hs)
    #plt.plot(be[0:100],hb)
    #plt.show()
    
    lr = np.empty((10000,2))
    xb = 0
    xs = 0
       
    for k in range(0, 10000):
        xb=hb[k]
        xs=hs[k]
      
        if(xs==0):
            lr[k][0]=k
            lr[k][1]=0.0
      
        else:
            lr[k][0]=k
            lr[k][1]=(xb/(float)(xs+xb))
        
    for p in range(0, 10000):
        for q in range(0, 10000-1):
            if(lr[q][1]>lr[q+1][1]):
                bint=lr[q+1][0]
                lt=lr[q+1][1]
                lr[q+1][0]=lr[q][0]
                lr[q+1][1]=lr[q][1]
                lr[q][0]=bint
                lr[q][1]=lt
          
    l = np.empty((10000))
    ls = np.empty((10000))
    lb = np.empty((10000))
      
    for k in range(0, 10000):
        
        l[k] = lr[k][1]
        bin = (int)(lr[k][0])
        ls[k] = hs[bin]
        lb[k] = hb[bin]
      
    #ls.Scale(1.0/ls.Integral())
    #lb.Scale(1.0/lb.Integral())
    
    #print ls, lb
    
    ls = ls/np.sum(ls)
    lb = lb/np.sum(lb)
    
    #print ls
    #print lb
    
    sums=0
    sumb=0
    
    x = np.empty((10000))
    #x = np.arange(0,10000)
    y = np.empty((10000))
    #y = np.arange(0,10000)
    
    cut=0
    for k in range(0, 10000):
        sums+=ls[k]
        sumb+=lb[k]
        x[k]=sumb
        y[k]=sums
        
        #if k%==0:
        #    print k, sumb, sums, x[k], y[k]
        
    area = sklearn.metrics.auc(x,y)
    
    return (area)