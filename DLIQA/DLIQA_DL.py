import numpy as np
import sys
import os
import pandas as pd


def dl_train():
    
    from keras.models import Sequential
    from keras.layers.convolutional import Conv2D
    from keras.layers import Conv2D,AveragePooling2D,Dropout,Flatten, Conv1D,AveragePooling1D, Dense, Merge, Activation
    from keras.layers.normalization import BatchNormalization
    from keras.optimizers import Adam
    
    import keras.backend as K #For compile


#------------------------------------------------------------------------------------

    model2D = Sequential()


    model2D.add(Conv2D(filters=64, kernel_size=(3,3), strides=(1,1), activation='relu',padding='valid', data_format='channels_first', input_shape=(8, 120, 128)))

    model2D.add(Conv2D(filters=64, kernel_size=(3,3), strides=(1,1), activation='relu',padding='valid', data_format='channels_first'))

    model2D.add(AveragePooling2D(pool_size=(2,2),strides=(2,2),data_format='channels_first'))

    model2D.add(Dropout(0.25)) 

    model2D.add(Conv2D(filters=128, kernel_size=(3,3), strides=(1,1), activation='relu',padding='valid', data_format='channels_first'))

    model2D.add(Conv2D(filters=128, kernel_size=(3,3), strides=(1,1), activation='relu',padding='valid', data_format='channels_first'))

    model2D.add(AveragePooling2D(pool_size=(2,2),strides=(2,2),data_format='channels_first'))

    model2D.add(Dropout(0.25)) 

    model2D.add(Flatten())

#------------------------------------------------------------------------------------

    model1D_rips = Sequential()


    model1D_rips.add(Conv1D(filters=128, kernel_size=3, strides=1, activation='relu',padding='valid', input_shape=(200, 36)))

    model1D_rips.add(Conv1D(filters=128, kernel_size=3, strides=1, activation='relu',padding='valid'))

    model1D_rips.add(AveragePooling1D(pool_size=2,strides=2))

    model1D_rips.add(Dropout(0.25)) 

    model1D_rips.add(Conv1D(filters=256, kernel_size=3, strides=1, activation='relu',padding='valid'))

    model1D_rips.add(Conv1D(filters=256, kernel_size=3, strides=1, activation='relu',padding='valid'))

    model1D_rips.add(AveragePooling1D(pool_size=2,strides=2))

    model1D_rips.add(Dropout(0.25)) 

    model1D_rips.add(Flatten())

#------------------------------------------------------------------------------------

    model1D_ripscharge = Sequential()

    model1D_ripscharge.add(Conv1D(filters=128, kernel_size=3, strides=1, activation='relu', input_shape=(100, 50)))
    model1D_ripscharge.add(Conv1D(filters=128, kernel_size=3, strides=1, activation='relu'))
    
    #model1D_ripscharge.add(Conv1D(filters=128, kernel_size=3, strides=1, activation='relu',padding='valid', input_shape=(100, 50)))

    #model1D_ripscharge.add(Conv1D(filters=128, kernel_size=3, strides=1, activation='relu',padding='valid'))

    model1D_ripscharge.add(AveragePooling1D(pool_size=2,strides=2))

    model1D_ripscharge.add(Dropout(0.25)) 
    
    model1D_ripscharge.add(Conv1D(filters=256, kernel_size=3, strides=1, activation='relu'))
    model1D_ripscharge.add(Conv1D(filters=256, kernel_size=3, strides=1, activation='relu'))

    #model1D_ripscharge.add(Conv1D(filters=256, kernel_size=3, strides=1, activation='relu',padding='valid'))

    #model1D_ripscharge.add(Conv1D(filters=256, kernel_size=3, strides=1, activation='relu',padding='valid'))

    model1D_ripscharge.add(AveragePooling1D(pool_size=2,strides=2))

    model1D_ripscharge.add(Dropout(0.25)) 

    model1D_ripscharge.add(Flatten())

#------------------------------------------------------------------------------------

    def f(M):
        m1, m2, m3 = M
        Y = m1 * K.concatenate([m2, m3])
        Y = K.reshape(Y, (-1, 2, 2))
        Y = K.sum(Y, 2)
        return Y

    #merged = Merge([model2D, model1D_rips, model1D_ripscharge], mode=f, output_shape=lambda _: (None, 2))
    merged = Merge([model2D, model1D_rips, model1D_ripscharge],mode='concat')
    
    
    
    #merged_expert=Merge([model1D_rips,model1D_ripscharge],mode='concat')
    #merged = Merge([model2D, merged_expert], mode='mul')

    tot_model= Sequential()
    tot_model.add(merged)

    tot_model.add(Dense(4096,activation='tanh'))
    
    tot_model.add(Dense(4096,activation='tanh'))
    
    tot_model.add(Dense(4096,activation='relu'))
    
    tot_model.add(Dense(4096,activation='relu'))

    tot_model.add(Dropout(0.5)) 
    
    tot_model.add(Dense(1, activation='softmax'))

    model2D.summary()
    model1D_rips.summary()
    model1D_ripscharge.summary()
    

    tot_model.summary()
    print('ready')
    def mean_pred(y_true, y_pred):
        return K.mean(y_pred)
    
    train_2D=np.random.rand(12,8,120,128)
    train_1D_rips=np.random.rand(12,200,36)
    train_1D_ripscharge=np.random.rand(12,100,50)
    target=np.random.rand(12,1)
    

    adam = Adam(lr=0.0001, decay=0.01)
    tot_model.compile(loss='mean_squared_error',
            optimizer=adam,
              metrics=['accuracy', mean_pred])
    
    print('Start training')
    
    #tot_model.fit([train_2D,train_1D_rips,train_1D_ripscharge],target,
     #     epochs=150,
      #    batch_size=16)
    
    print('Training done')





#------------------------------------------------------------------------------------

def main():
    
    dl_train()
    exit()
    valuefile= './CnM.featuresNPqDNZ'
    value_df=pd.read_csv(valuefile,delim_whitespace=1)
    data={}
    cross_val=open('./cross_val_sets/cross_val_set_1')
    
    lines = cross_val.read().splitlines()
    for protein in lines:    
        for file in os.listdir('./data/' + protein):
            if file.endswith("_feature_rips_charge_1D.npz"):
                    name=file[:-27]
                    
                    data_1D_charge = np.load(file)
                    data_1D = np.load(name + '_feature_rips_1D.npy')
                    data_2D = np.load(name + '_feature_alpha_2D.npy')
                    
                    for key,array in data_1D_charge.items():
                        data.setdefault('train_1D_ripscharge', []).append(data_1D_charge[key])
                    for key,array in data_1D.items():
                        data.setdefault('train_1D_rips', []).append(data_1D[key])
                    for key,array in data_2D.items():
                        data.setdefault('train_2D', []).append(data_2D[key])
                        
                    data.setdefault('target', []).append(value_df.loc[((value_df['#'] == name),'CPscore')].values[0])
                    

        

if __name__ == '__main__':
    main()
