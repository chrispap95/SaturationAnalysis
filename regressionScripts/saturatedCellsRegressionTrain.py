import keras
from keras import backend
import uproot, numpy, h5py, pandas, root_numpy
import sklearn.model_selection
import sklearn.preprocessing
from ROOT import TCanvas, TH2D, TLine, gStyle

# Define custom Root Mean Square Error metric
def rmse(y_true, y_pred):
	return backend.sqrt(backend.mean(backend.square(y_pred - y_true), axis=-1))

def me(y_true, y_pred):
	return backend.mean((y_pred - y_true), axis=-1)

treename = 't1'
filename = {}
upfile = {}
params = {}
df = {}

filename['file1'] = 'TrainingSamples.nosync/out_E500to1500Eta1p62Phi0p0_converted.root'
outtag = 'Keras_saturated_E500to1500Eta1p62Phi0p0_3hl_25nd'
VARS = ['layer','n1','n2','n3','n4','n5','n6',
        'un1','un2','un3','un4','un5','un6',
        'dn1','dn2','dn3','dn4','dn5','dn6',
        'nup','ndown']
TARGET = ['rechit']

upfile['file1'] = uproot.open(filename['file1'])
params['input'] = upfile['file1'][treename].arrays(VARS)
params['target'] = upfile['file1'][treename].arrays(TARGET)

df['input'] = pandas.DataFrame(params['input'],columns=VARS)
df['target'] = pandas.DataFrame(params['target'],columns=TARGET)

NDIM = len(VARS)
inputs = keras.layers.Input(shape=(NDIM,), name = 'input')
hLayer = keras.layers.Dense(25, name = 'hidden1', kernel_initializer='he_uniform', activation='relu')(inputs)
hLayer = keras.layers.Dense(25, name = 'hidden2', kernel_initializer='he_uniform', activation='relu')(hLayer)
hLayer = keras.layers.Dense(25, name = 'hidden3', kernel_initializer='he_uniform', activation='relu')(hLayer)
outputs = keras.layers.Dense(1, name = 'output', activation='linear')(hLayer)

# creae the model
model = keras.models.Model(inputs=inputs, outputs=outputs)
# compile the model
model.compile(optimizer='adam', loss='mean_squared_error',metrics=[rmse,me,'mape'])
# print the model summary
model.summary()

df_all = pandas.concat([df['input'],df['target']],axis=1)
dataset = df_all.values
X = dataset[:,0:NDIM]
Y = dataset[:,NDIM]

X_train_val, X_test, Y_train_val, Y_test = sklearn.model_selection.train_test_split(X, Y, test_size=0.2, random_state=7)

# preprocessing: standard scalar
scaler = sklearn.preprocessing.StandardScaler().fit(X_train_val)
X_train_val = scaler.transform(X_train_val)
X_test = scaler.transform(X_test)

# early stopping callback
early_stopping = keras.callbacks.EarlyStopping(monitor='val_loss', patience=10)

# model checkpoint callback
# this saves our model architecture + parameters into dense_model.h5
model_checkpoint = keras.callbacks.ModelCheckpoint(outtag+'.h5', monitor='val_loss',
                                   verbose=0, save_best_only=True,
                                   save_weights_only=False, mode='auto',
                                   period=1)

# Train
print('\n'+80*'='+'\n'+'Fit model on training data\n'+80*'='+'\n')
history = model.fit(X_train_val,
                    Y_train_val,
                    epochs=1000,
                    batch_size=1024,
                    verbose=1, # switch to 1 for more verbosity
                    callbacks=[early_stopping, model_checkpoint],
                    validation_split=0.25)

# Evaluate the model on the test data using `evaluate`
print('\n\n'+120*'='+'\n'+'Training has finished. Evaluating on test data...\n'+120*'='+'\n')
results = model.evaluate(X_test, Y_test,verbose=0)
print('test_loss = %.5f, test_rmse = %.5f, test_me = %.5f, test_mape = %.5f'%(results[0],results[1],results[2],results[3]))

df['dense'] = model.predict(X_test) # add prediction to array

# Make a few plots
hist_data = numpy.concatenate((Y_test.reshape(Y_test.shape[0],1),df['dense'],df['dense']-Y_test.reshape(Y_test.shape[0],1)),axis=1)

gStyle.SetOptStat(0)

canvas1 = TCanvas('canvas1','canvas1',1)
hist1 = TH2D('hist1',';rechit_{true}[GeV];rechit_{regression}[GeV]',100,25,100,100,25,100)
root_numpy.fill_hist(hist1,hist_data[:,[0,1]])
hist1.Draw('colz')
ln1 = TLine(25,25,100,100)
ln1.SetLineStyle(2)
ln1.Draw()

canvas2 = TCanvas('canvas2','canvas2',1)
hist2 = TH2D('hist2',';rechit_{true}[GeV];rechit_{regression}-rechit_{true}[GeV]',100,25,100,100,-20,20)
root_numpy.fill_hist(hist2,hist_data[:,[0,2]])
hist2.Draw('colz')
ln2 = TLine(25,0,100,0)
ln2.SetLineStyle(2)
ln2.Draw()

# Store output to a root file
def get_features_from_file(filename='', treename='', branches=[]):
    t = root_numpy.root2array(filename, treename=treename, branches=branches)
    t = t.view(numpy.float32).reshape(t.shape + (-1,))
    return t

def write_prediction_to_file(features, scaler, model, filename='',treename='',branch=''):
    data_out = features
    features = features[:,0:NDIM]
    features = scaler.transform(features)
    y_predict_all = model.predict(features) # normal numpy array
    data_out = numpy.concatenate((data_out,y_predict_all),axis=1)
    dtype = numpy.dtype([('layer',numpy.float32),
                        ('n1',numpy.float32),('n2',numpy.float32),('n3',numpy.float32),
                        ('n4',numpy.float32),('n5',numpy.float32),('n6',numpy.float32),
                        ('un1',numpy.float32),('un2',numpy.float32),('un3',numpy.float32),
                        ('un4',numpy.float32),('un5',numpy.float32),('un6',numpy.float32),
                        ('dn1',numpy.float32),('dn2',numpy.float32),('dn3',numpy.float32),
                        ('dn4',numpy.float32),('dn5',numpy.float32),('dn6',numpy.float32),
                        ('nup',numpy.float32),('ndown',numpy.float32),
                        ('event',numpy.float32),('rechitsum',numpy.float32),
                        ('rechit',numpy.float32),('rechit_dnn', numpy.float32)])

    data_out = numpy.core.records.fromarrays(data_out.transpose(), dtype=dtype) # structured numpy array
    root_numpy.array2root(data_out, filename, treename=treename, mode='recreate')

OUT_VARS = ['layer',
            'n1','n2','n3','n4','n5','n6',
            'un1','un2','un3','un4','un5','un6',
            'dn1','dn2','dn3','dn4','dn5','dn6',
            'nup','ndown',
            'event',
            'rechitsum',
            'rechit']

X_all = get_features_from_file(filename['file1'],
                               treename=treename,
                               branches=OUT_VARS)

write_prediction_to_file(X_all,
                         scaler,
                         model,
                         filename='keras_results/'+outtag+'.root',
                         treename=treename,
                         branch='dense')
