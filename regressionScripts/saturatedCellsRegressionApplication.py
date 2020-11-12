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

filename['file1'] = 'EvaluationSamples2/out_E2000Eta1p62Phi0p0_0_converted.root'
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
df_all = pandas.concat([df['input'],df['target']],axis=1)
dataset = df_all.values
X = dataset[:,0:NDIM]
Y = dataset[:,NDIM]

# preprocessing: standard scalar
scaler = sklearn.preprocessing.StandardScaler().fit(X)
X = scaler.transform(X)

# load json and create model
dependencies = {
    'rmse': rmse,
    'me': me
}
model = keras.models.load_model("dense_model.h5", custom_objects=dependencies)
model.summary()
print("Loaded model from disk")

#df['dense'] = model.predict(X) # add prediction to array

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
                         filename='keras_results/flatRegressionResult_E2000Eta1p62Phi0p0_0_lr12.root',
                         treename=treename,
                         branch='dense')
