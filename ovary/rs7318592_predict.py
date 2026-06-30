import pandas as pd 
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
import pickle
import os

curr_dir=os.path.dirname(os.path.abspath(__file__))

df=pd.read_csv(curr_dir+'/input_files/input_expression_ovary.txt',sep='\t')
num_samples=df.shape[1]-1
df.t_name=df.t_name.str.split('.').str[0]

feat_df=pd.read_csv(curr_dir+'/selected_transcripts/rs7318592.txt',sep='\t')

df = feat_df.merge(df, on='t_name', how='left')
df.loc[len(df)] = ['SNP']+[1]*num_samples

X_train=df.iloc[:, 1:].values
X_train=X_train.T

stdscaler=pickle.load(open(curr_dir+'/models/rs7318592/rs7318592_stdscaler.sav', 'rb'))
model=pickle.load(open(curr_dir+'/models/rs7318592/rs7318592_LR.sav', 'rb'))

X_train=stdscaler.transform(X_train)
pred_prob=np.round(model.predict_proba(X_train)[:,1]*100,2)

final_df=pd.DataFrame(columns=['sample_id','pred_proba'])
final_df['sample_id']=['rs7318592']
final_df['pred_proba']=pred_prob
final_df.to_csv(curr_dir+'/output_files/rs7318592_pred_ouput.txt',index=False,sep='\t')

