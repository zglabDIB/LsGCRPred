import pandas as pd 
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
import pickle
import os

curr_dir=os.path.dirname(os.path.abspath(__file__))

df=pd.read_csv(curr_dir+'/input_files/input_expression_breast.txt',sep='\t')
num_samples=df.shape[1]-1
df.t_name=df.t_name.str.split('.').str[0]

feat_df=pd.read_csv(curr_dir+'/selected_transcripts/rs2366152.txt',sep='\t')

df = feat_df.merge(df, on='t_name', how='left')
df.loc[len(df)] = ['SNP']+[1]*num_samples


X_train=df.iloc[:, 1:].values
X_train=X_train.T

stdscaler=pickle.load(open(curr_dir+'/models/rs2366152/rs2366152_stdscaler.sav', 'rb'))
model=pickle.load(open(curr_dir+'/models/rs2366152/rs2366152_LR.sav', 'rb'))


X_train=stdscaler.transform(X_train)
pred_prob=np.round(model.predict_proba(X_train)[:,1]*100,2)

final_df=pd.DataFrame(columns=['sample_id','pred_proba'])
final_df['sample_id']=['rs2366152']
final_df['pred_proba']=pred_prob
final_df.to_csv(curr_dir+'/output_files/rs2366152_pred_ouput.txt',index=False,sep='\t')
