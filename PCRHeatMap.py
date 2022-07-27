import pandas as pd
import io
import numpy as np
from numpy.ma.core import copy
import seaborn as sns
from matplotlib import pyplot as plt

class RunFile:

    def __init__(self, fn):
        self.fn = fn
        self.df = pd.read_csv(io.BytesIO(uploaded[fn]), header=12)
    
    def last_index(self, col):
        endIndex = 0
        for i, x in enumerate(col):
            if pd.isna(x) == True:
                #print(i)
                endIndex = i
                break
        return(endIndex)

    def new_df(self, df, df_type):
        type_dict = {0:'MS2 Ct', 1:'ORF1ab Ct', 2:'S gene Ct', 3:'N gene Ct', 4:'Interpretive Result'}
        x = df[['Well', 'Well Number', type_dict[df_type]]].copy().fillna(0)
        filler_df = pd.DataFrame({"Well": ["M23", "M24", "N23"], "Well Number":[310, 311, 334], type_dict[df_type]: [0, 0, 0]})
        x = x.append(filler_df, ignore_index=True)
        return(x)

    def reshape_df(self, df, df_type):
        type_dict = {0:'MS2 Ct', 1:'ORF1ab Ct', 2:'S gene Ct', 3:'N gene Ct', 4:'Interpretive Result'}
        rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P']
        df = df.sort_values(by=['Well Number'], ignore_index=True)
        npdf = df[type_dict[df_type]].to_numpy().reshape(16,24)
        xdf = pd.DataFrame(npdf, columns = list(range(1,25)))
        xdf['Row'] = rows
        xdf = xdf.set_index('Row')
        return(xdf)

    def well_dict(self):
        wellnumberlist = list(range(0,384))
        wellvaluelist = []
        for wellnumber in wellnumberlist:
            wellvaluelist.append(0)

        welldict = dict(zip(wellnumberlist, wellvaluelist))
        return(welldict)

    def mask_dict(self):
        wellnumberlist = list(range(0,384))
        wellvaluelist = []
        for wellnumber in wellnumberlist:
            wellvaluelist.append(True)

        maskdict = dict(zip(wellnumberlist, wellvaluelist))
        return(maskdict)

    def df_to_dict(self, df, df_type):
        type_dict = {0:'MS2 Ct', 1:'ORF1ab Ct', 2:'S gene Ct', 3:'N gene Ct', 4:'Interpretive Result'}
        x = df[type_dict[df_type]].fillna(0)
        y = df['Well Number'].fillna(0)
        vals = []
        wells = []
        for i in x:
            vals.append(i)
        for j in y:
            wells.append(j)
        z = dict(zip(wells, vals))
        return(z)

    def sample_dict(self):
        firstColumn = self.df['Batch ID']
        df1 = self.df[0:self.last_index(firstColumn)]
        x = df1['Sample Name'].fillna(0)
        y = df1['Well Number'].fillna(0)
        vals = []
        wells = []
        for i in x:
            vals.append(i)
        for j in y:
            wells.append(j)
        z = dict(zip(wells, vals))
        return(z)

    def dict_to_df(self, dfdict, df_type):
        type_dict = {0:'MS2 Ct', 1:'ORF1ab Ct', 2:'S gene Ct', 3:'N gene Ct', 4:'Interpretive Result'}
        refdict = self.well_dict()
        welllist = []
        vallist = []
        
        for i in dfdict:
            x = i
            y = dfdict[i]
            refdict[x] = y

        for i in refdict:
            x = i
            y = refdict[i]
            welllist.append(x)
            vallist.append(y)
        
        df_data = {'Well Number': welllist, type_dict[df_type]: vallist}
        newdf = pd.DataFrame(data=df_data)
        return(newdf)

    def mask_df(self, dftype):
        type_dict = {0:'MS2 Ct', 1:'ORF1ab Ct', 2:'S gene Ct', 3:'N gene Ct', 4:'Interpretive Result'}
        value_dict = self.sample_dict()
        ref_mask = self.mask_dict()
        welllist = []
        vallist = []

        for i in value_dict:
            x = i
            y = value_dict[i]
            if type(y) is str:
                ref_mask[x] = False
        for i in ref_mask:
            x = i
            y = ref_mask[x]
            welllist.append(x)
            vallist.append(y)

        df_data = {'Well Number': welllist, type_dict[dftype]: vallist}
        newdf = pd.DataFrame(data=df_data)
        newmaskdf = self.reshape_df(newdf, dftype)
        return(newmaskdf)

    def find_min(self, df):
        l = df.to_numpy().flatten()
        dfvals = []
        for i in l:
            if i > 0:
                dfvals.append(i)
        return(min(dfvals))
    
    def find_max(self, df):
        l = df.to_numpy().flatten()
        dfvals = []
        for i in l:
            if i > 0:
                dfvals.append(i)
        return(max(dfvals))


    def find_mean(self, df):
        l = df.to_numpy().flatten()
        dfvals = []
        for i in l:
            if i > 0:
                dfvals.append(i)
        ave = sum(dfvals)/len(dfvals)
        return(ave)

    def d_frame(self, dftype):
        firstColumn = self.df['Batch ID']
        df1 = self.df[0:self.last_index(firstColumn)]
        if self.last_index(firstColumn) == 381:
            df2 = self.new_df(df1, dftype)
            df3 = self.reshape_df(df2, dftype)
            return(df3)
        elif self.last_index(firstColumn) < 381:
            new_dict = self.df_to_dict(df1, dftype)
            df2 = self.dict_to_df(new_dict, dftype)
            df3 = self.reshape_df(df2, dftype)
            return(df3)
   
    def plot_it(self, dftype, minval=18, addmask=False):
        type_dict = {0:'MS2 Ct', 1:'ORF1ab Ct', 2:'S gene Ct', 3:'N gene Ct', 4:'Interpretive Result'}
        x = self.d_frame(dftype)
        mask = self.mask_df(dftype)
        plt.figure(figsize=(15,8))
        plt.title(self.df['Batch ID'][0] + ": " + type_dict[dftype], fontsize=20)
        minimum_value = self.find_min(x)
        maximum_value = self.find_max(x)
        mean_value = self.find_mean(x)
        print(f'Min: {minimum_value}, Max: {maximum_value}, Mean: {mean_value}')
        if type(minval) is str:
            minimum = self.find_min(x)
        else:
            minimum = minval
        if addmask == True:
            sns.heatmap(x, vmin=minimum, mask=mask, annot=True, linewidths=1, cmap="icefire")
        else:
            sns.heatmap(x, vmin=minimum, annot=True, linewidths=1, cmap="icefire")
            

# y = RunFile('C3037FHC_analysis_data_1654193379100 - Sheet1.csv')
# y.plot_it(0, minval='auto', addmask=True)
            
