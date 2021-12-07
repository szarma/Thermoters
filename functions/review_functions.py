import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def import_Zhou(
    shinedalgarno = "ACAGGAAACA",
    upElement = "AATGAGCTG",
    ShineDelgarno_toll = 3,
    seq_length_toll = 3,
):
    filename = "review_datasets/Data_for_Srdjan.xlsx"
    isheet = 4
    print(f"importing sheet no.{isheet} from {filename}")
    data = pd.read_excel(filename, sheet_name=isheet)
    
    data.columns = [c.strip().lower() for c in data.columns]
    
    print (f"adding the up element ({upElement}) to all the sequences")
    data["sequence"] = [upElement+s for s in data['sequence']]
    data["seq_length"] = data.sequence.apply(len)
    print ("here is a tally of the sequence lengths now")
    print (data["seq_length"].value_counts().sort_index())
    seqLength = data["seq_length"].value_counts().index[0]
    print (f"considering a consensus seqeunce length of {seqLength}")
    
    
    
    print (f"Considering Shine-Delgarno: '{shinedalgarno}'")
    data["ShineDelgarnoPos"] = data['sequence'].apply(lambda xi: xi.find( shinedalgarno) )
    print ("Here is a tally of Shine-Delgarno positions (-1 means SD is not being found, i.e. is being mutated)")
    print(data["ShineDelgarnoPos"].value_counts())
    ShineDelgarnoPos = data["ShineDelgarnoPos"].value_counts().index[0]
    print (f"Considering {ShineDelgarnoPos} as the consensus SD position.")
    fltr = (np.abs(data["ShineDelgarnoPos"]-ShineDelgarnoPos)<=ShineDelgarno_toll) & (np.abs(data["seq_length"]-seqLength)<=seq_length_toll)
    print (f"Only accepting sequences of length {seqLength} +/- {seq_length_toll}, and with Shine Delgarno at position {ShineDelgarnoPos} +/- {ShineDelgarno_toll}.")
    print (f"The number of sequences went from {data.shape[0]} to {fltr.sum()}.")
    data = data[fltr]
    
    
    print ("aligning sequences wrt to ShineDelgarno position...")
    
    reseq = []
    for ir,row in data.iterrows():
        delta = ShineDelgarnoPos-row.ShineDelgarnoPos
        if delta>=0:
            seq1 = "A"*delta + row.sequence
        else:
            seq1 = row.sequence[-delta:]

        deltaL = seqLength-len(seq1)
        if deltaL>=0:
            seq1 += "A"*deltaL
        else:
            seq1 = seq1[:deltaL]
        if len(seq1)!=seqLength:
            break
        reseq += [seq1]
    assert len( pd.Series.value_counts(list(map(len,reseq)))) == 1
    
    data["sequence_to_use"] = reseq
    data["loglums"] = np.log10(data["sfgfp/od600"])
    
    return {
        "dataset":data,
        "ShineDelgarnoPos":ShineDelgarnoPos
    }

def import_Hossain():
    # shinedelgarno = "aggag"
    data = pd.read_excel("review_datasets/Data_for_Srdjan.xlsx", sheet_name=0)
    data.columns = [c.lower() for c in data.columns]
    # data["promoter sequence"].apply(lambda xi: xi.find(shinedalgarno)).value_counts()
    # So, no Shine-Delgarno anywhere.
    data['promoter length'] = data["promoter sequence"].apply(len)
    print ("Here is a tally of sequence lengths")
    print(data["promoter length"].value_counts())
    print("The short one is the last one, and is used for normalization, so we exclude it from the dataset.")
    
    out = {"normRow": data.iloc[-1].copy()}
    data = data.iloc[:-1]

    # plt.hexbin(
    #     data['normalized-tx-rate-mean'].values,
    #     data['normalized-tx-rate-stdev'].values,
    #     cmap="hot",
    #     mincnt=1,
    #     xscale="log",
    #     yscale="log",
    #     gridsize=70
    # );
    # mr = np.geomspace(.01, 1e4)
    # k = (data['normalized-tx-rate-stdev']/data['normalized-tx-rate-mean']).mean()
    # plt.plot(mr, mr*k)
    # plt.gca().set_aspect("equal")
    # plt.grid()

    data["loglums"] = np.log10(data["normalized-tx-rate-mean"])

    data["loglums-stdev"] = data["normalized-tx-rate-stdev"]/data["normalized-tx-rate-mean"]/np.log(10)
    out["dataset"] = data
    return out


def import_Utrecho():
    
    print ("Importing sheet 3 from review_datasets/Data_for_Srdjan.xlsx. That's where the sequence are.")
    dataseq = pd.read_excel("review_datasets/Data_for_Srdjan.xlsx", sheet_name=3)
    print ("Dropping duplicates.")
    dataseq = dataseq.drop_duplicates()

    print ("Importing expressions from review_datasets/GSE108535_sigma70_variant_data.txt")
    dataexp = pd.read_csv("review_datasets/GSE108535_sigma70_variant_data.txt", sep=" ")
    dataexp.name = [name.replace("_","-") for name in dataexp.name]
    print ("Merging the two.")
    data = dataexp.merge(dataseq, on="name")
    data.columns = [c.strip() for c in data.columns] 
    data.columns = [c.lower() for c in data.columns]
    # data['sequence'].apply(lambda xi: xi.find(shinedalgarno)).value_counts()
    # So, no Shine-Delgarno anywhere.
    # data['promoter length'] = data["sequence"].apply(len)
    # data["promoter length"].value_counts()

    data['loglums'] = np.log10(data['rna_exp_average'])

    return {"dataset":data}


def import_Johns():
    
    print ("Importing sheet 1 from review_datasets/Data_for_Srdjan.xlsx.")
    data = pd.read_excel("review_datasets/Data_for_Srdjan.xlsx", sheet_name=1)
    data.columns = [c.strip() for c in data.columns] 
    data.columns = [c.lower() for c in data.columns]
    select = np.isfinite(data['expression in m9'])
    print (f"There is {select.shape[0]} sequences in the library, but {(~select).sum()} of them are not finite. I'll remove them.")
    data = data[select]
    select = data['expression in m9']>0
    print (f"In addition, there are {(~select).sum()} sequences with 0 expression. I'll remove them too.")
    data = data[select]
    data['loglums'] = np.log10(data['expression in m9'])

    return {"dataset":data}

def import_36N():
    data = pd.read_csv("36N_seqences/36N_constitutive.csv")
    data["loglums"] = data["estimate_facs"]
    return {"dataset":data}