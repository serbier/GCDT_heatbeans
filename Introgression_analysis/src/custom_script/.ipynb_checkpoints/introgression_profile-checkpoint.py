import argparse
import time 
import numpy as np
import pandas as pd
import csv
from collections import defaultdict
import random 

import csv
import numpy as np
import matplotlib.pyplot as plt


def getFreqs(pathFile):
    nsamples = defaultdict(list)
    MAF = defaultdict(dict)
    with open(pathFile, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for i, line in enumerate(reader):
            if i != 0:
                nsamples[line[0]].append(line[3])
                genotypes = line[4:]
                freqs = [[x.split(':')[0],float(x.split(':')[1])] for x in genotypes]
                MAF[line[0]][line[1]] = dict()
                for i in freqs:
                    if str(i[1]) == 'nan':
                        try:
                            del MAF[line[0]][line[1]]
                        except:
                            pass
                    else:
                        MAF[line[0]][line[1]][i[0]] = i[1]
    return MAF

class AlleleP:
    def __init__(self,allele, freqs, priorP=None):
        self.allele = allele
        self.freqs = [freq[self.allele] for freq in freqs]
        
        if priorP != None:
            self.priorP = priorP
        else:
            self.priorP = [1.0/len(freqs)] * len(freqs)
        self.P = self.getAlleleP()
        self.backgroundProbabilities = list()
        self.getBackgroundProbabilities()

        
    def getAlleleP(self):
        P = 0
        for freq, priorP in zip(self.freqs, self.priorP):
            P += freq*priorP    
        return P
    
    def getGroupPgivenAllele(self, freq, priorP):
        try:
            P = freq*priorP/self.P
        except ZeroDivisionError:
            P = 0
        return P
    def getBackgroundProbabilities(self):
        for freq, priorP in zip(self.freqs, self.priorP):
            Pbackground = self.getGroupPgivenAllele(freq, priorP)
            self.backgroundProbabilities.append(Pbackground)
        
        
def joinSnpsPos(keys):
    joinpositions = list()
    for positions in keys:
        joinpositions.extend(positions)
    out = list(set(joinpositions))
    return out

def initializeBackgroundP(freqs):
    groups = len(freqs)
    probs = defaultdict(dict)
    #chromosomes
    for chrom in freqs[0].keys():
        positions = [x[chrom].keys() for x in freqs]
        joinPositions = joinSnpsPos(positions)
        #snps
        for snp in joinPositions:            
            alleles = freqs[0][chrom][snp].keys()
            probs[chrom][snp] = dict()
            for allele in alleles:
                try:
                    ifrqs = [x[chrom][snp] for x in freqs]
                    if len(ifrqs) == groups:
                        probs[chrom][snp][allele] = AlleleP(allele, ifrqs)
                except:
                    try:
                        del probs[chrom][snp]
                    except:
                        pass
    return probs

def prettyPrint(data):
    alleles = list(data.keys())
    freqs = [data[a].freqs for a in alleles]
    backgroundProbabilities = [data[a].backgroundProbabilities for a in alleles]
    marginal = [data[a].P for a in alleles]
    probs = pd.DataFrame()
    for n, allele in enumerate(alleles):
        probs['f '+allele] = freqs[n]
        probs['p '+allele+'| Pi'] = backgroundProbabilities[n]
    marginal.extend(['-','-'])
    probs.loc[-1] = marginal
    probs.index = ['Vulgaris', 'Acutifolius', 'P(Ai)']
    return probs

def getSimillarityMatrix(df, sample, backgrounds):
    similarity = pd.DataFrame(df[['CHROM', 'POS']])
    for haplotype in backgrounds:
        if haplotype != sample:
            similarity[haplotype] = (df[sample]==df[haplotype])
            similarity[haplotype].loc[df[haplotype].isnull()] = np.nan
            similarity[haplotype].loc[df[sample].isnull()] = np.nan
    similarity = similarity.replace(0,-1)
    similarity['GENPOS'] = df['GENPOS']
    return similarity


def changeScale(x):
    return 3*(x-(1/3.0))
   

def getWindowBackgroundInformativeness(window, sample):
    score = 0 
    for i, row in window.iterrows():
        if not pd.isnull(row[sample]):
            allele = row[sample].split('/')[0]
            try:
                backgroundProbs = outP[row['CHROM']][str(row['POS'])][allele].backgroundProbabilities
                for backgroundP in backgroundProbs[1:]:
                    score += backgroundP
            except:
                pass

    return score

def getWindowScores(window, sample):
    scores = list()
    unshared = 0
    for i, row in window.iterrows():
        if not pd.isnull(row[sample]):
            allele = row[sample].split('/')[0]
            try:
                backgroundProbs = np.array(outP[row['CHROM']][str(row['POS'])][allele].backgroundProbabilities)
                
                scores.append(backgroundProbs)
                
                
            except Exception as e:
                unshared += 1
    
    backgroundScores = np.array(scores).sum(axis=0)
    n = sum(backgroundScores)
    return backgroundScores/(n-unshared)


def sumMaximumSlidingWindow(df, step, threshold,sample):
    nSNP = df.shape[0]
    name = sample
    chrm = df.iloc[0].CHROM
    dfWindows = pd.DataFrame(columns=['CHROM','n','nSNP','p0','p1','Md_p','GENPOS','POS','pi', 'pf','pr1', 'pr2'])
    p0 = 0
    p1 = threshold 
    n = 0
    flag = False
    #print("%s: %s"%(chrm,nSNP))
    while  p1 < nSNP :
        p = 0
        while p <= threshold:
            if p1 < nSNP:
                wSampleAlleles =  df.iloc[p0:p1]
                try:
                    p = getWindowBackgroundInformativeness(wSampleAlleles, sample)
                
                except Exception  as e: 
                    #print(e)
                    pass
                p1 += step
            else:
                flag = True
                break
        pr1,pr2 = getWindowScores(wSampleAlleles, sample)
        md_p= (p0+ p1)/2.0
        avGenPos =  wSampleAlleles['GENPOS'].mean()
        avPos =  wSampleAlleles['POS'].mean()
        nsnp = p1 - p0
        dfWindows.loc[len(dfWindows)] =  [chrm, n,nsnp,p0,p1,md_p,avGenPos,avPos,wSampleAlleles['POS'].min(),wSampleAlleles['POS'].max(), pr1, pr2]
        p0 += step
        p1 = p0 + threshold
        n += 1 #id ventana
        if flag:
            break
        
    return dfWindows
    
def getData(samples, vcfFile, genposFile):
    fields = ['CHROM', 'POS']
    fields.extend(samples)
    df = pd.read_csv(vcfFile, compression='gzip', usecols=fields, sep='\t', skiprows=1, na_values=['./.'], index_col=False)
    genPos = pd.read_csv(genposFile)
#    genPos.columns = ['CHROM', 'POS', 'GENPOS']
    genPos['id'] = df.index
    df['id'] = df.index
    df = pd.merge(df, genPos[['GENPOS', 'id']], on='id')
    return df



def MaximumSlidingWindows(samples, df, threshold, step):
    columns = samples
    itime = time.time()
    similarities = list()
    for chrom, chData in df.groupby('CHROM'):
        bychmom = pd.DataFrame()
        print(chrom)
        chData = chData.dropna()  
        simChom = sumMaximumSlidingWindow(chData,step,threshold,samples[0])
        similarities.append(simChom)
    print(time.time() - itime)
    return(pd.concat(similarities))

def simulateWindows(df,sample,p0, p1):
    nSNP = df.shape[0]
    name = sample
    chrm = df.iloc[0].CHROM
    dfWindows = pd.DataFrame(columns=['CHROM','n','informativeness','p0','p1','Md_p','GENPOS','POS','pi', 'pf','pr1', 'pr2', 'pr3'])
    wSampleAlleles =  df.iloc[p0:p1]
    dontShared = 0
    try:
        p = getWindowBackgroundInformativeness(wSampleAlleles, sample)

    except Exception  as e: 
        print(e)
        dontShared += 1
        pass
    
    pr1,pr2,pr3 = getWindowScores(wSampleAlleles, sample)
    md_p= (p0+ p1-1)/2.0
    avGenPos =  wSampleAlleles['GENPOS'].mean()
    avPos =  wSampleAlleles['POS'].mean()
    nsnp =int(pr1 + pr2 + pr3)
    dfWindows.loc[len(dfWindows)] =  [chrm, (p1-p0)-dontShared,p,p0,p1-1,md_p,avGenPos,avPos,wSampleAlleles['POS'].min(),wSampleAlleles['POS'].max(), pr1, pr2, pr3]
    return dfWindows


WD='./../..'
file='GCDT_plates_fix_Vulgaris_v2_Bi_Q40_Dp3_imiss83_NS99_MAF2_He2_annotated'

acutFrq = getFreqs('%s/data/custom_script/%s_pop_parents_acut.freq' % (WD, file))
vulFrq = getFreqs('%s/data/custom_script/%s_pop_parents_vul.freq' % (WD, file))
outP = initializeBackgroundP([vulFrq,acutFrq])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Get the unknown origin regions of a child genotype based in the two parents genotype.')
    parser.add_argument('--s', nargs='+')
    parser.add_argument("--vcf", help="path to vcf")
    parser.add_argument("--genpos", help="path to genpos file")
    parser.add_argument("--out", help="path to genpos file")
    args = parser.parse_args()
    print(args)
    samples = args.s
    df = getData(samples,args.vcf, args.genpos)
    maxW = MaximumSlidingWindows(samples, df, 10, 10)
    maxW.to_csv(args.out+'_'.join(samples)+'.csv', sep='\t', index=False)
