import argparse
import time 
import numpy as np
import pandas as pd
import csv
from collections import defaultdict
from functools import reduce
import random 

import csv
import numpy as np
import matplotlib.pyplot as plt

def getFreqs(pathFile):
    nsamples = defaultdict(list)
    #Dict for store for each snp present in the population its frequencies
    MAF = defaultdict(dict)
    with open(pathFile, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        #iterate throught each SNP in frequencies file
        for i, line in enumerate(reader):
            if i != 0:
                nsamples[line[0]].append(line[3])
                genotypes = line[4:]
                freqs = [[x.split(':')[0],float(x.split(':')[1])] for x in genotypes]
                MAF[line[0]][int(line[1])] = dict()
                for i in freqs:
                    if str(i[1]) == 'nan':
                        try:
                            del MAF[line[0]][int(line[1])]
                        except:
                            pass
                    else:
                        MAF[line[0]][int(line[1])][i[0]] = i[1]
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
        
        

def initializeBackgroundP(freqs):
    groups = len(freqs)
    probs = defaultdict(dict)
    #chromosomes
    print("Computing probabilities for common SNPs between populations...")
    for chrom in freqs[0].keys():
        #get intersection of snps for each chrom between populations
        snpList = [list(freqPop[chrom].keys()) for freqPop in freqs]
        commonSNPs = reduce(np.intersect1d,snpList)
        #sizes by pop 
        sizes = [len(freqPop[chrom].keys()) for freqPop in freqs]
        print("%s: %s common SNPs %s" %(chrom, sizes, len(commonSNPs)))
        #snps
        for snp in commonSNPs:            
            alleles = freqs[0][chrom][snp].keys()
            probs[chrom][snp] = dict()
            for allele in alleles:
                ifrqs = [x[chrom][snp] for x in freqs]
                if len(ifrqs) == groups:
                    probs[chrom][snp][allele] = AlleleP(allele, ifrqs)
                else:
                    print(chrom, snp)
                
    return probs




def getWindowBackgroundInformativeness(window, sample):
    score = 0 
    for i, row in window.iterrows():
        if not pd.isnull(row[sample]):
            allele = row[sample].split('/')[0]
            try:
                backgroundProbs = outP[row['CHROM']][row['POS']][allele].backgroundProbabilities
                for backgroundP in backgroundProbs[1:]:
                    score += backgroundP
            except KeyError:
                pass
                #multiallelic snps caution
    return score

def getWindowScores(window, sample):
    scores = list()
    unshared = 0
    for i, row in window.iterrows():
        allele = row[sample].split('/')[0]
        try:
            backgroundProbs = np.array(outP[row['CHROM']][row['POS']][allele].backgroundProbabilities)
            scores.append(backgroundProbs)
        except KeyError:
            pass
            #multiallelic snps caution
            #print(row[['CHROM', 'POS']])
    
    backgroundScores = np.array(scores).sum(axis=0)
    return backgroundScores/backgroundScores.sum()


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
        nsnp = wSampleAlleles.shape[0]
        dfWindows.loc[len(dfWindows)] =  [chrm, n,nsnp,p0,p1,md_p,avGenPos,avPos,wSampleAlleles['POS'].min(),wSampleAlleles['POS'].max(), pr1, pr2]
        p0 += step
        p1 = p0 + threshold
        n += 1 #id ventana
        if flag:
            break
        
    return dfWindows
    
def getData(targetSamples,samplesFile, vcfFile, genposFile):
    #get sample index befor read vcf for load just certain columns
    samples = pd.read_csv(samplesFile, header=None)
    samples.rename(columns={0:'Genotype'}, inplace = True)
    samples['indexPos'] = samples.index + 2
    samplesIndex = list(samples[samples['Genotype'].isin(targetSamples)].indexPos)
    fields = [0, 1]
    fields.extend(samplesIndex)
    df = pd.read_csv(vcfFile, usecols=fields, sep='\t', na_values=['./.'], index_col=False, header = None, comment='#')
    colnames = ['CHROM', 'POS']
    colnames.extend(targetSamples)
    df.columns = colnames
    genPos = pd.read_csv(genposFile)
    genPos.rename(columns = ({'Chromosome':'CHROM', 'Position_bp': 'POS', 'Predicted_cM':'GENPOS'}), inplace=True)
    #    genPos.columns = ['CHROM', 'POS', 'GENPOS']
    genPos['id'] = df.index
    df['id'] = df.index
    df = pd.merge(df, genPos[['GENPOS', 'id']], on='id')
    
    return df



def MaximumSlidingWindows(samples, df, outP, threshold, step):
    columns = samples
    itime = time.time()
    similarities = list()
    for chrom, chData in df.groupby('CHROM'):
        bychmom = pd.DataFrame()
        comparableDf = chData[chData.POS.isin(list(outP[chrom].keys()))]
        initialSize = comparableDf.shape[0]
        comparableDf = comparableDf.dropna()  
        print("%s %s Not Genotyped SNPs: %s of %s" % (samples[0],chrom, (initialSize - comparableDf.shape[0]), initialSize))
        simChom = sumMaximumSlidingWindow(comparableDf,step,threshold,samples[0])
        similarities.append(simChom)
    print(time.time() - itime)
    return(pd.concat(similarities))


WD='./../..'
file='GCDT_plates_fix_Vulgaris_v2_Bi_Q40_Dp3_imiss83_NS99_MAF2_He2_annotated'

acutFrq = getFreqs('%s/data/custom_script/%s_pop_parents_acut.frq' % (WD, file))
vulFrq = getFreqs('%s/data/custom_script/%s_pop_parents_vul.frq' % (WD, file))
outP = initializeBackgroundP([vulFrq,acutFrq])
samplesFile = WD+'/data/custom_script/'+'samples.txt'


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Get the unknown origin regions of a child genotype based in the two parents genotype.')
    parser.add_argument('--s', nargs='+')
    parser.add_argument("--vcf", help="path to vcf")
    parser.add_argument("--genpos", help="path to genpos file")
    parser.add_argument("--out", help="path to genpos file")
    args = parser.parse_args()
    print(args)
    samples = args.s
    df = getData(samples,samplesFile, args.vcf, args.genpos)
    maxW = MaximumSlidingWindows(samples, df, outP, 10, 10)
    maxW.to_csv(args.out+"Intro_"+samples[0]+'.csv', sep='\t', index=False)
