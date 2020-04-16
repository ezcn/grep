#!/usr/bin/python3
import pandas as pd
import glob, argparse, re , random 


def selectFiles(mydir, listID):
    """ choose files if ID in list id is in the filename""" 
    fileList=[]
    for ind in listID: 
        fileList += [f for f in glob.glob(mydir) if re.search(ind, f) ]
    return fileList

def mergeThat(fileList):
    df = pd.DataFrame()
    for f in fileList:
        tmp = pd.read_table(f, index_col='index_x')
        df = df.append(tmp)
    return df


def main(): 
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="path to  input dir  ",type=str,required=True)
    parser.add_argument("-l", help="list of id ",type=str,required=True)
    parser.add_argument("-n", help="number of individual to sample ",type=int,required=True)
    parser.add_argument("-i", help="number of iterations ",type=int,required=True) 
    args = parser.parse_args()

    listID = [line.rstrip('\n') for line in open(args.l)]
    sampleToConsider=random.sample(listID, args.n)
    #print (sampleToConsider)
    
    fileToConsider=selectFiles(args.f, sampleToConsider)
    myd=mergeThat(fileToConsider)
    #myd.to_csv("cicci", sep="\t", index=True)
    cycle=0
    while cycle < args.i:
        cycle+=1

        ####### 1. remove genes with too many variants 
        #number of variants/impact per gene  TO DO: plot for stats 
        #myd[['index_x' , 'impact', 'gene_symbol']].drop_duplicates().groupby(['gene_symbol', 'impact']).count()
    
        #number of variants per gene 
        #count number of variants per gene, make adictionary at the first cycle, update dictionary at other cycles 
        if cycle==1: 
            variantsPerGene=myd[['index_x' , 'gene_symbol']].drop_duplicates().groupby(['gene_symbol']).count()#.to_dict()
        else
            tmp=myd[['index_x' , 'gene_symbol']].drop_duplicates().groupby(['gene_symbol']).count()
            variantsPerGene.merge(tmp, 'outer')



    #genes per sample 
    #myd[[ 'gene_symbol', 'sample']].drop_duplicates().groupby(['gene_symbol','sample', 'impact']).count()

  # una volta fatto per una volta mettiamo ciclo che fa per args.i volte, ma questo lo vediamo dopo 

if __name__ == "__main__": 
    main() 
