#!/usr/bin/python3
import pandas as pd
import glob, argparse, re , random 

""" python3 scr/grep.control.py -f 'testdata/hgdp/*' -l testdata/hgdp/id.list  -n 4  -i 2  """


def selectFiles(mydir, listID):
    """ choose files if ID in list id is in the filename""" 
    fileList=[]
    for ind in listID: 
        fileList += [f for f in glob.glob(mydir) if re.search(ind, f) ]
    return fileList

def mergeThat(fileList):
    df = pd.DataFrame()
    for f in fileList:
        tmp = pd.read_table(f) 
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
    
    cycle=0
    while cycle < args.i: 
        cycle+=1

        #choose a random sample     
        sampleToConsider=random.sample(listID, args.n)
        print (sampleToConsider)
    
        # build myd data frame with data from  sample choosen
        fileToConsider=selectFiles(args.f, sampleToConsider)
        myd=mergeThat(fileToConsider)
        #myd.to_csv("cicci", sep="\t", index=True)

        ####### 1. remove genes with too many variants 
        #number of variants/impact per gene  TO DO: plot for stats 
        #myd[['index_x' , 'impact', 'gene_symbol']].drop_duplicates().groupby(['gene_symbol', 'impact']).count()
    
        #count number of variants per gene and genes per sample , make a data frame at the first cycle, update data frame  at other cycles 
        if cycle==1: 
            #number of variant per gene 
            variantsPerGene=myd[['index_x' , 'gene_symbol']].drop_duplicates().groupby(['gene_symbol']).count()#.to_dict()
            #genes per sample 
            genesPerSample=myd[[ 'gene_symbol', 'sample']].drop_duplicates().groupby(['gene_symbol']).count()

        else: 
            tmpv=myd[['index_x' , 'gene_symbol']].drop_duplicates().groupby(['gene_symbol'], as_index=False).count()
            variantsPerGene=variantsPerGene.join(tmpv.set_index('gene_symbol'), on='gene_symbol', how='outer',  lsuffix='_variantsPerGene', rsuffix='_tmpv')
            
            tmpg=myd[[ 'gene_symbol', 'sample']].drop_duplicates().groupby(['gene_symbol'], as_index=False).count()
            genesPerSample=genesPerSample.join(tmpg.set_index('gene_symbol'), on='gene_symbol', how='outer', lsuffix='_genesPerSample', rsuffix='_tmpg')
        
    variantsPerGene.to_csv('ciccivar', sep='\t')    
    genesPerSample.to_csv('ciccigene', sep='\t')

if __name__ == "__main__": 
    main() 
