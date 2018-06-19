import pandas as pd
from collections import OrderedDict
import pickle
from difflib import SequenceMatcher
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

class Test:
    
    path = '/Users/advaitbalaji/'
    downloads = 'Downloads/'
    desktop = 'Desktop/'
    
    def createBackbonePrecursor(self, fname):
        
        with open(self.path+self.desktop+fname,'r') as inpf, open(self.path+self.desktop+'MappingPrecursor.txt','w+') as of:
            processes = []
            lines = inpf.readlines()
            for line in lines[1:]:
                if line == '\n':
                    continue
                processes.append(line.strip('-').split('\t')[0]+'\t\n')
            of.write('MBCO\tGO\n')
            of.writelines(processes)
    
    def copyPreviousMappings(self,fname):
        df = pd.read_csv(self.path+self.desktop+fname, sep ='\t', header = 0)
        mbco = df['MBCO'].tolist()
        go = df['GO'].tolist()
        data = pd.read_csv(self.path+self.desktop+'MappingPrecursor.txt', sep = '\t', header = 0)
        data_p = data['MBCO'].tolist()
        for p in mbco:
            if p in data_p:
                data.loc[data_p.index(p),'GO'] = go[mbco.index(p)]
        data.to_csv(self.path+self.desktop+'MappingPrecursor.txt', sep = '\t', header = True, index  = False)
    
    def read(self,destination,fname,delimiter,header):
        if destination.lower() == 'downloads':
            return pd.read_csv(self.path+self.downloads+fname, sep = delimiter , header = header)
        return pd.read_csv(self.path+self.desktop+fname, sep = delimiter , header = header)
    
    def openFileLineWise(self,fname):
        with open(self.path+self.downloads+fname,"r") as f:
            return f.readlines()
    
    def writeDictToFileLineWise(self, d, name):
        with open(self.path+self.desktop+name+'.txt','w') as outf:
            for i, j in d.items():
                outf.write(i+'\t'+' '.join(j)+'\n')
        print("WRITTEN!...")    
    
    def prepareOrderedDict(self,information):
        return{info.split('\t')[0].split('(GO')[0].strip(): [ a.split(',')[0] for a in info.split('\t')[2:-1] ] for info in information}
        
                 
    
    def orderDictionaryWithKey(self,d):
        return OrderedDict(sorted(d.items(), key = lambda t: t[0].lower())) 


    def pickleDictToBinary(self,destination,name,d): 
        print("PICKLED!...")
        if destination.lower() == 'desktop':
            with open(self.path+self.desktop+name+'.pickle',"wb") as wf:
                pickle.dump(d, wf, protocol=pickle.HIGHEST_PROTOCOL)
        else:
            with open(self.path+self.downloads+name+'.pickle',"wb") as wf:
                pickle.dump(d, wf, protocol=pickle.HIGHEST_PROTOCOL)
        
            
    def readPickledDictFromBinary(self,destination,name):
        print("UNPICKLED!...")
        if destination.lower() == 'desktop':
            with open(self.path+self.desktop+name+'.pickle',"rb") as rf:
                return pickle.load(rf)
        else:
            with open(self.path+self.downloads+name+'.pickle',"rb") as rf:
                return pickle.load(rf)
        

    def similar(self, a, b):
        print(SequenceMatcher(None, a, b).ratio())
        
    def optionFoRIU(self,option):
        if option.lower() == 'every':
            return False,True
        elif option.lower() == 'atleastone':
            return True,False
        else:
            return True,True
    
    def createMappingFromBackbone(self,directory,inputfile):
        with open(self.path+directory+inputfile,'r') as inpf,open(self.path+self.desktop+"Mapping_preliminary.txt",'w') as of:
            of.write("MBCO\tGO\n")
            lines = inpf.readlines()
            for line in lines:
                line = line.strip('-')
                try:
                    a,b = line.split('\t')[0],line.split('\t')[1]
                except IndexError:
                    continue
                of.write(a+'\t'+b)    
            print("PRELIMINARY MAPPING CREATED SUCCESSFULLY! CHECK FOR 3 COLUMNS, MAY NOT BE SEPARATED!...")
            print("USE CLEAN NaN TO AVOID Nan ERROR AND GENERATE FINAL MAPPING!...")
            
    def createFinalMappingWithoutNaN(self, directory, inputfile):
        df = pd.read_csv(self.path+directory+inputfile, sep = '\t', header = 0)
        df.dropna(axis = 0,inplace  = True)
        df = df.reset_index(drop=True)
        df.to_csv(self.path+self.desktop+'Mapping.txt', sep = '\t', header = True, index = False)
    
    def prepareGODictFromJensHansenFormat(self,directory,inputfile):
        df = pd.read_csv(self.path+directory+inputfile, sep = '\t', header = 0)
        df = df[['Source_name','Target']]
        process_names = df['Source_name'].unique().tolist()
        with open(self.path+directory+"GO_dict.txt","w") as of:
            of.write("Process\tGenes\n")
            for process in process_names:
                genes = " ".join(df['Target'][df['Source_name'] == process].unique().tolist())
                of.write(process+"\t"+genes+"\n")
                print("WRITTEN PROCESS "+process+"!...")
        print("CREATED GO_DICT FROM HANSEN FORMAT!...")
    
    def jaccardAnalysis(self,t):
        df = t.read('desktop','Mapping.txt','\t',0)
        mbcP = df.MBCO.tolist() 
        goP = df.GO.tolist()
        ref_df = t.read('desktop','GO_dict.txt','\t',0)
        data = t.read('downloads','Validation_results_summary_combined.jns','\t',0)
        data.drop(['Description','ReadWrite_identified_terms','ReadWrite_mbc_versions'], axis = 1, inplace = True)
    
        true_positive = ['True_positive','Accepted_positive']
        false_positive = ['False_positive','Denied_positive','Misinterpreted_term','Sibling']
        
        with open('/Users/advaitbalaji/Desktop/JaccardMap.txt','w') as of:
            of.write('MBCGO\tJI\n')
            for process in mbcP:
                genes = data['Symbol'][data['Subcellular_process'] == process].tolist()
                validations = data['Validation'][data['Subcellular_process'] == process].tolist()
                mbcstatus = []
                for v in validations:
                    if v in true_positive:
                        mbcstatus.append('True positive')
                    else:
                        mbcstatus.append('False postitive')
                genes = [gene for gene in genes if mbcstatus[genes.index(gene)] == 'True positive']
                check_process = goP[mbcP.index(process)].split(';')
                for p in check_process:
                    gene_list = "".join(ref_df['Genes'][ref_df.Process == p].tolist()).split(' ')
                    gene_list[-1] = gene_list[-1].split('\n')[0]
                    of.write(process+'::'+p+'\t'+str(len(list(set(genes) & set(gene_list)))/ len(list(set(genes) | set(gene_list))))+'\n')
        
        
        #t.prepareLevel21correlation()
        #t.prepareLevel32correlation()
        
    def modJaccardAnalysis(Self,t):
        df = t.read('desktop','Mapping.txt','\t',0)
        mbcP = df.MBCO.tolist() 
        goP = df.GO.tolist()
        ref_df = t.read('desktop','GO_dict.txt','\t',0)
        data = t.read('desktop','Supplementary Table S32 - gene-SCP associations_Gene_symbol.txt','\t',0)
        data = data[['ProcessName','Symbol']]
        process_flag = data['ProcessName'].unique().tolist()
        #print(data.head(10))
        with open('/Users/advaitbalaji/Desktop/JaccardMap.txt','w') as of:
            of.write('MBCGO\tJI\n')
            for process in mbcP:
                if process in process_flag:
                    genes = data[data['ProcessName'] == process]['Symbol'].unique().tolist()
                    check_process = goP[mbcP.index(process)].split(';')
                    for p in check_process:
                        gene_list = "".join(ref_df['Genes'][ref_df.Process == p].tolist()).split(' ')
                        gene_list[-1] = gene_list[-1].split('\n')[0]
                        of.write(process+'::'+p+'\t'+str(len(list(set(genes) & set(gene_list)))/ len(list(set(genes) | set(gene_list))))+'\n')
        
    def splitMaptoChildren(self,t,inpf):
        df_32 = t.read('desktop','Level32.csv','\t',0)
        process_curr = df_32['Process'].tolist()
        parent_p = df_32['Parent_process'].tolist()
        map2 = []
        map3 = []
        with open('/Users/advaitbalaji/Desktop/'+inpf,'r') as f, open('/Users/advaitbalaji/Desktop/JaccardMap2.txt','w+') as of2, open('/Users/advaitbalaji/Desktop/JaccardMap3.txt','w+') as of3:
            lines = f.readlines()
            for line in lines[1:]:
                a,b = line.split('\t')
                process = a.split('::')[0]
                if process in process_curr:
                    map3.append(line)
                elif process in parent_p:
                    map2.append(line)
                else:
                    continue
            of2.writelines(map2)
            of3.writelines(map3)
                    
                    
                
        
    def createLevel2ProcessMappedJI(self,t):
        df_32 = t.read('desktop','Level32.csv','\t',0)
        process_curr = df_32['Process'].tolist()
        parent_p = df_32['Parent_process'].tolist()
        with open('/Users/advaitbalaji/Desktop/JaccardMap.txt','r') as inpf,open('/Users/advaitbalaji/Desktop/JaccardMapParentL2.txt','w') as outpf:
            outpf.write('Level2ProcessMBCO\tJI\n')
            lines = inpf.readlines()
            for line in lines[1:]:
                a,b = line.split('\t')
                process = a.split('::')[0]
                if process in process_curr:
                    level2p = parent_p[process_curr.index(process)]
                elif process in parent_p:
                    level2p = process
                else:
                    continue
                outpf.write(level2p+'\t'+b)
            
            
            
    def createLevel1ProcessMappedJI(self,t,inp):
        df_32 = t.read('desktop','Level32.csv','\t',0)
        df_21 = t.read('desktop','Level21.csv','\t',0)
        number = inp[-5]
        print(number)    
        with open('/Users/advaitbalaji/Desktop/'+inp,'r') as inpf,open('/Users/advaitbalaji/Desktop/JaccardMapParent'+number+'.txt','w') as outpf:
            outpf.write('Level1ProcessMBCO\tJI\n')
            lines = inpf.readlines()
            for line in lines[1:]:
                a,b = line.split('\t')
                process = a.split('::')[0]
                #print(process)
                level1p = "".join(df_21['Parent_process'][df_21['Process'] == "".join(df_32['Parent_process'][df_32['Process'] == process].tolist()) ].tolist())
                if level1p == "":
                    level1p = "".join(df_21['Parent_process'][df_21['Process'] == process].tolist())
            
                #print(str(lines.index(line))+". "+level1p)
                #print("")
                if not level1p == "":
                    outpf.write(level1p+'\t'+b)
                else:
                    continue
                
        #data_f = t.read('desktop','JaccardMap.txt','\t',0)
        #data_f = data_f[["MBCGO","JI"]][data_f["JI"] != 0]
        #ax = sns.boxplot(x="MBCGO", y="JI", data=data_f, whis=np.inf, orient = 'v')
        #ax = sns.swarmplot(x="MBCGO", y="JI", data=data_f, color=".2", orient = 'v')
        #plt.show()
    
    def createMultipleBoxSwarmPlot(self,inputf1,inputf2):
        data1 = pd.read_csv(self.path+self.desktop+inputf1, sep = '\t', header = 0)
        data2 = pd.read_csv(self.path+self.desktop+inputf2, sep = '\t', header = 0)
        fig, ax = plt.subplots(figsize=(9,9), ncols=1, nrows=2, sharex = True, squeeze = False)
        #plt.subplots_adjust(hspace = 0.7, wspace = 0.7)
        plt.suptitle("MBCO processes mapped against GO biological processes")
        sns.boxplot(x="Level1ProcessMBCO", y="JI", data=data2, whis=np.inf,ax =  ax[0][0])
        sns.swarmplot(x="Level1ProcessMBCO", y="JI", data=data2, color=".2", size=3, ax = ax[0][0])
        ax[0][0].set_xlabel('')
        ax[0][0].set_title('Level3 Children')
        ax[1][0].set_title("Level2 Children")
        sns.boxplot(x="Level1ProcessMBCO", y="JI", data=data1, whis=np.inf,ax =  ax[1][0])
        sns.swarmplot(x="Level1ProcessMBCO", y="JI", data=data1, color=".2", size=3, ax = ax[1][0])
        plt.xticks(rotation=90)
        plt.tick_params(axis='x', which='major', labelsize=7)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95], w_pad = 4, h_pad = 5)
        plt.show()
    
    def createBoxSwarmPlot(self,inputf):
        data = pd.read_csv(self.path+self.desktop+inputf, sep = '\t', header = 0)
        plt.figure(figsize=(13,13))
        plt.title("MBCO processes mapped against GO biological processes")
        ax = sns.boxplot(x="Level2ProcessMBCO", y="JI", data=data, whis=np.inf)
        ax = sns.swarmplot(x="Level2ProcessMBCO", y="JI", data=data, color=".2", size=3)
        #plt.figure(figsize=(20,10))
        plt.xlabel("Level 2 Processes in MBC Ontology")
        plt.ylabel("Jaccard Index")
        plt.xticks(rotation=90)
        plt.tick_params(axis='x', which='major', labelsize=7)
        plt.tight_layout()
        plt.show()
        
    def prepareLevel21correlation(self):
        data = pd.read_csv('/Users/advaitbalaji/Desktop/Mbc_backbone_regular.csv', sep = '\t', header = 0)
        data = data[['Process','Parent_process']][data['Level'] == 2]
        data = data.reset_index(drop=True)
        data.to_csv('/Users/advaitbalaji/Desktop/Level21.csv', sep = '\t', header = True, index  = False)
    
    def prepareLevel32correlation(self):
        data = pd.read_csv('/Users/advaitbalaji/Desktop/Mbc_backbone_regular.csv', sep = '\t', header = 0)
        data = data[['Process','Parent_process']][data['Level'] == 3]
        data = data.reset_index(drop=True)
        data.to_csv('/Users/advaitbalaji/Desktop/Level32.csv', sep = '\t', header = True, index  = False)
            
                
    def createSummary(self):
        df = t.read('desktop','Mapping.txt','\t',0)
        mbcP = df.MBCO.tolist()
        goP = df.GO.tolist()
        ref_df = t.read('desktop','GO_dict.txt','\t',0)
        data = t.read('downloads','Validation_results_summary_combined.jns','\t',0)
        data.drop(['Description','ReadWrite_identified_terms','ReadWrite_mbc_versions'], axis = 1, inplace = True)
    
        true_positive = ['True_positive','Accepted_positive']
        false_positive = ['False_positive','Denied_positive','Misinterpreted_term','Sibling']
    
        #a,b = t.optionForRIU('both')
    
        with open('/Users/advaitbalaji/Desktop/missing_genes_per_MBCprocess.txt','w') as metaf:
            metaf.write('MBCProcess\tMissingGenes\n')
    

    
        with open('/Users/advaitbalaji/Desktop/Summary.txt',"w") as sf:
            sf.write("Gene\tMBCO Process\tMBCStatus\tGOStatus\tGO Process description\tIntersection\tUnion\n")
        
            for process in mbcP:
                genes = data['Symbol'][data['Subcellular_process'] == process].tolist()
                validations = data['Validation'][data['Subcellular_process'] == process].tolist()
                #print(goP[mbcP.index(process)])
                check_process = goP[mbcP.index(process)].split(';')
                        
                description = ['' for i in range(len(check_process))]
                GO_megagenelist = []
                #print("This if the first description ",description)
                #print(genes)
                for gene in genes:
                    #print(gene,type(gene))
                    for p in check_process:
                        gene_list = "".join(ref_df['Genes'][ref_df.Process == p].tolist()).split(' ')
                        gene_list[-1] = gene_list[-1].split('\n')[0]
                        #print(gene in gene_list)
                        GO_megagenelist.extend(gene_list)
                        #if not gene in gene_list:
                            #description[check_process.index(p)] == 'False Positive'
                        if gene in gene_list:
                            description[check_process.index(p)] = 'True Positive'
                        else:
                            description[check_process.index(p)] = 'False Positive'
                    #print(description)
                    #print('')
                    if 'False Positive' in description:
                        if not 'True Positive' in description:
                            gostatus = "Not Present in any mapped GO Process"
                            intersection,union = 'False','False'
                        else:
                            gostatus = "Present in "+str(description.count('True Positive'))+"/"+str(len(description))+" mapped GO Process"
                            intersection,union = 'False','True'
                    else:
                        gostatus = "Present in every mapped GO Process"
                        intersection,union = 'True','True'
                    #print(process,gostatus)              
                    process_descriptions = []
                    for i,j in zip(check_process,description):
                        process_descriptions.append(i+"("+{'False Positive':'Not Present','True Positive': 'Present','': 'False Positive'}[j]+")")
                    process_descriptions = ",".join(process_descriptions)
                    mbcstatus = 'False Positive' if validations[genes.index(gene)] in false_positive else 'True positive'    
                    sf.write(gene+'\t'+process+'\t'+mbcstatus+'\t'+gostatus+'\t'+process_descriptions+'\t'+intersection+'\t'+union+'\n')
            
                with open('/Users/advaitbalaji/Desktop/missing_genes_per_MBCprocess.txt','a+') as metaf:
                   metaf.write(process+'\t'+','.join(list(set([g for g in GO_megagenelist if not g in genes])))+"\n")    


        



if __name__ == '__main__':
    t = Test()
    #lines = t.openFileLineWise('GO_Biological_Process_2017.txt')
    #d = t.prepareOrderedDict(lines)  
    #d = t.orderDictionaryWithKey(d)
    #t.pickleDictToBinary('desktop','pickled_GO',d)
    #d = t.readPickledDictFromBinary('desktop','pickled_GO')
    #df = t.read('downloads','validation_results_summary_combined.jns','\t',header = 0)
    #t.writeDictToFileLineWise(d,'GO_dict')
    
    #t.similar('activation of MAPKK activity', 'Acetylcholine-mediated control of postsynaptic potential')
    
    #t.createMappingFromBackbone('Downloads/','Mbc_backbone_go_processes_m.txt')
    #t.prepareGODictFromJensHansenFormat('Downloads/','gene_association_upgraded_Homo_sapiens_2017June17.jns')
    
    #t.jaccardAnalysis(t)
    #t.createLevel1ProcessMappedJI(t,'JaccardMap2.txt')
    #t.createLevel1ProcessMappedJI(t,'JaccardMap3.txt')
    #t.createLevel2ProcessMappedJI(t)
    #t.createBoxSwarmPlot('JaccardmapParentL2.txt')
    #t.createBackbonePrecursor('Mbc_backbone.txt')
    #t.prepareLevel21correlation()
    #t.prepareLevel32correlation()
    #t.copyPreviousMappings('Mapping.txt')
    #t.createFinalMappingWithoutNaN('Desktop/','MappingPrecursor.txt')
    #t.modJaccardAnalysis(t)
    #t.createMultipleBoxSwarmPlot('JaccardmapParent.txt')
    #t.splitMaptoChildren(t,'JaccardMap.txt')
    t.createMultipleBoxSwarmPlot('JaccardMapParent2.txt','JaccardMapParent3.txt')
                 
'''
#Frequent Words Function
def FrequentWordsWithMismatches( s, k, d ):
    counts = {}
    for i in range(len(s)-k+1):
        for neighbor in neighbors(s[i:i+k],d):
            print(neighbor)
            if neighbor not in counts:
                counts[neighbor] = 0
            counts[neighbor] += 1
    m = max(counts.values())
    return [kmer for kmer in counts if counts[kmer] == m]

#Finding neighbouring kmers
def neighbors( s, d ):
    #print(s, d)
    if d == 0:
        return [s]
    if len(s) == 1:
        return ['A','C','G','T']
    out = []
    for neighbor in neighbors(s[1:],d):
        if hamming(s[1:],neighbor) < d:
            out.extend(['A'+neighbor,'C'+neighbor,'G'+neighbor,'T'+neighbor])
        else:
            out.append(s[0] + neighbor)
    print(out)    
    return out

#Hamming Distance between kmers
def hamming( s, t ):
    return sum([s[i] != t[i] for i in range(len(s))])

#s = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
s = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
k = 4
d = 1
print(' '.join(FrequentWordsWithMismatches(s,k,d)))

'''



'''
sequence = 'CATTCCAGTACTTCATGATGGCGTGAAGA'

skew = []
prev_val = 0
skew.append(prev_val)
for i in sequence:
    if i == 'C':
        prev_val -= 1
        skew.append(prev_val)
    elif i == 'G':
        prev_val += 1
        skew.append(prev_val)
    else:
        skew.append(prev_val)
#print(skew)
indexes = [i for i,j in enumerate(skew) if j == max(skew)]
print(indexes)
'''

'''
prev_val = 0
min_val = 0
indexes = []
indexes.append(0)
for i in range(len(Genome)):
    if Genome[i] == 'C':
        prev_val -= 1
        if prev_val < min_val:
            indexes = []
            indexes.append(i+1)
            min_val = prev_val
        elif prev_val == min_val:
            indexes.append(i+1)
        else:
            pass            
    elif Genome[i] == 'G':
        prev_val += 1
    else:
        pass
print(indexes)
'''







'''
with open("/Users/advaitbalaji/Desktop/EcGenome.txt","r") as f:
    genome = f.read()
    
k, Length , t = 9,500,3
overall_count = []
for i in range(len(genome) - Length + 1):
    count_dict ={}
    mod_genome = genome[i:i+Length]
    for j in range(len(mod_genome) - k + 1):
        if mod_genome[j:j+k] in count_dict:
            count_dict[mod_genome[j:j+k]] += 1
            if count_dict[mod_genome[j:j+k]] >= t:
                overall_count.append(mod_genome[j:j+k])
        else:
            count_dict[mod_genome[j:j+k]] = 1
#print("YES")
print(" ".join(list(set(overall_count))))
'''            


'''
with open("/Users/advaitbalaji/Desktop/VBGenome.txt","r") as f:
    Genome = f.read()

Pattern = "CTTGATCAT"



pat_index = []
for i in range(len(Genome) - len(Pattern) + 1):
    if Genome[i:i+len(Pattern)] == Pattern:
        pat_index.append(str(i))
print(" ".join(pat_index))

'''

'''
Text = 'CGGAGGACTCTAGGTAACGCTTATCAGGTCCATAGGACATTCA'
for k in range(3,4):
    count_dict = {}
    max_count = 0
    for i in range(len(Text) - k + 1):
        if Text[i:i+k] in count_dict:
            count_dict[Text[i:i+k]] += 1
            if max_count < count_dict[Text[i:i+k]]:
                max_count = count_dict[Text[i:i+k]]
        else:
            count_dict[Text[i:i+k]] = 1

    max_kmers = []
    for i, j in count_dict.items():
        if j == max_count:
            max_kmers.append(i)
    print(str(k)+": "+" ".join(max_kmers)+ " "+ str(max_count))
'''