# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 14:04:21 2021

@author: Alexander Kratz
"""

from Bio import SeqIO
from Bio import Seq
import pandas as pd 
import os
import subprocess

os.chdir("..") #Change directory to location of config csv files. Reads should be in a subfolder named "reads"
#%% Definition of subprocesses

def Fastq2Reads(FileName):
    #Read Fastq, return a list of all reads in string format
    Reads=[]
    for record in SeqIO.parse("reads//"+FileName,"fastq"):
        Reads.append(str(record.seq))
    return Reads


def CatCountUniqueReads(SegPrimers,Reads):
    #Categorize and count unique reads
    # This function takes a dataframe of primers paired with their segment numbers
    # and returns a list of dictionaries, where the Nth dictionary in the list
    # takes a read as the key and contains the count for how many times that read appeared

    UniqueReads=[]
    for n in range(0,len(SegPrimers)):
        UniqueReads.append(dict()) # Create a dictionary for each segment
    
    #count=0 #Initialize counter for progress bar
    
    for Read in Reads:
        
        #count=count+1 
        #if(count*10/len(Reads)>1):#If we have progressed 10% of the way through the job
        #    print(".",end="")     #Then print a .
        #    count=0               #And reset out count to 0
        
        for seg in range(0,len(SegPrimers)):                           #For each segment
            if(Read.find( SegPrimers.at[seg,"Primer"].upper() ) == 0): #If the read starts with the primer designated for that segment
                SegmentName=SegPrimers.at[seg,"Segment"]               #Get the name for this segment
                
                if(Read not in UniqueReads[SegmentName].keys() ):             #If this read has not been seen in this Segment before
                    UniqueReads[SegmentName][Read]=0                          #Initialize the count for this read to 0
                
                UniqueReads[SegmentName][Read]=UniqueReads[SegmentName][Read]+1   #Then, increase count by one
                
                break #Since we have found which primer this read goes to, we can break out of the segment-primer loop and go to the next read
    
    return(UniqueReads) #Return the dictionaries of counts of unique reads in each 


def CatCount2Fasta(UniqueReads): #Takes a list of dictionaries of categorized reads and puts them into a fasta

    try:os.mkdir("output//"+Condition+"//Categorized Reads")
    except:pass    
    
    for Segment in range(0,len(UniqueReads)): #For each segment in the list of dictionaries
        SegReads = UniqueReads[Segment]       #Extract th
        
        SortedReads = dict(sorted(SegReads.items(),key=lambda item:item[1],reverse=True))

        FileHandle=open("output/"+Condition+"/Categorized Reads/Segment_"+str(Segment)+".fasta","x")
        SeqRecList=[]
        for key in SortedReads.keys():
            Sequence=Seq.Seq(key)
            SeqRec=SeqIO.SeqRecord(Sequence,"Count: "+str(SortedReads[key]),description="")
            SeqRecList.append(SeqRec)
        SeqIO.write(SeqRecList,FileHandle,"fasta")
        FileHandle.close()

def CallWater(verbose=False):
    #Calls Water to allign categorized reads to the reference
    #Water will output the reads as a fasta which contains paired alignments between the reference sequence and the NGS reads   
    #Call with Verbose=True to print each command as its called          
    try:os.mkdir("output//"+Condition+"//Aligned Reads")
    except:pass
    
    for Segment in range(0,len(SegPrimers)):
        WaterCom='water -asequence=ref3CL.fasta -bsequence="output/'+Condition+'/Categorized Reads/Segment_'+str(Segment)+'.fasta" -gapopen=10 -gapextend=1 -aformat3="fasta" outfile="output//'+Condition+"/Aligned Reads/Segment_"+str(Segment)+'_Aligned.fasta"'
        if(verbose):print(WaterCom)
        subprocess.run(WaterCom,capture_output=True)
        #print("\tSegment "+str(Segment)+" completed") 

def AnalyzeResidueVariants():#Reads aligned reads of each segment sequentially, then puts them into frame  
    for Segment in range(0,len(SegPrimers)):
        Records = list(SeqIO.parse("output//"+Condition+"//Aligned Reads/Segment_"+str(Segment)+"_Aligned.fasta","fasta"))
        RecIter = iter(Records) # create an iterator of the alligned reads
        
        ReferenceAlign=next(RecIter,'')  #Get the first aligned reference from the water output file
        ReadAlign=next(RecIter,'')       #Get the first aligned NGS read from the water output file
        
        
        while(ReferenceAlign):#While we are still getting alignment pairs from the alignment file
            ReadCount=int(ReadAlign.description.split(" ")[1])  #Get the count of how many times this read appeared
            if("-" in ReferenceAlign or "-" in ReadAlign):            # If we have an insertion or deletion
                SumStats["Indel"]=SumStats["Indel"]+ReadCount         # then discard the read, noting it as a indel for summary stats
            
            if not("-" in ReferenceAlign or "-" in ReadAlign):
                
                
                Frame=str(RefSeq.seq).find(str(ReferenceAlign.seq))%3                   #Determine what frame the reference starts at
                StartNucleotide=[0,2,1][Frame]                                          #This converts the Frame calculated above into what nucleotide is the first in frame codon
                EndNucleotide=StartNucleotide+3*int((len(ReadAlign.seq)-Frame)/3)       #This converts the Frame into the location of the last in-frame codon
                
                TranslatedRead=ReadAlign[StartNucleotide:EndNucleotide].translate() #Translate the In Frame nucleotides to an AA string
                
                FirstResidue=RefAA.find(str(ReferenceAlign[StartNucleotide:EndNucleotide].translate().seq)) #Translate the aligned reference sequence and locate it in the reference AA sequence, finding what AA positions our read covers
                
                if("X" in TranslatedRead): #If any amino acids are missing due to sequencing errors, discard the read
                    SumStats["X-containing"]=SumStats["X-containing"]+ReadCount #note down for summary statistics
                    break
                
                for Residue in range(0,len(TranslatedRead)):                                                                                                      #For each residue in the translated, in frame read
                    TranslatedCountDic[Residue+FirstResidue][TranslatedRead[Residue]]=TranslatedCountDic[Residue+FirstResidue][TranslatedRead[Residue]]+ReadCount #Add the number of times this read was seen to the count of how many times the amino acid was seen at that residue
                    if(TranslatedRead[Residue]=='*'):break #If you see a stop codon, do not continue counting further residues
                    
            ReferenceAlign=next(RecIter,'') #Get the next aligned reference from the water output file
            ReadAlign=next(RecIter,'')      #Get the next aligned NGS read from the water output file
    return

#%% Main loop
try:os.mkdir("output")
except:pass

RefSeq = SeqIO.read("ref3CL.fasta","fasta")
RefAA=str(RefSeq.translate().seq)
AllAA=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','*']

Conditions=pd.read_csv("config\Conditions.csv") #Get a Dataframe containing the Conditions from a CSV in the folder
SegPrimers=pd.read_csv("config\SegmentPrimers.csv") #Get a Dataframe containing the primers for each segment 

SummaryStatsDF=pd.DataFrame(columns=("Condition","Total Reads","Indel","X-containing")) #Create a dataframe for summary stats for each run

for Condition in Conditions.Conditions:
    SumStats={}#Create a dictionary to hold temporary summary stats before conversion to dataframe
    SumStats["Condition"]=Condition 
    SumStats["Indel"]=0
    SumStats["X-containing"]=0
        
    try:os.mkdir("output//"+Condition)
    except:pass
    
    #print("\n"+Condition)
    #print("Reading")
    AllReads=Fastq2Reads(Condition+".fastq")
    
    SumStats["Total Reads"]=len(AllReads)
    
    #print("Categorizing",end='')
    UniqueReads=CatCountUniqueReads(SegPrimers, AllReads)
    
    #print("\nSaving To Fasta......")
    CatCount2Fasta(UniqueReads)
    
    #print("Aligning...")
    CallWater()
    
    
    TranslatedCountDic=[]             #Create a list of dictionaries that takes a amino acid and returns a count 
    for i in range(0,len(RefAA)):     #For each residue in the reference AA sequence
        TranslatedCountDic.append({}) #Create the dictionary that will take amino acid and return the count at that residue
        TranslatedCountDic[i]['Residue']=i+1 #Since python numbers from 0, set each residue equal to i+1
        TranslatedCountDic[i]['Reference']=RefAA[i] #Set each reference amino acid equal to the amino acid from the reference sequence
        for AA in AllAA:                #For every Amino acid
            TranslatedCountDic[i][AA]=0 #Initialize the count of each amino acid 'AA' at each residue 'i' to zero
            
    AnalyzeResidueVariants() #Fill the 
    
    VariantsDF=pd.DataFrame(data=TranslatedCountDic,columns=("Residue","Reference")+tuple(AllAA)) #Create a dataframe with the columns Residue, Reference, and then one column for each Amino acid, and the counts for each amino acid at each residue from the list of dictionaries 
    
    SummaryStatsDF=SummaryStatsDF.append(SumStats,ignore_index=True) #Append the summary statistics for this condition to the dataframe that is holding the sumary stats
    
    VariantsDF.to_csv("output\\"+Condition+".csv",index=False) #Output the counts of each amino acid at each residue to an output csv
    
SummaryStatsDF.to_csv("output\Summary_Stats.csv",index=False) #Output the total summary stats for all conditions to an output csv

