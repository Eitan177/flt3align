from Bio import Align
from Bio.Seq import Seq
import requests
import pandas as pd
import numpy as np
import re
import streamlit as st 
# Creating sample sequences
with st.form(key='parameters'):
    reference=st.text_input('coordinates of breakpoint breakpoints of ITD', 'chr13:28034050-28034134',help='Example: chr13:28034050-28034134')
    seq_var=st.text_input('sequence of read to align','ATTTGGCACATTCCATTCTTACCAAACTCTAAATTTTCTCTTGGAAACTCCCATTTGAGATCATATTCATATTCTTGGCCGTGGTGCAGAAACATTTGGCACATTCCATTCTTACCAAACTCTAAATTTTCTCTTGGAAACTCCCATTTG')
    sequence_to_see_flank =st.number_input('flanking sequence to display in alignments', min_value=5, max_value=500, value=200, step=1)
    submit_button = st.form_submit_button(label='Submit')
if submit_button:
    chrom=reference.split(':')[0]
    start=reference.split(':')[1].split('-')[0]
    end=reference.split(':')[1].split('-')[1]
    left_flank_start=int(start)-sequence_to_see_flank
    right_flank_end=int(end)+sequence_to_see_flank
    ucsc_variant_seq=requests.get("https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom="+chrom+";start="+start+";end="+end).json()['dna']
    itd=ucsc_variant_seq*2
    left_flank = requests.get("https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom="+chrom+";start="+str(left_flank_start)+";end="+start).json()['dna']
    right_flank = requests.get("https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom="+chrom+";start="+end+";end="+str(right_flank_end)).json()['dna']
    seq_ref=left_flank+itd+right_flank
    st.write('sequence retrieved')

    # Finding similarities
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -2
    alignments = aligner.align(seq_ref, seq_var)
    st.write('pairwise alignment finished')
    
    #seqshow1= [a for a in str(alignments[0].sequences[0])] 
    #seqshow2= [a for a in str(alignments[0].sequences[1])]
    ref_a=re.search('-',alignments[0]._get_row(0))
    var_a=re.search('-',alignments[0]._get_row(1))

    if ref_a != None and var_a != None:
        rr=ref_a.span()+var_a.span()
        min_val, max_val = np.min(rr), np.max(rr)
        max_val=max_val+np.max((len(re.findall('-',alignments[0]._get_row(0))),len(re.findall('-',alignments[0]._get_row(1)))))
    elif ref_a != None:
        rr=ref_a.span()
        min_val, max_val = np.min(rr), np.max(rr)
        max_val=max_val+len(re.findall('-',alignments[0]._get_row(0)))
    else:
        rr=var_a.span()   
        min_val, max_val = np.min(rr), np.max(rr)
        max_val=max_val+len(re.findall('-',alignments[0]._get_row(1))) 
    array1=[]
    array2=[]

    def rev(it):
        "Reverses an interable and returns it as a string"
        return ''.join(reversed(it))
    
    for ch2 in alignments[0]._get_row(1):
        array2.append(ch2)
    for ch1 in alignments[0]._get_row(0):
        array1.append(ch1)    


    towrite1=pd.concat((pd.DataFrame(array1).T,pd.DataFrame(array2).T))

    st.write(towrite1)
    chart_data1 = towrite1.apply(lambda x: (True,x.iloc[0]==x.iloc[1]))
    st.write(chart_data1,width=towrite1.shape[1])
    chart_data1=chart_data1.T
    chart_data1.columns=['ref','match']
    chart_data1['minmatch']=0
    chart_data1['minmatch'][np.min(np.where(chart_data1['match'])[0]):]=1
    chart_data1['maxmatch']=0
    chart_data1['maxmatch'][np.max(np.where(chart_data1['match'])[0]):]=1
    st.write(chart_data1.apply(lambda  x: x['minmatch']== 1 and x['maxmatch']==0 and x['match']==0))
    st.write(chart_data1)
    itdlen=np.max(np.where(chart_data1['match'])) -np.min(np.where(chart_data1['match']))  - len(ucsc_variant_seq) 
    st.write(itdlen)
    st.write(len(ucsc_variant_seq))
    st.bar_chart(chart_data1,y='ref',color='match',width=towrite1.shape[1])
    
    st.write(alignments[0]._get_row(0))
    st.write(alignments[0]._get_row(1))
