from Bio import Align
from Bio.Seq import Seq
import requests
import pandas as pd
import numpy as np
import re
import streamlit as st 


def do_analysis(alignments, titlewords,withitd):
    array1=[]
    array2=[]    
    for ch2 in alignments[0]._get_row(1):
        array2.append(ch2)
    for ch1 in alignments[0]._get_row(0):
        array1.append(ch1)    


    towrite1=pd.concat((pd.DataFrame(array1).T,pd.DataFrame(array2).T))
    st.write('Alignment '+titlewords)
    st.write(towrite1)
    chart_data1 = towrite1.apply(lambda x: (True,x.iloc[0]==x.iloc[1]))
    #st.write(chart_data1,width=towrite1.shape[1])
    chart_data1=chart_data1.T
    chart_data1.columns=['ref','match']
    chart_data1['minmatch']=0
    chart_data1['minmatch'][np.min(np.where(chart_data1['match'])[0]):]=1
    chart_data1['maxmatch']=0
    chart_data1['maxmatch'][np.max(np.where(chart_data1['match'])[0]):]=1
    chart_data1['insert']=chart_data1.apply(lambda  x: x[2]== 1 and x[3]==0 and x[1]==0,axis=1)
    chart_data1['refseq']=towrite1.apply(lambda x:x.iloc[0])
    chart_data1['varseq']=towrite1.apply(lambda x:x.iloc[1])  

    st.bar_chart(chart_data1,y='ref',color='match',width=towrite1.shape[1])
    st.write('Matching bases between sequence and reference')
    st.write(np.sum(chart_data1['match']))
    if withitd:
        st.write('inserted base number')
        st.write(np.sum(chart_data1['insert']))    
    return chart_data1



# Creating sample sequences
with st.form(key='parameters'):
    reference=st.text_input('coordinates of breakpoint breakpoints of ITD', 'chr13:28034050-28034134',help='Example: chr13:28034050-28034134')
    seq_var=st.text_input('sequence of read to align','ATTTGGCACATTCCATTCTTACCAAACTCTAAATTTTCTCTTGGAAACTCCCATTTGAGATCATATTCATATTCTTGGCCGTGGTGCAGAAACATTTGGCACATTCCATTCTTACCAAACTCTAAATTTTCTCTTGGAAACTCCCATTTG')
    sequence_to_see_flank =st.number_input('flanking sequence to display in alignments', min_value=2, max_value=500, value=200, step=1)
    revc=st.radio('reverse and/or complement the hypothetical duplication',['no','reverse','complement','reversecomplement'])
    ality=alignmenttypec=st.radio('alignment strategy',['local','global'])
    reconstructed=st.checkbox('include alignment of the reconstructed ITD')
    submit_button = st.form_submit_button(label='Submit')
if submit_button:
    chrom=reference.split(':')[0]
    start=reference.split(':')[1].split('-')[0]
    end=reference.split(':')[1].split('-')[1]
    left_flank_start=int(start)-sequence_to_see_flank
    right_flank_end=int(end)+sequence_to_see_flank
    ucsc_variant_seq=requests.get("https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom="+chrom+";start="+start+";end="+end).json()['dna']
    if revc=='reversecomplement:
        itd=ucsc_variant_seq+str(Seq(ucsc_variant_seq).reverse_complement())
    elif revc=='reverse':
        itd=ucsc_variant_seq+str(Seq(ucsc_variant_seq).reverse())
    elif revc=='complement':
        itd=ucsc_variant_seq+str(Seq(ucsc_variant_seq).complement())
    else:    
        itd=ucsc_variant_seq*2
    normal=ucsc_variant_seq
    left_flank = requests.get("https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom="+chrom+";start="+str(left_flank_start)+";end="+str(int(start)-1)).json()['dna']
    right_flank = requests.get("https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom="+chrom+";start="+str(int(end)+1)+";end="+str(right_flank_end)).json()['dna']
    seq_ref=left_flank+itd+right_flank
    seq_normal=left_flank+normal+right_flank

    # Finding similarities
    aligner = Align.PairwiseAligner()
    aligner.mode = ality
    aligner.open_gap_score = -2
    #aligner.extend_gap_score = -0.1
    #aligner.target_end_gap_score = 1.0
    #aligner.query_end_gap_score = -1.0
    alignments = aligner.align(seq_ref, seq_var)
    
    chart_data1=do_analysis(alignments,'query with input sequence',True)

    alignments2=aligner.align(seq_normal, seq_var)
    chart_data2=do_analysis(alignments2,'query without insert',False)
    
    itdbase=len(ucsc_variant_seq)

    itdseq=ucsc_variant_seq+''.join(chart_data1['varseq'][chart_data1['insert']])

    itdseqalignments = aligner.align(seq_ref, itdseq)
    if reconstructed:
        do_analysis(itdseqalignments,'query with ITD',False)
    st.write('ITD')
    st.write(itdseq)
    st.write('ITD Length')
    st.write(itdbase+np.sum(chart_data1['insert']))    
    #seqshow1= [a for a in str(alignments[0].sequences[0])] 
    #seqshow2= [a for a in str(alignments[0].sequences[1])]


