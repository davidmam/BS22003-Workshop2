# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 10:03:50 2016

@author: dmamartin
"""
import random
startpos=117479963
cdsstart=133
exons=[x.split('..') for x in "1..185,24291..24401,29072..29180,50937..51152,                     54314..54403,55286..55449,56586..56711,60138..60384,                     62054..62146,68679..68861,79502..79693,107777..107871,                     110391..110477,111972..112695,114968..115096,                     122864..122901,123570..123820,126712..126791,                     130557..130707,131619..131846,134651..134751,                     147560..147808,162476..162631,172880..172969,                    184726..184898,185497..185602,186946..188703".replace(' ','').split(',')]

exonlist=[]
for e in exons:
    exonlist.append([int(x) for x in e])
exons=[]
ecount=0
lensofar=0
for e in exonlist:
    ecount=ecount+1
    exon={'exon': ecount, 'Start':e[0], 'End':e[1], 'GenomeStart':e[0]+startpos-1, 'length':e[1]+1-e[0], 'lengthsofar':lensofar}
    #if lensofar==0:
    #    lensofar=-132
    lensofar+=e[1]+1-e[0]
    exons.append(exon)
    
efh=open('exons.txt','w')
efh.write('exon_frame1\tff0000\nexon_frame2\t00ff00\nexon_frame3\t0000ff\n')

for e in exons:
    cds=e['lengthsofar']-cdsstart+1
    frame= ((3-cds%3)%3)+1
    efh.write('\t'.join(['exon %s'%e['exon'], 'CFTR_reference', '-1','%s'%(e['lengthsofar']+1),'%s'%(e['lengthsofar']+e['length']), 'exon_frame%s'%frame ])+'\n')
efh.close()    

synmut=[]
fh=open('cftr_snp_syn.vcf')
for p in fh.readlines():
    if p[0]=='#':
        continue
    synmut.append(p.strip().split())

synonymous=[]
for m in synmut:
    syn={'Genomic':int(m[1]), 'reference':m[3], 'type':'SYN','variants':m[4].split(',')}
    ex=0
    while int(m[1]) >exons[ex]['GenomeStart']+exons[ex]['length']:
        ex+=1
    if int(m[1])<exons[ex]['GenomeStart']:
        continue
    if int(m[1]) >exons[ex]['GenomeStart']+exons[ex]['length']:
        continue
    syn['mRNA']=int(m[1])+1-exons[ex]['GenomeStart']+exons[ex]['lengthsofar']
    synonymous.append(syn)

mrna=''.join(open('cftr_mrna.fasta').readlines()[1:]).replace('\n','')

def normalmrna(prob=100):
    seq=mrna
    for p in synonymous:
       for v in p['variants']:
           if random.randint(0,prob)==int(prob/2):
               seq=seq[:p['mRNA']-1]+v+seq[p['mRNA']:]
    return seq
    
def formatseq(acc, seq, seqwidth=60):
    text='>%s\n'%acc
    pos=0
    while pos<len(seq)-seqwidth:
        text=text+seq[pos:pos+seqwidth]+'\n'
        pos+=seqwidth
    text=text+seq[pos:]+'\n'
    return text

nsm=[]
fh=open('cftr_snp_nsm.vcf')
for p in fh.readlines():
    if p[0]=='#':
        continue
    nsm.append(p.strip().split())
nonsynonymous=[]
for m in nsm:
    mut={'Genomic':int(m[1]), 'reference':m[3], 'type':'NSM','variants':m[4].split(',')}
    ex=0
    while int(m[1]) >exons[ex]['GenomeStart']+exons[ex]['length']:
        ex+=1
    if int(m[1])<exons[ex]['GenomeStart']:
        continue
    if int(m[1]) >exons[ex]['GenomeStart']+exons[ex]['length']:
        continue
    mut['mRNA']=int(m[1])+1+exons[ex]['lengthsofar']-exons[ex]['GenomeStart']
    nonsynonymous.append(mut)

def patientmrna(syn=2, nonsyn=1, outfile=None):
    muts=random.sample(synonymous,syn)+random.sample(nonsynonymous,nonsyn)
    seq=mrna
    for m in muts:
        var=''
        for v in m['variants']:
            seq=seq[:m['mRNA']-1]+v+seq[m['mRNA']:]
            var='%s%s%s%s'%(m['type'],m['reference'],m['mRNA'],v)
        if outfile:
            outfile.write('\t%s'%var)
    return seq

# now read in a list of matriculation numbers and generate patient data.
fh=open('pateintlist.txt')
fh.readline()
pfh=open('mutations.txt','w')
for p in fh.readlines():
    try:
        paid=int(p.split('\t')[3])
        ofh=open('data/patient_%s.fasta'%paid,'w')
        random.seed(paid)
        pfh.write(str(paid))
        ofh.write(formatseq('patient_%s'%paid,patientmrna(outfile=pfh)))
        pfh.write('\n')
        ofh.close()
    except:
        print('error on line: %s'%p)
pfh.close()