# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
updated to match CFTR exercise
"""

import random

errmsg='''An error occurred in your code. 
Please check the previous cells have all been run correctly and in order'''
msg=' then run the test again'
def test1(v ):
    try:    
        assert v['enzyme']=='CTTAAG'
        print('variable enzyme created correctly')
        assert v['cftr_mrna'].find(v['enzyme'])==992
        print('Enzyme cut site found in mRNA sequence')
    except:
        print(errmsg)

def test1a(v):
    if 'cut1' not in v :
        print("Try using 'cut1 = ' to set the value of cut1")
        return
    if v['cut1']==2225:
        print('cut1 value set correctly')
        if 'cut2' not in v:
            print("Try using 'cut2 = ' to set the value of cut2")
            return    
        if v['cut2']==3713:
            print('Correct answers. cut1 and cut2 set correctly')
        else:
            print('cut2 value incorrect')
    else:
        print('cut1 value incorrect')
        
def test2(v):
    try:
        if v['frag1']==2225:
            print("Fragment 1 correct" )
        else:
            print("Fragment 1 incorrect")
        if v['frag2']==1488 :
            print("Fragment 2 correct" )
        else:
            print("Fragment 2 incorrect")
        if v['frag3']== 2419 :
            print("Fragment 3 correct" )
        else:
            print("Fragment 3 incorrect")
    except:
        print(errmsg)
        
def test4(v):
    ps='aattggaagcaaatgacatcacagcaggtcagagaaaaagggttgagcggcaggcacccagagtagtaggtctttggcattaggagcttgagcccagacggccctagcagggaccccagcgcccgagagaccatgcagaggtcgcctctggaaaaggccagcgttgtctccaaactttttttcagctggaccagaccaattttgaggaaaggatacagacagcgcctggaattgtcagacatataccaaatcccttctgttgattctgctgacaatctatctgaaaaattggaaagagaatgggatagagagctggcttcaaagaaaaatcctaaactcattaatgcccttcggcgatgttttttctggagatttatgttctatggaatctttttatatttaggggaagtcaccaaagcagtacagcctctcttactgggaagaatcatagcttcctatgacccggataacaaggaggaacgctctatcgcgatttatctaggcataggcttatgccttctctttattgtgaggacactgctcctacacccagccatttttggccttcatcacattggaatgcagatgagaatagctatgtttagtttgatttataagaagactttaaagctgtcaagccgtgttctagataaaataagtattggacaacttgttagtctcctttccaacaacctgaacaaatttgatgaaggacttgcattggcacatttcgtgtggatcgctcctttgcaagtggcactcctcatggggctaatctgggagttgttacaggcgtctgccttctgtggacttggtttcctgatagtccttgccctttttcaggctgggctagggagaatgatgatgaagtacagagatcagagagctgggaagatcagtgaaagacttgtgattacctcagaaatgattgaaaatatccaatctgttaaggcatactgctgggaagaagcaatggaaaaaatgattgaaaacttaagacaaacagaactgaaactgactcggaaggcagcctatgtgagatacttcaatagctcagccttcttcttctcagggttctttgtggtgtttttatctgtgcttccctatgcactaatcaaaggaatcatcctccggaaaatattcaccaccatctcattctgcattgttctgcgcatggcggtcactcggcaatttccctgggctgtacaaacatggtatgactctcttggagcaataaacaaaatacaggatttcttacaaaagcaagaatataagacattggaatataacttaacgactacagaagtagtgatggagaatgtaacagccttctgggaggagggatttggggaattatttgagaaagcaaaacaaaacaataacaatagaaaaacttctaatggtgatgacagcctcttcttcagtaatttctcacttcttggtactcctgtcctgaaagatattaatttcaagatagaaagaggacagttgttggcggttgctggatccactggagcaggcaagacttcacttctaatggtgattatgggagaactggagccttcagagggtaaaattaagcacagtggaagaatttcattctgttctcagttttcctggattatgcctggcaccattaaagaaaatatcatctttggtgtttcctatgatgaatatagatacagaagcgtcatcaaagcatgccaactagaagaggacatctccaagtttgcagagaaagacaatatagttcttggagaaggtggaatcacactgagtggaggtcaacgagcaagaatttctttagcaagagcagtatacaaagatgctgatttgtatttattagactctccttttggatacctagatgttttaacagaaaaagaaatatttgaaagctgtgtctgtaaactgatggctaacaaaactaggattttggtcacttctaaaatggaacatttaaagaaagctgacaaaatattaattttgcatgaaggtagcagctatttttatgggacattttcagaactccaaaatctacagccagactttagctcaaaactcatgggatgtgattctttcgaccaatttagtgcagaaagaagaaattcaatcctaactgagaccttacaccgtttctcattagaaggagatgctcctgtctcctggacagaaacaaaaaaacaatcttttaaacagactggagagtttggggaaaaaaggaaGAATTCtattctcaatccaatcaactctatacgaaaattttccattgtgcaaaagactcccttacaaatgaatggcatcgaagaggattctgatgagcctttagagagaaggctgtccttagtaccagattctgagcagggagaggcgatactgcctcgcatcagcgtgatcagcactggccccacgcttcaggcacgaaggaggcagtctgtcctgaacctgatgacacactcagttaaccaaggtcagaacattcaccgaaagacaacagcatccacacgaaaagtgtcactggcccctcaggcaaacttgactgaactggatatatattcaagaaggttatctcaagaaactggcttggaaataagtgaagaaattaacgaagaagacttaaaggagtgcttttttgatgatatggagagcataccagcagtgactacatggaacacataccttcgatatattactgtccacaagagcttaatttttgtgctaatttggtgcttagtaatttttctggcagaggtggctgcttctttggttgtgctgtggctccttggaaacactcctcttcaagacaaagggaatagtactcatagtagaaataacagctatgcagtgattatcaccagcaccagttcgtattatgtgttttacatttacgtgggagtagccgacactttgcttgctatgggattcttcagaggtctaccactggtgcatactctaatcacagtgtcgaaaattttacaccacaaaatgttacattctgttcttcaagcacctatgtcaaccctcaacacgttgaaagcaggtgggattcttaatagattctccaaagatatagcaattttggatgaccttctgcctcttaccatatttgacttcatccagttgttattaattgtgattggagctatagcagttgtcgcagttttacaaccctacatctttgttgcaacagtgccagtgatagtggcttttattatgttgagagcatatttcctccaaacctcacagcaactcaaacaactggaatctgaaggcaggagtccaattttcactcatcttgttacaagcttaaaaggactatggacacttcgtgccttcggacggcagccttactttgaaactctgttccacaaagctctgaatttacatactgccaactggttcttgtacctgtcaacactgcgctggttccaaatgagaatagaaatgatttttgtcatcttcttcattgctgttaccttcatttccattttaacaacaggagaaggagaaggaagagttggtattatcctgactttagccatgaatatcatgagtacattgcagtgggctgtaaactccagcatagatgtggatagcttgatgcgatctgtgagccgagtctttaagttcattgacatgccaacagaaggtaaacctaccaagtcaaccaaaccatacaagaatggccaactctcgaaagttatgattattgaGAATTCacacgtgaagaaagatgacatctggccctcagggggccaaatgactgtcaaagatctcacagcaaaatacacagaaggtggaaatgccatattagagaacatttccttctcaataagtcctggccagagggtgggcctcttgggaagaactggatcagggaagagtactttgttatcagcttttttgagactactgaacactgaaggagaaatccagatcgatggtgtgtcttgggattcaataactttgcaacagtggaggaaagcctttggagtgataccacagaaagtatttattttttctggaacatttagaaaaaacttggatccctatgaacagtggagtgatcaagaaatatggaaagttgcagatgaggttgggctcagatctgtgatagaacagtttcctgggaagcttgactttgtccttgtggatgggggctgtgtcctaagccatggccacaagcagttgatgtgcttggctagatctgttctcagtaaggcgaagatcttgctgcttgatgaacccagtgctcatttggatccagtaacataccaaataattagaagaactctaaaacaagcatttgctgattgcacagtaattctctgtgaacacaggatagaagcaatgctggaatgccaacaatttttggtcatagaagagaacaaagtgcggcagtacgattccatccagaaactgctgaacgagaggagcctcttccggcaagccatcagcccctccgacagggtgaagctctttccccaccggaactcaagcaagtgcaagtctaagccccagattgctgctctgaaagaggagacagaagaagaggtgcaagatacaaggctttagagagcagcataaatgttgacatgggacatttgctcatggaattggagctcgtgggacagtcacctcatggaattggagctcgtggaacagttacctctgcctcagaaaacaaggatgaattaagtttttttttaaaaaagaaacatttggtaaggggaattgaggacactgatatgggtcttgataaatggcttcctggcaatagtcaaattgtgtgaaaggtacttcaaatccttgaagatttaccacttgtgttttgcaagccagattttcctgaaaacccttgccatgtgctagtaattggaaaggcagctctaaatgtcaatcagcctagttgatcagcttattgtctagtgaaactcgttaatttgtagtgttggagaagaactgaaatcatacttcttagggttatgattaagtaatgataactggaaacttcagcggtttatataagcttgtattcctttttctctcctctccccatgatgtttagaaacacaactatattgtttgctaagcattccaactatctcatttccaagcaagtattagaataccacaggaaccacaagactgcacatcaaaatatgccccattcaacatctagtgagcagtcaggaaagagaacttccagatcctggaaatcagggttagtattgtccaggtctaccaaaaatctcaatatttcagataatcacaatacatcccttacctgggaaagggctgttataatctttcacaggggacaggatggttcccttgatgaagaagttgatatgccttttcccaactccagaaagtgacaagctcacagacctttgaactagagtttagctggaaaagtatgttagtgcaaattgtcacaggacagcccttctttccacagaagctccaggtagagggtgtgtaagtagataggccatgggcactgtgggtagacacacatgaagtccaagcatttagatgtataggttgatggtggtatgttttcaggctagatgtatgtacttcatgctgtctacactaagagagaatgagagacacactgaagaagcaccaatcatgaattagttttatatgcttctgttttataattttgtgaagcaaaattttttctctaggaaatatttattttaataatgtttcaaacatatataacaatgctgtattttaaaagaatgattatgaattacatttgtataaaataatttttatatttgaaatattgactttttatggcactagtatttctatgaaatattatgttaaaactgggacaggggagaacctagggtgatattaaccaggggccatgaatcaccttttggtctggagggaagccttggggctgatgcagttgttgcccacagctgtatgattcccagccagcacagcctcttagatgcagttctgaagaagatggtaccaccagtctgactgtttccatcaagggtacactgccttctcaactccaaactgactcttaagaagactgcattatatttattactgtaagaaaatatcacttgtcaataaaatccatacatttgtgtgaaa'
    try:
        prettyseq=v['prettyseq']
    except:
        print(errmsg)
        return
    if prettyseq==ps:
        print('congratulations. Your sequence is suitably pretty')
        return
    if prettyseq[cut1+6:cut1+12]==ecori.lower() or prettyseq[cut2+6:cut2+12]==ecori.lower():
        print("Fail. Don't forget to skip the restriction site when adding the fragment" )
    else:
        print('Wrong answer - get each part with text[start:end+1] and use .lower() as required')
        
def test5(v):
    msg=' then run the test again'
    if 'noncutter' not in v or 'largest_frag' not in v or 'largest_cutter' not in v or 'smallest_frag' not in v or 'smallest_cutter' not in v:
        print('Make sure you have properly defined all five variables'+msg)
        return
    if v['noncutter'].upper()=='GGTCTC' or v['noncutter'].upper()=='BSAI' \
    or v['noncutter'].upper()=='CCCGGG' or v['noncutter'].upper()=='SMAI':
        print('Correct. BsaI and SmaI do not cut.')
    else:
        print('check the results for your non-cutting enzyme'+msg)
        return
        
    if v['largest_frag'] == 6132:
        print('Choose an enzyme that does cut'+msg)
        return
    if v['largest_frag'] ==5644 and v['largest_cutter'].lower()=='mvni':
        print('Correct largest fragment')
    else:
        print('check your results for the largest fragment carefully'+msg)
    if v['smallest_frag'] ==488 and v['smallest_cutter'].lower()=='mvni':
        print('Correct smallest fragment')
    else:
        print('check your results for the smallest fragment carefully'+msg)
        
def test6a(v):
    if v['exon'][:10]=='ATGCAGAGGT' and v['exon'][-10:]=='TTTTTTTCAG':
        print('Exon 1 found')
    elif v['exon'][:10]=='CTGGACCAGA' and v['exon'][-10:]=='AATTGGAAAG':
        print('Exon 2 found')
    elif v['exon'][:10]=='AGAATGGGAT' and v['exon'][-10:]=='ATATTTAGGG':
        print('Exon 3 found')
    elif v['exon'][:10]=='GAAGTCACCA' and v['exon'][-10:]=='TTATAAGAAG':
        print('Exon 4 found')
    elif v['exon'][:10]=='ACTTTAAAGC' and v['exon'][-10:]=='ATTTGATGAA':
        print('Exon 5 found')
    elif v['exon'][:10]=='GGACTTGCAT' and v['exon'][-10:]=='TGAAGTACAG':
        print('Exon 6 found')
    elif v['exon'][:10]=='AGATCAGAGA' and v['exon'][-10:]=='ACTTAAGACA':
        print('Exon 7 found')
    elif v['exon'][:10]=='AACAGAACTG' and v['exon'][-10:]=='CAAAATACAG':
        print('Exon 8 found')
    elif v['exon'][:10]=='GATTTCTTAC' and v['exon'][-10:]=='CTGGGAGGAG':
        print('Exon 9 found')
    elif v['exon'][:10]=='GGATTTGGGG' and v['exon'][-10:]=='AGCAGGCAAG':
        print('Exon 10 found')
    elif v['exon'][:10]=='ACTTCACTTC' and v['exon'][-10:]=='ACTAGAAGAG':
        print('Exon 11 found')
    elif v['exon'][:10]=='GACATCTCCA' and v['exon'][-10:]=='CTTTAGCAAG':
        print('Exon 12 found')
    elif v['exon'][:10]=='AGCAGTATAC' and v['exon'][-10:]=='TATTTGAAAG':
        print('Exon 13 found')
    elif v['exon'][:10]=='CTGTGTCTGT' and v['exon'][-10:]=='AGACTTAAAG':
        print('Exon 14 found')
    elif v['exon'][:10]=='GAGTGCTTTT' and v['exon'][-10:]=='TCTGGCAGAG':
        print('Exon 15 found')
    elif v['exon'][:10]=='GTGGCTGCTT' and v['exon'][-10:]=='TCCTTGGAAA':
        print('Exon 16 found')
    elif v['exon'][:10]=='CACTCCTCTT' and v['exon'][-10:]=='TTGAAAGCAG':
        print('Exon 17 found')
    elif v['exon'][:10]=='GTGGGATTCT' and v['exon'][-10:]=='CTTCATCCAG':
        print('Exon 18 found')
    elif v['exon'][:10]=='TTGTTATTAA' and v['exon'][-10:]=='GAATCTGAAG':
        print('Exon 19 found')
    elif v['exon'][:10]=='GCAGGAGTCC' and v['exon'][-10:]=='TTAACAACAG':
        print('Exon 20 found')
    elif v['exon'][:10]=='GAGAAGGAGA' and v['exon'][-10:]=='GGATAGCTTG':
        print('Exon 21 found')
    elif v['exon'][:10]=='ATGCGATCTG' and v['exon'][-10:]=='TGGCCAGAGG':
        print('Exon 22 found')
    elif v['exon'][:10]=='GTGGGCCTCT' and v['exon'][-10:]=='GATACCACAG':
        print('Exon 23 found')
    elif v['exon'][:10]=='AAAGTATTTA' and v['exon'][-10:]=='TGCAGATGAG':
        print('Exon 24 found')
    elif v['exon'][:10]=='GTTGGGCTCA' and v['exon'][-10:]=='TGGATCCAGT':
        print('Exon 25 found')
    elif v['exon'][:10]=='AACATACCAA' and v['exon'][-10:]=='ACAATTTTTG':
        print('Exon 26 found')
    elif v['exon'][:10]=='GTCATAGAAG' and v['exon'][-10:]=='AAGGCTTTAG':
        print('Exon 27 found')
    else:
        print('limits incorrectly specified. %s'%errmsg)
        
def test6(v):
    if 'variance' not in v:
        print('ensure that your code runs without an error'+msg)
        return
    if v['variance']==0:
        print('ensure your code adds to the sumofsquares variable'+msg)
        return
    if v['variance']==437112:
        print('you should add to the sumofsquares, not replace it. Check your code'+msg)
        return
    if v['variance'] != 873290:
        print(errmsg)
        return
    else:
        print('Correct. The variance has been calculated.')
        
def test7(v):
    if 'lastcut' not in v or 'lastcutter' not in v:
        print('you must run the previous cell without error first. Please fix your code'+msg)
        return
    if v['lastcut']==6068:
        print('Correct. The longest fragment has been successfully found')
    else:
        print(errmsg)

def test8(v):
    if 'singlecutters' not in v or not v['singlecutters']:
        print('check your code to make sure you append something to the list of single cutters'+msg)
        return
    if len(v['singlecutters'])==67:
        print('Correct. There are 92 enzymes in the list that cut exactly once')
        return
    else:
        print('check your code. Especially check that you have indented it properly and the if statement is correct'+msg)
        
def test9(v):
    if '' in v['overlaps']:
        print('Check you have got the smaller value in the first position when findign your overlap in the pattern'+msg)
        return
    if len(v['overlaps'])==42:
        print('congratulations, you have found the overlaps')
        return
    else:
        print(errmsg)
        
def test10(v):        
    if not 'stickies' in v:
        print('Ensure your code runs without errors'+msg)
        return
    try:
        assert len(v['stickies'])==14
        assert len(v['stickies']['Blunt'])==4
    except:
        print(errmsg)
        return
    print('All tests passed - stickies are good to go!')
    
def test11(v):
    gcutters=[]
    for enz in elist:
        if sequence.count(enz['pattern'])==1:
            gcutters.append(enz)
    try:
        assert len(gcutters) == len(v['gcutters'])
        print('Correctly identified the single cutting enzymes')
        assert v['gene']==sequence
        assert v['matricnum']
        print('correctly retrieved your gene sequence with your matriculation number')
    except:
        print('Are you sure you are checking the correct value and adding to the list? '+msg)
        
def test12(v):
    firstcutsite=len(gene)
    lastcutsite=0
    firstcutter='' #name of our first cutting enzyme
    lastcutter='' #name of our last cutting enzyme
    for gc in gcutters:
        if gc['overlap'] in stickies:
            # yes, we found a compatible cutter
            ecut=gene.find(gc['pattern']) # add something to make this the correct postion - where does the enzyme actually cut?
            if ecut < firstcutsite:
                firstcutsite=ecut
                firstcutter=gc
            #add an if statement and block to ckeck if it is the last cutter.
            if ecut > lastcutsite:
                lastcutsite=ecut
                lastcutter=gc
    try:
        assert v['firstcutsite']==firstcutsite
        assert v['firstcutter']['name']==firstcutter['name']
        print('first enzyme found: %s cutting at %s'%(firstcutter['name'],firstcutsite) )
        assert v['firstcutsite']==firstcutsite
        assert v['firstcutsite']==firstcutsite
        print('second enzyme found: %s cutting at %s'%(lastcutter['name'],lastcutsite))
    except:
        print(errmsg)        
        
def test13(v):
    try:
        assert v['firstcutter']
        assert v['lastcutter']
    except:
        print('Ensure you have run the previous code successfully before running this test '+msg)
        return
    try:
        assert v['plasmidcutsite1']
        assert v['plasmidcutter1']
        assert v['plasmidcutsite2']
        assert v['plasmidcutter2']
    except:
        print('ensure you have successfully defined the plasmid cut sites '+msg)
        return
    ps1=[x['pcut'] for x in stickies[firstcutter['overlap']]]    
    pn1=[x['name'] for x in stickies[firstcutter['overlap']]]
    ps2=[x['pcut'] for x in stickies[lastcutter['overlap']]]
    pn2=[x['name'] for x in stickies[lastcutter['overlap']]]
    try:
        assert v['plasmidcutsite1'] in ps1        
        assert v['plasmidcutsite2'] in ps2
        assert v['plasmidcutter1'] in pn1
        assert v['plasmidcutter1'] in pn2
    except:
        print('plasmid cut sites incorrectly defined. '+msg)

def test14(v):
    try:
        assert v['revcomp']
    except:
        print('Ensure you have run the cell to define the revcomp() method if you need it.')
    fragstart=0
    fragend=0
    frag=''
    cloneseq=''
    try:
        if v['plasmidcutsite1']<v['plasmidcutsite2']:
            fragstart=v['firstcutsite']+v['firstcutter']['forcut']-1
            fragend=v['lastcutsite']+v['lastcutter']['forcut']-1
            frag=sequence[fragstart:fragend]
            pstart=v['plasmidcutsite1']+v['plasmidcutter1']['forcut']-1
            pend=v['plasmidcutsite2']+v['plasmidcutter2']['forcut']-1
            cloneseq=pseq[:pstart]+frag+pseq[pend:]
        else:
            fragstart=v['firstcutsite']+v['firstcutter']['revcut']-1
            fragend=v['lastcutsite']+v['lastcutter']['revcut']-1
            frag=revcomp(sequence[fragstart:fragend])
            pend=v['plasmidcutsite1']+v['plasmidcutter1']['forcut']-1
            pstart=v['plasmidcutsite2']+v['plasmidcutter2']['forcut']-1
            cloneseq=pseq[:pstart]+frag+pseq[pend:]
        assert cloneseq==v['clonedsequence']
        print('Congratulations, you have the correctly cloned sequence')
    except:
        feedback='''Not quite there. Things to check:
            1. Have you got your sequence the right way round? Check which sticky end comes first in the plasmid.
            2. Did you take into account the sticky ends? Where does the enzyme actually cut?
            3. If your gene goes in in reverse orientation did you reverse and complement the sequence?
            4. Think about the cut sites - if it is reversed thene the cut sites will be on the reverse strand
            5. Are the sites where you think they are? Check to see the fragments start and stop with the correct sequence.
            '''
        print(feedback)
            
            
def getmysequence(seed):
    globals().update({'matricnum':int(seed)})
    #pass
   
    #seed=12345
    snum=int(seed)
    random.seed(snum)
    global elist
    elist=[]
    global olist
    olist={}
    for e in open('enzymes.txt').readlines():
        elist.append(dict(zip(['name', 'pattern','forcut','revcut'], e.strip().split('\t'))))
        elist[-1]['forcut']=int(elist[-1]['forcut'])
        elist[-1]['revcut']=int(elist[-1]['revcut'])    
        overlap='External'
        if elist[-1]['forcut'] > len(elist[-1]['pattern']) or elist[-1]['revcut'] > len(elist[-1]['pattern']):
            # enzyme cuts outwith the pattern
            overlap='External'
        elif elist[-1]['forcut']==elist[-1]['revcut']:
            #fill this in 
            overlap='Blunt'
        elif elist[-1]['forcut'] < elist[-1]['revcut']: #fill this in and don't forget the : at the end - this should get the 5' overhangs
            overlap=elist[-1]['pattern'][elist[-1]['forcut']:elist[-1]['revcut']].lower() 
        else:
            # this should be the upper case overlap so set overlap to the overlap in upper case. 
            #remember to put the smaller value first when getting the overlap
            overlap=elist[-1]['pattern'][elist[-1]['revcut']:elist[-1]['forcut']]
        elist[-1]['overlap']=overlap      
        try:
            olist[overlap].append(elist[-1])
        except:
            olist[overlap]=elist[-1:]
    global pseq
    pseq=open('plasmid.txt').read()
    global cutters
    cutters=[]
    for e in elist:
        if pseq.count(e['pattern'])==1 and pseq.find(e['pattern'])>650 and pseq.find(e['pattern'])<750 and len(e['pattern'])>=4:
            cutters.append(e)
            e['pcut']=pseq.find(e['pattern'])
    #print(cutters)
    global stickies
    stickies={}
    for x in cutters:
        try:
            stickies[x['overlap']].append(x)
        except:
            stickies[x['overlap']]=[x]    
    if 'External' in stickies:
        del stickies['External']
    if 'Blunt' in stickies:
        del stickies['Blunt']
    #print(stickies)
    seqlen=snum%10000 +10000
    global sequence
    sequence=''
    for s in range(seqlen):
        sequence+=random.sample(['A','C','T','G'],1)[0]
    minpos=len(sequence)
    startcut=False
    startpat=''
    endcut=False
    endpat=''
    usable=0
    for e in elist:
        if sequence.count(e['pattern'])==1 and e['overlap'] in stickies:
            if sequence.find(e['pattern'])<=400:
                startcut=True
                startpat=e['pattern']
            if sequence.find(e['pattern'])>=seqlen-400:
                endcut=True
                endpat=e['pattern']
        if sequence.count(e['pattern'])==0 and e['overlap'] in stickies:
            usable=usable+1
    random.seed(snum)
    pats=random.sample(stickies, 2)    
    if usable==0 and (endcut or startcut):
        pats=pats[:1]
    elif usable>1 or (endcut or startcut):
        pats=[]
    for p in pats:
        sequence.replace(random.sample(olist[p],1)[0]['pattern'],'')
        
        # remove random patterns from the list of potential cloning enzymes.
        # TODO 
    if not startcut:
        if endpat:
            sequence=addsticky(sequence, elist, stickies).replace(endpat,'')+sequence
        else:            
            sequence=addsticky(sequence, elist, stickies)+sequence
    if not endcut:
        if startpat:       
            sequence=sequence+addsticky(sequence,elist,stickies).replace(startpat,'')
        else:            
            sequence=sequence+addsticky(sequence,elist,stickies)
    global gcutters
    gcutters=[]
    for enz in elist:
        if sequence.count(enz['pattern'])==1:
            gcutters.append(enz)
    firstcutsite=len(sequence)
    firstcutter=''
    lastcutsite=0
    lastcutter=''
    for gc in gcutters:
        if gc['overlap'] in stickies:
            # yes, we found a compatible cutter
            ecut=sequence.find(gc['pattern'])
            if ecut < firstcutsite:
                firstcutsite=ecut
                firstcutter=gc
            #add an if statement and block to ckeck if it is the last cutter.
            if ecut > lastcutsite:
                lastcutsite=ecut
                lastcutter=gc
    while abs(stickies[firstcutter['overlap']][0]['pcut']-stickies[lastcutter['overlap']][0]['pcut']) <10:
         newcutter=random.sample(cutters, 1)[0]
         if abs(newcutter['pcut']-stickies[firstcutter['overlap']][0]['pcut']) > 10:
             newenz=random.sample(olist[newcutter['overlap']],1)[0]
             sequence=sequence.replace(newenz['pattern'],'').replace(firstcutter['pattern'], newenz['pattern'])
             firstcutter=newenz    
    return sequence    
        
        
        
        
def addsticky(sequence,elist,stickies):
    random.seed(len(sequence))
    newseq=''
    newsite=False
    while not newsite:
        for n in range(200):
            newseq=newseq+random.sample(['A','C','G','T'],1)[0]
        for e in elist:
            if e['overlap'] in stickies and sequence.count(e['pattern'])==0:
                newsite=True
                if newseq.count(e['pattern'])>1:
                    newseq=newseq.replace(e['pattern'],'',newseq.count(e['pattern'])-1)
                elif newseq.count(e['pattern']) == 0:
                    splice=random.sample(range(200),1)[0]
                    newseq=newseq[:splice]+e['pattern']+newseq[splice:]
    return newseq
            
        
        
        
        
        