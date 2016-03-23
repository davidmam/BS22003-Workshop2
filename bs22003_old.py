


import random 

class bi():

	geneticCode = {
"AAA":"K",
"AAC":"N",
"AAG":"K",
"AAT":"N",
"ACA":"T",
"ACC":"T",
"ACG":"T",
"ACT":"T",
"AGA":"R",
"AGC":"S",
"AGG":"R",
"AGT":"S",
"ATA":"I",
"ATC":"I",
"ATG":"M",
"ATT":"I",
"CAA":"Q",
"CAC":"H",
"CAG":"Q",
"CAT":"H",
"CCA":"P",
"CCC":"P",
"CCG":"P",
"CCT":"P",
"CGA":"R",
"CGC":"R",
"CGG":"R",
"CGT":"R",
"CTA":"L",
"CTC":"L",
"CTG":"L",
"CTT":"L",
"GAA":"E",
"GAC":"D",
"GAG":"E",
"GAT":"D",
"GCA":"A",
"GCC":"A",
"GCG":"A",
"GCT":"A",
"GGA":"G",
"GGC":"G",
"GGG":"G",
"GGT":"G",
"GTA":"V",
"GTC":"V",
"GTG":"V",
"GTT":"V",
"TAA":"*",
"TAC":"Y",
"TAG":"*",
"TAT":"Y",
"TCA":"S",
"TCC":"S",
"TCG":"S",
"TCT":"S",
"TGA":"*",
"TGC":"C",
"TGG":"W",
"TGT":"C",
"TTA":"L",
"TTC":"F",
"TTG":"L",
"TTT":"F" }

	codons={}
	
	cfmutationlist=((31,"L",'CTC'),(117,"P","CCC"),(178,"R","AGA"),(225,"R","CGT"),(455,"E","GAG"),
                (549,"N","AAT"), (1066,"H","CAT")
              )
	nonsynlist = ((44,"V","GTT"),
(75,"Q","CAA"),
(138,"P","CCA"),
(170,"H","CAT"),
(182,"G","GGT"),
(322,"M","ATG"),
(351,"S","AGT"),
(351,"S","TCT"),
(353,"H","CAT"),
(353,"H","CAC"),
(467,"F","TTT"),
(470,"M","ATG"),
(506,"M","ATG"),
(506,"V","GTC"),
(507,"V","GTC"),
(508,"C","TGT"),
(532,"E","GAG"),
(562,"I","ATA"),
(576,"A","GCA"),
(605,"F","TTT"),
(654,"G","GGT"),
(668,"C","TGT"),
(693,"L","CTT"),
(693,"L","TTA"),
(693,"L","TTG"),
(903,"H","CAT"),
(909,"I","ATC"),
(967,"S","TCG"),
(1067,"V","GTC"),
(1162,"L","CTA"),
(1220,"I","ATA"),
(1453,"W","TGG") )

	asingle ={"A":"Ala", "C":"Cys", "D":"Asp", "E":"Glu", 
          "F":"Phe", "G":"Gly", "H":"His", "I":"Ile", 
          "K":"Lys", "L":"Leu", "M":"Met", "N":"Asn",
          "P":"Pro", "Q":"Gln", "R":"Arg", "S":"Ser", 
          "T":"Thr", "V":"Val", "W":"Trp", "Y":"Tyr","*":"STOP"}
	atriple ={}

	mrna=''
	
	def __init__(self):
	
		for c in self.geneticCode.keys():
			try:
				self.codons[self.geneticCode[c]].append(c)
			except:
				self.codons[self.geneticCode[c]]=[c]

		for x in self.asingle.keys():
			self.atriple[self.asingle[x]]=x
		try:
			fh=open(


	def createPatients(number):
		#sample of 4 diseases
		muts=random.sample(cfmutationlist,3)
		muts.append((551,"D","GAT"))
		nsmuts=random.sample(nonsynlist, 15)
		patients=[]
		normals=[]
		for p in range(number):
			patmrna=mrna
			dm=random.sample(muts, 1)
			nsm=random.sample(nsmuts, 6)
			syn=random.sample(range(len(mrna)/3),10)
			for m in syn:
				cs=3*(m-1)
				ce=3*m
				aa=geneticCode[mrna[cs:ce]]
				patmrna=patmrna[:cs]+random.sample(codons[aa],1)[0]+patmrna[ce:]
			for m in nsm:
				cs=3*(m[0]-1)
				ce=3*m[0]
				aa=geneticCode[mrna[cs:ce]]
				patmrna=patmrna[:cs]+m[2]+patmrna[ce:]
			for m in dm:
				cs=3*(m[0]-1)
				ce=3*m[0]
				aa=geneticCode[mrna[cs:ce]]
				patmrna=patmrna[:cs]+m[2]+patmrna[ce:]
			patients.append(("patient", p,patmrna))
		for p in range(number):
			patmrna=mrna
			dm=random.sample(muts, 1)
			nsm=random.sample(nsmuts, 6)
			syn=random.sample(range(len(mrna)/3),10)
			for m in syn:
				cs=3*(m-1)
				ce=3*m
				aa=geneticCode[mrna[cs:ce]]
				patmrna=patmrna[:cs]+random.sample(codons[aa],1)[0]+patmrna[ce:]
			for m in nsm:
				cs=3*(m[0]-1)
				ce=3*m[0]
				aa=geneticCode[mrna[cs:ce]]
				patmrna=patmrna[:cs]+m[2]+patmrna[ce:]        
			patients.append(('normal', p+number,patmrna))
		return patients
