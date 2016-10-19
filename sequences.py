from string import *
from copy import copy
from copy import deepcopy
import re
import exceptions
import commands
dna_strict ='atgcATGC'
dna_lax ='atgcnATGCN'
dna_degen ='atgcrywsmkhdbvnATGCRYWSMKHDBVN'
alpha='abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
complement_bases = {'a':'t','c':'g', 'g':'c', 't':'a', 'u':'a', 'r':'y', 'y':'r','w':'w','s':'s','m':'k','k':'m','h':'d','d':'h','b':'v','v':'b','n':'n'}


genetic_code = {'ttt': 'F','tct': 'S','tat': 'Y','tgt': 'C','ttc': 'F','tcc': 'S','tac': 'Y','tgc': 'C','tta': 'L','tca': 'S','taa': '*', 'tga': '*','ttg': 'L','tcg': 'S','tag': '*','tgg': 'W','ctt': 'L','cct': 'P','cat': 'H','cgt': 'R','ctc': 'L','ccc': 'P','cac': 'H','cgc': 'R','cta': 'L','cca': 'P','caa': 'Q','cga': 'R','ctg': 'L','ccg': 'P','cag': 'Q','cgg': 'R','att': 'I','act': 'T','aat': 'N','agt': 'S','atc': 'I','acc': 'T','aac': 'N','agc': 'S','ata': 'I','aca': 'T','aaa': 'K','aga': 'R','atg': 'M','acg': 'T','aag': 'K','agg': 'R','gtt': 'V','gct': 'A','gat': 'D','ggt': 'G','gtc': 'V','gcc': 'A','gac': 'D','ggc': 'G','gta': 'V','gca': 'A','gaa': 'E','gga': 'G','gtg': 'V','gcg': 'A','gag': 'E','ggg': 'G', 'nnn':'X','---':'-'}

degenerated_bases={'r':['g','a'], 'y':['t','c'],'m':['a','c'],'k':['g','t'],'s':['g','c'],'w':['a','t'],'h':['a','t','c'], 'b':['g','t','c'],'v':['g','c','a'],'d':['g','a','t']}

KD_scale = {'A':1.8,'C':0.5, 'D':-3.5, 'E':-3.5, 'F':2.8, 'G':-0.4, 'H':-3.2, 'I':4.5, 'K':-3.9, 'L':3.8, 'M':1.9, 'N':-3.5, 'P':-1.6, 'Q':-3.5, 'R':-4.5, 'S':-0.8, 'T':-0.7, 'V':4.2, 'W':-0.9, 'Y':-1.3}

GVH_scale= {'A':0.267, 'C':1.806, 'D':-2.303, 'E':-2.442, 'F':0.427, 'G':0.16, 'H':-2.189, 'I':0.971, 'K':-2.996, 'L':0.623, 'M':0.136, 'N':-1.988, 'P':-0.451, 'Q':-1.814, 'R':-2.749, 'S':-0.119, 'T':-0.083, 'V':0.721, 'W':-0.875, 'Y':-0.386}
GES_scale = {'A':1.6, 'C':2.0, 'D':-9.2, 'E':-8.2, 'F':3.7, 'G':1.0, 'H':-3.0, 'I':3.1, 'K':-8.8, 'L':2.8, 'M':3.4, 'N':-4.8, 'P':-0.2, 'Q':-4.1, 'R':-12.3, 'S':0.6, 'T':1.2, 'V':2.6, 'W':1.9, 'Y':-0.7}

def guess_type (seq):
    type = 'dna'
    for element in seq:
        if element in alpha:
            if element in 'efilpqEFILPQ':
                type = 'protein'
                break
    return type

def make_sequence (name, sequence, annotation=[]):
    type = guess_type(sequence)
    if type == 'protein':
        p = Protein (name, sequence, annotation)
        return p
    elif type == 'dna':
        dna = Dna (name, sequence,annotation)
        return dna
    else:
        raise TypeError

def motif_parser (infile):
        motif = []
        line = infile.readline()
        while line:
            list = split(line)
            m= Motif(list[0], list[1])
            motif.append(m)
            line = infile.readline()

        return motif

def get_seq (id, db='sp'):
    cmd="golden " + db + ":" + id
    status, output = commands.getstatusoutput(cmd)
    if status != 0:
        return status, "A problem occured: " + output + '\nstatus: ' + str(status) + '\n'
    return None, output


def no_stop_protein(protein):
    if '*' not in protein:
	return protein
    no_stop = replace(protein, '*', 'X')
    return no_stop

def reverse_all(prot_seq, old_dna, check = True):
    new_dna = {}
    for protein in prot_seq.keys():
	#print protein
	for dna in old_dna.keys():
	    #print dna
	    if protein == dna:
		#print 'hola'
		if check == True:
		    code = old_dna[dna].coding()
		    print 'check'
		else:
		    code = old_dna[dna].sequence
		ndna= prot_seq[protein].reverse_translation(code)
		new_dna[dna] = ndna
		#print protein
		break
    return new_dna


def max_length(prot_str):
    #print prot_str
    max = 0
    start = prot_str.find('M')
    while start != -1:
	stop = prot_str.find('*', start)
	if stop == -1:
	    stop = len(prot_str) -1
	orf = stop - start
	if orf > max:
	    max = orf
	start = prot_str.find('M', start+1)
    #print 'max', max
    return max


def degen_synon(codon):
    queu =[]
    p = []
    for base in codon:
	p+=[base]
    queu += [p]
    for elem in queu:
	for char in elem:
	    if char not in 'atgc':
		i = index(elem, char)
		opt = degenerated_bases[char]
		temp = copy(elem)
		for alt in opt:
		    temp[i] = alt

		    queu += deepcopy([temp])

    new_list=[]
    for entry in queu:
	resolved = True
	for base in entry:
	    if base not in 'atgc':
		resolved = False
		break
	if resolved:
	    new_codon = join(entry)
	    new_codon = replace(new_codon, ' ','')
	    new_list+= [new_codon]
    list_aa=[]
    for c in new_list:
	aa = genetic_code[c]
	list_aa+=[aa]
    for i in range(1, len(list_aa)):
	if list_aa[0] != list_aa[i]:
	    return 'nnn'
    return c

class Sequence:
    alphabet = 'atgc'

    """Creates objects, either proteins or dna that contain the name and the sequence"""
    def __init__ (self, name=None, sequence=None, annotation=''):
        self.name= name
        self.sequence= sequence
        self.annotation=annotation
        #self.annotation.append(annotation)
        self.type = None
	self.compart = None
	self.subtype = None
	self.Hsubtype = None

    def __str__(self):
        text = self.name +'\n'+ str(self.annotation) +'\n' + self.sequence
        return text

    def to_fasta (self, include_annot = False):
        fasta = '>'+ self.name

	if include_annot == True:
	    for elem in self.annotation:
		fasta = fasta +' '+ elem
	fasta+='\n'+self.sequence+'\n'
        return fasta

    def annotate (self, text):
        self.annotation.append(text)

    def print_nice (self, n):
        """Prints a sequence in lines of n characters"""
        seq= self.sequence.clean()
        for base in range(0, len(seq),n):
            print seq[base:base+n]


    def __getitem__(self,element):
        item = self.sequence[element]
        return item


    def __delitem__ (self,position):
        seq = self.sequence[0:position] + self.sequence[position+1:]
        return seq

    def remove (self,position):
        seq = self.sequence[0:position] + self.sequence[position+1:]
        return seq

    def __setitem__ (self,position,value):
        seq = self.sequence[0:postion] + value +self.sequence[position+1:]
        return seq

    def __add__ (self, other):
        if self.__class__ == other.__class__:
            return self.sequence + other.sequence
        else:
            raise SeqError, 'Cannot concatenate sequences from different classes'


    def __eq__ (self,other):
        """ compares if the two sequences are identical"""
        if self.sequence == other.sequence:
            return True
        else:
            return False

    def __len__ (self):
        """returns the length of the sequence"""
        lenght = len(self.sequence)
	#n = count(self.sequence,'-')
        return lenght

    def __contains__ (self,element):
        """checks if the sequence contains the element"""

        if element in self.sequence:
            return True
        else:
            return False

    def id(self):
        seq_id = self.name
        return seq_id

    def clean(self):
        """Clean a sequence, keeping only alphabetic characters, and then cheking if ther are undefined bases in the sequence. It takes as parameters a sequence, the character to keep (by default the alphabet)"""
	limpia = ''
	for base in self.sequence:
		if base in alpha:
			if base in dna_strict:
				limpia = limpia + base
			else:
				limpia = limpia + 'n'
	return lower(limpia)

    def set_compart(self, list_compart):
	for compart in list_compart:
	    if compart in self.name:
		self.compart = compart
		break

    def set_type (self, types = {}):
	for type in types.keys():
	    if self.sequence == types[type]:
		self.type = type
		break
	if self.type == None:
	    self.type = self.name
	    types[self.name]= self.sequence
	return types


class Dna(Sequence):


    alphabet ='atgcATGC'

    def gc_percent(self):
        """ returns the percentage of gc in a sequence"""
        seq= self.sequence.clean_seq()
        gc = 0
        for char in seq:
            if (char == 'c') or (char == 'g'):
                gc = gc +1
        gc_per= float(gc)/len(seq) *100
        return gc_per


    def complement (self):
        """Returns the complementary sequence of a DNA strand"""
        comp= ''
        for base in self.sequence:
            comp= comp + complement_bases[base]

        return comp


    def reverse_complement (self):
        """ Returns the reverse of a DNA sequence """
        seq = self.sequence.complement()
        rev= ''
        pos2= -(len(seq)+1)
        for charac in range (-1,pos2,-1):
            rev = rev + seq[charac]
        return rev

    def seq_checker (self, code=dna_lax):
        """Cheks if the sequence contains degenetated bases"""

        for base in self:
            if base in dna_x:
                seq_type = True
            else:
                seq_type = False
        return seq_type


    def base_count(self):
        """Returns a list of tuples with the numer of a, t c, g, n or other characteres in a sequence"""
        seq= self.clean(dna_lax)
        count_a = count_t = count_g = count_c = count_n = count_o = 0
        for base in seq:
            if base == 'a':
                count_a = count_a + 1
            elif base == 't':
                count_t = count_t + 1
            elif base == 'g':
                count_g = count_g + 1
            elif base == 'c':
                count_c = count_c + 1
            elif base == 'n':
                count_n = count_n + 1
            else:
                count_o = count_o + 1
	return [('a',count_a),('t', count_t),('g', count_g),('c', count_c),('n', count_n),('o', count_o)]



    def translation (self, phase=None):
	#print self.sequence
	if phase == None:
	    phase = self.find_orf()
        #seq = self.clean()
	seq = lower(self.sequence)
        prot = ''
        i = phase -1
        for i in range (i,len(seq),3):
            codon = seq[i:i+3]
	    #print codon
            if len(codon) == 3:
		new_codon = ''
                if "-" in codon:
                    codon = 'nnn'
		if 'n' in codon:
		    codon ='nnn'
                for i in range(3):
		    if codon[i] in 'yrmkswhbvd':
			new_codon = degen_synon(codon)
			break
		if new_codon:
		    codon = new_codon
                prot = prot + genetic_code[codon]
        return prot

    def degenerate_stops(self, phase = 1):
		seq = self.clean()
		no_stop_seq = ''
		i = phase -1
		for i in range (i, len(seq), 3):
			codon = seq[i:i+3]
			if len(codon) < 3:
				break
			if codon == 'taa' or codon == 'tag' or codon == 'tga':
				codon = 'nnn'
			no_stop_seq = no_stop_seq + codon
		return no_stop_seq


    def trans_6_phases (self):

        """Finds the 3 phases of tanslation in a sequence"""
        trans={}
        for phase in range (1,4):
            trans[phase]= self.translation()


        for phase in range(4,7):
            trans[phase]= (self.reverse.complement).translation()

    def find_orf(self):
	#print self.name
	phase =1
	max = 0
	for i in range(1,4):
	    orf= self.translation(i)
	    n = max_length(orf)
	    #print n
	    if n > max:
		max = n
		phase = i
	#print max
	return phase

    def coding(self):
	phase = self.find_orf()
	prot = self.translation(phase)
	shift = phase -1
	start = 0
	R_start = start
	max= 0
	while start != -1:
	    stop = prot.find('*', start +1)
	    if stop == -1:
		stop = len(prot) -1
	    n = stop - start
	    if n > max:
		max = n
		R_start = start
		R_stop = stop
	    start = prot.find('*', start+1)
	    if start == -1:
		 break
	    else:
		start+=1
	nstart = prot.find('M', R_start)
	nstart = nstart*3 + shift
	nstop = R_stop*3+1 + shift
	#print nstart, nstop
	coding = self.sequence[nstart:nstop]
	return coding

    def ispresent (self, enz):
        """ finds if a motif is present in a sequence"""
        enz = lower (enz)
        seq = self.sequence.clean_seq()
        if enz in seq:
            return True
        else:
            return False


    def isunique (self,enz):
        """ Finds if a restriction site is unique in a sequence"""
        enz= lower(enz)
        seq= self.sequence.clan_seq()
        if count (seq,enz) == 1:
            return True
        else:
            return False


    def distance (self,enz1,enz2):
        """ Gives the distance between two restrictiion sites"""
        enz1= lower (enz1)
        enz2= lower (enz2)
        seq= self.sequence.clean_seq()
        if enz1 in seq and enz2 in seq:
            enz1_site= find (seq, enz1)
            enz2_site= find (seq, enz2)
            if enz1_site > enz2_site :
                distance= enz1_site - enz2_site
                return distance
            else:
                distance= enz2_site -enz1_site
                return distance
        else:
            return 0





class Protein (Sequence) :

    alphabet = 'ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy'


    def reverse_translation(self, old_dna):
	new_dna =''
	i = 0
	for aa in self.sequence:
	    if aa == '-':
		new_dna = new_dna + '---'
	    else:
		new_dna = new_dna + old_dna[i:i+3]
		i = i+3
	#redo_old = replace(new_dna, '-','')
	#print len(redo_old) - len(old_dna)
	#print len(old_dna), i
	#print len(self.sequence)*3, count(self.sequence, '-')*3
	dna = make_sequence(self.name, new_dna,[])
	return dna


    def hydrophobe (self, scale, length=21):
        """Calculates the hydrophobicity in a window of length n, in a sequence using an specific scale"""
        seq= self.sequence.clean()
        j = 1
        i = 0
        result = {}
        while i <len(seq):
            interval= seq[i:i+length]
            hphobe = 0 + scale[seq[i]]
            result[j] = hphobe
            i = i + length
            j = j +1
        return result


class FastaIterator:
    def __init__ (self,fh, accession = False):
        self.fh = fh
        self.line= None
	self.accession = accession

    def next (self):
        seq = ''
	name = ''
        a=None
        if self.line == None:
            self.line = self.fh.readline()
        if self.line != '' :
	    #print self.line
            if self.line[0:1] == '>':
                heading = self.line.split()
                name = heading[0][1:]
		if '\n' in name:
		    name = replace(name, '\n', '')
		if self.accession == True:
		    name = name.split('|')[3]
		if len(heading)>2:
		    annotation = heading[1:]
		else:
		    annotation = []
		#print name, annotation
                self.line = self.fh.readline()

            while self.line and self.line[0:1] != '>' :
                seq = seq + self.line

                self.line = self.fh.readline()
            seq = replace(seq,'\n','')
            seq = replace(seq,'\r','')
	    seq = replace(seq, ' ','')
            a = make_sequence(name,seq, annotation)
            return a
        else:
            return None

class readInterleaved:
    def __init__(self, fh):
	self.fh = fh
	self.seq_dico = {}
    def read(self):
	names = []
	self.line = self.fh.readline()
	self.line = self.line.replace('\n','')
	a = self.line.split()
	nseq = int(a[0])
	# print nseq
	length = int(a[1])
	self.line = self.fh.readline()
	i = 0
	while i< nseq:
	    self.line = self.line.replace('\n','')
	    a = self.line.split()
	    names.append(a[0])
	    seq = ''
	    for elem in a[1:]:
		seq+= elem
	    self.seq_dico[a[0]] = make_sequence(a[0],seq,[])
	    self.line = self.fh.readline()
	    i+= 1
	while len(self.seq_dico[names[-1]].sequence)<length:
	    self.line= self.line.replace('\n','')
	    for elem in names:
		seq = self.line.strip()
		seq = seq.replace(' ','')
		self.seq_dico[elem].sequence = self.seq_dico[elem].sequence + seq
		self.line = self.fh.readline()
		self.line = self.line.replace('\n', '')
	return self.seq_dico





class ProtectedSeq:
    def __init__ (self,sequence):
        self.sequence=sequence

    def __getattr__(self,name):
        return getattr(self.sequence,name)

    def __setitem__ (self,position,value):
        if self.sequence[position] not in 'atgcATGC-':
            self.sequence= self.sequence[:position] + value + self.sequence[position+1:]

        else:
            raise SeqError, 'You are not allowed to modify defined bases'

    def __delitem__ (self,position):
        if self.sequence[position] == '-':
            self.sequence= self.sequence[:position] + self.sequence[position+1:]

        else:
            raise SeqError, 'You are only allowed to remove gaps'






class Motif :
    def __init__ (self, name=None, expression=None):
        self.name = name
        self.expression= expression
        self.pattern = re.compile(expression)

    def __str__ (self):
        text = self.name+' '+ self.expression
        return text

    def find (self, protein,pos=0):
        pos = []
        match = self.pattern.search(protein.sequence)
        if match == None:
            return Match(self,protein,None)
        while match !=None:
            pos.append(match.span())

            match = self.pattern.search(protein.sequence, match.start()+1)

        return Match(self,protein,pos)



class Match :
    def __init__ (self, motif=None, protein= None, pos=None):
        self.motif = motif
        self.protein = protein
        self.position= pos
        if self.position != None:
            self.found = True
        else:
            self.found = False

    def show (self):
        motif_seq = []
        for element in self.position:
            m= self.protein.sequence[element[0]:element[1]]
            motif_seq.append(m)
        return motif_seq

    def __str__ (self):
        matches = ''
        for element in self.position:
            matches = matches + str(element)

        text = self.motif.expression +' '+ 'protein'+ ' '+ self.protein.name+' '+'position'+ matches
        return text



class SeqError(exceptions.Exception):
    def __init__ (self, args=None):
        self.args = args
