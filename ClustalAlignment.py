import re
class ClustalAlignment:
      fliename=None
      f=None
      sequence_dictionary=None
      Conservation=None
      align_length=None
      average_length=None
      total_count=None
      format_id=None #length from start of seqID until start of sequence, used to calculate num of white spaces for formatting reasons
      

      def __init__(self,filename):
            self.filename=filename
            self.f= open(self.filename,'r')
            self.sequence_dictionary={}
            self.Conservation=''
            self.align_length=0
            self.average_length=0
            self.total_count=0
            self.format_id=0             
                       
            
      def load_alignment_data(self): #creates a dictionary of id and sequences, calls the set_alignment_length method to init the length, average length and count values of the class
            self.f.readline()#eliminates first line of file
            self.init_sequence_dictionary()
            line = self.f.readline()#moves to next line after block of sequences
            while line!='':#ends at the end of the file
                  block_info=self.get_block_data(line)
                  seq_id=block_info[0]
                  sequences=block_info[1]
                  for i in range(len(sequences)):
                        original_sequence=self.sequence_dictionary[seq_id[i]]
                        update_seq=original_sequence+sequences[i]
                        self.sequence_dictionary[seq_id[i]]=update_seq
                  line = self.f.readline()#moves to next line after block of sequences
            self.set_alignment_length()
            
            
      
      def get_block_data(self,line):
            seq_id=[]
            sequences=[]
            while line=='\n': #removes initial blanklines before start of alignment
                  line=self.f.readline()
            while line[0]!=' ':
                  firstspace=line.index(' ')
                  
                            
                  #determinmation of the start of the sequence
                  i=firstspace
                  while line[i]==' ':
                        i=i+1
                                  
                  # removes \n from end of line and spilts by white space if the cumulative count of residues is present
                  seq_id.append(line[0:firstspace])
                  seq_data= ((line[i:]).rstrip()).split(' ')
                  sequences.append(seq_data[0])
                  end=len(seq_data[0])#calculates length of sequence, critical to extract correct conservation pattern including whitespaces at end of sequence
                  line=self.f.readline()#moves to next line
            self.format_id=i           
            con_pattern=(line[i:]).rstrip()
            missing_white_spaces = ' '*(end-len(con_pattern))
            self.Conservation = self.Conservation+ con_pattern+ missing_white_spaces
            return (seq_id,sequences)
            
            
            
      def init_sequence_dictionary(self):      
            line = self.f.readline()#Moves to first blank line in file
            dic_info=self.get_block_data(line)
            seq_id=dic_info[0]
            sequences=dic_info[1]
            for i in range(len(sequences)):
                  self.sequence_dictionary.update({seq_id[i]:sequences[i]})
                  
    
      def set_alignment_length(self): #sets the length of alignment and average length of sequences, also sets total_count of sequences
            
            #Since the sequences are all aligned with gaps(-), all the string sequences in the sequence_dictionary are the same length, thus the alignment length can be given by any sequence.
            
            sequences=self.sequence_dictionary.values()
            self.align_length=len(sequences[0])
            
            #to determine the average length, all gaps must be excluded. The conservation sequence must also be excluded
            self.total_count = len(sequences)
            total_length=0
            for i in range(self.total_count):
                  for j in sequences[i]:
                        if j=='-':
                              continue
                        else:
                              total_length+=1
            self.average_length=total_length/self.total_count
            
      
      def seq_ids(self):
            seq_ids = self.sequence_dictionary.keys()
            return seq_ids
      
      def match_percentage(self): #returns match percentage
            count_matches=0.0
            for i in self.Conservation:
                  if i =='*':
                        count_matches+=1
            match_percentage= (count_matches/self.align_length)*100.0
            return "%.2f"%(match_percentage)+"%"
      
      def x_matches(self):
            
            #creates a dictionary with X matches as keys, from 2 matches up to total number of sequences. Sets all matches to zero.
            match_dictionary={}
            for i in range(2,self.total_count+1):
                  match_dictionary.update({i:0})
            
            sequences=self.sequence_dictionary.values()      
            already_checked=[] #list of residues that have already been checked
            
            for res in range(self.align_length):
                  residueList=[]
                  already_checked=[] #list of residues that have already been checked
                  for seq in sequences:
                        residueList.append(seq[res])#creates a list of all residues at the same index from each sequence
                  
                  for j in range(len(residueList)-1):
                        matches=1
                        if residueList[j] in already_checked or residueList[j]=='-':
                              continue
                        else:
                              for k in range(j+1,len(residueList)):
                                    if residueList[j]==residueList[k]:
                                          matches+=1
                        
                              already_checked.append(residueList[j])       
                              if matches>=2:
                                    match_dictionary[matches]+=1
                              
            return match_dictionary
      
      def segment(self,start,end):
            segments=''
            s=''
            seq_ids=self.sequence_dictionary.keys()
            len_id=len(seq_ids[0])
            num_blocks = ((end-start)+1)/60
            block=0
            
            
            for block in range(num_blocks):
                  s=''
                  con_pattern=' '*self.format_id
                  for ID in seq_ids:
                        seq=self.sequence_dictionary[ID]
                        s=s+ID+' '*(self.format_id-len(ID))
                        for base in range((start-1)+(60*block),(start-1)+60+(60*block)):
                              s=s+seq[base]
                        s=s+'\n'
                  con_pattern=con_pattern+self.Conservation[start-1+60*block:base]+'\n'
                  s+=con_pattern
                  segments+=s
                  
            if ((end-start)+1)%60!=0:
                  s=''
                  con_pattern=' '*self.format_id
                  for ID in seq_ids:
                        seq=self.sequence_dictionary[ID]
                        s=s+ID+' '*(self.format_id-len(ID))                       
                        for base in range((start-1)+(60*block),end):
                              s=s+seq[base]
                        s=s+'\n'
                  con_pattern=con_pattern+self.Conservation[start-1+60*block:end]+'\n'
                  s+=con_pattern                  
                  segments+=s
            
            return segments
      
      def seq_info(self,sID):
            seq=self.sequence_dictionary[sID]
            
            seq_ex_gaps=self.remove_gaps(seq)
            
            length=len(seq_ex_gaps)
            
            count=[0,0,0,0]#holds counts for bases: count[0]=A,count[1]=C,count[2]=G,count3=[T]
            
            for base in seq_ex_gaps:
                  if base=='A' or base=='a':
                        count[0]+=1
                  else:
                        if base=='C' or base=='c':
                              count[1]+=1
                        else:
                              if base=='G' or base=='g':
                                    count[2]+=1
                              else:
                                    if base=='T' or base=='t':
                                          count[3]+=1
                                    
            seq=''
            for i in range(len(seq_ex_gaps)):
                  seq+=seq_ex_gaps[i]
                  if (i+1)%60==0:
                        seq+='\n'
            
            return (sID,length,count,seq)
      
      def remove_gaps(self,seq): #removes '-' from the sequences
            seq_ex_gaps=''
            for base in seq:
                  if base!='-':
                        seq_ex_gaps=seq_ex_gaps+base
            return seq_ex_gaps
      
      
      def glycosylate(self): #returns seqID, glycosyated seq, and a flag (false if no signatures are found)
            seqIDs=self.sequence_dictionary.keys()
            protein_sequences={}
            glycosylated_sequences=[]

            for ID in seqIDs:
                  
                  seq=self.remove_gaps(self.sequence_dictionary[ID])
                  pro = self.translate(seq)
                  protein_sequences.update({ID:pro})
            
           
            gly_exp = re.compile("[Nn][^Pp][SsTt]")
            
            for ID in seqIDs:
                  seq=protein_sequences[ID]
                  has_glycosylation=False
                  index=0
                  found=gly_exp.search(seq,index)
                  while found!=None:
                        index=found.start()
                        pattern= seq[index:index+3].upper()
                        seq=seq.replace(seq[index:index+3],pattern,1)
                        has_glycosylation=True
                        found=gly_exp.search(seq,index+1)
                  glycosylated_sequences.append((ID,seq,has_glycosylation))
            
            return glycosylated_sequences
      
      
      def translate(self,seq):
            
            amino_acids = {'ATG':'m','TGA':'stop','TAA':'stop','TAG':'stop','TTT':'f','TTC':'f','CTT':'l','CTC':'l','CTG':'l','CTA':'l','ATT':'i','ATC':'i','ATA':'i','GTA':'v','GTC':'v','GTG':'v','GTT':'v','TCA':'v','TCT':'v','TCG':'s','TCC':'s','CCA':'p','CCC':'p','CCT':'p','CCG':'p','ACT':'t','ACG':'t','ACC':'t','ACA':'t','GCG':'a','GCC':'a','GCA':'a','GCT':'a','TAT':'y','TAC':'y','CAT':'h','CAC':'h','CAA':'q','CAG':'q','AAT':'n','AAC':'n','AAA':'k','AAG':'k','GAT':'d','GAC':'d','GAA':'e','GAG':'e','GGA':'g','GGG':'g','GGC':'g','GGT':'g','CGA':'r','CGC':'r','CGG':'r','CGT':'r','TGT':'c','TGC':'c','TGG':'w','AGT':'s','AGC':'s'}
            seq=seq.upper()
            start=seq.index('ATG')
            codon=seq[start:start+3]
            protein=""
            while amino_acids[codon]!='stop':
                  protein=protein+amino_acids[codon]
                  start+=3
                  codon=seq[start:start+3]
                  if codon=='':#in the case of a stop codon not occuring
                        break            
            return protein          
      
      
      def export_to_fasta(self,fname,file_opt,seq_opt,case):                
            file_opt=file_opt.upper()
            seq_opt=seq_opt.upper()
            seq_ids=[]
            sequences=[]
            
            if not case:
                  if seq_opt=='P':
                        proteins = self.glycosylate()
                        for pro in proteins:
                              seq_ids.append(pro[0])
                              sequences.append(pro[1].lower())
                        
                  else:
                        if seq_opt=='D':
                              seq_ids= self.sequence_dictionary.keys()
                              for ID in seq_ids:
                                    seq=self.remove_gaps(self.sequence_dictionary[ID])
                                    sequences.append(seq.lower())
                                                     
                                    
            else:
                  if seq_opt=='P':
                        proteins = self.glycosylate()
                        for pro in proteins:
                              seq_ids.append(pro[0])
                              sequences.append(pro[1])
                                                         
                  else:
                        if seq_opt=='D':
                              seq_ids=self.sequence_dictionary.keys()
                              for ID in seq_ids:
                                    seq=self.remove_gaps(self.sequence_dictionary[ID])
                                    seq=self.glycosylate_DNA(seq)
                                    sequences.append(seq)
                                            
                                    
            if file_opt=='S':
                  w=open(fname,'w')
                  fasta=''
                  for i in range(len(sequences)):
                        fasta=fasta+">"+seq_ids[i]+'\n'
                        seq=sequences[i]
                        for j in range(len(seq)):
                              fasta+=seq[j]
                              if (j+1)%60==0:
                                    fasta+='\n'
                        fasta+='\n\n'
                  w.write(fasta)
                  w.close()
            else:
                  if file_opt=='M': #names files according to ID
                        for i in range(len(sequences)):
                              f=open(seq_ids[i]+'.fasta','w')
                              f.write('>'+seq_ids[i]+'\n')
                              f.write(sequences[i])
                              f.close()
                              
                              
      
      def glycosylate_DNA(self,seq): #returns DNA sequence with bases of gylosylation in capitals
            seq=seq.lower()
            asparginine=['aat','aac']
            proline=['cca','ccc','cct','ccg']
            serine=['tca','tcc','tcg','tct','agt','agc']
            threonine=['aca','acc','acg','act']
            for i in range(0,len(seq),3):
                  codon1= seq[i:i+3]
                  codon2= seq[i+3:i+6]
                  codon3= seq[i+6:i+9]
                  if (codon1 in asparginine) and (not(codon2 in proline)) and ((codon3 in threonine) or (codon3 in serine)):
                        glycos=seq[i:i+9]
                        seq=seq.replace(glycos,glycos.upper(),1)
            return seq
      
      def percentage_identity(self):# returns a list of identity calulations for all sequences against each other
            seqIDs=self.sequence_dictionary.keys()
            perc_identities=[]
            for i in range(len(seqIDs)-1):
                  for j in range(i+1,len(seqIDs)):
                        countMatches=0.0;
                        seq1=self.sequence_dictionary[seqIDs[i]]
                        seq2=self.sequence_dictionary[seqIDs[j]]
                        for k in range(self.align_length):
                              if seq1[k]==seq2[k] and seq1[k]!='-' and seq2[k]!='-':
                                    countMatches+=1
                        perc_ID="%.2f"%((countMatches/self.align_length)*100)+"%"
                        perc_identities.append((seqIDs[i],seqIDs[j],perc_ID)) #returns IDs of each sequence and % Match
            return perc_identities
                                    
                        
                  
                       
            
                  
                  
                        
                                                   
                                    
                              
                        
            
                  
                  
                  
            
            
      
