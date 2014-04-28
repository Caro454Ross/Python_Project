import ClustalAlignment

#Display menu:

filename = None
file_loaded=False #becomes true once file has been successfully loaded
num_of_file=0# counts different number of alignment files

option = raw_input('*'*59+'\n* MULTIPLE ALIGNMENT ANALYZER'+' '*29+'*\n'+'*'*59+'\n* Select an option from below:'+' '*28+'*\n*'+' '*57+'*\n\t1)Open a Multiple Alignment File\t\t(O)\t*\n*\t2)Alignment Information\t\t\t(A)\t*\n*\t3)Alignment Explorer\t\t\t(E)\t*\n*\t4)Information per sequence\t\t(I)\t*\n*\t5)Analysis of Glycosylation Signatures\t(S)\t*\n*\t6)Export to fasta\t\t\t(X)\t*\n*\t7)Percentage Identity Calculator\t\t(P)\t*\n\t8)Exit\t\t\t\t\t(Q)\t*\n*'+' '*(59-len('File: '+str(filename)+' *'))+'File: '+str(filename)+' *'+'\n*'+'*'*59+'\n* Option: ').upper()

      

while option!='Q' and option!='8':
        if option=='1' or option=='O':
                filename=raw_input('Enter a valid path for Mutliple Alignment File: ')
                try:
                        alignment=ClustalAlignment.ClustalAlignment(filename)
                        try:
                                alignment.load_alignment_data()
                                file_loaded=True
                                print "File successfully loaded"
                                menu=raw_input("Press[Enter] to display menu: ")
                        except IndexError:
                                print "ERROR!!!!!: FILE FORMAT IS CORRUPTTED"
                                menu=raw_input("Press[Enter] to display menu: ")
                except IOError:
                        print "ERROR!!!!!: FILE NOT FOUND"
                        menu=raw_input("Press[Enter] to display menu: ")
        else:
                if file_loaded:
                        if option=='2' or option=='A':
                                print "\nFile name\t"+filename
                                seqIDs=alignment.sequence_dictionary.keys()
                                s=''
                                for ID in seqIDs:
                                        s+=ID+','
                                print 'Sequences\t'+s
                                print 'Length\t\t'+str(alignment.align_length)
                                print 'Average length\t'+str(alignment.average_length)
                                print '% of matches\t'+alignment.match_percentage()
                                print 'Number of X-matches:'
                                
                                match_dic=alignment.x_matches()
                                keys = match_dic.keys()
                                for match in keys:
                                        print '\t\t['+str(match)+']\t'+str(match_dic[match])
                                
                                menu=raw_input("Press[Enter] to display menu: ")
                        else:
                                if option=='3' or option =='E':
                                        option='E'
                                        while option=='E':
                                                try:
                                                        start= int(raw_input("Enter start of segment: "))
                                                        end=int(raw_input("Enter end of segment: "))
                                                        print '\n'+alignment.segment(start,end)
                                                        option=raw_input("Press[Enter] to display menu or E to for another segment: ").upper()
                                                except IndexError:
                                                        print "Segment not found in sequence:\nSegment may be too large"
                                                        option=raw_input("Press[Enter] to display menu or E to for another segment: ").upper()
                                else:
                                        if option=='4' or option=='I':
                                                option = 'I'
                                                while option=='I':
                                                        try:
                                                                ID=raw_input("Please enter sequence ID: ")
                                                                seq_info=alignment.seq_info(ID)
                                                                print "ID:\t"+ID+'\nLength:\t'+str(seq_info[1])+'\nBase Frequency:'
                                                                base_freq=seq_info[2]
                                                                print '\t\t[A]: '+str(base_freq[0])
                                                                print '\t\t[C]: '+str(base_freq[1])
                                                                print '\t\t[G]: '+str(base_freq[2])
                                                                print '\t\t[T]: '+str(base_freq[3])
                                                                print 'Sequence:\n'+seq_info[3]
                                                                option=raw_input("Press[Enter] to display menu or I for another sequence: ").upper()
                                                        except KeyError:
                                                                print "Sequence Not found: ID incorrect"
                                                                option=raw_input("Press[Enter] to display menu or I for another sequence: ").upper()
                                                        
                                                        
                                                
                                        else:
                                                if option=='5' or option=='S':
                                                        try:
                                                                not_found=''#represents IDs for seq without any signatures
                                                                signature_info=alignment.glycosylate()
                                                                for seq in signature_info:
                                                                        if seq[2]:
                                                                                print 'Glycosylation signatures found in '+seq[0]+':\n'+seq[1]
                                                                        else:
                                                                                not_found+=seq[0]+','
                                                                print 'There are no signatures were in: '+not_found
                                                                menu=raw_input("Press[Enter] to display menu: ")
                                                        except ValueError:
                                                                print 'ERROR! Sequences in file are out of frame or are not DNA sequence'
                                                                menu=raw_input("Press[Enter] to display menu: ")
                                                
                                                else:
                                                        if option=='6' or option=='X':
                                                                fname=''
                                                                case=False
                                                                
                                                                file_opt=raw_input("Press M to create Multiple Files(one per sequence)\nOR press S to create a Single file with all sequences: ").upper()
                                                                
                                                                
                                                                seq_opt=raw_input("Do you want to save: DNA (Input: D) or Proteins (Input: P): ")
                                                                
                                                                case_opt=raw_input("Do you want to indicate glycosylation signatures by case [Y/N]: ").upper()
                                                                
                                                                if case_opt=='Y':
                                                                        case=True
                                                                if file_opt=='S':
                                                                        fname=raw_input("Enter the name of the file to save sequences: ")
                                                                alignment.export_to_fasta(fname,file_opt,seq_opt,case)
                                                                menu=raw_input("Press[Enter] to display menu: ")
                                                        else:
                                                                if option=='7' or option=='P':
                                                                        print "PERCENTAGE IDENTITY CALCULATOR:\nPercentage identities per sequence-sequence match:\n"
                                                                        
                                                                        perc_id_info=alignment.percentage_identity()
                                                                        for match in perc_id_info:
                                                                                print match[0]+"\n"+match[1]+"\nPercentage Identity = "+match[2]+'\n'
                                                                        menu=raw_input("Press[Enter] to display menu: ")
                                                                else:
                                                                        print "OPTION DOES NOT EXIST"
                                                                        menu=raw_input("Press[Enter] to display menu: ")

                else:
                        print "NO FILE HAS BEEN SUCCESSFULLY LOADED!"
                        menu=raw_input("Press [Enter] to display menu: ")
                        
        option = raw_input('*'*59+'\n* MULTIPLE ALIGNMENT ANALYZER'+' '*29+'*\n'+'*'*59+'\n* Select an option from below:'+' '*28+'*\n*'+' '*57+'*\n\t1)Open a Multiple Alignment File\t\t(O)\t*\n*\t2)Alignment Information\t\t\t(A)\t*\n*\t3)Alignment Explorer\t\t\t(E)\t*\n*\t4)Information per sequence\t\t(I)\t*\n*\t5)Analysis of Glycosylation Signatures\t(S)\t*\n*\t6)Export to fasta\t\t\t(X)\t*\n*\t7)Percentage Identity Calculator\t\t(P)\t*\n\t8)Exit\t\t\t\t\t(Q)\t*\n*'+' '*(59-len('File: '+str(filename)+' *'))+'File: '+str(filename)+' *'+'\n*'+'*'*59+'\n* Option: ').upper()
 
          
             
            
            

