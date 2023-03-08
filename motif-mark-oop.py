import cairo
import argparse
import re
import math


# Argparse code to pass in paths to files

def get_args():
    parser = argparse.ArgumentParser(description="A program to find and plot binding motifs within a gene.")
    parser.add_argument("-f", "--fasta", help = "Path to fasta file with gene sequence", required = True, type = str)
    parser.add_argument("-m","--motif", help = "Path to file containing list of motifs", required = True, type = str)
    return parser.parse_args()

args =  get_args()

fasta_path = args.fasta
motif_path = args.motif

class Gene:

    def __init__(self,gene_name, sequence):
        self.gene = gene_name #gene name from fasta header
        self.sequence = sequence #gene sequence from fasta
        self.exon_list = [] #list to hold exon objects for use in gene methods
        self.motif_list = [] #list to hold motif objects for use in gene methods

    def find_exon(self):
        exon_start = [] #list to store the start positions of each exon
        exon_end = [] #list to store the end position of each exon

        for index in range(len(self.sequence)-1):

            #If lower case into upper case (start of exon)
            if self.sequence[index].islower() and self.sequence[index+1].isupper():
                exon_start.append(index+1)

            #If upper case into lower case (end of exon)
            elif self.sequence[index].isupper() and self.sequence[index+1].islower():
                exon_end.append(index+1)

        #check if sequence starts with exon
        if self.sequence[0].isupper():
            exon_start.append(0)
        #check if sequence ends with exon 
        if self.sequence[-1].isupper():
            exon_end.append(len(self.sequence))

        
        for index in range(len(exon_start)):

            start = exon_start[index]
            end = exon_end[index]
            exon_seq = self.sequence[start:end] #uses start and end positions to grab sequence from whole gene

            self.exon_list.append(Exon(exon_seq, start,self)) #creates and stores exon object within gene list




    def find_motif(self):
        #loop through motif list that was populated when motif objects were created
        for motifs in self.motif_list:

            #grabs motif with regex options to look for ambiguous bases
            motif = motifs.corr_seq

            #finds sequences that match patterns from regex expressions
            for found_motif in re.findall(motif, self.sequence, re.IGNORECASE):

                #finds positions of matches within whole sequence and stores 
                for pos in re.finditer('(?={0})'.format(re.escape(found_motif)),self.sequence, re.IGNORECASE):
                    motifs.position.append(pos.start())

        

class Exon:

    def __init__(self, sequence, position,read_obj):
        self.sequence = sequence #Exon sequence
        self.position = position #Exon starting position
        read_obj.exon_list.append(self) #Populate gene exon list

    def plot_exon(self, x, y, context):

        context.rectangle(x + self.position, y - 12.5, len(self.sequence), 25)
        context.fill()

class Motif:

    def __init__(self, sequence,gene_obj):
        self.sequence = sequence #Motif sequence from file
        self.corr_seq = sequence #Motif sequence that will be changed to include regex expressions
        self.position = [] #List to hold all starting positions for motif in a given gene
        gene_obj.motif_list.append(self) # Stores motif object in a list in the gene 
        self.gene_seq = gene_obj.sequence # Stores gene seq to be used for correct 
    
    def motif_correction(self):
        seq = self.gene_seq.lower()

        #If no u in gene, will replace u's in motif with [ut], else will replace t with [ut]
        #Corrects base pairs for ambiguous pairs

        if 'u' not in seq and 'U' not in seq:
            IUPAC_pairs = {'u': '[ut]', 'r':'[ag]', 'y': '[ctu]','k':'[gtu]','m':'[ac]','s':'[cg]', 'w': '[atu]','b': '[cgtu]','d':'[agtu]','h':'[actu]','v':'[acg]','n':'[acgtu]'}
        else:
            IUPAC_pairs = {'t': '[ut]', 'r':'[ag]', 'y': '[ctu]','k':'[gtu]','m':'[ac]','s':'[cg]', 'w': '[atu]','b': '[cgtu]','d':'[agtu]','h':'[actu]','v':'[acg]','n':'[acgtu]'}

        for key, value in IUPAC_pairs.items():
            self.corr_seq = self.corr_seq.lower().replace(key,value)

    def plot_motif(self,x,y, context,num_motif, total_motif):


        length = len(self.sequence)

        for pos in self.position:

            context.rectangle(x + pos, y - (total_motif - num_motif)*2, length,  (total_motif - num_motif)*4)

            context.fill()

def read_fasta_sequence(file_path):

    # Dictionary with gene name as key and gene object (eventually) as value
    genes = {}

    with open(file_path, 'r') as fasta:

        for line in fasta:
            line = line.strip('\n')

            #Checks for header line
            if line[0] == '>':

                header = line.replace(' ','_').strip('_').strip('>')
                genes[header] = ''
                

            else:
                #concatenate sequence lines from multiple lines in fasta

                genes[header] += line

    for gene, sequence in genes.items():
        genes[gene] = Gene(gene,sequence) #Initialize gene object

    return genes

def read_motifs(file_path, gene_obj):

    # Dictionary where motif is the key and motif object is the value
    motif_dict = {}

    with open(file_path, 'r') as motifs:

        for line in motifs:
            line = line.strip('\n')
            motif_dict[line] = Motif(line,gene_obj)
            motif_dict[line].motif_correction() #Correct sequence for ambiguous bases

    return motif_dict

def plot_genes(gene_dict,fasta_path, bevel = 10):


    #list of colors for motifs
    colors = [(255,0,0), (0,255,0), (0,0,255), (255,255,0), (0,255,255)]

    num_genes = len(gene_dict)

    #Parse fasta name for output image name
    fasta_name = fasta_path.split('/')[-1].strip('.fasta').strip('.fa')


    longest_seq = 0
    longest_motif = 0

    #Find longest sequence and motif for scaling surface
    for genes in gene_dict.values():
        num_motif = len(genes.motif_list)
        if len(genes.sequence) > longest_seq:
            longest_seq = len(genes.sequence)

        for motif in genes.motif_list:
            if len(motif.corr_seq) > longest_motif:
                longest_motif = len(motif.corr_seq)


    # Set size of surface
    width = longest_seq + longest_motif + 50 
    height = (num_genes+1) * 50

    #Minimum surface size
    if width < 100:
        width = 100
    if height < 100:
        height = 100

    surface = cairo.SVGSurface('test.svg',width, height)

    context = cairo.Context(surface)


    #Paints background white
    context.set_source_rgb(1,1,1)
    context.paint()

    #Creates rectangle for legend
    context.set_line_width(5)
    context.rectangle(width-longest_motif - 50, -5, width,  25 + 25*num_motif)


    context.set_line_width(1)

    for index, genes in enumerate(gene_dict.values()):
        context.set_source_rgb(0,0,0)

        # Move to location of read
        x = (longest_seq - len(genes.sequence))/2 + 10
        y = index + 50* (1 +  index)
        

        # Draws line for intron
        context.move_to(x, y)
        context.line_to(x + len(genes.sequence), y)
        context.stroke()
        
        #Writes read name from fasta
        context.move_to(10, y-15)
        context.show_text(genes.gene)


        for exon in genes.exon_list:
            #Draw exon for gene
            exon.plot_exon(x, y, context)

        for num, motif in enumerate(genes.motif_list):

            #Specify colors for a motif
            color1 = colors[num][0]/255
            color2 = colors[num][1]/255
            color3 = colors[num][2]/255

            context.set_source_rgb(color1,color2,color3)
            
            #plots motif
            motif.plot_motif(x, y, context, num, len(genes.motif_list))


            #plots motif legend entry
            leg_x = width - longest_motif - 20
            leg_y = 25 + 25*num

            context.rectangle(leg_x-15,leg_y - 10, 10, 10)
            context.fill()

            context.set_source_rgb(0,0,0)
            context.move_to(leg_x, leg_y)
            context.show_text(motif.sequence)
        
    surface.write_to_png(f"{fasta_name}.png") # Output to PNG

    surface.finish()




#Call function to read fasta
gene_list = read_fasta_sequence(fasta_path)

#Create motif objects in each gene object
for genes in gene_list.values():
    read_motifs(motif_path, genes)

    #Find exons and motifs within a gene
    genes.find_exon()
    genes.find_motif()

#plot everything
plot_genes(gene_list, fasta_path)