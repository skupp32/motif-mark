import cairo
import argparse
import re
import math

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
        self.gene = gene_name
        self.sequence = sequence 
        self.exon_list = []
        self.motif_list = []

    def find_exon(self):
        exon_start = []
        exon_end = []

        for index in range(len(self.sequence)-1):


            if self.sequence[index].islower() and self.sequence[index+1].isupper():
                exon_start.append(index+1)

            elif self.sequence[index].isupper() and self.sequence[index+1].islower():
                exon_end.append(index+1)


        if self.sequence[0].isupper():
            exon_start.append(0)
        if self.sequence[-1].isupper():
            exon_end.append(len(self.sequence))

        
        for index in range(len(exon_start)):

            start = exon_start[index]
            end = exon_end[index]
            exon_seq = self.sequence[start:end]

            self.exon_list.append(Exon(exon_seq, start,self))




    def find_motif(self):

        for motifs in self.motif_list:
            motif = motifs.corr_seq

            for found_motif in re.findall(motif, self.sequence, re.IGNORECASE):

                for pos in re.finditer('(?={0})'.format(re.escape(found_motif)),self.sequence, re.IGNORECASE):
                    motifs.position.append(pos.start())

        

class Exon:

    def __init__(self, sequence, position,read_obj):
        self.sequence = sequence
        self.position = position
        read_obj.exon_list.append(self)

    def plot_exon(self, x, y, context):

        context.rectangle(x + self.position, y - 12.5, len(self.sequence), 25)
        context.fill()

class Motif:

    def __init__(self, sequence,gene_obj):
        self.sequence = sequence
        self.corr_seq = sequence
        self.position = []
        gene_obj.motif_list.append(self) 
        self.gene_seq = gene_obj.sequence 
    
    def motif_correction(self):
        seq = self.gene_seq.lower()
        #IUPAC_pairs = {'u': '[ut]','t': '[ut]', 'r':'[ag]', 'y': '[ctu]','k':'[gtu]','m':'[ac]','s':'[cg]', 'w': '[atu]','b': '[cgtu]','d':'[agtu]','h':'[actu]','v':'[acg]','n':'[acgtu]'}

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
    genes = {}
    with open(file_path, 'r') as fasta:

        for line in fasta:
            line = line.strip('\n')
            if line[0] == '>':

                header = line.replace(' ','_').strip('_').strip('>')
                genes[header] = ''
                

            else:

                genes[header] += line

    for gene, sequence in genes.items():
        genes[gene] = Gene(gene,sequence)

    return genes

def read_motifs(file_path, gene_obj):

    motif_dict = {}

    with open(file_path, 'r') as motifs:

        for line in motifs:
            line = line.strip('\n')
            motif_dict[line] = Motif(line,gene_obj)
            motif_dict[line].motif_correction()

    return motif_dict

def plot_genes(gene_dict,fasta_path, bevel = 10):

    colors = [(255,0,0), (0,255,0), (0,0,255), (255,255,0), (0,255,255)]
    num_genes = len(gene_dict)

    fasta_name = fasta_path.split('/')[-1].strip('.fasta').strip('.fa')


    longest_seq = 0
    longest_motif = 0

    for genes in gene_dict.values():
        num_motif = len(genes.motif_list)
        if len(genes.sequence) > longest_seq:
            longest_seq = len(genes.sequence)

        for motif in genes.motif_list:
            if len(motif.corr_seq) > longest_motif:
                longest_motif = len(motif.corr_seq)


    
    width = longest_seq + longest_motif + 50 
    height = (num_genes+1) * 50

    if width < 100:
        width = 100
    if height < 100:
        height = 100

    surface = cairo.SVGSurface('test.svg',width, height)

    context = cairo.Context(surface)

    context.set_source_rgb(1,1,1)
    context.paint()

    context.set_line_width(5)
    context.rectangle(width-longest_motif - 50, -5, width,  25 + 25*num_motif)


    context.set_line_width(1)

    for index, genes in enumerate(gene_dict.values()):
        context.set_source_rgb(0,0,0)

        x = (longest_seq - len(genes.sequence))/2 + 10
        y = index + 50* (1 +  index)
        
        context.move_to(x, y)
        context.line_to(x + len(genes.sequence), y)
        context.stroke()
        
        context.move_to(10, y-15)
        context.show_text(genes.gene)


        for exon in genes.exon_list:
            exon.plot_exon(x, y, context)

        for num, motif in enumerate(genes.motif_list):
            color1 = colors[num][0]/255
            color2 = colors[num][1]/255
            color3 = colors[num][2]/255

            context.set_source_rgb(color1,color2,color3)
            
            motif.plot_motif(x, y, context, num, len(genes.motif_list))

            leg_x = width - longest_motif - 20
            leg_y = 25 + 25*num

            context.rectangle(leg_x-15,leg_y - 10, 10, 10)
            context.fill()

            context.set_source_rgb(0,0,0)
            context.move_to(leg_x, leg_y)
            context.show_text(motif.sequence)
        
    surface.write_to_png(f"{fasta_name}.png") # Output to PNG

    surface.finish()





gene_list = read_fasta_sequence(fasta_path)

for genes in gene_list.values():
    read_motifs(motif_path, genes)
    genes.find_exon()
    genes.find_motif()

plot_genes(gene_list, fasta_path)