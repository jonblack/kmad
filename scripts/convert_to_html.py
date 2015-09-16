#!/usr/bin/python
from itertools import cycle, islice
from sys import argv
# usage: ./convertEncodedToColoured.py infile outfile codon_length mode
# modes:
#   0 - colouring without features
#   1 - motifs colouring
#   2 - domains colouring
#   3 - phosphorylations colouring
#   4 - create multiple files in different colouring modes


def unwrap(alignment):
    new = []
    for i in alignment:
        if i.startswith('>'):
            new.append(i)
            new.append("")
        elif not i.startswith("## PROB"):
            new[-1] += i
        else:
            break
    return new

def srepeat(string, n):
   return ''.join(islice(cycle(string), n))


def beginHTMLdocument(filename):
    template = open("/home/joanna/Documents/scripts/regcol_template.html").read()
    filename.write(template)


def endHTMLdocument(filename):
    htmlCode='</doc>\n</pre>\n</body>\n</html>\n'
    filename.write(htmlCode)


def createHTMLsequence(sequence, codon_length, mode):
    htmlSequence = []
    for i,charI in enumerate(sequence):
        if i%codon_length == 0:
            if "-" != charI and " " != charI:
                if mode == 0: # 0 mode means colouring without features
                        htmlSequence+=["<span class=feat%s>%s</span>" % (charI.upper(),charI)]
                elif codon_length == 7:
                    if mode == 2 and sequence[i+2:i+4] != "AA":
                        htmlSequence+=["<span class=featDOMAIN%s>%s</span>" % (sequence[i+2:i+4],charI)]
                    elif mode == 1 and sequence[i+5:i+7] != "AA" :
                        htmlSequence+=["<span class=featMOTIF%s>%s</span>" % (sequence[i+5:i+7],charI)]
                    elif mode == 3 and sequence[i+4] != 'A':
                        featclass = ptm_dict[sequence[i+4]]
                        htmlSequence+=["<span class=feat{}>{}</span>".format(featclass, charI) ]
                    else:
                        #htmlSequence+=charI
                        htmlSequence+=["<span class=feat%s>%s</span>" % (charI.upper(),charI)]

                else:
                    htmlSequence+=["<span class=feat%s>%s</span>" % (charI.upper(),charI)]
            elif charI == " ":
                htmlSequence+=[" "]
            else:
                htmlSequence+=["-"]
    return "".join(htmlSequence)+"\n"
def numbering(seqLength, codon_length):
    numbers = range(1,seqLength/codon_length)
    result = srepeat(" ",15) + '0'
    for i in numbers:
        if i % 5 == 0:
            result += srepeat(" ", 5-len(str(i))) + str(i)
    return result+"\n"

cod_length = int(argv[3])
mode = int(argv[4])
infile = unwrap(open(argv[1]).read().splitlines())
newline = numbering(len(infile[1]),cod_length)
if mode != 4:
    outfile = open(argv[2],'w')
    beginHTMLdocument(outfile)
    outfile.write(newline)
else:
    outfile_dom = open(argv[2]+"_dom",'w')
    outfile_phosph = open(argv[2]+"_phosph",'w')
    outfile_motifs = open(argv[2]+"_motifs",'w')
    outfile_reg = open(argv[2]+"_reg",'w')
    beginHTMLdocument(outfile_dom)
    beginHTMLdocument(outfile_phosph)
    beginHTMLdocument(outfile_motifs)
    beginHTMLdocument(outfile_reg)
    outfile_dom.write(newline)
    outfile_phosph.write(newline)
    outfile_motifs.write(newline)
    outfile_reg.write(newline)
for i in range(len(infile)):
    if "## PROBABILITIES" in infile[i]:
        break
    if i%2 == 0:
        seqID = str((i/2) + 1)+"."+infile[i].split()[0].split('|')[-1]
        if len(seqID) > 13:
            seqID = seqID[0:12]
        newline=seqID+srepeat(" ",15-len(seqID))
        if mode != 4:
            outfile.write(newline)
        else:
            outfile_dom.write(newline)
            outfile_phosph.write(newline)
            outfile_motifs.write(newline)
            outfile_reg.write(newline)
    else:
        if mode != 4:
            outfile.write(createHTMLsequence(infile[i].rstrip('\n'), cod_length, mode))
        else:
            outfile_reg.write(createHTMLsequence(infile[i].rstrip('\n'), cod_length, 0))
            outfile_motifs.write(createHTMLsequence(infile[i].rstrip('\n'), cod_length, 1))
            outfile_dom.write(createHTMLsequence(infile[i].rstrip('\n'), cod_length, 2))
            outfile_phosph.write(createHTMLsequence(infile[i].rstrip('\n'), cod_length, 3))
if mode != 4:
    endHTMLdocument(outfile)
else:
    endHTMLdocument(outfile_dom)
    endHTMLdocument(outfile_reg)
    endHTMLdocument(outfile_phosph)
    endHTMLdocument(outfile_motifs)
if __name__ == "__main__":

