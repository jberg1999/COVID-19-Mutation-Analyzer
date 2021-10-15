import sys
import pandas
import numpy
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Align

# AQUIRES RANGE OF SEQUENCES TO WORK ON BASED ON ARGUMENTS
pID = sys.argv[1]
start = int(sys.argv[2])
end = int(sys.argv[3])

print(pID, "started")

# AQUIRES THE REFERENCE SEQUENCE
seqparse = SeqIO.parse("refseq.fasta", "fasta")
refseq = next(seqparse)
refseq.id = refseq.id.split(".")[0]
#print(refseq.id)

# CREATES THE MUTATIONS TABLE
mutations=pandas.DataFrame(columns=['Position','Mutation_Id','Type','New,Seq_Id',
'Date','Location','Amino_Acid','Gene','Protein','Week','Month','Year',
'Effect_Type','Effect,Region','Country'])


#GETS THE SEQUENCES TABLE
sequences= pandas.read_csv("seq"+pID+".csv")
pandas.set_option("precision",0)





# FUNCTION THAT SETS UP THE ALIGNER OBJECT TO USE NCBI BLAST SCORING
def alignerSetUp():
    aligner = Align.PairwiseAligner()
    aligner.match_score = 1
    aligner.mismatch_score = -2
    aligner.gap_score = -2.5
    aligner.mode = 'local'
    return aligner

#SET UP ALIGNER
aligner = alignerSetUp()



'''DATA CLEAN IS A FUNCTION THAT IS USED TO MODIFY DATA IN THE NEW DATAFRAME SO THAT IT MATCHES
THE FORMAT OF THE EXISTING SEQUENCES DATAFRAME THAT IT WILL BE ADDED TO LATER. IT CREATES NEW
COLUMNS SUCH AS TOTAL_MUTATIONS AND N_COUNT THAT WILL BE FILLED LATER AND PARSES THE LOCATION
AND DATE COLUMNS TO BE MORE USEFULL.
'''

def data_clean(newSequences):
    newSequences["date"] = pandas.to_datetime(newSequences["Collection_Date"])
    newSequences["Month"] = pandas.DatetimeIndex(newSequences["date"]).month
    newSequences["Year"] = pandas.DatetimeIndex(newSequences["date"]).year
    newSequences["DayOfWeek"] = pandas.DatetimeIndex(newSequences["date"]).weekday
    newSequences["Week"] = pandas.DatetimeIndex(newSequences["date"]).weekofyear
    newSequences["Week"] = newSequences.apply(lambda row: weekGetter(row), axis=1)

    #NEW COLUMNS TO BE FILLED AFTER SEQUENCE IS ANALYSED
    newSequences["Insertion_Count"] = numpy.nan
    newSequences["Deletion_Count"] = numpy.nan
    newSequences["Substitution_Count"] = numpy.nan
    newSequences["Total_Mutations"] = numpy.nan
    newSequences["Insertions"] = numpy.nan
    newSequences["Deletions"] = numpy.nan
    newSequences["Substitutions"] = numpy.nan
    newSequences["Aligned"] = numpy.nan
    newSequences["N_Count"] = numpy.nan  # ambiguous bases
    return newSequences



"""
THIS FUNCTION IS RESPONSIBLE FOR ALIGNING EACH NEW SEQUENCE AND LOCATING ALL OF
THE MUTATIONS IN EACH. IT USES THE ALIGNER TO PERFORM LOCAL ALIGNMENT, AND RUNS THROUGH
EACH ALIGNMENT. IT LOCATES MUTATIONS AND CALLS GENE_FINEDER, GENE.PEPTIDEFINDER,
AND PEPTIDE.MUTATIONEFFECT FOR DETAILED MUTATION INFORMATION. IT CALLS MUTADDER TO ADD MUTATIONS
TO THE MUTATION TABLE, AND ONCE THE ALIGNMENT IS DONE, IT TALLIES SEQUENCE LEVEL INFORMATION AND
ADDS THE SEQUENCE INTO THE MASTER SEQUENCES TABLE.
"""
def alignmentHandler(ref, seqdata):
    global mutations
    # THESE LINES GET ALL POSSIBLE EQUALLY SCORING ALIGNMENTS AND CHOSE THE FIRST ONE
    alignments = aligner.align(seqdata["sequence"], ref)
    alignment = alignments[0]

    # VARIABLES TALLY THE COUNT OF EACH OF THESE TYPES OF MUTATIONS
    insertions = []
    substitutions = []
    deletions = []
    ref = alignment.query
    subject = alignment.target
    '''
    BIOPYTHON ALIGNMENTS ARE STORED AS HITZONES LIKE THIS:
    [[SUBJECT HITS],[REFERENCE HITS]]

    [[(1-5),(5-8)],[(1-4)(5-9)]]

    EACH TUPLE REPRESENTS A RANGE WITHIN EITHER THE SUBJECT OR REFERENCE THAT IS
    UNINTERUPTED BY INDELS BUT CAN STILL HAVE SUBSTITUTIONS. A INSERTION IN A SEQUENCE WILL CAUSE
    THAT BASE TO NOT BE MATCHED BY THE OTHER AND WILL BREAK THE HITZONES OF BOTH SEQUENCES.
    THE APPROACH IS TO SCAN EACH PAIR OF HITZONES FOR SUBSTITUTIONS AND THEN TO LOOK AT THE GAP
    BETWEEN BOUNDARIES OF HITZONES TO SEE WHICH SEQUENCE HAS THE UNMATCHED BASE THEREFRE WETHER
    THE SUBJECT HAS AN INSERTION OR DELETION RELATIVE TO THE REFERENCE. IN THE ABOVE CASE, THE
    REFENCE HAS A GAP OF 1 SO THERE IS A DELETION IN THE SUBJECT SEQUENCE.
    '''
    # ASIGNS THE SUBJECT AND REFENCE HITZONE FROM THE ALIGNMENT OBJECT
    subHits = alignment.aligned[0]
    refHits = alignment.aligned[1]
    N = 0 # NUMBER OF AMBIGUOUS BASES IN A SEQUENCE

    # THIS OUTER FOR LOOP RUNS FOR ALL THE PAIRS OF HITZONES
    for zone in range(0, len(subHits)):
        subZoneStart = subHits[zone][0] # START OF ZONE IN SUBJECT
        refZoneStart = refHits[zone][0] # START OF ZONE IN REFERENCE


        # THIS LOOPS OVER EACH BASE IN THE HITZONE AND CHECKS FOR MISMATCHES
        for b in range(0, (subHits[zone][1] - subHits[zone][0])):
            position = refZoneStart + b
            newBase = subject[subZoneStart + b]
            refBase = ref[position]

            if newBase != refBase:  # SUBSTITUTION
                position = position + 1 # BIOLOGY STARTS AT 1

                #USED FOR MUT_ID SORTABILITY IN WEBSITE
                string_position = str(position)
                string_position = (5-len(string_position))*"0" + string_position



                if (newBase in "GATC"):  # VALID SUBSTITUTION
                    Mutation = "Sub: " + string_position + " " + refBase + " -> " + newBase #FORMATING ID COL
                    substitutions.append([position, newBase])

                    # CHECK FOR SPECIAL CASE THAT MUTATION IS IN MULIPLE GENES
                    if position in range(27756, 27760): # NUCS 27756-59 ARE IN BOTH ORF7a AND 7b
                        mutAdder(position, Mutation, "Substitution", newBase, seqdata, "ORF7a & ORF7b", "muliple",
                                 "multiple", "multiple","multiple")

                    # ELSE, PROCEED AS NORMAL BY GETTING INFORMATION ABOUT SUBSTITUTION
                    else:
                        gene = geneFinder(position) # FINDS WHICH GENE MUTATION IS IN
                        if gene != "Non-coding":
                            # GENE IS CODING
                            # GETS PROTEIN AND EFFECT OF THE MUTATION
                            protein = gene.peptideFinder(position)
                            details = protein.mutationEffect(position, "Substitution", newBase)
                            AA = ((position - protein.start) // 3) + 1

                            #ADDS MUTATION TO THE MASTER TABLE
                            mutAdder(position, Mutation, "Substitution", newBase, seqdata, gene.name, protein.name, AA, details[0],
                                     details[1])
                        else:
                            #MUTATION IS NON-CODING
                            mutAdder(position, Mutation, "Substitution", newBase, seqdata, gene=gene)

                else:
                    N += 1

        # EXAMINES GAP BETWEEN ZONES TO IDENTIFY WHICH SEQUENCE HAS AN EXTRA BASE
        # WE EXCLUDE AREAS OUTSIDE OF THE EXTREMES OF THE HITZONES

        if zone < len(subHits) - 1:  #IF THE ZONE IS NOT THE LAST ONE
            # NOTE THAT SUBGAP > 1 IS ACTUALLY A MISSING BASE IN THE REFERENCE AND VIS VERSA
            subGap = subHits[zone + 1][0] - subHits[zone][1]
            refGap = refHits[zone + 1][0] - refHits[zone][1]

            if subGap != 0:  # INSERTION DUE TO UNMATCHED SUBJECT BASES
                new = subject[subHits[zone][1]:subHits[zone + 1][0]] # ALL OF THE BASES BETWEEN ZONES
                position = refHits[zone][1] + 1
                insertions.append([position, new])
                string_position = str(position)
                string_position = (5-len(string_position))*"0" + string_position
                Mutation = "Ins: " + string_position + " + " + new
                t = "Insertion"


            elif refGap != 0:  # DELETION DUE TO MISSING SUBJECT BASES
                new = "-" * refGap # LENGTH OF DELETION
                position = refHits[zone][1] + 1
                deletions.append(position)
                string_position = str(position)
                string_position = (5-len(string_position))*"0" + string_position
                Mutation = "Del: " + str(string_position) + " - " + str(position + (len(new)-1))
                t = "Deletion"


            # CHECK FOR SPECIAL CASE OF DEL IN MULTIPLE proteins
            if position in range(27756, 27760): # NUCS 27756-59 ARE IN BOTH ORF7a AND 7b
                        mutAdder(position, Mutation, t, new, seqdata, "ORF7a & ORF7b", "muliple",
                                 "multiple", "multiple","multiple")

            # PROCEED AS NORMAL
            else:
                end_pos = position + len(new) - 1
                gene = geneFinder(position)
                gene_end = geneFinder(end_pos)

                # SOME RARE POSSIBILITIES NOT ACCOUNTED FOR

                #option 1: non coding into coding. solution: new start for effect will be start of new gene? promoter concerns?
                # translate after deletion. if first base of first aug is out of old frame with oud aug then loss. Otherwise calculate aa loss
                #option 2: coding into non coding. solution: nonstop mutation effect?
                #option 3: different protiens from same gene. solution: frameshift or clevage site alteration?
                #option 4: different proteins from different genes. solution: nonstop and promoter concerns?


                if gene != "Non-coding":
                    protein = gene.peptideFinder(position)
                    details = protein.mutationEffect(position, t, new)
                    AA = ((position - protein.start) // 3) + 1
                    mutAdder(position,Mutation, t, new, seqdata, gene.name, protein.name, AA, details[0],
                             details[1])
                else:
                    mutAdder(position, Mutation, t, new, seqdata, gene=gene)


    # AFTER ALIGNMENT HAS BEEN ANALYSED, DATA IS WRITTEN INTO THE ROW AND IT IS ADDED TO THE TABLE
    global sequences
    line = seqdata
    line["Insertion_Count"] = len(insertions)
    line["Deletion_Count"] = len(deletions)
    line["Substitution_Count"] = len(substitutions)
    line["Total_Mutations"] = line["Insertion_Count"] + line["Deletion_Count"] + line["Substitution_Count"]
    line["Insertions"] = insertions
    line["Deletions"] = deletions
    line["Substitutions"] = substitutions
    try:
        line["Aligned"] = len(alignments)
    except OverflowError:
        line["Aligned"] = 9223372036854775807

    line["N_Count"] = N
    line["Start"] = refHits[0][0]
    line["End"] = refHits[-1][1]
    line["Alignment_Length"] = line["End"] - line["Start"]

    sequences.loc[line['id']] = line  # OVERRIDE EXISTING


'''
TAKES INPUT DATA FROM ALIGNMENT HANDLER AND ADDS IT TO THE MUTATION DATAFRAME
'''
def mutAdder(pos, Mutation, type, new, data, gene=None, protein=None, AA=None, effectType=None, effect=None):
    global mutations
    mutations = mutations.append({"Position": pos, "Mutation_Id": Mutation, "Type": type, "New": new,"Seq_Id": data["id"], "Date": data["date"],
                                  "Location": data["Location"], "Month": data["Month"], "Year": data["Year"],
                                  "Week": data["Week"], "Gene": gene, "Protein": protein, "Amino_Acid": AA,
                                  "Effect_Type": effectType, "Effect": effect}, ignore_index=True)


"""
IS USED TO DELETE ANY EXISTING MUTATIONS FROM A GIVEN SEQUENCE. THIS OCCURS WHEN
A SEQUENCE THAT HAS ALREADY BEEN PROCESSED IS PROCESSED AGAIN DUE TO IT BEING UPDATED
ON NCBI OR IF THE PROGRAM IS RERUN
"""
def mutationCleanse(id):
    global mutations
    mutations = mutations[mutations.Seq_Id != id]

"""
GETS THE WEEK AS A NUMBER RELATIVE TO WHEN THE OUTBREAK STARTED IN DECEMBER
"""
def weekGetter(row):
    if row["Year"] == 2019:
        row["Week"] = 1
    else:
        row["Week"] = 52 * (row["Year"] - 2020) + (row["Week"] + 1)
    return row["Week"]


# MUTATION CLASSIFICATIONS

'''OPEN READING FRAME CLASS IS USED TO  IDENTIFY WHICH ORF A MUTATION IS IN.
IT HAS A FOR IDENTIFYING WHICH PROTEIN A MUTATION WOULD BE IN
'''
class ORF:
    def __init__(self, name, start, end):
        self.start = start
        self.end = end
        self.name = name
        self.peptides = []

    def setPeptides(self, peptides):
        self.peptides = peptides

    #CALLED IN ALIGNMENT HANDLER TO GET THE PEPTIDE OF A MUTATION
    def peptideFinder(self, pos):
        for pep in self.peptides:
            if pos >= pep.start and pos <= pep.end:
                result = pep
        return result

# USED TO DETERMINE THE AFFECT A MUTAION WOULD HAVE ON A PEPTIDES AMINO ACID SEQUENCE
class Protein:
    def __init__(self, name, start, end):
        self.start = start
        self.end = end
        self.name = name

    '''
    USED TO GET THE MUTATION EFFECT IN ALIGNMENT HANDLER. BASED ON WHAT TYPE A
    MUTATION IS KNOWN TO BE AND ITS LOCATION, IT WILL CALCULATE ANY CODON CHANGES
    OR FRAME SHIFTS
    INPUTS: POSITION, TYPE OF MUTATION, NEW BASES
    '''
    def mutationEffect(self, pos, type, new):
        '''WITH INDELS, THE MIAN ISSUE IS WHETHER OR NOT AN INDEL CAUSES A FRAMSHIFT
        WE CALCULATE THE FRAME SHIFT BY THE REMAINDER OF THE LENGTH OF THE INDEL / 3
        '''
        if type == "Insertion":
            offset = len(new) % 3
            if offset == 0:
                #NO FRAME SHIFT
                effectType = "In-frame Insertion"
                x = Seq(new)
                x = x.translate()
                effect = "insertion of " + x # SHOULD BE A SERIES OF AMINO ACIDS
            else:
                # FRAME SHIFT
                effectType = "Frameshift"
                effect = "+" + str(offset) + "frameshift"

        elif type == "Deletion":
            offset = len(new) % 3
            if offset == 0:
                effectType = "In-frame deletion"
                effect = str(len(new)//3) + " AA deletion"
            else:
                effectType = "Frameshift"
                effect = "-" + str(offset) + "frameshift"
        else:
            # SUBSTITUTION
            '''DETEMINES THE POSITION WITHEN A CODON BY TAKING REMAINDER OF THE STARTING
            POSITION OF THE PEPTIDE AND THE MUTATION POSITION. (0,1,OR,2)
            '''
            slot = (pos-self.start)%3
            codonStart = pos - slot # STARTING POSITION OF CODON IN REFERENCE SEQUENCE

            #GETS THE CODON BASES FROM REFERENCE SEQUENCE
            bases = refseq[codonStart - 1:codonStart + 2] #position of refseq starts at 0 so is off by 1
            old = bases.translate() # OLD AMINO ACID
            bases = bases.seq.tomutable() # NESSISARY TO REPLACE THE BASE WITH THE CHANGED BASE
            bases[slot] = new # CHANGES THE BASE IN THE APPROPRIATE SLOT INTO THE NEW BASE
            bases = bases.toseq()
            changed = str(bases.translate()) # NEW AMINO ACID
            old = str(old.seq) # NEEDED FOR COMPARISON

            if changed == old: # NO EFFECT
                effectType = "Silent"
                effect = "None"
            elif changed == "*": # SYMBOL FOR STOP CODON
                effectType = "Nonsense"
                effect = "Premature stop"
            elif old == "*": # STOP CODON IS OVERWRITTEN
                effectType = "Nonstop"
                effect = "Extended peptide"
            else: # AMINO ACID IS CHANGED
                effectType = "Missense"
                effect = old + " -> " + changed

        return [effectType, effect]


# USED TO FIND THE GENE A MUTATION IS IN. PROBABLY COULD BE REPLACED BY SOMETHING BETTER
def geneFinder(pos):
    gene = "Non-coding"
    if pos <= 21555:
        if pos >= 266:
            gene = ORF1ab
    else:
        if pos <= 25384:
            if pos >= 21563:
                gene = S
        elif pos <= 26220:
            if pos >= 25393:
                gene = ORF3a
        elif pos <= 26472:
            if pos >= 26245:
                gene = E
        elif pos <= 27191:
            if pos >= 26523:
                gene = M
        elif pos <= 27387:
            if pos >= 27202:
                gene = ORF6
        elif pos <= 27759:
            if pos >= 27756:
                gene = "Both!!!"  # need to figure out how to handle this...
            elif pos >= 27394:
                gene = ORF7a
        elif pos <= 27887:
            if pos >= 27759:
                gene = ORF7b
        elif pos <= 28259:
            if pos >= 27894:
                gene = ORF8
        elif pos <= 29533:
            if pos >= 28274:
                gene = N
        elif pos <= 29674:
            if pos >= 29558:
                gene = ORF10
    return gene


#CODE BELOW CREATES OBJECTS FOR PROTEINS AND ORFS BASED ON NCBI REFERENCE GENOME DATA

# GENES
ORF1ab = ORF("ORF1ab", 266, 21555)
S = ORF("S", 21563, 25384)
ORF3a = ORF("ORF3a", 25393, 26220)
E = ORF("E", 26245, 26472)
M = ORF("M", 26523, 27191)
ORF6 = ORF("ORF6", 27202, 27387)
ORF7a = ORF("ORF7a", 27394, 27759)
ORF7b = ORF("ORF7b", 27756, 27887)
ORF8 = ORF("ORF8", 27894, 28259)
N = ORF("N", 28274, 29533)
ORF10 = ORF("ORF10", 29558, 29674)


#PROTEINS

# CREATES AND ASSIGNS PROTEINS TO APPROPRIATE GENES
a = Protein("YP_009725297.1", 266, 805)
b = Protein("YP_009725298.1", 806, 2719)
c = Protein("YP_009725299.1",2720,8554)
d = Protein("YP_009725300.1", 8555, 10054)
e = Protein("YP_009725301.1", 10055, 10972)
f = Protein("YP_009725302.1", 10973, 11842)
g = Protein("YP_009725303.1", 11843, 12091)
h = Protein("YP_009725304.1", 12092, 12685)
i = Protein("YP_009725305.1", 12686, 13024)
j = Protein("YP_009725306.1", 13025, 13441) # last protein in both ab and a
k = Protein("YP_009725307.1", 13442, 16236)
l = Protein("YP_009725308.1", 16237, 18039)
m = Protein("YP_009725309.1", 18040, 19620)
n = Protein("YP_009725310.1", 19621, 20658)
o = Protein("YP_009725311.1", 20659, 21552)
ORF1ab.setPeptides([a,b,c,d,e,f,g,h,i,j,k,l,m,n,o])

a = Protein("YP_009724390.1", 21563, 25384)
S.setPeptides([a])

a = Protein("YP_009724391.1", 25393, 26220)
ORF3a.setPeptides([a])

a = Protein("YP_009724392.1", 26245, 26472)
E.setPeptides([a])

a = Protein("YP_009724393.1", 26523, 27191)
M.setPeptides([a])

a = Protein("YP_009724394.1", 27202, 27387)
ORF6.setPeptides([a])

a = Protein("YP_009724395.1", 27394, 27759)
ORF7a.setPeptides([a])

a = Protein("YP_009725318.1", 27756, 27887)
ORF7b.setPeptides([a])

a = Protein("YP_009724396.1", 27894, 28259)
ORF8.setPeptides([a])

a = Protein("YP_009724397.2", 28274, 29533)
N.setPeptides([a])

a = Protein("YP_009725255.1", 29558, 29674)
ORF10.setPeptides([a])



# PREPS SEQUENCE TABLE FOR ANALYSIS
sequences = data_clean(sequences)

#RUNS ANALYSIS
processed = 0
print("process", pID, "about to start aligning!")
for index, row in sequences.iterrows():
    alignmentHandler(refseq.seq, row)

    #AUTOSAVE USED FOR PROCESSING LARGE BATCHES
    if processed % 1000 == 0:
        print("process", pID,': processed', processed)
        sequences.to_csv("seq"+pID+".csv")
        mutations.to_csv("mut"+pID+".csv", index=False)
        with open("log.txt", "a") as file:
            file.write(pID+" | ")
            file.write(pandas.datetime.now().strftime("%m/%d/%Y, %H:%M:%S"))
            file.write(" | ")
            file.write(str(processed))
            file.write("\n")
    processed += 1



# FINAL OUTPUT IS WRITTEN TO TWO FILES
sequences.drop(columns=["sequence"]).to_csv("seq"+pID+".csv")
mutations.to_csv("mut"+pID+".csv", index=False)
