import os,sys
from Bio import SeqIO
from Bio import Seq
import matplotlib.pyplot as plt

read1 = sys.argv[1]
read2 = sys.argv[2]
condaDir = sys.argv[3]
cutoff = int(sys.argv[4])
outputFolder = sys.argv[5]

#Creating jellyfish Database
os.system("mkdir "+outputFolder)
#print("\n"+condaDir+"/bin/jellyfish count -m 17 -s 100M -C -t 6 -o "+outputFolder+"/mer_counts.jf "+read1+" "+read2)
os.system(condaDir+"/bin/jellyfish count -m 17 -s 100M -C -t 6 -o "+outputFolder+"/mer_counts.jf "+read1+" "+read2)

#Collecting longest sequence per genotype
longestSeq = {}
for seq_record in SeqIO.parse("geneClassification.fasta","fasta"):
    idList = str(seq_record.id).split("_")
    if not (idList[0],idList[-1]) in longestSeq:
        longestSeq[(idList[0],idList[-1])] = ""
    if len(str(seq_record.seq))> len(longestSeq[(idList[0],idList[-1])]):
        longestSeq[(idList[0],idList[-1])] = str(seq_record.seq)
    


hyperVariableGenes = ['rl12','ul9','ul1','ul120','rl5a','ul74','ul73','ul146','ul20','rl6','ul139','ul11','rl13']

for gene in hyperVariableGenes:
    alignmentLength = 0
    kmerCountFile = open(outputFolder+"/"+gene+"_kmerCounts.txt","w")
    
    totReadsPerGenotype = {}
    kmerDict = {}
    dataToPlot = {}
    foundGenotypes = set()

    genotypeFasta = open(outputFolder+"/genotypeFasta.fasta","w")

    #Load in memory all kmers of each genotype of the gene
    infile = open("mainDB_seqs_filtered.txt")
    while True:
        line = infile.readline().rstrip()
        if not line:
            break
        fields = line.split("\t")
        if fields[0] == gene:
            if not fields[1] in kmerDict:
                kmerDict[fields[1]] = fields[2]
    infile.close()

    #find all the kmer position in the longest sequence for a genotype
    kmerPos = {}
    kmerPosList = {}
    for genotype in kmerDict:
        
        print("positioning kmers for genotype",genotype)
        if not genotype in totReadsPerGenotype:
            totReadsPerGenotype[genotype] = 0
        if not genotype in kmerPos:
            kmerPos[genotype] = {}
        if not genotype in kmerPosList:
            kmerPosList[genotype]=[]
        for kmer in kmerDict[genotype].split(","):
            if not kmer in kmerPos[genotype]:
                kmerPos[genotype][kmer] = longestSeq[(gene,genotype)].find(kmer)
            kmerPosList[genotype].append(longestSeq[(gene,genotype)].find(kmer))

    
        
    #Search all the kmers of each genotype within the reads with jellyfish
    for genotype in kmerDict:
        kmerCountFile.write(genotype+":\n")
        print("Analyzing combination %s / %s combination" %(gene,genotype))
        os.system(condaDir+"/bin/jellyfish query "+outputFolder+"/mer_counts.jf "+kmerDict[genotype].replace(","," ")+" >"+outputFolder+"/jellyfishOutput.txt")

        infile = open(outputFolder+"/jellyfishOutput.txt")
        while True:
            line = infile.readline().rstrip()
            if not line:
                break
            fields = line.split(" ")
            
            #if a kmer is found with occurances > cutoff the corresponding sequence is reported in fasta file
            if int(fields[1])>cutoff:
                kmerCountFile.write(fields[0]+"\t"+fields[1]+"\n")
                foundGenotypes.add(genotype)
                if not genotype in dataToPlot:
                    dataToPlot[genotype] = []
                if fields[0] in kmerPos[genotype]:
                    dataToPlot[genotype].append((kmerPos[genotype][fields[0]],int(fields[1])))
                    totReadsPerGenotype[genotype]+=int(fields[1])
                else:
                    dataToPlot[genotype].append((kmerPos[genotype][Seq.reverse_complement(fields[0])],int(fields[1])))
                    totReadsPerGenotype[genotype]+=int(fields[1])
        infile.close()
    
    kmerCountFile.close()




    #Write the longest sequence for each found genotype
    for genotype in foundGenotypes:
        genotypeFasta.write(">"+genotype+"\n"+longestSeq[(gene,genotype)]+"\n")
    genotypeFasta.close()

    #Align the found sequences with mafft
    os.system(condaDir+"bin/mafft --auto "+outputFolder+"/genotypeFasta.fasta > "+outputFolder+"/alignedGenotypeFasta.fasta")

    #Calculate the new coordinates of the kmers in the aligned sequences
    kmerNewPosDict = {}
    
    for seq_record in SeqIO.parse(outputFolder+"/alignedGenotypeFasta.fasta","fasta"):
        alignedGenotype = str(seq_record.id)
        alignmentLength = len(str(seq_record.seq))
        if not alignedGenotype in kmerNewPosDict:
            kmerNewPosDict[alignedGenotype] = {}
        oldPos = 0
        sequence = str(seq_record.seq)
        for a in range(len(sequence)):
            if not sequence[a]=='-':
                if not oldPos in kmerNewPosDict[alignedGenotype]:
                    kmerNewPosDict[alignedGenotype][oldPos]=a
                    oldPos+=1
    

    

    #Create data for genotypes ideogram
    kmerPositionIdeogram = {}
    observedGenotype = len(dataToPlot)+1
    
    for genotype in dataToPlot:
        observedGenotype-=1
        if not genotype in kmerPositionIdeogram:
            kmerPositionIdeogram[genotype] = []
        for position in kmerPosList[genotype]:
            kmerPositionIdeogram[genotype].append((kmerNewPosDict[genotype][position],observedGenotype))


    observedGenotype = len(dataToPlot)+1



    #Plot the kmer coverage per position
    fig,(ax1,ax2,ax3) = plt.subplots(3,1,figsize=(10,10),gridspec_kw={'height_ratios': [10, 1,0.5]})
    for genotype in dataToPlot:
        x = []
        y = []
        for coordinates in dataToPlot[genotype]:
            x.append(kmerNewPosDict[genotype][coordinates[0]])
            y.append(coordinates[1])
        ax1.bar(x,y,alpha=0.3,label=genotype)
    
    ax1.legend()
    ax1.set_xlabel("Sequence position",fontsize=12,fontweight='bold')
    ax1.set_ylabel("Number of reads",fontsize=12,fontweight='bold')
    ax1.set_xlim(0,alignmentLength)


    for genotype in dataToPlot:
        x = []
        y = []
        for coordinates in kmerPositionIdeogram[genotype]:
            x.append(coordinates[0])
            y.append(coordinates[1])

        ax2.scatter(x,y,label=genotype,marker="s")
        ax2.plot((0,alignmentLength),(y[0],y[0]))
    ax2.set_xlim(0,alignmentLength)

        
    
    ax2.axis('off')
    #ax2.legend()
    
    ax2.set_ylim(0,observedGenotype+1)
    ax1.set_title(gene,fontsize=12,fontweight='bold')

    freqString = "Average kmer coverage:\n "
    ax3.text(0,1.1,freqString,fontsize=10,fontweight='bold')
    yPos = 1.1
    for genotype in dataToPlot:
        yPos -= 0.5
        ax3.text(0,yPos,genotype+":",fontsize=10,fontweight='bold')
        ax3.text(0.2,yPos,(str(float(totReadsPerGenotype[genotype])/float(len(kmerPosList[genotype])))[:5]),fontsize=10 )


    ax3.axis('off')
    fig.savefig(outputFolder+"/"+gene+"_kmerCov.png",dpi=300)



                
