import conan_constants, os
from conan_node import Node
from conan_bool_helpers import isHMM
from conan_aln_helpers import readFasta, readStockholm, maxIdCDhit
from conan_aln_filters import hmmFiltering, residueFiltering, seqbasedFiltering
from conan_graph_helpers import getNodes, getNodeFrequency, writeSeqmodeNetwork, nodesFiltering, genCorrelationNetwork, genCorrelationNetworkWithAntiCorr, filterRedundancy, correctedNresidues, writeCommunities, writeBackbone, getFullJacardCoefficient


def printHelp():
	print("###########################################################################\n"
		"CONAN - Co-Variation Network Analyzer\n"
		"###########################################################################\n\n"
		"This software analyzes a residue co-variation network with the aim to emphasize local evolutionary constraints and "
		"consequently detect functional and specificity determinant sites. Its only mandatory input is a multiple sequence "
		"alignment in one of the following formats: Pfam, Stockholm or FASTA; and the output consists of several networks and "
		"communities files from several cuts of the network.\n\n"
		"The mandatory inputs consits of:\n"
		"-i <filename> - A multiple sequence alignment file\n"
		"-o <directory> - An output directory path\n\n"
		"Optional parameters:\n"
		"-p <value> - Minimum correlation p-value (in -log(x)) [DEFAULT: 15]\n"
		"-O <value> - Minimum Occupancy value. It is used to remove fragments of sequences. (0-1) [DEFAULT: 0.8]\n"
		"-I <value> - Maximum Identity value. It is used to remove high identity sequences. (0-1) [DEFAULT: 0.8]\n"
		"-f <value> - Minimum node frequency. It removes nodes (residues) that rarely occurs in the alignment. (0-1) [DEFAULT: 0.05]\n"
		"-F <value> - Maximum node frequency. It removes higly conserved nodes (residues). (0-1) [DEFAULT: 0.95]\n"
		"-m <0 or 1> - Method to Statistically validate the nodes.\n\t0 - Tumminello (Network based validation - DEFAULT)\n\t1 - DRCN (Frequency based validation)\n"
		"-e <0 or 1> - Include marginally conservation properties.\n"
		"-a <0, 1, 2>- Anticorrelations calculation modes:\n0: Do not calculate\n1: Only for Communities\n2: Calculate all\n[DEFAULT: 0]\n"
		"-seq <header>- Sequence correlation mode.\nProgram dumps escores for each residue pair in a sequence.\nE.g. Sequence AEIOU: A1-E2, A1-I3, A1-O4, A1-U5 ... O4-U5 .\nHeader must not contain > for fasta msa.\n "
		"\t0 - Consider only co-variation between amino acids.\n"
		"\t1 - Also include stereochemical and structural amino acids properties. [DEFAULT]\n"
		"\nMore information at http://www.biocomp.icb.ufmg.br/conan")

######GENERAL VARIABLES##########
hmm = False

######PARSE PARAMETERS##########
if len(sys.argv) < 2:
   printHelp()
   sys.exit()

argv = sys.argv[1:]

inputfile = ''
outputdir = ''
minocc = 0.8
maxid = 0.8
minfreq = 0.05
maxfreq = 0.95
method = 0
netthr = 5
marginal = 1
min_pv = 15
anticorr=0
seqmode=False

try:
	opts, args = getopt.getopt(argv,"i:o:p:O:I:f:F:m:e:h")
	#print(opts)
	#print(args)
except getopt.GetoptError:
	print("###########################################################################\n"
		"CONAN - Co-Variation Network Analyzer\n"
		"###########################################################################\n\n"
		"This software analyzes a residue co-variation network with the aim to emphasize local evolutionary constraints and "
		"consequently detect functional and specificity determinant sites. Its only mandatory input is a multiple sequence "
		"alignment in one of the following formats: Pfam, Stockholm or FASTA; and the output consists of several networks and "
		"communities files from several cuts of the network.\n\n"
		"This software was still not pubblished, but further information about this methodology can be accessed at:\n\n"
		"da Fonseca, N. J., Afonso, M. Q. L., de Oliveira, L. C., & Bleicher, L. (2018). A new method bridging graph "
		"theory and residue co-evolutionary networks for specificity determinant positions detection. Bioinformatics.\n\n"
		"The mandatory inputs consits of:\n"
		"-i <filename> - A multiple sequence alignment file\n"
		"-o <directory> - An output directory path\n"
		"Use -h for access the optional parameters.\n"
		"You can also find more information at http://www.biocomp.icb.ufmg.br/conan")
	sys.exit(2)
for opt, arg in opts:
	if opt == '-h':
		printHelp()
		sys.exit()
	elif opt in ("-i"):
		inputfile = arg
	elif opt in ("-o"):
		outputdir = arg
	elif opt in ("-p"):
		min_pv = float(arg)
	elif opt in ("-O"):
		minocc = float(arg)
	elif opt in ("-I"):
		maxid = float(arg)
	elif opt in ("-f"):
		minfreq = float(arg)
	elif opt in ("-F"):
		maxfreq = float(arg)
	elif opt in ("-m"):
		method = int(arg)
	elif opt in ("-e"):
		marginal = int(arg)
    elif opt in ("-a"):
        anticorr = int(arg)
    elif opt in ("-seq"):
        seqmode = str(arg)
if inputfile == "":
	print("The input file is a mandatory parameter.")
	sys.exit()
if outputdir == "":
	print("The output directory is a mandatory parameter.")
	sys.exit()

###########READ ALIGNMENT#############
##ATUALIZAR PARA LER FASTA
msa = {}
nodes = set()
nodes_frequencies = []

###Check MSA type
fr = open(inputfile,'r')

line = fr.readline()
fr.close()

if line[0] == '>':
	msa = readFasta(inputfile)
else:
	msa = readStockholm(inputfile)

N_msa = len(msa)

#########CREATE FOLDER############
if outputdir[len(outputdir)-1] == '/':#Remove /
	outputdir = outputdir[:-1]
if not os.path.exists(outputdir):
	os.makedirs(outputdir)

autoskip = False
if autoskip == False or os.path.exists(outputdir + "/filtered.txt") == False:
	######FILTER MSA BY OCC############
	print("Filtering the Alignment...")
	hmm = isHMM()
	if hmm and (minocc > 0.0):
		msa = hmmFiltering()
	N_msaOCC = len(msa)
	#print(N_msaOCC)

	######FILTER BY MAXID#############
	if maxid > 0.0:
		msa = maxIdCDhit()
	N_msaOCC_ID = len(msa)

	######FILTER RESIDUES###############
	if marginal == 0:
		res_rem,nodes_frequencies = residueFiltering()
    if seqmode:
        print(" Trimming alignment based on seqheader: " + seqmode)
        res_remx = seqbasedFiltering(seqmode)

	#######WRITE FILTERED MSA#########
	fw = open(outputdir + "/filtered.txt",'w')
	for name,sequence in msa.items():
		fw.write(name + "\t" + sequence + "\n")
	fw.close()
else:
	print("Reloading MSA from previous filtered state...")

######GENERATE NETWORK############
print("Calculating the correlations...")
nodes = getNodes()
antinetwork = False

if seqmode:
    #all is reset
    netthr = 0
    minfreq = 0.0
    maxfreq = 1.0
    min_pv = 0
if marginal == 1:
	nodes,res_rem,nodes_frequencies = nodesFiltering(nodes)
if anticorr == 0:
    network = genCorrelationNetwork(method,netthr,nodes_frequencies)
else:
    network, antinetwork = genCorrelationNetworkWithAntiCorr(method,netthr,nodes_frequencies)

if seqmode:
    print("Dumping Res1, Res2, Freq1, Freq2, Freq2|1, Freq1|2, P-Value for: " + seqmode)
    seqPath =  outputdir + "/" + header + "_graph.txt"
    writeSeqmodeNetwork(seqPath, network, antinetwork)
    quit()

#####COMMUNITY DETECTION##########
print("Detecting communities...")
#outputcomm = outputdir + "/comms.txt"

comPath = outputdir + "/communities/"
if not os.path.exists(comPath):
	os.makedirs(comPath)
backPath = outputdir + "/backbones/"
if not os.path.exists(backPath):
	os.makedirs(backPath)

fw = open(outputdir + "/cutoff.txt",'w')

print(" Correlation filtering and parsing...")
print("  1. Minimum p-value filtering...")
#Filter by p-value
network = [(ni,nj,w,jc) for (ni,nj,w,jc) in network if w >= min_pv]
print("  1. Done!")
print("  2. Minimum Jacard Coeficient filtering...")
#Filter by Min Jacard Coeficient Value (Reduce complexity)
network = [(ni,nj,w,jc) for (ni,nj,w,jc) in network if jc >= conan_constants.MIN_JC]
print("  2. Done!")
print("  3. Jacard Coeficient sorting...")
#Sort by Jacard Coeficient
network = sorted(network, key=lambda tup: tup[3])
print("  3. Done!")

if anticorr > 0:
    print(" Anticorrelation filtering and parsing...")
    #Filter by p-value
    print("  1. Minimum p-value filtering...")
    antinetwork = [(ni,nj,w,jc) for (ni,nj,w,jc) in antinetwork if w >= min_pv]
    print("  1. Done!")
    #Filter by Min Jacard Coeficient Value (Reduce complexity)
    # print("  2. Minimum Jacard Coeficient filtering...")
    # antinetwork = [(ni,nj,w,jc) for (ni,nj,w,jc) in antinetwork if jc >= conan_constants.MIN_JC]
    #Sort by Jacard Coeficient
    print("  3. Jacard Coeficient sorting...")
    antinetwork = sorted(antinetwork, key=lambda tup: tup[3])
    print("  3. Done!")
    if anticorr == 1:
        #Filter by residue in network
        print("  3b. Filtering to Communities residues only...")
        uniques_network = set([network[0] for a in network]).union(set([network[1] for a in network]))
        antinetwork = [(ni,nj,w,jc) for (ni,nj,w,jc) in antinetwork if ni in uniques_network or nj in uniques_network]
        print("  3b. Done!")
    #Filter by minimum freq in MSA
    print("  4. Filtering by min freq...")
    antinetwork = [(ni,nj,w,jc) for (ni,nj,w,jc) in antinetwork if getNodeFrequency(ni) >= conan_constants.MIN_CORR_FREQ and getNodeFrequency(nj) >= conan_constants.MIN_CORR_FREQ]
    print("  4. Done!")
    

communities = []
detected_nodes = {}
Gnx = nx.Graph()
extra_edgesJC = {}

print(" Community assignment start...")
###Last State Variables###
lsv_communities = []
lsv_N_residues = 0

#First Iteration. 
if len(network) > 0:
    #A single edge is extracted from a network and added to a community
	ni,nj,pv,jc = network.pop()
    #node to comm index detected_nodes dictionary is initialized for the two Nodes with 0 as start point
	comm = set([ni,nj])
	detected_nodes[ni] = 0
	detected_nodes[nj] = 0
    #first community is pushed to communities and this composes the last state variable communities
	communities.append(comm)
	lsv_communities = communities
    #last variable residues is equal to 2 and Gnx variable receives edge
	Gnx.add_edge(ni,nj,weight=jc,pvalue=pv)
	lsv_N_residues = 2
print( str(len(network)) + " remaining...")
#Post iterations. This lasts until all edges are removed from community
while(len(network) > 0):
    #a new edge is extracted from the network
	ni,nj,pv,jc = network.pop() 
    #and added to a new network variable
	Gnx.add_edge(ni,nj,weight=jc,pvalue=pv)

    #if both edges are present in detected_nodes dictionary keys
	if ni in detected_nodes and nj in detected_nodes:
        #communities of I and J are defined as the result of dictionary (index for each community)
		commI = detected_nodes[ni]
		commJ = detected_nodes[nj]
        #if these communities differ
		if commI != commJ:#Merge communities
            #if I index is bigger than J index, invert their order
			if commI > commJ:
                #I becomes J and J becomes I
				temp = commI
				commI = commJ
				commJ = temp
            #since they already exist, community of smallest index receives all Nodes from community of biggest index
			communities[commI] = communities[commI].union(communities[commJ])
            #second, now redundant, community is deleted
			del communities[commJ]
			for res,comm in detected_nodes.items():
                #and detected_nodes node to comm index dictionary is refactored for previous residues to I
				if comm == commJ:
					detected_nodes[res] = commI
                #and for other residues with other indexa bigger than J as index-1
				elif comm > commJ:
					detected_nodes[res] = comm-1
        #if index I does not differ from index J, simply go on. you are in the same community with the same nodes, weirdly
		else:
			continue #Both nodes were already detected
    #if only the first node is inside a community
	if ni in detected_nodes:#Add Nj to Ni community
		commI = detected_nodes[ni]
		communities[commI].add(nj)
		detected_nodes[nj] = commI
    #if only the second node is inside a community
	elif nj in detected_nodes:#Add Ni to Nj community
		commJ = detected_nodes[nj]
		communities[commJ].add(ni)
		detected_nodes[ni] = commJ
    #if both nodes are not found inside any previously defined community
	else:#New Cluster
		comm = set([ni,nj])
		i = len(communities)
		detected_nodes[ni] = i
		detected_nodes[nj] = i
		communities.append(comm)

	#Validate and Write State
	comms = []
    #for each merged, new or updated(addition) community
	for i,comm in enumerate(communities):
        #If both Residues and Properties are found in a same community, they are merged.
        # Properties with a higher number of gruped aminoacids are preferentialy kept in relation to others
		fixed_comm = filterRedundancy(comm)
		#N_residues+=len(fixed_comm)
		comms.append(fixed_comm)
    
    #the number of unique Nodes in communities with 2+ residues is then calculated
	N_residues = correctedNresidues(comms)

    #Possible Storage and Writing in lsv with new edge extraction happens if number of unique Nodes in communities with 2+ residues is BIGGER (New Nodes added to Comm)
	if N_residues > lsv_N_residues:#Write lsv
		#Gnx_filtered = filterNetwork(Gnx.copy(),lsv_communities)
		Gnx_filtered = Gnx.copy()
		
		writeCommunities(comPath + str(lsv_N_residues),lsv_communities)
		writeBackbone(backPath + str(lsv_N_residues),Gnx_filtered)
		Ncoms = len(lsv_communities)

		jc_avg,minJC = getFullJacardCoefficient(Gnx,lsv_communities,extra_edgesJC)
		if minJC < conan_constants.MIN_JC:
            #if addition of new Node causes decrease in Avg Jaccard Coefficient below the constant threshold (0.5), 
            #it is ignored and a new iteration begin
			break
        
        #if jaccard restriction is not violated, new lsv is written to file
		fw.write(str(lsv_N_residues) + " " + str(Ncoms) + " " + str(jc_avg) + "\n")
		print("N. Residues: " + str(lsv_N_residues) + "\tN. communities: " + str(Ncoms) + " JC avg: " + str(jc_avg) + " Min JC: " + str(minJC))
        print( str(len(network)) + " remaining...")
        
        #and both communities and lsv_residues are saved for the next iteration
		lsv_N_residues = N_residues
		lsv_communities = comms
	else:
        #Storage of communities but not lsv_residues and a new round of edge extraction happens
        # if number of unique Nodes in communities with 2+ residues is UNCHANGED or SMALLER than previous
		lsv_communities = comms
fw.close()