#!/usr/bin/python

##GenotypingAnalysis.py by Rohan Maddamsetti
##This script does a preliminary analysis of the Illumina VeraCode GoldenGate
##SNP assay that me and Jeff put together.

##The first thing to do is analyze SNPs from the standard clonal mixtures.
##This is written in Enthought's Python 2.6 distribution to use their nice libraries.

from __future__ import print_function
import csv
import re
import pprint
import os
import matplotlib.pyplot as plt
import sys ##this is just for incremental testing.

###############################################################################
#CLASS DEFINITIONS
###############################################################################

##For reference:
##illumina_fields = ['SNP Name', 'Position', 'Well', 'Allele1 - Top', 'Allele2 - Top', 'GC Score', 'Theta', 'R', 'X', 'Y', 'X Raw', 'Y Raw', 'B Allele Freq', 'Log R Ratio']
class SNPdata:
    """ class that store the data for a single SNP (a line from the Illumina data file) """
    def __init__(self,name,position,well,allele1,allele2,gc_score,theta,r,x,y,x_raw,y_raw,b_allele_freq,log_r_ratio):
        self.name = name
        self.position = position
        self.well = well
        self.allele1 = allele1
        self.allele2 = allele2
        self.gc_score = gc_score
        self.theta = theta
        self.r = r
        self.x = x
        self.y = y
        self.x_raw = x_raw
        self.y_raw = y_raw
        self.b_allele_freq = b_allele_freq
        self.log_r_ratio = log_r_ratio

    def printme(self):
        print(self.name,self.position,self.allele1,self.allele2,self.gc_score,self.theta,\
              self.r,self.x,self.y,self.x_raw,self.y_raw,self.b_allele_freq,self.log_r_ratio,sep='\t')

##For reference:
#plate_fields = ['Well', 'Strain', 'Project', 'Pop', 'Gen', 'M or C']
class Welldata:
    """Class that stores the plate description metadata (a line from the Plate data file)"""
    def __init__(self,well,strain,project,pop,gen,m_or_c):
        self.well = well
        self.strain = strain
        self.project = project
        self.pop = pop
        self.gen = gen
        self.m_or_c = m_or_c

    def printme(self):
        print(self.well,self.strain,self.project,self.pop,self.gen,self.m_or_c,sep='\t')

class Sample:
    """This class contains all SNPs relevant to a single sample, and the relevant Welldata."""
    def __init__(self,well_data,snp_list=[]):
        self.well_data = well_data
        self.well_id = well_data.well
        self.snp_list = []
        ##check snpdata to make sure it comes from a single well.
        for candidate in snp_list:
            if candidate.well == self.well_id:
                self.snp_list.append(candidate)
            else:
                raise TypeError("Sample can only contain data from one well")

    def addSNP(self,snp):
        """add a SNP to the Sample's list of SNPs."""
        if snp.well == self.well_id:                
            self.snp_list.append(snp)
        else:
            raise TypeError("Sample can only contain data from one well")

class WellSNPpair:
    """This class contains a paired Welldata and SNPdata."""
    def __init__(self,welldata,snpdata):
        if welldata.well == snpdata.well:
            self.welldata = welldata
            self.snpdata = snpdata
        else:
            raise TypeError("WellSNPpair must be from the same well")

class SNP:
    """This class contains WellSNPpairs for a SNP."""
    def __init__(self,pairs=[]):
        self.pairs = []
        self.name = ''
        self.clone_pairs = []
        self.mixedclone_pairs = []
        self.mixedpop_pairs = []
        if pairs:
            self.name = pairs[0].snpdata.name
        for pair in pairs:
            if self.name == pair.snpdata.name:
                self.pairs.append(pair)
                the_string = pair.welldata.m_or_c
                if re.match('^M$',the_string):
                    self.mixedpop_pairs.append(pair)
                elif re.match('^MC$',the_string):
                    self.mixedclone_pairs.append(pair)
                elif re.match('^C',the_string):
                    self.clone_pairs.append(pair)
            else:
                raise TypeError("Bad SNP was almost made")

    def add_pair(self,pair):
        """add a SNP to the list of SNPS."""
        if self.pairs:
            if self.name == pair.snpdata.name:
                self.pairs.append(pair)
                the_string = pair.welldata.m_or_c
                if re.match('^M$',the_string):
                    self.mixedpop_pairs.append(pair)
                elif re.match('^MC$',the_string):
                    self.mixedclone_pairs.append(pair)
                elif re.match('^C',the_string):
                    self.clone_pairs.append(pair)
            else:
                raise TypeError("Bad WellSNPpair used!")
        else:
            self.name = pair.snpdata.name
            self.pairs.append(pair)
            the_string = pair.welldata.m_or_c
            if re.match('^M$',the_string):
                self.mixedpop_pairs.append(pair)
            elif re.match('^MC$',the_string):
                self.mixedclone_pairs.append(pair)
            elif re.match('^C',the_string):
                self.clone_pairs.append(pair)

#####################################################################
def makeWells(csv_plate_filename):
    """This function makes a dict of welldata indexed by the well name, i.e. 'A4'"""
    """Input:the filename of a csv file. Output:A dict of welldata."""
    wellReader = csv.reader(open(csv_plate_filename,"rU"),dialect='excel')
##skip first row, and take 96 data points for the plate.
    line_counter = 0
    ze_wells = {}
    for row in wellReader:
        rowdata = row[:6]
        line_counter += 1
        if line_counter == 1:
            continue
    ##This next condition ignores the unused samples in the plate.
        if line_counter <= 97:
            my_well = Welldata(*rowdata)
            well_name = my_well.well
            ze_wells[well_name] = my_well
        else:
            break
    return ze_wells

def wellSNPpairs_sort(SNPpairs_list,is_mixedclones=False):
    """This sorts a list of wellSNPpairs by generation (timepoint)."""
    if is_mixedclones:
    ##grab first match before %.
        deco = []
        for i, pair in enumerate(SNPpairs_list):
            my_regex = re.compile('^(\d+)%')
            the_match = my_regex.match(pair.welldata.strain)
            mix_ratio = the_match.group(1)
            deco.append((int(mix_ratio),i,pair))
    else:
    ##Schwartzian Transform!
        deco = [ (int(pair.welldata.gen),i,pair) for i, pair in enumerate(SNPpairs_list)]
    deco.sort()
    sorted_pairs = [pair for _, _, pair in deco]
    return sorted_pairs

def SNPgraph(snpname,sorted_pairs,type,dir=""):
    """This function graphs a SNP over time."""
##Skip SNPs whose name doesn't end in '-snp'.
    goodsnp = re.compile('-snp$')
    goodsnp_match = goodsnp.search(snpname)
    if not goodsnp_match:
        print(snpname + " was not matched.")
    else:
        if plt.figure():
            plt.close('all')
        xvals = []
        yvals = []
        mycwd = os.getcwd()
        if dir:
            if not os.path.exists(dir):
                os.makedirs(dir)
                os.chdir(dir)
            else:
                os.chdir(dir)
        if type == 'Mixed Clones':
            my_regex = re.compile('^(\d+)%')
            for pair in sorted_pairs:
                the_match = my_regex.match(pair.welldata.strain)
                mix_ratio = the_match.group(1)
                xvals.append(float(mix_ratio))
                yvals.append(float(pair.snpdata.b_allele_freq))
            myGraph = plt.plot(xvals,yvals,color='red')
            #myGraph = plt.bar(xvals,yvals,width=5,bottom=0,color='red')
            my_xlabel = plt.xlabel('% REL 607 of mixture')
            my_axis = plt.axis([0,100,0,1.5])
            my_ylabel = plt.ylabel('B Allele Frequency')
            my_title = "Graph for " + snpname
            the_title = plt.title(my_title)
        ##save graphs to file.
            filename = "SNP_" + snpname +".jpeg"
            plt.savefig(filename)
            print(filename, "has been saved.")
            #print "xvals are: ",xvals," yvals are: ",yvals
            os.chdir(mycwd)
        else:
            xvals = range(1,len(sorted_pairs)+1)
            for pair in sorted_pairs:
                yvals.append(float(pair.snpdata.b_allele_freq))
#        print "xvals are: ",xvals," yvals are: ",yvals
            #theGraph = plt.bar(xvals,yvals,width=1,bottom=0,color='red')
            theGraph = plt.plot(xvals,yvals,color='red')
            my_xlabel = plt.xlabel('Ranked by Generations')
            my_axis = plt.axis([1,len(sorted_pairs)+1,0,1.5])
            my_ylabel = plt.ylabel('B Allele Frequency')
            my_title = "Graph for " + snpname + ' ' + type 
            the_title = plt.title(my_title)
        ##save graphs to file.
            filename = "SNP_" + snpname +".jpeg"
            plt.savefig(filename)
            plt.close('all')
            print(filename, "has been saved.")

            os.chdir(mycwd)
        
def SampleGraph(sample,my_type,dir=""):
    """This function graphs a set of snps from a sample."""
    if plt.figure():
        plt.close('all')
    xvals = range(1,97)
    snp_names = []
    yvals = []
    mycwd = os.getcwd()
    dirstatus = 0
    if dir:
        dirstatus = 1
        if not os.path.exists(dir):
            os.makedirs(dir)
            os.chdir(dir)
        else:
            os.chdir(dir)
    my_generation = sample.well_data.gen
    
    ##sort SNPs alphanumerically.
###    sorted_snp_list = SNPlist_sort(sample.snp_list)
    for snp in sample.snp_list:
        allele_freq = snp.b_allele_freq
        my_name = snp.name
        snp_names.append(my_name)
        yvals.append(float(allele_freq))
    myGraph = plt.bar(xvals,yvals,width=1,bottom=0)
    my_xlabel = plt.xlabel('SNPS')
    my_axis = plt.axis([0,96,0,1])
    my_title = "Histogram for " + my_type + " at generation " + my_generation
    filename = my_type + "_" + my_generation + ".jpeg"
    if my_type == 'Mixed Clones':
        ratio = sample.well_data.strain
        my_title = "Histogram " "at " + ratio
        p = re.compile('/')
        fixed_ratio = p.sub('-',ratio)
        filename = my_type + "_" + fixed_ratio + ".jpeg"
    graph_title = plt.title(my_title)
    plt.savefig(filename)
    print(filename, "has been saved.")
    plt.close('all')
    if dirstatus:
        os.chdir(mycwd)

def sortSamples(sample_list,is_mixedclones=False):
    """takes a list of samples and sorts them by timepoint, or mixture."""
    if is_mixedclones:
        deco = []
        for i, sample in enumerate(sample_list):
            my_regex = re.compile('^(\d+)%')
            the_match = my_regex.match(sample.well_data.strain)
            mix_ratio = the_match.group(1)
            deco.append((int(mix_ratio),i,sample))
    else:
        ##Schwartzian Transform!
        deco = [ (int(sample.well_data.gen),i,sample) for i,sample in enumerate(sample_list)]
    deco.sort()
    sorted_samples = [sample for _, _, sample in deco]
    return sorted_samples

def printSamples(sample_list,type):
    """Takes a list of clones, mixed pops, etc., and prints to stdout."""
    print(90*'*')
    print(type, 'samples')
    print(90*'*')
    for sample in sample_list:
        print('Well','Strain','Project','Pop','Gen','M or C',sep='\t')
        sample.well_data.printme()
        print('SNP Name','Position','Sample ID','Allele1 - Top','Allele2 - Top',\
              'GC Score', 'Theta', 'R', 'X', 'Y', 'X Raw', 'Y Raw', 'B Allele Freq',\
              'Log R Ratio', sep='\t')
        for snp in sample.snp_list:
            snp.printme()
        print()
        print(90*'#')
        
wells = makeWells("plate1.csv")

##plate is a dict of Sample objects, keyed by well name. No SNP data yet.
plate = dict([(k,Sample(v)) for k,v in wells.iteritems()])

##snp_plate is a dict of SNP objects, keyed by SNP name.
##It is currently empty.
snp_plate = {}

##For formatting the sample data
row_dict = {1:'A',2:'B',3:'C',4:'D',5:'E',6:'F',7:'G',8:'H'}

plate_data = open("Plate 1_FinalReport.txt")
line_counter  = 0
for line in plate_data:
    line_counter += 1
    if line_counter < 11:
        continue
    line.strip()
    ##This ugly code is for formatting the well name from the data file.
    line_data = line.split('\t')
    sampleid = line_data[2]
    sampleid = sampleid[15:]
    row,column = sampleid.split('_')
    row = row[1:]
    column = column[1:]
    row = int(row)
    column = int(column)
    my_well = str(row_dict[row]) + str(column)
    line_data[2] = my_well
    ##Make the new SNPdata object.
    new_snp = SNPdata(*line_data)
##add new_snp to the correct Sample object in plate.
    plate[new_snp.well].addSNP(new_snp)
##make WellSNPpair object.
    my_well_obj = wells[new_snp.well]
    new_pair = WellSNPpair(my_well_obj,new_snp)
##check snp_plate.
    snp_name = new_snp.name
    if snp_name in snp_plate: ##add the pair.
        snp_plate[snp_name].add_pair(new_pair)
    else: ##make a new SNP entry.
        ##because the argument takes a list of pairs.
        snp_plate[snp_name] = SNP([new_pair])

##Now, plate contains all the Sample objects, which each contain data for the SNPs for
##that well; and snp_plate contains all the SNP objects, containing sample-wide data for a
##single SNP.

##This data should now be split into Clones, Mixed Clones, and MixedPop data.

clone_samples = []
mixedclone_samples = []
mixedpop_samples = []
for k,sample in plate.iteritems():
    the_string = sample.well_data.m_or_c
    if re.match('^M$',the_string):
        mixedpop_samples.append(sample)
    elif re.match('^MC$',the_string):
        mixedclone_samples.append(sample)
    elif re.match('^C',the_string):
        clone_samples.append(sample)

sorted_clone_samples = sortSamples(clone_samples)
sorted_mixedclone_samples = sortSamples(mixedclone_samples,True)
sorted_mixedpop_samples = sortSamples(mixedpop_samples)

printSamples(sorted_clone_samples,'CLONES')
printSamples(sorted_mixedclone_samples,'MIXED CLONES')
printSamples(sorted_mixedpop_samples,'MIXED POPULATIONS')

##To check the quality of the data, for each snp I should separate the data into
##clone_pairs, mixedclone_pairs, and mixedpop_pairs, sort each group by timepoint,
##and graph allele frequency per SNP using matplotlib.
##good data will not have random peaks in the middle.
    
##take the clone_pairs from the SNP plate.
#for k,snp in snp_plate.iteritems():
#    sorted_clones_pairs = wellSNPpairs_sort(snp.clone_pairs)
#    SNPgraph(k,sorted_clones_pairs,'Clones','IlluminaGenotyping/CloneSNPs')
#    sorted_mixed_pop_pairs = wellSNPpairs_sort(snp.mixedpop_pairs)
#    SNPgraph(k,sorted_mixed_pop_pairs,'Mixed Pop','IlluminaGenotyping/MixedPopSNPs')
#    sorted_mixed_clone_pairs = wellSNPpairs_sort(snp.mixedclone_pairs,True)
#    SNPgraph(k,sorted_mixed_clone_pairs,'Mixed Clones','IlluminaGenotyping/MixedCloneSNPs')
    

##To examine the data in a nice way, I could graph allele frequency for each SNP in each
##sample, with SNP order on the x-axis and SNP frequency on the y-axis, 
##and each separate graph would be a different snapshot of allele frequencies over time.

#for my_sample in sorted_clone_samples:
#    SampleGraph(my_sample,'Clone','IlluminaGenotyping/CloneSamples')
#for my_sample in sorted_mixedclone_samples:
#    SampleGraph(my_sample,'Mixed Clones','IlluminaGenotyping/MixedCloneSamples')
#for my_sample in sorted_mixedpop_samples:
#    SampleGraph(my_sample,'Mixed Pop','IlluminaGenotyping/MixedPopSamples')

