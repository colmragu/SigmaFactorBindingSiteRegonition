#Copyright 2014 Colm Herbert 


#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import print_function
from Bio import SeqIO
import difflib
import pdb
import sys
import re
from collections import Counter
#Globals
global minlen

def sortbyNcdna(GeneId2Ncdna, GeneId2Note, gb_file):
  scan_pos = 0
  f = open(gb_file,"r")
  ncdna_last = "S"
  count =0
  try:
    gb_record = SeqIO.read(f, "genbank")
  except:
    print ("Skipping %s , cannot read file"%gb_file)
    return (GeneId2Ncdna, GeneId2Note) 
  f.close()
  for feature in gb_record.features[1:]: #Select all but the first entry in gbk file
    feature = pad_gb_features(feature)
    ncdna = str(gb_record.seq[scan_pos:feature.location.start.position+1]) #ncdna = non coding dna
    if ncdna == "":
      GeneId2Ncdna, GeneId2Note = addfeaturetodict(ncdna_last, GeneId2Ncdna, GeneId2Note, feature)
    else:
      GeneId2Ncdna, GeneId2Note = addfeaturetodict(ncdna, GeneId2Ncdna, GeneId2Note, feature)
      ncdna_last= ncdna
    scan_pos = feature.location.end.position
  print (len(GeneId2Ncdna))
  return(GeneId2Ncdna, GeneId2Note)

def addfeaturetodict(ncdna,GeneId2Ncdna, GeneId2Note, feature):
#  
  GeneId = feature.qualifiers["db_xref"][0]
  if not "CDS" in feature.type:
    return(GeneId2Ncdna, GeneId2Note)

  if feature.qualifiers.has_key("product"):
    feature.qualifiers["note"]  = str(feature.qualifiers["note"]) +  feature.qualifiers["product"][0]

  if not GeneId2Ncdna.has_key(GeneId):
    GeneId2Ncdna[GeneId]=[ncdna]
    GeneId2Note[GeneId]=formatestring(feature.qualifiers["note"]  ,int (feature.location.start) ,int (feature.location.end))
  else:
    GeneId2Ncdna[GeneId].append(ncdna)

  return(GeneId2Ncdna, GeneId2Note) 

def formatestring(note, start, end ):
  formatedString = note ,start ,end
  return(formatedString)

def pad_gb_features(feature):
  ### Make sure all gb objects conatain a dic entry for note, db_xref
  if not feature.qualifiers.has_key("note"):
    feature.qualifiers["note"] = [""]
  if not feature.qualifiers.has_key("db_xref"):
    feature.qualifiers["db_xref"]=["nodb_xref"]
  return feature

def dictbyCommonsubstring(GeneId2Ncdna, GeneId2Note,commonsubstrings ,outputfile):
  genidByCommonsubstring = {}
  for gene in GeneId2Ncdna.keys():
    for element in GeneId2Ncdna[gene]:
      substring_match = [substring1 for substring1 in commonsubstrings if re.compile(substring1).search(element) is not None]    
      for substring in substring_match:
        if gene in genidByCommonsubstring.keys():
          genidByCommonsubstring[gene].append(substring)
        else:
          genidByCommonsubstring[gene]=[]
          genidByCommonsubstring[gene].append(substring)
  print (len(genidByCommonsubstring.keys()))
  return(genidByCommonsubstring)

def dictby2Commonsubstring(GeneId2Ncdna, GeneId2Note,commonsubstrings ,outputfile):
  GeneId2commonsubstring={}
  dictbycommonsubstring={}
  inv_ncdna_sort={} 
  percent_total=len(GeneId2Ncdna.keys())
  
  print ("\nBuilding dict")
  for y, gene in enumerate(GeneId2Ncdna.keys()):   
    print ("\r".join([str(y*100/percent_total) ,"%"]),end="")
#    print ("\r" + str(float(y*100/len(GeneId2Ncdna.keys()))) + "%",end='')
    substring_pairs=[substring1 for substring1 in GeneId2Ncdna[gene]]
    GeneId2commonsubstring.update({substring_pair:{} for substring_pair in substring_pairs if substring_pair not in GeneId2commonsubstring.keys()})
    for sharedSubstrings in substring_pairs:
        GeneId2commonsubstring[sharedSubstrings].update({gene:GeneId2Note[gene]})
  return(GeneId2commonsubstring)

def printbyCommonsubstring(GeneId2Ncdna, GeneId2Note,outputfile):
  f = open(outputfile, "w")
  ncdna= [x for y in GeneId2Ncdna.values() for x in y]
  commonsubstrings = findcommonsubstrings(ncdna, 5, 5)
  genidByCommonsubstring = dictbyCommonsubstring(GeneId2Ncdna, GeneId2Note,commonsubstrings ,outputfile)   
  dictbycommonsubstring = dictby2Commonsubstring(genidByCommonsubstring, GeneId2Note,commonsubstrings ,outputfile)
  for ncdna in sorted(dictbycommonsubstring.items(), key= lambda x:len(x[1]), reverse=True):
    print(ncdna[0], file = f)
    print_element(dictbycommonsubstring,ncdna[0],f)

def print_sorted(ncdna_sort, outputfile):
  f = open(outputfile, "w")
  commonsubstrings = findcommonsubstrings(ncdna_sort.keys(), 5, 10)
  print (commonsubstrings,file = f)
  for ncdna in reversed(sorted(ncdna_sort.keys(), key= lambda x: sortbycommonsubstring(x, commonsubstrings))):
    print (ncdna, file=f)
    print_element(ncdna_sort,ncdna,f)
  f.close()

def findcommonsubstrings(ncdna_keys, minlen2, minoccurances):
  global minlen 
  minlen = minlen2 
  mostcommonsubstrings=[] 
  ncdna_key= sorted(ncdna_keys,reverse=True)
  ncdna_pairs = [(ncdna1, ncdna2) 
                for x, ncdna1 in enumerate(ncdna_key) 
                for ncdna2 in ncdna_key[x+1:-1] 
                if len(ncdna1) > 10 and len(ncdna2)>10]
  print ("Finding common substrings")
  commonsubstrings = [element 
                      for ncdna_pair in ncdna_pairs
                      for element in ncdna_commonsubstrings(ncdna_pair)]
  commonsubstrings = Counter(commonsubstrings)
  for k in reversed(sorted(commonsubstrings.keys(), key=commonsubstrings.__getitem__)):
    if commonsubstrings[k] >= minoccurances:
      mostcommonsubstrings.append(k)
    else:
      break
  mostcommonsubstrings = sorted(mostcommonsubstrings,key=len, reverse=True)
  return (mostcommonsubstrings)

def ncdna_commonsubstrings(key):
  matches =  [x for x in difflib.SequenceMatcher(None, key[0], key[1]).get_matching_blocks() if x[2]>minlen]

  all_seq_paired = [sigmapinrowjoin(key[0][test_seq[0]:test_seq[0] + test_seq[2]] ,abs(test_seq[0]+test_seq[2]-test_seq2[0]) , key[0][test_seq2[0]:test_seq2[0] + test_seq2[2]]) 
                    for x, test_seq in enumerate(matches) 
                    for test_seq2 in matches[x+1:-1]]
  return(all_seq_paired)

def sigmapinrowjoin(sigma,  distance, pinrow):
#  return (".{%s}".join([sigma, pinrow]))%(distance)
  return (".+".join([sigma, pinrow]))
def sortbycommonsubstring(ncdna, commonsubstrings):
  contains_substring=[0]*len(commonsubstrings)
  for i, substring in enumerate(commonsubstrings):
    if substring in ncdna:
      contains_substring[i]=1
  return (tuple(contains_substring))

def print_element(ncdna_sort,ncdna, f):
  for gene,info in sorted(ncdna_sort[ncdna].items(), key=lambda x:x[1][1]):
      print (info, file=f)
  print ("", file=f)

def sortall(gb_files): 
  GeneId2Ncdna={}
  GeneId2Note={}
  for gb_file in gb_files:
    print(gb_file)
    GeneId2Ncdna, GeneId2Note = sortbyNcdna(GeneId2Ncdna, GeneId2Note, gb_file)
  return(GeneId2Ncdna, GeneId2Note)

def parse_commandline():
  help_message = "Incorrect command you need to supply a gbk file \n python sortbyNcdna.py -sort data/NC_009930.gbk outputfile \n python sortbyNcdna.py -sort data/Azospirillum_B510_uid46085/NC_013860.gbk Azospirillum_B510_uid46085.txt\n python sortbyNcdna.py -search protein data/Azospirillum_B510_uid46085/NC_013860.gbk Azospirillum_B510_uid46085.txt" 
  if len(sys.argv) > 1:
    if sys.argv[1] == "-sort":
      gb_files = sys.argv[2:-1]
      outputfile = sys.argv[-1]
      GeneId2Ncdna, GeneId2Note = sortall(gb_files)
      printbyCommonsubstring(GeneId2Ncdna, GeneId2Note,outputfile)
#      print_sorted(ncdna_sorted,outputfile)
    else: 
      print(help_message)
      return(0)
  else:
    print(help_message)
    return(0)



parse_commandline()






