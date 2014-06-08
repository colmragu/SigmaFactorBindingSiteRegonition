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

def sortbyNcdna(sortbyNCDNA, gb_file):
  scan_pos = 0
  f = open(gb_file,"r")
  ncdna_last = "S"
  count =0
  try:
    gb_record = SeqIO.read(f, "genbank")
  except:
    print ("Skipping %s , cannot read file"%gb_file)
    return (sortbyNCDNA) 
  f.close()
  for feature in gb_record.features[1:]: #Select all but the first entry in gbk file
    feature = pad_gb_features(feature)
    ncdna = str(gb_record.seq[scan_pos:feature.location.start.position+1]) #ncdna = non coding dna
    if ncdna == "":
      sortbyNCDNA = addfeaturetodict(ncdna_last,sortbyNCDNA, feature)
    else:
      sortbyNCDNA = addfeaturetodict(ncdna,sortbyNCDNA, feature)
      ncdna_last= ncdna
    scan_pos = feature.location.end.position
  return(sortbyNCDNA)

def addfeaturetodict(ncdna,sortbyNCDNA, feature):
#  pdb.set_trace()
  if not "CDS" in feature.type:
    return(sortbyNCDNA)

  if feature.qualifiers.has_key("product"):
    feature.qualifiers["note"]  = str(feature.qualifiers["note"]) +  feature.qualifiers["product"][0]

  if not sortbyNCDNA.has_key(ncdna):
    sortbyNCDNA[ncdna]={}
  if not sortbyNCDNA[ncdna].has_key(feature.qualifiers["db_xref"][0]):
    sortbyNCDNA[ncdna][feature.qualifiers["db_xref"][0]] ={}
  sortbyNCDNA[ncdna][feature.qualifiers["db_xref"][0]] =formatestring(feature.qualifiers["note"]  ,int (feature.location.start) ,int (feature.location.end)) #add more information
  return(sortbyNCDNA) 

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

def dictbyCommonsubstring(ncdna_sort,commonsubstrings ,outputfile):
  dictbycommonsubstring = dict((substring,{}) for substring in commonsubstrings)
  for substring in commonsubstrings:
    for element in ncdna_sort.keys():
      if substring in element:
        dictbycommonsubstring[substring].update(ncdna_sort[element])
  return(dictbycommonsubstring)    

def printbyCommonsubstring(ncdna_sort, outputfile):
  f = open(outputfile, "w")
  commonsubstrings = findcommonsubstrings(ncdna_sort.keys(), 5, 5)
  dictbycommonsubstring = dictbyCommonsubstring(ncdna_sort,commonsubstrings ,outputfile)  
  for ncdna in sorted(dictbycommonsubstring.items(), key= lambda x:len(x[1]), reverse=True):
    print(ncdna[0], file = f)
    print_element(dictbycommonsubstring,ncdna[0],f)

def print_sorted(ncdna_sort, outputfile):
  f = open(outputfile, "w")
  commonsubstrings = findcommonsubstrings(ncdna_sort.keys(), 5, 5)
  print (commonsubstrings,file = f)
  for ncdna in reversed(sorted(ncdna_sort.keys(), key= lambda x: sortbycommonsubstring(x, commonsubstrings))):
    print (ncdna, file=f)
    print_element(ncdna_sort,ncdna,f)
  f.close()

def findcommonsubstrings(ncdna_keys, minlen, minoccurances):
  commonsubstrings = {}
  mostcommonsubstrings=[]
  ncdna_key= sorted(ncdna_keys,reverse=True)
  for x, key in enumerate(ncdna_keys):
    print (len(ncdna_keys))
    print (x)
    for key2 in ncdna_keys[x:-1]:
      if key == key2: 
        continue
      seq_matcher = difflib.SequenceMatcher(None, key, key2)
      all_seq = seq_matcher.get_matching_blocks()
#      pdb.set_trace()
      for test_seq in all_seq:
        if test_seq[2] > minlen:
          seq = key[test_seq[0]:test_seq[0] + test_seq[2]]
          if seq in commonsubstrings.keys():
            commonsubstrings[seq]+=1
          else:
            commonsubstrings[seq]=1
  for k in reversed(sorted(commonsubstrings.keys(), key=commonsubstrings.__getitem__)):
    if commonsubstrings[k] >= minoccurances: 
      mostcommonsubstrings.append(k)
    else:
      break      
  return (mostcommonsubstrings)

def sortbycommonsubstring(ncdna, commonsubstrings):
  contains_substring=[0]*len(commonsubstrings)
  for i, substring in enumerate(commonsubstrings):
    if substring in ncdna:
      contains_substring[i]=1
  return (tuple(contains_substring))

def print_sorted_contains_note(ncdna_sort, note , outputfile):
# ornts all genes for sorted under a ncdna key if any gene contains a given note
  f = open(outputfile, "w")
  matches=[]
  for ncdna in sorted(ncdna_sort.keys(), key=len):
    key = search_string_in_notes(note, ncdna_sort, ncdna)
    if key != 0:
      matches.append(ncdna)
    else:
      continue
  for match in sorted(matches, key=len):
    print (match, file=f)
    print_element(ncdna_sort,match, f)
  f.close()

def search_string_in_notes(note, ncdna_sort, ncdna):
  for organizm in sorted(ncdna_sort[ncdna]):
    for gene in ncdna_sort[ncdna][organizm].keys():
      if note.lower() in str(ncdna_sort[ncdna][organizm][gene][0]).lower():
        return 1
  return 0

def print_element(ncdna_sort,ncdna, f):
  for gene,info in sorted(ncdna_sort[ncdna].items(), key=lambda x:x[1][1]):
      print (info, file=f)
  print ("", file=f)

def sortall(gb_files): 
  ncdna_sort = {}
  for gb_file in gb_files:
    print(gb_file)
    ncdna_sort = sortbyNcdna(ncdna_sort, gb_file)
  return(ncdna_sort)

def parse_commandline():
  help_message = "Incorrect command you need to supply a gbk file \n python sortbyNcdna.py -sort data/NC_009930.gbk outputfile \n python sortbyNcdna.py -sort data/Azospirillum_B510_uid46085/NC_013860.gbk Azospirillum_B510_uid46085.txt\n python sortbyNcdna.py -search protein data/Azospirillum_B510_uid46085/NC_013860.gbk Azospirillum_B510_uid46085.txt" 
  if len(sys.argv) > 1:
    if sys.argv[1] == "-search":
      gb_files = sys.argv[3:-1]
      outputfile = sys.argv[-1]
      search_string = sys.argv[2]
      ncdna_sorted = sortall(gb_files)
      print_sorted_contains_note(ncdna_sorted,search_string, outputfile)  
    elif sys.argv[1] == "-sort":
      gb_files = sys.argv[2:-1]
      outputfile = sys.argv[-1]
      ncdna_sorted = sortall(gb_files)
      printbyCommonsubstring(ncdna_sorted,outputfile)
#      print_sorted(ncdna_sorted,outputfile)
    else: 
      print(help_message)
      return(0)
  else:
    print(help_message)
    return(0)
parse_commandline()
