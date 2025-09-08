import csv
import pickle
import os,sys
import random
import importlib
import subprocess
import pandas as pd 
import json, re, time
import urllib.request, urllib.parse, urllib.error
import traceback

from datetime import datetime
from operator import itemgetter

try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET


def print_log(string, end="\n"):
    print(string, end=end)

### GDS interface
def getGDSSamples(date_region=False, path = "gdsSamples.txt"):
    """Will run the predefined query and return a list of GDS ids
    NOTE: this returns ALL GDS samples which are of SRA type i.e.
    ALL CHIP-SEQ, RNA-SEQ, etc.
    """
    ret = []
    species_information = "%20AND%20(homo%20sapiens[Organism]%20OR%20mus%20musculus[Organism])"
    species_information = ""
    #REAL URL
    URL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" \
    "db=gds&term=SRA[Sample%20Type]%20AND%20gse[Entry%20Type]{}&retmax=1000000&usehistory=y".format(species_information)
    if date_region:
        maxTime = date_region.split('-')[1]
        minTime = date_region.split('-')[0]
        print(date_region)
        URL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" \
        "db=gds&term=SRA[Sample%20Type]%20AND%20gse[Entry%20Type]{}&" \
        "mindate={}&maxdate={}&datetype=pdat&retmax=1000000&usehistory=y".format(species_information, minTime, maxTime)
    try:
        print_log("getGDSSample: %s" % URL) # output record
        f = urllib.request.urlopen(URL)
        root = ET.fromstring(f.read())
        f.close()
        #Get the IDList
        tmp = root.findall("IdList/Id")
        ret = [i.text for i in tmp]
        #write to disk
        f = open(path, "w")
        for l in ret:
            f.write("%s\n" % l)
        f.close()
        print("Refresh %s"%path)
    except:
        print("Exception in user code:")
        print('-' * 60)
        traceback.print_exc(file=sys.stdout)
        print('-' * 60)
    return ret


def isXML(doc):
    """TEST if it is a valid geo XML record
    NOTE: first lines are-
    <?xml version="1.0" encoding="UTF-8" standalone="no"?>
    """
    f = doc.split("\n")
    return f[0].strip() == """<?xml version="1.0" encoding="UTF-8" standalone="no"?>"""


def gse_idToAcc(gdsId):
    """Given a GDS id, e.g. 300982523, tries to give a GDS accession, e.g.
    GSM982523

    NOTE: there is an algorithm: acc = "GSM"+gdsId[1:] (strip leading 0s)
    """
    # Cut = dropping of the "3" (which indicates sample) and removal of leading
    # leading 0s
    cut = gdsId[1:].lstrip("0")
    return "GSE%s" % cut

def getSyncLog(infoStr):
    """ouput the record to DoneGsmXml.log file
    """
    os.system('echo "[%s] %s"' % (time.strftime('%H:%M:%S'), infoStr))


def getGeoXML(accession, error_file_name_path,path='geo_gse'):
    """HANDLES GEO XML records--i.e. our GEO XML librarian!
    Given a GEO ACCESSION ID, return the xml record for it
    (making the urllib call)"""
    # path pattern: EXAMPLE-GSE1126513 geo/GSE112/GSE112651
    # path = os.path.join(_ppath, ddir)
    if not os.path.exists(path):
        os.mkdir(path)
    path_trunc = os.path.join(path, accession[0:6])
    if not os.path.exists(path_trunc):
        os.mkdir(path_trunc)
    subdir = os.path.join(path_trunc, accession)
    if not os.path.exists(subdir):
        os.mkdir(subdir)
    path = os.path.join(subdir, "%s.xml" % accession)
    if os.path.exists(path):
        f = open(path)
        docString = f.read()
        f.close() 
    else:
        # print accession
        URL = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s&view=quick&form=xml&targ=self" % accession
        try:
            with urllib.request.urlopen(URL) as response:
        # Read the content
                content = response.read()
                # Decode bytes to string (assuming UTF-8 encoding)
                docString = content.decode('utf-8')
            if not isXML(docString):  # try again
                # signal.alarm(180)
                with urllib.request.urlopen(URL) as response:
        # Read the content
                    content = response.read()
                    # Decode bytes to string (assuming UTF-8 encoding)
                    docString = content.decode('utf-8')
            if isXML(docString):
                # write to file
                f = open(path, "w")
                f.write(docString)
                f.close()
                # getSyncLog(proxy.values()[0]+'\t'+accession + '\n')# output record
            else:
                print(accession)
                print(
                    "ERROR: accession is NOT xml. (The accession may be deleted from GEO repository)")
                f1 = open('gsm_notXML.txt', 'a')
                f1.write(accession + '\n')
                f1.close()
        except:
            try:
                print("Network fluctuations. Wait for 5s.")
                time.sleep(5)
                with urllib.request.urlopen(URL) as response:
            # Read the content
                    content = response.read()
                    # Decode bytes to string (assuming UTF-8 encoding)
                    docString = content.decode('utf-8')
                if not isXML(docString):  # try again
                    # signal.alarm(180)
                    with urllib.request.urlopen(URL) as response:
            # Read the content
                        content = response.read()
                        # Decode bytes to string (assuming UTF-8 encoding)
                        docString = content.decode('utf-8')
                if isXML(docString):
                    # write to file
                    f = open(path, "w")
                    f.write(docString)
                    f.close()
                    # getSyncLog(proxy.values()[0]+'\t'+accession + '\n')# output record
                else:
                    print(accession)
                    print(
                        "ERROR: accession is NOT xml. (The accession may be deleted from GEO repository)")
                    f1 = open('gsm_notXML.txt', 'a')
                    f1.write(accession + '\n')
                    f1.close()
            except Exception as e:
                print("Exception in user code: GEOXML")
                print(f"Error parsing {accession}: {e}")
                print('-' * 60)
                traceback.print_exc(file=sys.stdout)
                print('-' * 60)
                docString = None
                print(f"Error parsing {accession}")
                with open(error_file_name_path, 'a', newline='', encoding='utf-8') as csvfile:
                                fieldnames = ['gse','gsm']
                                writer = csv.DictWriter(csvfile,fieldnames=fieldnames, delimiter='\t')
                                writer.writerow({'gse': accession,'gsm':'no GSE'})
    return docString, path

def getGSMXML(gseid,accession, error_file_name_path, path):
    """HANDLES GSM XML records--i.e. our GSM XML librarian!
    Given a GSM ACCESSION ID, return the xml record for it
    (making the urllib call)"""
    subdir = os.path.join(path, "%s.xml" % accession)
    if os.path.exists(subdir):
        f = open(subdir)
        docString = f.read()
        f.close()
    else:
        # print accession
        URL = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s&view=quick&form=xml&targ=self" % accession
        try:
            with urllib.request.urlopen(URL) as response:
        # Read the content
                content = response.read()
                # Decode bytes to string (assuming UTF-8 encoding)
                docString = content.decode('utf-8')
            if not isXML(docString):  # try again
                # signal.alarm(180)
                with urllib.request.urlopen(URL) as response:
        # Read the content
                    content = response.read()
                    # Decode bytes to string (assuming UTF-8 encoding)
                    docString = content.decode('utf-8')
            if isXML(docString):
                # write to file
                f = open(subdir, "w")
                f.write(docString)
                f.close()
                # getSyncLog(proxy.values()[0]+'\t'+accession + '\n')# output record
            else:
                print(accession)
                print(
                    "ERROR: accession is NOT xml. (The accession may be deleted from GEO repository)")
                f1 = open('gsm_notXML.txt', 'a')
                f1.write(accession + '\n')
                f1.close()
        except:
            try:
                print("Network fluctuations. Wait for 5s.")
                time.sleep(5)
                with urllib.request.urlopen(URL) as response:
            # Read the content
                    content = response.read()
                    # Decode bytes to string (assuming UTF-8 encoding)
                    docString = content.decode('utf-8')
                if not isXML(docString):  # try again
                    # signal.alarm(180)
                    with urllib.request.urlopen(URL) as response:
            # Read the content
                        content = response.read()
                        # Decode bytes to string (assuming UTF-8 encoding)
                        docString = content.decode('utf-8')
                if isXML(docString):
                    # write to file
                    f = open(subdir, "w")
                    f.write(docString)
                    f.close()
                    # getSyncLog(proxy.values()[0]+'\t'+accession + '\n')# output record
                else:
                    print(accession)
                    print(
                        "ERROR: accession is NOT xml. (The accession may be deleted from GEO repository)")
                    f1 = open('gsm_notXML.txt', 'a')
                    f1.write(accession + '\n')
                    f1.close()
            except:
                print("Exception in user code: GSMXML")
                print(f"{accession} is error")
                print('-' * 60)
                traceback.print_exc(file=sys.stdout)
                print('-' * 60)
                docString = None
                with open(error_file_name_path, 'a', newline='', encoding='utf-8') as csvfile:
                    fieldnames = ['gse','gsm']
                    writer = csv.DictWriter(csvfile,fieldnames=fieldnames, delimiter='\t')
                    writer.writerow({'gse': gseid,'gsm':accession})
    return docString, subdir

def read_xmlfile(xml_path):
    ''' 
        Read XML file's content.
    '''
    tree = ET.parse(xml_path)
    root = tree.getroot()
    return root

def get_element_data(xml_root, element_name, attribute_name):
    '''
        Retrieve the content of the characteristic names from the XML file by given element.
    '''
    elements = xml_root.findall(f".//{{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}}{element_name}")
    element_dict = {}
    element_list = []
    for element in elements:
        attribute_value = element.text.strip() if element.text else None
        if attribute_value:
            element_list.append(attribute_value)
    if element_list:
        element_dict[attribute_name] = element_list
    return element_dict


def get_GSE_dict(xml_root):
    '''
        Get the content of GSE by given element and save them in a dictionary. 
    '''
    GSE_dict = {}
    GSE_dict.update(get_element_data(xml_root, 'Summary', 'GSE_Summary'))
    GSE_dict.update(get_element_data(xml_root, 'Overall-Design', 'GSE_Overall_Design'))
    GSE_dict.update(get_element_data(xml_root, 'Title', 'GSE_Title'))
    GSE_dict.update(get_element_data(xml_root, 'Supplementary-Data', 'GSE_URL'))
    GSE_dict.update(get_element_data(xml_root, 'Last-Update-Date', 'GSE_Last_Update_Date'))
    GSE_dict.update(get_element_data(xml_root, 'Submission-Date', 'GSE_Submission_Date'))
    GSE_dict.update(get_element_data(xml_root, 'Release-Date', 'GSE_Release_Date'))
    GSE_dict.update(get_element_data(xml_root, 'Pubmed-ID', 'GSE_Pubmed_ID'))
    GSE_dict.update(get_element_data(xml_root, 'Accession', 'Accession'))
    GSE_dict.update(get_element_data(xml_root, 'Type', 'Type'))
    return GSE_dict

def get_GSM_dict(xml_root):
    '''
        Get the content of GSM by given element and save them in a dictionary. 
    '''
    GSM_dict = {}
    GSM_dict.update(get_element_data(xml_root, 'Title', 'GSM_Title'))
    GSM_dict.update(get_element_data(xml_root, 'Accession', 'GSM_Accession'))
    GSM_dict.update(get_element_data(xml_root, 'Characteristics', 'GSM_Characteristics'))
    GSM_dict.update(get_element_data(xml_root, 'Organism', 'GSM_Organism'))
    GSM_dict.update(get_element_data(xml_root, 'Supplementary-Data', 'GSM_GSM_URL'))
    GSM_dict.update(get_element_data(xml_root, 'Last-Update-Date', 'GSM_Last_Update_Date'))
    GSM_dict.update(get_element_data(xml_root, 'Submission-Date', 'GSM_Submission_Date'))
    GSM_dict.update(get_element_data(xml_root, 'Source', 'GSM_Source'))
    GSM_dict.update(get_element_data(xml_root, 'Release-Date', 'GSM_Release_Date'))
    GSM_dict.update(get_element_data(xml_root, 'Treatment-Protocol', 'GSM_Treatment_Protocol'))
    GSM_dict.update(get_element_data(xml_root, 'Extract-Protocol', 'GSM_Extract_Protocol'))
    GSM_dict.update(get_element_data(xml_root, 'Data-Processing', 'GSM_Data_Processing'))
    return GSM_dict

def match_KeyWord(xmlContent, key):
    '''
    Extract fields containing a specific keyword from the content of an XML file 
    and store these fields along with their matched content in a dictionary.
    '''
    res = {}
    for field in xmlContent.keys():
        if field in list(xmlContent.keys()):
            tmp = list(set(re.findall(r'%s' % key, xmlContent[field][0].replace(
                '-', ' ').replace('_', ' '), re.I)))
            if tmp:
                res[field] = tmp  # a list in dict
    return res


def match_scRNAseq(xmlContent):
    """
    match key words, like single cell, single cell RNA-seq, and sequencing platform,etc
    return a dict, contains geo items and matched keys
    """
    if not xmlContent:
        os.system('echo "No XML content"')
        return {}
    match_res = {}
    # 1. match with single cell words, remove special characters, like '-'
    for key1 in ['scrna sequencing','single cell', 'scrnaseq', 'single cells', 'scrna seq', 'singlecell rnaseq', 'singlecell transcriptome']:
        tmp = match_KeyWord(xmlContent, key1)
        if tmp:
            for i in tmp.keys():
                match_res[i].extend(tmp[i]) if i in match_res.keys(
                ) else match_res.update({i: tmp[i]})
    
    return match_res


def match_scATACseq(xmlContent):
    """
    match scATACseq key words,
    return a dict, contains geo items and matched keys
    """
    if not xmlContent:
        os.system('echo "No XML content"')
        return {}
    match_res = {}
    singleCellATACSeqKeywords = ['scatacseq', 'sc atacseq', 'singlecell atacseq', 'singlecell chromatin accessibility',
                                 'single cell atacseq',  'single cell chromatin accessibility', 'single nucleus atacseq',
                                 'singlecell assay for transposase accessible chromatin', 'sciatacseq', 'dsciatacseq',
                                 'singlenucleus atacseq', '(sci)atacseq', 'scatac', 'single cell profiling of chromatin accessibility',
                                 'scthsseq', 'sndropseq','atac','atac_seq','atacseq', 'atac seq']
    # 1. match with single cell words, remove special characters, like '-'
    for key1 in singleCellATACSeqKeywords:
        tmp = match_KeyWord(xmlContent, key1)
        if tmp:
            for i in tmp.keys():
                match_res[i].extend(tmp[i]) if i in match_res.keys(
                ) else match_res.update({i: tmp[i]})
    # print(match_res)
    return match_res

def match_ST(xmlContent):
    """
    match ST key words, 
    return a dict, contains geo items and matched keys
    """
    if not xmlContent:
        os.system('echo "No XML content"')
        return {}
    match_res = {}
    spatial_keywords = ['visium','stereoseq','stereo-seq','seqscope','spatial transcriptomics',
                        'spatial transcriptome','spatial gene expression','spatial rna sequencing']
    for key1 in spatial_keywords:
        tmp = match_KeyWord(xmlContent, key1)
        if tmp:
            for i in tmp.keys():
                match_res[i].extend(tmp[i]) if i in match_res.keys(
                ) else match_res.update({i: tmp[i]})
    return match_res

def match_other(xmlContent):
    """
    match other key words, 
    return a dict, contains geo items and matched keys
    """
    if not xmlContent:
        os.system('echo "No XML content"')
        return {}
    match_res = {}
    singleCellATACSeqKeywords = ["stereo-seq"," coxmx 6000", "visium hd", "xenium",'stereo seq',]
    # 1. match with single cell words, remove special characters, like '-'
    for key1 in singleCellATACSeqKeywords:
        tmp = match_KeyWord(xmlContent, key1)
        if tmp:
            for i in tmp.keys():
                match_res[i].extend(tmp[i]) if i in match_res.keys(
                ) else match_res.update({i: tmp[i]})
    return match_res

def match_clip(xmlContent):
    """
    match clip key words, 
    return a dict, contains geo items and matched keys
    """
    if not xmlContent:
        os.system('echo "No XML content"')
        return {}
    match_res = {}
    ClipKeywords = ["clip-seq", "clip seq", "hits-clip", "hits clip", 
                    "par-clip", "par clip", 'iclip', 'eclip', 'irclip',
                    'clip', 'crosslinking immunoprecipitation', 'crosslinking and immunoprecipitation',
                    'crosslinking-immunoprecipitation', 'crosslinking-and immunoprecipitation',
                    'chirp-seq', 'chirp seq', 'chromatin isolation by rna purification', 
                    'chart-seq', 'chart seq', 'capture hybridization analysis of rna targets',
                    'rap-seq', 'rap seq', 'rna antisense purification', 'rap-ms', 'rap ms', 
                    'rnacompete', 'selex', 'systematic evolution of ligands by exponential enrichment',
                    'rna-interactome capture', 'rna interactome capture', 
                    'caric', 'rbdmap', 'enhanced rna interactome capture']
    # 1. match with single cell words, remove special characters, like '-'
    for key1 in ClipKeywords:
        tmp = match_KeyWord(xmlContent, key1)
        if tmp:
            for i in tmp.keys():
                match_res[i].extend(tmp[i]) if i in match_res.keys(
                ) else match_res.update({i: tmp[i]})
    return match_res

def get_meta_url(gseid, 
                 gse_output_csv,
                 gse_gsm_output_csv,
                 error_file_name_path,
                 ttype,
                 path = 'geo_gse',
                 only_GSE = False):
        '''
        download GSE XML and get GSE information,
        only_GSE determines whether to download the GSM when the GSE meets the conditions. 
        If only_GSE is Ture , it will only download and obtain the XML and meta information for GSE.
        '''
        _, gse_xml_path = getGeoXML(gseid, error_file_name_path, path)
        if os.path.isfile(gse_xml_path):
            try:
                xml_root = read_xmlfile(gse_xml_path)
                gse_dir_path = os.path.join(path, gseid)
                GSE_meta_dict = get_GSE_dict(xml_root)
                accession_values = GSE_meta_dict.get('Accession', None)
                GSE_URL_value = GSE_meta_dict.get('GSE_URL', None)
                GSE_Title_value = GSE_meta_dict.get('GSE_Title', None)
                GSE_Pubmed_ID_value = GSE_meta_dict.get('GSE_Pubmed_ID', None)
                GSE_Last_Update_Date_value = GSE_meta_dict.get('GSE_Last_Update_Date', None)
                GSE_Release_Date_value = GSE_meta_dict.get('GSE_Release_Date', None)
                ret = {}
                if 'sc-rna-seq' in ttype:
                    # single cell RNA-seq, Library-strategy must be RNA-Seq
                    res = match_scRNAseq(GSE_meta_dict)
                    if res:
                        ret['sc-rna-seq'] = res
                    res = None
                if 'sc-atac-seq' in ttype:
                    # single cell ATAC-seq, Library-strategy must be ATAC-Seq
                    #if 'Genome binding/occupancy profiling by high throughput sequencing' in GSE_meta_dict['Type']:
                    res = match_scATACseq(GSE_meta_dict)
                    if res:
                        ret['sc-atac-seq'] = res
                    res = None
                if 'spatial' in ttype:
                    res = match_ST(GSE_meta_dict)
                    if res:
                        ret['spatial'] = res
                    res = None
                if 'other' in ttype:
                    res = match_other(GSE_meta_dict)
                    if res:
                        ret['other'] = res
                    res = None
                if 'clip' in ttype:
                    res = match_clip(GSE_meta_dict)
                    if res:
                        ret['clip'] = res
                    res = None
                matched_key = ','.join(list(ret.keys()))
                matched_key_value = list(ret.values())
                # Write to GSE_output.csv
                if ret:
                    with open(gse_output_csv, 'a', newline='') as csvfile:
                        csvwriter = csv.writer(csvfile, delimiter='\t')
                        csvwriter.writerow([gseid, GSE_Title_value, GSE_Pubmed_ID_value, 
                                            GSE_Last_Update_Date_value, GSE_Release_Date_value, 
                                            GSE_URL_value, matched_key, matched_key_value])
                ## Decide whether to crawl GSM
                    if not only_GSE:
                        for gsm in accession_values:
                            if gsm.startswith('GSM'):
                                _, gsm_xml_path = getGSMXML(gseid, gsm, error_file_name_path, gse_dir_path)
                                if os.path.isfile(gsm_xml_path):
                                    gsm_xml_root = read_xmlfile(gsm_xml_path)
                                    gsm_URL = get_element_data(gsm_xml_root, 'Supplementary-Data', 'GSM_URL')
                                    GSM_Organism = get_element_data(gsm_xml_root, 'Organism', 'GSM_Organism')
                                    gsm_URL_value = gsm_URL.get("GSM_URL", None)
                                    GSM_Organism_value = GSM_Organism.get('GSM_Organism', None)
                                    # Write to GSE_GSM_output.csv
                                    with open(gse_gsm_output_csv, 'a', newline='') as gsm_csvfile:
                                        csvwriter_gsm = csv.writer(gsm_csvfile, delimiter='\t')
                                        csvwriter_gsm.writerow([gseid, gsm, GSM_Organism_value, gsm_URL_value])
                                else:
                                    print(f"Error parsing {gseid},{gsm}")
                                    with open(error_file_name_path, 'a', newline='', encoding='utf-8') as csvfile:
                                        fieldnames = ['gse','gsm']
                                        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')
                                        writer.writerow({'gse': gseid,'gsm':gsm})
            except:
                print(f"Error parsing {gseid}")
                with open(error_file_name_path, 'a', newline='', encoding='utf-8') as csvfile:
                    fieldnames = ['gse','gsm']
                    writer = csv.DictWriter(csvfile,fieldnames=fieldnames, delimiter='\t')
                    writer.writerow({'gse': gseid,'gsm':'GSE_XML_file error.'})
        else:
            print(f"Error parsing {gseid}")
            with open(error_file_name_path, 'a', newline='', encoding='utf-8') as csvfile:
                fieldnames = ['gse','gsm']
                writer = csv.DictWriter(csvfile,fieldnames=fieldnames, delimiter='\t')
                writer.writerow({'gse': gseid,'gsm':'no GSE'})
        