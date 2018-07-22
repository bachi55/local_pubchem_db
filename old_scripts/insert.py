import re
import sys
import os
import binascii
import commands
import pybel
import psycopg2
import psycopg2.pool

def insert_sdf(con):
    cur = con.cursor()
    a = 1
    b = 25000
    data_dir = "/share/kepaco/metabolites/pubchem_compound/"
    for i in range(2930): #pubchem has 2120 (2930) compound files 
        start = str(a).rjust(9,"0")
        end = str(b).rjust(9,"0")
        sdf = data_dir+"Compound_" + start + "_" + end + ".sdf"
        print "processing ",sdf,"now..."
        try:
            c_sdf = open(sdf)
        except IOError:
            print "When processing",sdf, "cause IOError"
            a = a + 25000
            b = b + 25000
            continue
        parse_sdf(c_sdf,cur)
        con.commit() # commit 
        c_sdf.close()
        #move to parse next file
        a = a + 25000
        b = b + 25000

def parse_sdf(sdf,cur):
    data = sdf.read()
    start = 0
    while True:
        end_pos = data.find("$$$$",start)
#        print end_pos
        if end_pos == -1:
            break
        sdf_s = data[start:end_pos-1]
        sdf_s = sdf_s.replace("'","")
        cid_s = re.findall("<PUBCHEM_COMPOUND_CID>\n([0-9]+)",sdf_s)[0]
        start = end_pos+4
#        insert_sql = "INSERT INTO pubchem_sdfs VALUES (%s, \'%s\');" % (cid_s, sdf_s)
#        print cid_s
#        cur.execute(insert_sql)


def insert_info(con,threadname):
    
    """ parse sdf and insert mass, molecular formula, inchi, fingerprint to a table"""
    cur = con.cursor()
    for cid in xrange(53150001,73250001): # now has 73250000 compounds
#        print 'thread',thread_name,':',cid
        sql1 = "SELECT sdf FROM pubchem_sdfs WHERE cid = %d;" % cid
        cur.execute(sql1)
        rows = cur.fetchall()
        if len(rows) == 0: # don't have this id in db
            continue
        sdf = rows[0][0]

        #######################
#        testsql = "SELECT cid FROM pubchem_info WHERE cid = %d;" % cid
#        cur.execute(testsql)
#        rows = cur.fetchall()
#        if len(rows) != 0: # already in the database
#            continue        
        #######################

        values = parse_sdf_into(threadname,cid,sdf)
        sql2 = "INSERT INTO pubchem_info VALUES " + str(values) + ";"
        cur.execute(sql2)
        if cid % 10000 == 0:
            con.commit()
    con.commit()

def parse_sdf_into(threadname,cid,sdf):
    babel_path = "/home/shenh1/Software/Babel/bin"

    MF = "NULL"
    mass = 0
    inchi = "NULL"
    smiles = "NULL"
    subkeys = "NULL"
    babelfp = "NULL"
    pubchemfp= "NULL"
    
    lines = sdf.split("\n")
    n_line = len(lines)
    i = 0
    while i < n_line:
        if lines[i].find("<PUBCHEM_EXACT_MASS>") != -1:
            mass = float(lines[i+1].strip())
            i = i+1
        if lines[i].find("PUBCHEM_CACTVS_SUBSKEYS") != -1:
            subkeys = lines[i+1].strip()
            i = i+1
        if lines[i].find("<PUBCHEM_IUPAC_INCHI>") != -1:
            inchi = lines[i+1].strip()
            i = i+1
        if lines[i].find("PUBCHEM_MOLECULAR_FORMULA") != -1:
            MF = lines[i+1].strip()
            i = i+1
        if lines[i].find("PUBCHEM_OPENEYE_CAN_SMILES") != -1:
            smiles = lines[i+1].strip()
            i = i+1
        i = i+1

#    print MF,mass,inchi,smiles,subkeys,babel_fp,pubchem_fp
    #pubchemfp = get_pubchem_fp(subkeys)
    #babelfp = get_528_babel_fp_fast(threadname,inchi,cid)
    return (cid,MF,mass,inchi,smiles,babelfp,pubchemfp)

def get_pubchem_fp(subkey):
    b = ""
    fp = subkey
    t = binascii.a2b_base64(fp) # t is ASCII charaters.     
    for s in t:
        b = b + bin(ord(s))[2:].rjust(8,'0')
    if len(b) != 920:
        print "something wrong with the conversion from HEX64 to binary data!"
        return "NULL"
    else:
        return "".join(b[32:913])


def get_528_babel_fp_fast(threadname, inchi, cid):
    try:
        mol = pybel.readstring("inchi",inchi)
    except IOError:
        print cid,"convert from inchi failed"
        return "NULL"
    fp3 = mol.calcfp(fptype='FP3') # -ofpt -xfFP3 >, no "-h"
    fp4 = mol.calcfp(fptype='FP4')
    fpm = mol.calcfp(fptype='MACCS')
    fp3b = ['0']*55
    fp4b = ['0']*307
    fpmb = ['0']*166
    for i in fp3.bits:
        fp3b[i-1] = '1'
    for i in fp4.bits:
        fp4b[i-1] = '1'
    for i in fpm.bits:
        fpmb[i-1] = '1'
    return "".join(fp3b)+"".join(fp4b)+"".join(fpmb)

def get_528_babel_fp(threadname,babel_path, inchi,cid):

    mid = "temp"+str(threadname)
    infn = mid+".inchi"; ofile = mid+".fpt"
    f = open(infn,"w")
    f.write(inchi)
    f.close()
    commands.getoutput("%s %s -ofpt -xfFP3 -h > %s" % (
                       babel_path+"/babel", infn, ofile))
    commands.getoutput("%s %s -ofpt -xfFP4 -h >> %s" % (
                       babel_path+"/babel", infn, ofile))
    commands.getoutput("%s %s -ofpt -xfMACCS -h >> %s" % (
                       babel_path+"/babel", infn, ofile))
    #fp3 fp4 maccs has 55, 307, 166 functional groups                      
    featcounts = [55, 307, 166]
    f = open(ofile)
    lines = f.readlines()
    lines = map(str.strip, lines)
    f.close()
    if len(lines) == 0:
        return "NULL"
    # the file contains three blocks on fingerprints:                      
    # FP3, FP4 and MACCS                                                   
    checksums = map(int, re.findall("([0-9]+) bits", " ".join(lines)))
    data = filter(lambda x: not x.startswith(">"), lines)
    fp3 = data[0]
    fp4 = data[1:4]
    maccs = data[4:6]
    binstr = ["","",""]
    # fp3                                                                  
    hexstr = fp3
    for block in hexstr.split():
        binstr[0] += bin(int(block, 16))[2:].rjust(32, "0")
    binstr[0] = binstr[0][-55:]

    hexstr = " ".join(fp4)
    for block in hexstr.split():
        binstr[1] += bin(int(block, 16))[2:].rjust(32, "0")
    binstr[1] = binstr[1][-307:]

    hexstr = " ".join(maccs)
    for block in hexstr.split():
        binstr[2] += bin(int(block, 16))[2:].rjust(32, "0")
    binstr[2] = binstr[2][-166:]

    binstr = binstr[0] + binstr[1] + binstr[2]
#        binstr = binstr[2] + binstr[1] + binstr[0]                           
#        binstr = binstr[::-1]                                                 
    try:
        assert sum(checksums) == binstr.count("1")
    except AssertionError:
        binstr = "NULL"
        print cid, "Assertion error!"
    commands.getoutput("rm -f %s %s" % (infn, ofile))
    return binstr

def debug(con):
    cur = con.cursor()
    sql = "SELECT sdf FROM pubchem_sdfs WHERE cid IN (1,2,3);"    
    cur.execute(sql)
    rows = cur.fetchall()
    for row in rows:
        print row    

def insert_mass(con,threadname):
    
    """ parse sdf and insert mass, molecular formula, inchi, fingerprint to a table"""
    cur = con.cursor()
    for cid in xrange(53150001,73250001): # now has 73250000 compounds
#        print 'thread',thread_name,':',cid
        sql1 = "SELECT sdf FROM pubchem_sdfs WHERE cid = %d;" % cid
        cur.execute(sql1)
        rows = cur.fetchall()
        if len(rows) == 0: # don't have this id in db
            continue
        sdf = rows[0][0]

        #######################
#        testsql = "SELECT cid FROM pubchem_info WHERE cid = %d;" % cid
#        cur.execute(testsql)
#        rows = cur.fetchall()
#        if len(rows) != 0: # already in the database
#            continue        
        #######################

        values = parse_sdf_into(threadname,cid,sdf)
        sql2 = "INSERT INTO pubchem_info VALUES " + str(values) + ";"
        cur.execute(sql2)
        if cid % 10000 == 0:
            con.commit()
    con.commit()

def update_mass(con,threadname):
    """ parse sdf and insert mass, molecular formula, inchi, fingerprint to a table"""
    cur = con.cursor()
    for cid in xrange(1,73250001): # now has 73250000 compounds
#        print 'thread',thread_name,':',cid
        sql1 = "SELECT sdf FROM pubchem_sdfs WHERE cid = %d;" % cid
        cur.execute(sql1)
        rows = cur.fetchall()
        if len(rows) == 0: # don't have this id in db
            continue
        sdf = rows[0][0]

        values = parse_sdf_into(threadname,cid,sdf)
        sql2 = "UPDATE pubchem_info SET mass='%s' where cid = %d;" % (values[2],cid)
        cur.execute(sql2)
        if cid % 10000 == 0:
            con.commit()
    con.commit()

def insert_sdf_tmp(con):
    data_dir = "/fs/project/kepaco/celine/Metabolites_identification_3000/data/data_20_ppm/metlin/download_sdf/"
    cur = con.cursor()
    for f in os.listdir(data_dir):
        cid = f[14:f.find(".")]
        data = open(data_dir+f).read()
        end_pos = data.find("$$$$",0)
        sdf_s = data[0:end_pos-1]
        sdf_s = sdf_s.replace("'","")
        cid_s = re.findall("<PUBCHEM_COMPOUND_CID>\n([0-9]+)",sdf_s)[0]
        insert_sql = "UPDATE pubchem_sdfs SET sdf='%s'WHERE cid='%s';" % (sdf_s,cid_s)
        #print cid_s, insert_sql
        cur.execute(insert_sql)        
    con.commit() # commit 


# deal with db connection
con = None
try:     
    con = psycopg2.connect(host='piki',database='shendb',user='shenh1',
                           port='5432', password='pubchem')
    #cur = con.cursor()
    #cur.execute('SELECT version()')          
    #ver = cur.fetchone()
    #print ver    
#    update_mass(con,'main')
#    debug(con)
    insert_sdf_tmp(con)
except psycopg2.DatabaseError, e:
    print 'Error %s' % e    
    sys.exit(1)        
finally:
    if con:
        con.close()

