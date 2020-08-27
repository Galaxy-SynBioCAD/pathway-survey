import json
import csv
import tempfile
import tarfile
import shutil
import os
import sys
import glob
import numpy as np
import random

sys.path.insert(0, '/home/mdulac/workspace/Galaxy-SynBioCAD/rpVisualiser/')
#import cli as rpVisualiser
import run as rpVisualiser


sys.path.insert(0, '/home/mdulac/workspace/Galaxy-SynBioCAD/rpReport/')
import rpToolServe as rpReport


#open the files
rp_name_conv = json.load(open('/home/mdulac/workspace/Galaxy-SynBioCAD/rpOptimise/pathway_survey_170620/rp_name_conv.json'))
fi = '/home/mdulac/Downloads/screened_pathways.csv'
screened_pathways = {}
with open(fi) as csv_fi:
    csv_reader = csv.reader(csv_fi)
    next(csv_reader)
    for row in csv_reader:
        screened_pathways[row[0]] = []
        screened_pathways[row[0]].append(row[1])
        screened_pathways[row[0]].append(row[2])
        screened_pathways[row[0]].append(row[3])
        screened_pathways[row[0]].append(row[4])
        screened_pathways[row[0]].append(row[5])


outpath = '/home/mdulac/workspace/Galaxy-SynBioCAD/rpOptimise/screened_pathway_survey_260820/'
#find the files and make the TAR
for measured_x in screened_pathways:
    os.mkdir(os.path.join(outpath, measured_x))
    with tempfile.TemporaryDirectory() as tmpOutputFolder:
        for rp_x in screened_pathways[measured_x]:
            #loop through all the files in the survey to extract the right SBML
            for uid in rp_name_conv[measured_x]:
                if rp_x in rp_name_conv[measured_x][uid]:
                    print('Found: '+str(measured_x)+' - '+str(rp_x)+' --> '+str(uid))
                    #open the corresponding tar file and move the right file in the output folder
                    with tempfile.TemporaryDirectory() as tmpInputFolder:
                        intar = os.path.join('/home/mdulac/workspace/Galaxy-SynBioCAD/rpOptimise/pathway_survey_170620/', 
                                             str(measured_x),
                                             str(uid),
                                             str(uid)+'.tar')
                        print('Trying to open the following tar: '+str(intar))
                        tar = tarfile.open(intar, mode='r')
                        tar.extractall(path=tmpInputFolder)
                        tar.close()
                        #move it
                        shutil.move(os.path.join(tmpInputFolder, 
                                                 rp_name_conv[measured_x][uid][rp_x]+'.sbml.xml'),
                                    os.path.join(tmpOutputFolder,
                                                 rp_x+'.sbml.xml'))
        with tarfile.open(os.path.join(outpath, measured_x, measured_x+'.tar'), mode='w:gz') as ot:
            for sbml_path in glob.glob(tmpOutputFolder+'/*'):
                file_name = str(sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', ''))
                file_name += '.sbml.xml'
                info = tarfile.TarInfo(file_name)
                info.size = os.path.getsize(sbml_path)
                ot.addfile(tarinfo=info, fileobj=open(sbml_path, 'rb'))


def chassisNames(rp2_report):
    print('+++++++++++++++++++')
    print(rp2_report)
    all_res = {}
    ###### report ########
    fh = open(rp2_report)
    for line in fh:
        if 'No Solutions' in line:
            continue
        if line.split(': ')[0]=='Molecule Name':
            all_res['target_name'] = line.split(': ')[1].replace('\n', '')
        if line.split(': ')[0]=='Molecule InChI':
            all_res['target_inchi'] = line.split(': ')[1].replace('\n', '')
        if line.split(': ')[0]=='GEM SBML':
            all_res['gem'] = line.split(': ')[1].replace('.sbml', '').replace('\n', '')
        if line.split(': ')[0]=='TopX':
            all_res['topx'] = int(line.split(': ')[1].replace('\n', ''))
        if line.split(': ')[0]=='Max Steps':
            all_res['max_steps'] = int(line.split(': ')[1].replace('\n', ''))
        if line.split(': ')[0]=='Rule Diameters':
            all_res['rr_diameters'] = line.split(': ')[1].replace('\n', '')
    fh.close()
    return all_res


#make the new tar files with the new names

rp_name_conv_new = {}
all_info = []
alt_path_names = ['A', 'B', 'C', 'D', 'E']
report_path = '/home/mdulac/workspace/Galaxy-SynBioCAD/rpOptimise/pathwayAnalysis_220620'
for measured_x in screened_pathways:
    rp_name_conv_new[measured_x] = {}
    with tempfile.TemporaryDirectory() as tmpInputFolder:
        with tempfile.TemporaryDirectory() as tmpOutputFolder:
            intar = os.path.join(outpath, 
                                 str(measured_x),
                                 str(measured_x)+'.tar')
            tar = tarfile.open(intar, mode='r')
            tar.extractall(path=tmpInputFolder)
            tar.close()
            num_files = len(glob.glob(tmpInputFolder+'/*'))
            hash_name = int(np.abs(hash(measured_x)))
            rp_name_conv_new[measured_x][str(hash_name)] = {}
            for i, y in zip(random.sample(range(0, num_files), num_files), glob.glob(tmpInputFolder+'/*')):
                #print('i: '+str(i)+' --> '+os.path.join(tmpOutputFolder, str(alt_path_names[i])+'.sbml.xml'))
                #print('y: '+str(y))
                shutil.move(y, os.path.join(tmpOutputFolder, str(alt_path_names[i])+'.sbml.xml'))
                rp_name_conv_new[measured_x][str(hash_name)][y.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', '')] = str(alt_path_names[i])
            if os.path.exists(os.path.join(outpath, measured_x, str(hash_name))):
                shutil.rmtree(os.path.join(outpath, measured_x, str(hash_name)))
            os.mkdir(os.path.join(outpath, measured_x, str(hash_name)))
            path_output_tar = os.path.join(outpath, measured_x, str(hash_name), str(hash_name)+'.tar')
            num_sol = len(glob.glob(tmpOutputFolder+'/*'))
            with tarfile.open(path_output_tar, mode='w:gz') as ot:
                for sbml_path in glob.glob(tmpOutputFolder+'/*'):
                    file_name = str(sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', ''))
                    file_name += '.sbml.xml'
                    info = tarfile.TarInfo(file_name)
                    info.size = os.path.getsize(sbml_path)
                    ot.addfile(tarinfo=info, fileobj=open(sbml_path, 'rb'))
            #now generate the html and report
            chassis_res = {}
            print(os.path.join(report_path, measured_x, 'rp2_report.txt'))
            if os.path.exists(os.path.join(report_path, measured_x, 'rp2_report.txt')):
                chassis_res = chassisNames(os.path.join(report_path, measured_x, 'rp2_report.txt'))
                chassis_res['gem'] = chassis_res['gem'].split('_')[0].upper()+'.'+chassis_res['gem'].split('_')[1]
                all_info.append([measured_x, chassis_res['gem'], chassis_res['target_name'], num_sol])
            else:
                print('No results for the following measured pathway: '+str(measured_x))
            rpVisualiser.main(path_output_tar,
                              'tar',
                              chassis_res['gem'],
                              chassis_res['target_name'],
                              str(hash_name),
                              os.path.join(outpath, measured_x, str(hash_name), str(hash_name)+'.html'))
            rpReport.runReport_hdd(path_output_tar,
                                   os.path.join(outpath, measured_x, str(hash_name), 'report.csv'))


with open(os.path.join(outpath, 'rp_name_conv.json'), 'w') as f:
    json.dump(rp_name_conv_new, f)
