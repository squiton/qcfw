import numpy as np
import json,multiprocessing

from pymatgen.io.qchem.inputs import QCInput
from pymatgen.io.qchem.outputs import QCOutput
from pymatgen.core import Molecule

from custodian.custodian import Custodian
from custodian.qchem.handlers import QChemErrorHandler
from custodian.qchem.jobs import QCJob

def QCOutput_to_encode(qcoutput,more_info=None):
    """ Converts a qcoutput object to an encoding
    from a .out file.
    """
    #Default recorded info
    requested_info = ['final_energy','errors','initial_geometry','last_geometry', 'species','charge','multiplicity']
    requested_info.append(more_info)

    data = {}
    # Assemble the compressed dictionary of results
    for info in qcoutput.as_dict()['data'].keys():
        if info in requested_info:
            data[info] = qcoutput.as_dict()['data'][info]
    # Return the reduced results in JSON compression
    return json.dumps(data)
    
    
def encode_to_QCInput(encode,rem,pcm=None,solvent=None):
    """Takes final geometry from encode, creates QCInput object
    """
    data = json.loads(encode, encoding='utf-8')
    
    try:
        opt_geom = Molecule(
            species=data.get('species'),
            coords=data.get('last_geometry'),
            charge=data.get('charge'),
            spin_multiplicity=data.get('multiplicity'))
    except KeyError:
         opt_geom = Molecule(
            species=data.get('species'),
            coords=data.get('initial_geometry'),
            charge=data.get('charge'),
            spin_multiplicity=data.get('multiplicity'))       
       
    NewInput = QCInput(molecule = opt_geom,rem=rem,pcm=pcm,solvent=solvent)
    return NewInput
    

def run_QChem(label,encode=None,rem=None,pcm=None,solvent=None):
    inname = label + '.inp'
    outname = label + '.out'
    logname = label + '.log'
    handlers = [QChemErrorHandler(input_file=inname,output_file=outname)]
    """if no encoding provided, assume first Firework in workflow and that input file is already written
    'label' is the name of the file without the extension (e.g. .inp, .out). otherwise, take encoding, 
    form new QCInput and write input file, then run.
    """   
    if encode!= None:
        qcin = encode_to_QCInput(encode=encode,rem=rem,pcm=pcm,solvent=solvent)
        qcin.write_file(inname)
    
    command='qchem'
    jobs = [
        QCJob(
            input_file=inname,
            output_file=outname,
            qchem_command = command,
            max_cores = multiprocessing.cpu_count(),
            qclog_file=logname
        )
    ]
    c = Custodian(handlers, jobs, max_errors=10)
    c.run()
    return parse_out_QChem(label=label)

def run_QChem_EDA(label,rem,encode=None,pcm=None,solvent=None):
    #custom function for determining partition for EDA
    inname = label + '.inp'
    outname = label + '.out'
    logname = label + '.log'
    handlers = [QChemErrorHandler(input_file=inname,output_file=outname)]
    """if no encoding provided, assume first Firework in workflow and that input file is already written
    'label' is the name of the file without the extension (e.g. .inp, .out). otherwise, take encoding, 
    form new QCInput and write input file, then run.
    """   
    """
    if encode!= None:
        if 'benzoate' in label:

            acid = False
            charge = int(benzoate_charges[i])
            multi = int(benzoate_multi[i])
        else:
            i_string = freq.split('-')[0][7:-1]
            i = int(i_string)
            acid = True
            charge = int(acid_charges[i])
            multi = int(acid_multi[i])
            
        qcin = encode_to_QCInput(encode=encode,rem=rem,pcm=pcm,solvent=solvent)
        frag = Fragmenter(molecule=qcin.molecule)
        all_frag_keys=frag.new_unique_frag_dict.keys()
        if acid:
            Benz = 'C7 H5 O2 E14'
            benz_charge = 0
        else:
            Benz = 'C7 H4 O2 E13'
            benz_charge = -1
        if Benz in all_frag_keys:
            k=0
            for n in list(all_frag_keys):
                if n == Benz:
                    if k % 2 == 1:
                        X_key = k - 1
                    else:
                        X_key = k + 1
                    break
                k = k+1
            X = list(all_frag_keys)[X_key]    
            X_mol = frag.new_unique_frag_dict[X][0].molecule
            Benz_mol = frag.new_unique_frag_dict[Benz][0].molecule
        else:
            print(freq + ' fragmentation failed, skipping')
            continue
        
        #Throws a fuss if we try to give it an eda-based charge. Give it a temporary false charge, +1 correct in qcinput
        Benz_mol.set_charge_and_spin(charge = benz_charge + 1, spin_multiplicity = 1)
        try:
            X_mol.set_charge_and_spin(charge = charge-benz_charge -1, spin_multiplicity = multi)
        except:
            print(freq+' has invalid spin_mult/charge for X, skipping')
            continue
        eda_frags = [X_mol, Benz_mol]    
        NewInput = QCInput(molecule = eda_frags, rem = NewRem) 

        #CHANGE FOR BENZOATE
        if acid:
            name_end = 7
        else:
            name_end = 8
        if int(isom[i]) == 1:
            file = freq[0:name_end] + i_string +'m-eda.inp'
            NewInput.write_file(file)
        elif int(isom[i]) == 2:
            file = freq[0:name_end] + i_string +'p-eda.inp'
            NewInput.write_file(file)        
    """
    command='qchem'
    jobs = [
        QCJob(
            input_file=inname,
            output_file=outname,
            qchem_command = command,
            max_cores = multiprocessing.cpu_count(),
            qclog_file=logname
        )
    ]
    c = Custodian(handlers, jobs, max_errors=10)
    c.run()
    return parse_out_QChem(label=label,keys=keys)

#Turns output into encoding
def parse_out_QChem(label):
    output = QCOutput(filename=label+'.out')
    return QCOutput_to_encode(output)