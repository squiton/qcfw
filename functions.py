import numpy as np
import json,multiprocessing

from ase.calculators.singlepoint import SinglePointCalculator as SPC
from ase.io import write
from ase import Atoms
from ase.constraints import dict2constraint
from ase.calculators.qchem import QChem

from pymatgen.io.qchem.inputs import QCInput
from pymatgen.io.qchem.outputs import QCOutput
from pymatgen.io.ase import AseAtomsAdaptor as PyASE
from pymatgen.core import Molecule

from custodian.custodian import Custodian
from custodian.qchem.handlers import QChemErrorHandler
from custodian.qchem.jobs import QCJob

def atoms_to_encode(images):
    """ Converts an list of atoms objects to an encoding
    from a .traj file.
    """

    if not isinstance(images, list):
        images = [images]

    # Convert all constraints into dictionary format
    constraints = [_.todict() for _ in images[0].constraints]
    for i, C in enumerate(constraints):

        # Turn any arrays in the kwargs into lists
        for k, v in list(C['kwargs'].items()):
            if isinstance(v, np.ndarray):
                constraints[i]['kwargs'][k] = v.tolist()

    # Convert any arrays from the parameter settings into lists
    keys = images[0].info
    for k, v in list(keys.items()):
        if isinstance(v, np.ndarray):
            keys[k] = v.tolist()

    data = {'trajectory': {}}
    # Assemble the compressed dictionary of results
    for i, atoms in enumerate(images):

        if i == 0:
            # For first images, collect cell and positions normally
            pos = atoms.get_positions()
            update_pos = pos

            #cell = atoms.get_cell()
            #update_cell = cell

            # Add the parameters which do not change
            data['numbers'] = images[0].get_atomic_numbers().tolist()
            data['pbc'] = images[0].get_pbc().tolist()
            data['constraints'] = constraints
            data['calculator_parameters'] = keys

        else:
            # For consecutive images, check for duplication
            # If duplicates are found, do not store it
            if np.array_equal(atoms.get_positions(), pos):
                update_pos = np.array([])
            else:
                pos = atoms.get_positions()
                update_pos = pos

            #if np.array_equal(atoms.get_cell(), cell):
                #update_cell = np.array([])
            #else:
                #cell = atoms.get_cell()
                #update_cell = cell

        if atoms._calc:
            '''
            nrg = atoms.get_potential_energy()
            #force = atoms.get_forces()
            #stress = atoms.get_stress()

            # Stage results and convert to lists in needed
            results = {
            'positions': update_pos,
            #'cell': update_cell,
            'energy': nrg,
            'forces': force,
            'stress': None}
            '''
            #results = {'positions': update_pos}.update(atoms._calc.results)
            results = atoms._calc.results
            results['positions'] = update_pos
        else:
            results = {
            'positions': update_pos,
            #'cell': update_cell
            }

        for k, v in list(results.items()):
            if isinstance(v, np.ndarray):
                results[k] = v.tolist()

        # Store trajectory, throwing out None values
        data['trajectory'][i] = {
            k: v for k, v in list(
                results.items()) if v is not None}

    # Return the reduced results in JSON compression
    
    #print(data)
    #return str(data).replace(''',''')
    return json.dumps(data)
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
    
def encode_to_atoms(encode):
    """ Dump the encoding to a local traj file.
    """

    # First, decode the trajectory
    data = json.loads(encode, encoding='utf-8')

    # Construct the initial atoms object
    atoms = Atoms(
    data['numbers'],
    data['trajectory']['0']['positions'],
    #cell=data['trajectory']['0']['cell'],
    pbc=data['pbc'])
    atoms.info = data['calculator_parameters']
    atoms.set_constraint([dict2constraint(_) for _ in data['constraints']])

    # Attach the calculator
    calc = SPC(
    atoms=atoms,
    energy=data['trajectory']['0'].get('energy'),
    forces=data['trajectory']['0'].get('forces'),
    stress=data['trajectory']['0'].get('stress'))
    atoms.set_calculator(calc)

    # Collect the rest of the trajectory information
    images = [atoms]
    for i in range(len(data['trajectory']))[1:]:
        atoms = atoms.copy()

        #if data['trajectory'][str(i)]['cell']:
            #atoms.set_cell(data['trajectory'][str(i)]['cell'])

        if data['trajectory'][str(i)]['positions']:
            atoms.set_positions(data['trajectory'][str(i)]['positions'])

        calc = SPC(
            atoms=atoms,
            energy=data['trajectory'][str(i)].get('energy'),
            forces=data['trajectory'][str(i)].get('forces'),
            stress=data['trajectory'][str(i)].get('stress'))
        atoms.set_calculator(calc)

        images += [atoms]

    # Write the traj file
    #if out_file:
        #write(out_file, images)

    return images[0]
    
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
    
def get_potential_energy(label,encoding=None,keys=None):
    #If no encoding provided, assume that input file is already written
    if encoding != None:
        #obtain atoms object and QChem keys
        mol = encode_to_atoms(encoding)
        #keys = mol.info
        mol.info = keys
    #Construct calculator and write input file, use max number of threads
    calc = QChem(label=label,**keys,nt=multiprocessing.cpu_count())
    if encoding != None:
        calc.write_input(atoms = mol, properties = keys)
        mol.set_calculator(calc)
        #Perform relaxation
        mol.get_potential_energy()
    else:
        calc.calculate()
    #Get the output encoding
    return parse_out(label=label,keys=keys)

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
def parse_out(label,keys):
    output = QCOutput(filename=label+'.out')
    if keys['jobtype']=='opt':
        mol=PyASE.get_atoms(output.data['molecule_from_last_geometry'],info=keys)
    else:
        mol=PyASE.get_atoms(output.data['initial_molecule'],info=keys)
    #if calc == None:
    calc = QChem(label=label,**keys)
    calc.read_results()
    mol.set_calculator(calc)
    encoding = atoms_to_encode(mol)
    return encoding
#Turns output into encoding
def parse_out_QChem(label):
    output = QCOutput(filename=label+'.out')
    return QCOutput_to_encode(output)