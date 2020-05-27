import json,multiprocessing,subprocess

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
    for n in range(len(qcoutput)):
        data['job_'+str(n)] = {}
        # Assemble the compressed dictionary of results
        for info in qcoutput[n].as_dict()['data'].keys():
            if info in requested_info:
                data['job_'+str(n)][info] = qcoutput[n].as_dict()['data'][info]
    for job in data:
        if data[job]['errors'] != None and data[job]['errors'] != []:
            raise Exception('Errors detected in one or more of the jobs: ' + str(data[job]['errors']))

    # Return the reduced results in JSON compression
    return json.dumps(data)
    
    
def encode_to_QCInput(encode,rem,pcm=None,solvent=None):
    """Takes final geometry from encode of last job, creates QCInput object
    """
    data = json.loads(encode, encoding='utf-8')
    data = data['job_'+str(len(data)-1)]
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
    

def run_QChem(label,encode=None,rem=None,pcm=None,solvent=None,more_info=None, self_correct=True):
    inname = label + '.inp'
    outname = label + '.out'
    logname = label + '.log'
    command='qchem'
    handlers = [QChemErrorHandler(input_file=inname,output_file=outname)]
    """If no encoding provided, assume this is the first Firework in workflow and that input file is already written.
    'label' is the name of the file without the extension (e.g. .inp, .out). otherwise, take encoding, 
    form new QCInput and write input file, then run.
    """   
    if encode!= None:
        qcin = encode_to_QCInput(encode=encode,rem=rem,pcm=pcm,solvent=solvent)
        qcin.write_file(inname)
    
    if self_correct:
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
    else: 
        job = QCJob(
                input_file=inname,
                output_file=outname,
                qchem_command = command,
                max_cores = multiprocessing.cpu_count(),
                qclog_file=logname
        )
        job.setup()
        p=job.run()
        p.wait()
        """
        qclog = open(logname, "w")
        current_command = ['qchem', '-nt', '20',inname]
        print(current_command)
        subprocess.run(current_command, stdout=qclog, shell=True)
        """

    try:
        output = [QCOutput(filename=outname)]
    except:
        output = QCOutput.multiple_outputs_from_file(QCOutput,filename)
    return QCOutput_to_encode(output,more_info=more_info)