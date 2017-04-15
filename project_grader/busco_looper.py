from os import walk, curdir
from subprocess import Popen, PIPE
from re import search


def walk_directory(directory_path):
    files = []
    for dirpath, dirnames, filenames in walk(directory_path):
        for file in filenames:
            if not file.endswith('.fasta'):
                continue
            files.append(file)
    return files


def call_process(command):
    p = Popen(command, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    stdout = stdout.decode()
    return stdout


def run_busco(fasta_file):
    match = search('(.*)[.]', fasta_file)
    name = match.group(1)
    if name == 'elizabeth_transfuse' or name == 'kaitlyn_transfuse':
        busco_cmd = 'python3 BUSCO.py -i {} -o {}_busco -l ../mammalia_odb9 -m tran --cpu 36'.format(fasta_file, name)
    elif name == 'nate_transfuse':
        busco_cmd = 'python3 BUSCO.py -i {} -o {}_busco -l ../metazoa_odb9 -m tran --cpu 36'.format(fasta_file, name)
    else:
        busco_cmd = 'python3 BUSCO.py -i {} -o {}_busco -l ../eukaryota_odb9 -m tran --cpu 36'.format(fasta_file, name)
    result = call_process(busco_cmd)
    mv_fasta_cmd = 'mv {} run_{}_busco/'.format(fasta_file, name)
    call_process(mv_fasta_cmd)
    with open('run_{}_busco/{}_busco_results.txt'.format(name, name), 'w') as f:
        f.write(result)


if __name__ == '__main__':
    fasta_files = walk_directory(curdir)
    n_files = len(fasta_files)
    print()
    print('-----------------------------------------')
    print('There are {} fasta files to run BUSCO on.'.format(n_files))
    print('-----------------------------------------')
    print()
    counter = 1
    for file in fasta_files:
        if counter == (n_files - 1):
            print('Running BUSCO on fasta file number {}. There is {} more fasta file to run BUSCO on after this one.'.format(counter, n_files - counter))
        else:
            print('Running BUSCO on fasta file number {}. There are {} more fasta files to run BUSCO on after this one.'.format(counter, n_files - counter))
        run_busco(file)
        counter += 1
