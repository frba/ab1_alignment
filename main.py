import os
from Bio import SeqIO, Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def group_2fasta_files(align_path_result_file, output_template_path, output_seq_fasta_path):
    result_file = open(align_path_result_file, 'w')
    print('Clone_name' + '\t' + 'Template_length' + '\t' + 'Score' + '\t' + '%')
    result_file.write('Clone_name' + 'Result' + '\n')
    for template_file in os.listdir(output_template_path):
        template = template_file[:-3]
        for sequence in os.listdir(output_seq_fasta_path):
            '''Matching sequence and template names'''
            if sequence.__contains__(template):
                template_read = SeqIO.read(open(os.path.join(output_template_path, template_file)), "fasta")
                sequence_read = SeqIO.read(open(os.path.join(output_seq_fasta_path, sequence)), "fasta")

                aligner = Align.PairwiseAligner()
                aligner.mode = 'local'
                aligner.match_score = 1
                aligner.mismatch_score = -1
                aligner.open_gap_score = -1
                aligner.extend_gap_score = -.1

                alignments = aligner.align(template_read.seq.upper(), sequence_read.seq.upper())

                try:
                    align = alignments[0]
                    # print(align.__format__('phyllip'))
                    perc = round((align.score/len(template_read.seq))*100,0)
                    print(str(sequence) + '\t' + str(len(template_read.seq)) + '\t' + str(round(align.score, 0)) + '\t' + str(perc))

                    if perc == 100:
                        result_file.write(sequence + ',' + 'Pass' + '\n')
                    else:
                        result_file.write(sequence + ',' + 'Fail' + '\n')
                except:
                    print(str(sequence) + '\t' + 'alignment_error' + '\t' + '0' + '\t' + '0')
                    result_file.write(sequence + ',' + 'Fail' + '\n')


def convert_ab12fasta(ab1_folder, output_seq_fasta_path):
    '''Convert ab1 to fasta using biopython'''
    for ab1_file in os.listdir(ab1_folder):
        try:
            records = SeqIO.parse(os.path.join(ab1_folder, ab1_file), "abi")
            count = SeqIO.write(records, os.path.join(output_seq_fasta_path, ab1_file[:-4] + '.fa'), "fasta")
            print("Converted %i records" % count)
        except:
            print(ab1_file + ' could not be converted.')


def create_fasta_from_csv(template_csv, output_path):
    '''This function converts a string listed in a csv file to a fasta file'''
    '''The fasta file created is used as a template sequence and used to align against the chromatogram .ab1'''
    header = template_csv.readline()
    for line in template_csv:
        sample_name, full_upstream, hU6, spacer, scaffold, prime_extension, directly_downstream, \
        full_downstream, concatenation, template_seq, insert_only, short_seq, full_seq, use_full_seq = line.strip(
            '\n').split(',')

        sequence_string = insert_only.upper()
        sequence_object = Seq(sequence_string)

        # Create a record
        record = SeqRecord(sequence_object,
                       id=sample_name,  # random accession number
                       name=sample_name,
                       description='')

        # Save as Fasta file
        file_path = output_path + str(sample_name) + ".fa"
        output_file = open(file_path, 'w')
        SeqIO.write(record, output_file, 'fasta')


align_path_result_file = '/path/to/result/file'
pathcsvfile = '/path/to/csv/file'
output_template_path = '/path/to/template/folder/'
ab1_folder = '/path/to/ab1/file'
output_seq_fasta_path = '/path/to/template/folder/'
template_csv = open(os.path.join(pathcsvfile), 'r')


if __name__ == '__main__':
    '''create template fasta file from csv'''
    create_fasta_from_csv(template_csv, output_template_path)

    ''' Convert ab1 to fasta files'''
    convert_ab12fasta(ab1_folder, output_seq_fasta_path)

    '''Combine fasta files and run the alignment'''
    group_2fasta_files(align_path_result_file, output_template_path, output_seq_fasta_path)

