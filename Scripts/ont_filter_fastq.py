
import argparse
import gzip

def filter_fastq(input_file, output_file, min_length, max_length):
    if input_file.endswith('.gz'):
        infile = gzip.open(input_file, 'rt')
    else:
        infile = open(input_file, 'rt')
    with infile, gzip.open(output_file, 'wt') as outfile:
        line_count = 0
        header = ''
        sequence = ''
        quality_header = ''
        quality = ''

        for line in infile:
            line_count += 1
            line = line.strip()
            if line_count % 4 == 1:  # Header line
                header = line
            elif line_count % 4 == 2:  # Sequence line
                sequence = line
            elif line_count % 4 == 3:  # Quality header line
                quality_header = line
            elif line_count % 4 == 0:  # Quality line
                quality = line

                seq_length = len(sequence)
                if min_length <= seq_length <= max_length:
                    # Write the FASTQ read to the output file
                    outfile.write(header + '\n')
                    outfile.write(sequence + '\n')
                    outfile.write(quality_header + '\n')
                    outfile.write(quality + '\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Filter FASTQ file based on read length')
    parser.add_argument('--input_file', dest = 'input_file', help='Path to the input FASTQ file')
    parser.add_argument('--output_file', dest = 'output_file', help='Path to the output FASTQ file')
    parser.add_argument('--min_length', dest = 'min_length', default = 1400, type=int, help='minimum read length')
    parser.add_argument('--max_length', dest = 'max_length', default = 1700, type=int, help='maximum read length')
    args = parser.parse_args()
    filter_fastq(args.input_file, args.output_file, args.min_length, args.max_length)
