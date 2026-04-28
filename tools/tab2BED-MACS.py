import pandas as pd
import sys
import os

def convert_tab_to_bed(tab_file, bed_file, ext_size=500):
    if not os.path.exists(tab_file):
        print(f"Error: could not read file {tab_file}")
        sys.exit(1)

    # Use error_bad_lines=False or on_bad_lines='skip' if your file has inconsistent tabs
    try:
        df = pd.read_csv(tab_file, sep='\t', low_memory=False)
    except Exception as e:
        print(f"Error processing file: {e}")
        sys.exit(1)

    # Convert column names to lowercase to match logic
    df.columns = map(str.lower, df.columns)

    with open(bed_file, 'w') as f:
        for i, row in df.iterrows():
            rname = str(row['rname'])
            
            # Use 'junction' column
            try:
                junction = int(row['junction'])
            except (ValueError, TypeError):
                continue # Skip rows with invalid junction data

            # Get the strand and determine if it counts as positive or negative
            raw_strand = str(row['strand']).strip()
            
            # Determine logic based on Perl script:
            # Perl '1' = Python '+' or '1'
            # Perl '-1' (or anything else) = Python '-' or '-1'
            is_plus_strand = raw_strand in ['1', '+']

            if is_plus_strand:
                # Logic for "+" strand (Perl strand == 1)
                rstart = junction - ext_size
                rend = junction
                strand_char = "-"  # Note: your Perl script flips these
            else:
                # Logic for "-" strand
                rstart = junction - 1
                rend = junction + ext_size - 1
                strand_char = "+"

            # BED format: chrom, start, end, name, score, strand
            line = f"{rname}\t{max(0, rstart)}\t{rend}\tTLX_{i+1}\t0\t{strand_char}\n"
            f.write(line)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python tab2BED-MACS.py <input.tab> <output.bed> [extsize]")
        sys.exit(1)

    tab_in = sys.argv[1]
    bed_out = sys.argv[2]
    ext = int(sys.argv[3]) if len(sys.argv) > 3 else 500

    convert_tab_to_bed(tab_in, bed_out, ext)