#!/usr/bin/env python3

from pathlib import Path
import sys

input_fasta=sys.argv[1]
output_dir=Path(sys.argv[2])

output_dir.mkdir(parents=True,exist_ok=True)

def sanitize_id(header):

    header=header.strip()

    parts=header.split("|")

    parts_no_label=[
        p for p in parts
        if not p.lower().startswith("label=")
    ]

    clean="|".join(parts_no_label)

    return (
        clean
        .replace("|","_")
        .replace("/","_")
        .replace("\\","_")
        .replace(" ","_")
    )


with open(input_fasta) as f:

    seq_name=None
    sequence=[]

    for line in f:

        if line.startswith(">"):

            if seq_name:

                outfile=output_dir/f"{seq_name}.fasta"

                with open(outfile,"w") as out:

                    out.write(f">{seq_name}\n")
                    out.write("".join(sequence))

            header=line[1:].strip()

            seq_name=sanitize_id(header)

            sequence=[]

        else:

            sequence.append(line)

    if seq_name:

        outfile=output_dir/f"{seq_name}.fasta"

        with open(outfile,"w") as out:

            out.write(f">{seq_name}\n")
            out.write("".join(sequence))

print("Query FASTA files written")