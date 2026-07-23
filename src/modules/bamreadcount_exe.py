import subprocess
import pandas as pd
import os


def bam_readcount_execution(bamreadcount, fasta, MAPQ, BAMfilePATH, coordinates, logger):
    try:
        logger.info(f"Running bam-readcount with: {[bamreadcount, '-f', fasta, '-q', MAPQ, BAMfilePATH, coordinates]}")
        bam_readcount_output = subprocess.check_output([bamreadcount, "-f", fasta, "-q", str(MAPQ), BAMfilePATH, coordinates], stderr=subprocess.DEVNULL).decode("utf-8") 
        output_lines = bam_readcount_output.splitlines()
        parsed_rows = []
        for line in output_lines:
            columns =line.split("\t")
            parsed_rows.append(columns)
        max_col = max(len(row) for row in parsed_rows)
        for i, row in enumerate(parsed_rows):
            if len(row) < max_col:
                parsed_rows[i] += [None] * (max_col - len(row))
        df = pd.DataFrame(parsed_rows, columns=[f"Column_{i+1}" for i in range(max_col)])
        return df, max_col
    except Exception as e:
        if not os.path.exists(fasta):
            logger.error(f"Error bam_readcount_execution(): Reference fasta file '{fasta}' not found.")
        if not os.path.exists(BAMfilePATH):
            logger.error(f"Error bam_readcount_execution(): BAM file '{BAMfilePATH}' not found.")
        if not os.path.exists(bamreadcount):
            logger.error(f"Error bam_readcount_execution(): bam-readcount executable file '{bamreadcount}' not found.")
        df = pd.DataFrame()
        return df, 0
        

def bam_readcount_data_extraction(df, individual, logger):
    try: 
        bases = ["A", "C", "G", "T", "N"]
        cols = [5, 6, 7, 8, 9] #A: C: G: T: N:
        counts = {}
        percentages = {}
        depth = 0 
        for base, col in zip(bases, cols):
            base_count_name = f"{base}_{individual}"
            counts[base_count_name] = int(df.iloc[0,col].split(":")[1])
            depth = depth + counts[base_count_name]
        for base in bases:
            base_count_name = f"{base}_{individual}"
            base_percentage_name = f"{base}_perc_{individual}"
            percentages[base_percentage_name] = round((counts[base_count_name]/depth)*100,2)
        return counts, percentages, depth
    except Exception as e:
        logger.error(f"Exception - function bam_readcount_data_extraction(): Error during data extraction")
        base_names = [f"{base}_{individual}" for base in ["A", "C", "G", "T", "N"]]
        counts = {name: 0 for name in base_names}
        percentages = {name.replace(individual, f"perc_{individual}"): 0.0 for name in base_names}
        return counts, percentages, 0
