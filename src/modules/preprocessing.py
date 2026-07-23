import pandas as pd

# load each trio information in a dataframe
def load_samples_info(file_path, logger):
    sample_info = pd.read_csv(file_path, delimiter=",")
    sample_info['trio'] = sample_info.index // 3
    return sample_info

def extract_variant_info(coordinates):
    chrom, posStart, var = coordinates.split(":")
    ref, alt = var.split("/")
    posEnd_ref = int(posStart) + len(ref) - 1
    posEnd_alt = int(posStart) + len(alt) - 1
    return chrom, posStart, posEnd_ref, posEnd_alt, ref, alt


def log_init(trio, bams_info, variant, logger):
    logger.info(f"\nANALISYS OF: {trio}")
    for key, value in bams_info.items():
        logger.info(f"{key}:{value}")
    logger.info(f"analyzed variant: {variant}")

def is_sheet_empty(sheet):
    for row in sheet.iter_rows(min_row=1, max_col=sheet.max_column, max_row=sheet.max_row):
        for cell in row:
            if cell.value is not None:
                return False
    return True

def remove_empty_sheets(workbook, logger):
    empty_sheets = [sheet for sheet in workbook.sheetnames if is_sheet_empty(workbook[sheet])]
    for sheet_name in empty_sheets:
        logger.warning(f"{sheet_name} sheet is empty: no data for this variants type")
        std = workbook[sheet_name]
        workbook.remove(std) 
    