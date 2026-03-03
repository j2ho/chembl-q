# chembl_curator/protein_filter.py

import os
import shutil
import subprocess
import logging
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Set
import requests
import numpy as np
from dataclasses import dataclass


# Non-biological ligands to exclude (from exclusion.py)
EXCLUDED_LIGANDS = set([
    # Sugars and glycans
    '045', '05L', '07E', '07Y', '08U', '09X', '0BD', '0H0', '0HX', '0LP', '0MK', '0NZ', '0UB', '0V4', '0WK', '0XY', '0YT', '10M', '12E', '145', '147', '149', '14T', '15L', '16F', '16G', '16O', '17T', '18D', '18O', '1CF', '1FT', '1GL', '1GN', '1LL', '1S3', '1S4',
    '1SD', '1X4', '20S', '20X', '22O', '22S', '23V', '24S', '25E', '26O', '27C', '289', '291', '293', '2DG', '2DR', '2F8', '2FG', '2FL', '2GL', '2GS', '2H5', '2HA', '2M4', '2M5', '2M8', '2OS', '2WP', '2WS', '32O', '34V', '38J', '3BU', '3DO', '3DY', '3FM',
    '3GR', '3HD', '3J3', '3J4', '3LJ', '3LR', '3MG', '3MK', '3R3', '3S6', '3SA', '3YW', '40J', '42D', '445', '44S', '46D', '46Z', '475', '48Z', '491', '49A', '49S', '49T', '49V', '4AM', '4CQ', '4GC', '4GL', '4GP', '4JA', '4N2', '4NN', '4QY', '4R1', '4RS',
    '4SG', '4UZ', '4V5', '50A', '51N', '56N', '57S', '5GF', '5GO', '5II', '5KQ', '5KS', '5KT', '5KV', '5L3', '5LS', '5LT', '5MM', '5N6', '5QP', '5SP', '5TH', '5TJ', '5TK', '5TM', '61J', '62I', '64K', '66O', '6BG', '6C2', '6DM', '6GB', '6GP', '6GR', '6K3',
    '6KH', '6KL', '6KS', '6KU', '6KW', '6LA', '6LS', '6LW', '6MJ', '6MN', '6PZ', '6S2', '6UD', '6YR', '6ZC', '73E', '79J', '7CV', '7D1', '7GP', '7JZ', '7K2', '7K3', '7NU', '83Y', '89Y', '8B7', '8B9', '8EX', '8GA', '8GG', '8GP', '8I4', '8LR', '8OQ', '8PK',
    '8S0', '8YV', '95Z', '96O', '98U', '9AM', '9C1', '9CD', '9GP', '9KJ', '9MR', '9OK', '9PG', '9QG', '9S7', '9SG', '9SJ', '9SM', '9SP', '9T1', '9T7', '9VP', '9WJ', '9WN', '9WZ', '9YW', 'A0K', 'A1Q', 'A2G', 'A5C', 'A6P', 'AAL', 'ABD', 'ABE', 'ABF',
    'ABL', 'AC1', 'ACR', 'ACX', 'ADA', 'AF1', 'AFD', 'AFO', 'AFP', 'AGL', 'AH2', 'AH8', 'AHG', 'AHM', 'AHR', 'AIG', 'ALL', 'ALX', 'AMG', 'AMN', 'AMU', 'AMV', 'ANA', 'AOG', 'AQA', 'ARA', 'ARB', 'ARI', 'ARW', 'ASC', 'ASG', 'ASO', 'AXP', 'AXR',
    'AY9', 'AZC', 'B0D', 'B16', 'B1H', 'B1N', 'B2G', 'B4G', 'B6D', 'B7G', 'B8D', 'B9D', 'BBK', 'BBV', 'BCD', 'BDF', 'BDG', 'BDP', 'BDR', 'BEM', 'BFN', 'BG6', 'BG8', 'BGC', 'BGL', 'BGN', 'BGP', 'BGS', 'BHG', 'BM3', 'BM7', 'BMA', 'BMX', 'BND',
    'BNG', 'BNX', 'BO1', 'BOG', 'BQY', 'BS7', 'BTG', 'BTU', 'BW3', 'BWG', 'BXF', 'BXP', 'BXX', 'BXY', 'BZD', 'C3B', 'C3G', 'C3X', 'C4B', 'C4W', 'C5X', 'CBF', 'CBI', 'CBK', 'CDR', 'CE5', 'CE6', 'CE8', 'CEG', 'CEZ', 'CGF', 'CJB', 'CKB', 'CKP',
    'CNP', 'CR1', 'CR6', 'CRA', 'CT3', 'CTO', 'CTR', 'CTT', 'D1M', 'D5E', 'D6G', 'DAF', 'DAG', 'DAN', 'DDA', 'DDL', 'DEG', 'DEL', 'DFR', 'DFX', 'DG0', 'DGO', 'DGS', 'DGU', 'DJB', 'DJE', 'DK4', 'DKX', 'DKZ', 'DL6', 'DLD', 'DLF', 'DLG', 'DNO',
    'DO8', 'DOM', 'DPC', 'DQR', 'DR2', 'DR3', 'DR5', 'DRI', 'DSR', 'DT6', 'DVC', 'DYM', 'E3M', 'E5G', 'EAG', 'EBG', 'EBQ', 'EEN', 'EEQ', 'EGA', 'EMP', 'EMZ', 'EPG', 'EQP', 'EQV', 'ERE', 'ERI', 'ETT', 'EUS', 'F1P', 'F1X', 'F55', 'F58', 'F6P',
    'F8X', 'FBP', 'FCA', 'FCB', 'FCT', 'FDP', 'FDQ', 'FFC', 'FFX', 'FIF', 'FK9', 'FKD', 'FMF', 'FMO', 'FNG', 'FNY', 'FRU', 'FSA', 'FSI', 'FSM', 'FSW', 'FUB', 'FUC', 'FUD', 'FUF', 'FUL', 'FUY', 'FVQ', 'FX1', 'FYJ', 'G0S', 'G16', 'G1P', 'G20', 'G28',
    'G2F', 'G3F', 'G3I', 'G4D', 'G4S', 'G6D', 'G6P', 'G6S', 'G7P', 'G8Z', 'GAA', 'GAC', 'GAD', 'GAF', 'GAL', 'GAT', 'GBH', 'GC1', 'GC4', 'GC9', 'GCB', 'GCD', 'GCN', 'GCO', 'GCS', 'GCT', 'GCU', 'GCV', 'GCW', 'GDA', 'GDL', 'GE1', 'GE3', 'GFP',
    'GIV', 'GL0', 'GL1', 'GL2', 'GL4', 'GL5', 'GL6', 'GL7', 'GL9', 'GLA', 'GLC', 'GLD', 'GLF', 'GLG', 'GLO', 'GLP', 'GLS', 'GLT', 'GM0', 'GMB', 'GMH', 'GMT', 'GMZ', 'GN1', 'GN4', 'GNS', 'GNX', 'GP0', 'GP1', 'GP4', 'GPH', 'GPK', 'GPM', 'GPO',
    'GPQ', 'GPU', 'GPV', 'GPW', 'GQ1', 'GRF', 'GRX', 'GS1', 'GS9', 'GTK', 'GTM', 'GTR', 'GU0', 'GU1', 'GU2', 'GU3', 'GU4', 'GU5', 'GU6', 'GU8', 'GU9', 'GUF', 'GUL', 'GUP', 'GUZ', 'GXL', 'GXV', 'GYE', 'GYG', 'GYP', 'GYU', 'GYV', 'GZL',
    'H1M', 'H1S', 'H2P', 'H3S', 'H53', 'H6Q', 'H6Z', 'HBZ', 'HD4', 'HNV', 'HNW', 'HSG', 'HSH', 'HSJ', 'HSQ', 'HSX', 'HSY', 'HTG', 'HTM', 'HVC', 'IAB', 'IDC', 'IDF', 'IDG', 'IDR', 'IDS', 'IDU', 'IDX', 'IDY', 'IEM', 'IN1', 'IPT', 'ISD', 'ISL', 'ISX', 'IXD',
    'J5B', 'JFZ', 'JHM', 'JLT', 'JRV', 'JSV', 'JV4', 'JVA', 'JVS', 'JZR', 'K5B', 'K99', 'KBA', 'KBG', 'KD5', 'KDA', 'KDB', 'KDD', 'KDE', 'KDF', 'KDM', 'KDN', 'KDO', 'KDR', 'KFN', 'KG1', 'KGM', 'KHP', 'KME', 'KO1', 'KO2', 'KOT', 'KTU', 'L0W', 'L1L',
    'L6S', 'L6T', 'LAG', 'LAH', 'LAI', 'LAK', 'LAO', 'LAT', 'LB2', 'LBS', 'LBT', 'LCN', 'LDY', 'LEC', 'LER', 'LFC', 'LFR', 'LGC', 'LGU', 'LKA', 'LKS', 'LM2', 'LMO', 'LNV', 'LOG', 'LOX', 'LRH', 'LTG', 'LVO', 'LVZ', 'LXB', 'LXC', 'LXZ', 'LZ0', 'M1F', 'M1P',
    'M2F', 'M3M', 'M3N', 'M55', 'M6D', 'M6P', 'M7B', 'M7P', 'M8C', 'MA1', 'MA2', 'MA3', 'MA8', 'MAB', 'MAF', 'MAG', 'MAL', 'MAN', 'MAT', 'MAV', 'MAW', 'MBE', 'MBF', 'MBG', 'MCU', 'MDA', 'MDP', 'MFB', 'MFU', 'MG5', 'MGC', 'MGL', 'MGS',
    'MJJ', 'MLB', 'MLR', 'MMA', 'MN0', 'MNA', 'MQG', 'MQT', 'MRH', 'MRP', 'MSX', 'MTT', 'MUB', 'MUR', 'MVP', 'MXY', 'MXZ', 'MYG', 'N1L', 'N3U', 'N9S', 'NA1', 'NAA', 'NAG', 'NBG', 'NBX', 'NBY', 'NDG', 'NFG', 'NG1', 'NG6', 'NGA', 'NGC',
    'NGE', 'NGK', 'NGR', 'NGS', 'NGY', 'NGZ', 'NHF', 'NLC', 'NM6', 'NM9', 'NNG', 'NPF', 'NSQ', 'NT1', 'NTF', 'NTO', 'NTP', 'NXD', 'NYT', 'OAK', 'OI7', 'OPM', 'OSU', 'OTG', 'OTN', 'OTU', 'OX2', 'P53', 'P6P', 'P8E', 'PA1', 'PAV', 'PDX', 'PH5',
    'PKM', 'PNA', 'PNG', 'PNJ', 'PNW', 'PPC', 'PRP', 'PSG', 'PSV', 'PTQ', 'PUF', 'PZU', 'QDK', 'QIF', 'QKH', 'QPS', 'QV4', 'R1P', 'R1X', 'R2B', 'R2G', 'RAE', 'RAF', 'RAM', 'RAO', 'RB5', 'RBL', 'RCD', 'RER', 'RF5', 'RG1', 'RGG', 'RHA', 'RHC',
    'RI2', 'RIB', 'RIP', 'RM4', 'RP3', 'RP5', 'RP6', 'RR7', 'RRJ', 'RRY', 'RST', 'RTG', 'RTV', 'RUG', 'RUU', 'RV7', 'RVG', 'RVM', 'RWI', 'RY7', 'RZM', 'S7P', 'S81', 'SA0', 'SCG', 'SCR', 'SDY', 'SEJ', 'SF6', 'SF9', 'SFU', 'SG4', 'SG5', 'SG6', 'SG7',
    'SGA', 'SGC', 'SGD', 'SGN', 'SHB', 'SHD', 'SHG', 'SIA', 'SID', 'SIO', 'SIZ', 'SLB', 'SLM', 'SLT', 'SMD', 'SN5', 'SNG', 'SOE', 'SOG', 'SOL', 'SOR', 'SR1', 'SSG', 'SSH', 'STW', 'STZ', 'SUC', 'SUP', 'SUS', 'SWE', 'SZZ', 'T68', 'T6D', 'T6P',
    'T6T', 'TA6', 'TAG', 'TCB', 'TDG', 'TEU', 'TF0', 'TFU', 'TGA', 'TGK', 'TGR', 'TGY', 'TH1', 'TM5', 'TM6', 'TMR', 'TMX', 'TNX', 'TOA', 'TOC', 'TQY', 'TRE', 'TRV', 'TS8', 'TT7', 'TTV', 'TU4', 'TUG', 'TUJ', 'TUP', 'TUR', 'TVD', 'TVG', 'TVM', 'TVS',
    'TVV', 'TVY', 'TW7', 'TWA', 'TWD', 'TWG', 'TWJ', 'TWY', 'TXB', 'TYV', 'U1Y', 'U2A', 'U2D', 'U63', 'U8V', 'U97', 'U9A', 'U9D', 'U9G', 'U9J', 'U9M', 'UAP', 'UBH', 'UBO', 'UDC', 'UEA', 'V3M', 'V3P', 'V71', 'VG1', 'VJ1', 'VJ4', 'VKN', 'VTB',
    'W9T', 'WIA', 'WOO', 'WUN', 'WZ1', 'WZ2', 'X0X', 'X1P', 'X1X', 'X2F', 'X2Y', 'X34', 'X6X', 'X6Y', 'XDX', 'XGP', 'XIL', 'XKJ', 'XLF', 'XLS', 'XMM', 'XS2', 'XXM', 'XXR', 'XXX', 'XYF', 'XYL', 'XYP', 'XYS', 'XYT', 'XYZ', 'YDR', 'YIO', 'YJM', 'YKR',
    'YO5', 'YX0', 'YX1', 'YYB', 'YYH', 'YYJ', 'YYK', 'YYM', 'YYQ', 'YZ0', 'Z0F', 'Z15', 'Z16', 'Z2D', 'Z2T', 'Z3K', 'Z3L', 'Z3Q', 'Z3U', 'Z4K', 'Z4R', 'Z4S', 'Z4U', 'Z4V', 'Z4W', 'Z4Y', 'Z57', 'Z5J', 'Z5L', 'Z61', 'Z6H', 'Z6J', 'Z6W', 'Z8H', 'Z8T', 'Z9D',
    'Z9E', 'Z9H', 'Z9K', 'Z9L', 'Z9M', 'Z9N', 'Z9W', 'ZB0', 'ZB1', 'ZB2', 'ZB3', 'ZCD', 'ZCZ', 'ZD0', 'ZDC', 'ZDO', 'ZEE', 'ZEL', 'ZGE', 'ZMR',
    # Ions
    '118', '119', '1AL', '1CU', '2FK', '2HP', '2OF', '3CO', '3MT', '3NI', '3OF', '4MO', '4PU', '4TI', '543', '6MO', 'AG', 'AL', 'ALF', 'AM', 'ATH', 'AU', 'AU3', 'AUC', 'BA', 'BEF', 'BF4', 'BO4', 'BR', 'BS3', 'BSY', 'CA', 'CAC', 'CD', 'CD1', 'CD3', 'CD5', 'CE',
    'CF', 'CHT', 'CO', 'CO5', 'CON', 'CR', 'CS', 'CSB', 'CU', 'CU1', 'CU2', 'CU3', 'CUA', 'CUZ', 'CYN', 'DME', 'DMI', 'DSC', 'DTI', 'DY', 'E4N', 'EDR', 'EMC', 'ER3', 'EU', 'EU3', 'F', 'FE', 'FE2', 'FPO', 'GA', 'GD3', 'GEP', 'HAI', 'HG', 'HGC', 'HO3',
    'IN', 'IR', 'IR3', 'IRI', 'IUM', 'K', 'KO4', 'LA', 'LCO', 'LCP', 'LI', 'LU', 'MAC', 'MG', 'MH2', 'MH3', 'MMC', 'MN', 'MN3', 'MN5', 'MN6', 'MO', 'MO1', 'MO2', 'MO3', 'MO4', 'MO5', 'MO6', 'MOO', 'MOS', 'MOW', 'MW1', 'MW2', 'MW3', 'NA2', 'NA5',
    'NA6', 'NAO', 'NAW', 'NET', 'NI', 'NI1', 'NI2', 'NI3', 'NO2', 'NRU', 'O4M', 'OAA', 'OC1', 'OC2', 'OC3', 'OC4', 'OC5', 'OC6', 'OC7', 'OC8', 'OCL', 'OCM', 'OCN', 'OCO', 'OF1', 'OF2', 'OF3', 'OH', 'OS', 'OS4', 'OXL', 'PB', 'PBM', 'PD', 'PER',
    'PI', 'PO3', 'PR', 'PT', 'PT4', 'PTN', 'RB', 'RH3', 'RHD', 'RU', 'SB', 'SE4', 'SEK', 'SM', 'SMO', 'SO3', 'T1A', 'TB', 'TBA', 'TCN', 'TEA', 'TH', 'THE', 'TL', 'TMA', 'TRA', 'V', 'VN3', 'VO4', 'W', 'WO5', 'Y1', 'YB', 'YB2', 'YH', 'YT3', 'ZCM', 'ZN', 'ZN2',
    'ZN3', 'ZNO', 'ZO3', 'ZR',
    # Buffers and crystallization aids
    '144', '15P', '1PE', '2F2', '2JC', '3HR', '3SY', '7N5', '7PE', '9JE', 'AAE', 'ABA', 'ACE', 'ACN', 'ACT', 'ACY', 'AZI', 'BAM', 'BCN', 'BCT', 'BDN', 'BEN', 'BME', 'BO3', 'BTB', 'BTC', 'BU1', 'C8E', 'CAD', 'CAQ', 'CBM', 'CCN', 'CIT', 'CL', 'CLR',
    'CM', 'CMO', 'CO3', 'CPT', 'CXS', 'D10', 'DEP', 'DIO', 'DMS', 'DN', 'DOD', 'DOX', 'EDO', 'EEE', 'EGL', 'EOH', 'EOX', 'EPE', 'ETF', 'FCY', 'FJO', 'FLC', 'FMT', 'FW5', 'GOL', 'GSH', 'GTT', 'GYF', 'HED', 'IHP', 'IHS', 'IMD', 'IOD', 'IPA', 'IPH',
    'LDA', 'MB3', 'MEG', 'MES', 'MLA', 'MLI', 'MOH', 'MPD', 'MRD', 'MSE', 'MYR', 'N', 'NA', 'NH2', 'NH4', 'NHE', 'NO3', 'O4B', 'OHE', 'OLA', 'OLC', 'OMB', 'OME', 'OXA', 'P6G', 'PE3', 'PE4', 'PEG', 'PEO', 'PEP', 'PG0', 'PG4', 'PGE', 'PGR',
    'PLM', 'PO4', 'POL', 'POP', 'PVO', 'SAR', 'SCN', 'SEO', 'SEP', 'SIN', 'SO4', 'SPD', 'SPM', 'SR', 'STE', 'STO', 'STU', 'TAR', 'TBU', 'TME', 'TPO', 'TRS', 'UNK', 'UNL', 'UNX', 'UPL', 'URE',
    # Water
    'HOH', 'H2O', 'WAT'
])


@dataclass
class PDBInfo:
    pdb_id: str
    method: str
    resolution: str
    chains: str  # e.g. "A=1-250,B=1-250"


@dataclass
class LigandInfo:
    ligand_name: str
    center: np.ndarray  # 3D coordinates
    chain: str


class ProteinFilter:
    """Filter protein structures based on PDB availability and binding site analysis."""

    def __init__(self, curated_dir: Path, log_level: str = "INFO"):
        """
        Args:
            curated_dir: Directory containing target subdirectories (uniprot IDs)
            log_level: Logging level
        """
        self.curated_dir = Path(curated_dir)
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(getattr(logging, log_level.upper()))

        # Setup console handler if not already configured
        if not self.logger.handlers:
            handler = logging.StreamHandler()
            handler.setLevel(getattr(logging, log_level.upper()))
            formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            handler.setFormatter(formatter)
            self.logger.addHandler(handler)

    def get_target_list(self) -> List[str]:
        """Get list of all target uniprot IDs from directory."""
        targets = []
        for item in self.curated_dir.iterdir():
            if item.is_dir():
                targets.append(item.name)
        return targets

    def fetch_uniprot_info(self, uniprot_id: str) -> Tuple[List[PDBInfo], Optional[str]]:
        """
        Fetch PDB structures and canonical sequence for a UniProt ID in one API call.

        Returns:
            (pdb_list, canonical_sequence) — sequence is None if unavailable
        """
        try:
            url = f'https://rest.uniprot.org/uniprotkb/{uniprot_id}'
            response = requests.get(url, headers={'Accept': 'application/json'})
            response.raise_for_status()
            data = response.json()

            pdb_list = []
            if 'uniProtKBCrossReferences' in data:
                for xref in data['uniProtKBCrossReferences']:
                    if xref['database'] == 'PDB':
                        pdb_id = xref['id']
                        properties = {p['key']: p['value'] for p in xref.get('properties', [])}
                        pdb_list.append(PDBInfo(
                            pdb_id=pdb_id,
                            method=properties.get('Method', '-'),
                            resolution=properties.get('Resolution', '-'),
                            chains=properties.get('Chains', '-')
                        ))

            sequence = None
            if 'sequence' in data and 'value' in data['sequence']:
                sequence = data['sequence']['value']

            return pdb_list, sequence

        except Exception as e:
            self.logger.error(f"Error fetching UniProt info for {uniprot_id}: {e}")
            return [], None

    def fetch_pdb_list(self, uniprot_id: str) -> List[PDBInfo]:
        """Fetch PDB structures for a given uniprot ID (backward-compatible wrapper)."""
        pdb_list, _ = self.fetch_uniprot_info(uniprot_id)
        return pdb_list

    def save_pdb_list(self, target_dir: Path, pdb_list: List[PDBInfo]):
        """Save PDB list to file."""
        pdb_list_file = target_dir / "pdbid.list"
        with open(pdb_list_file, 'w') as f:
            f.write("# PDBID method resolution chains\n")
            for pdb in pdb_list:
                f.write(f"{pdb.pdb_id}  {pdb.method}  {pdb.resolution}  {pdb.chains}\n")

    def download_pdb(self, pdb_id: str, output_dir: Path) -> bool:
        """
        Download PDB file using pdb_get command (if available) or RCSB web download.

        Args:
            pdb_id: PDB ID to download
            output_dir: Directory to save PDB file

        Returns:
            True if successful, False otherwise
        """
        pdb_file = output_dir / f"{pdb_id.lower()}.pdb"

        # Method 1: Try pdb_get (fast if available on cluster)
        if shutil.which('pdb_get'):
            try:
                result = subprocess.run(
                    ['pdb_get', pdb_id],
                    cwd=output_dir,
                    capture_output=True,
                    text=True,
                    timeout=60
                )

                # Check if PDB file was created
                if pdb_file.exists():
                    self.logger.debug(f"Downloaded {pdb_id} using pdb_get")
                    return True

                # pdb_get might create files with different naming
                pdb_files = list(output_dir.glob(f"*{pdb_id.lower()}*"))
                if pdb_files:
                    # Rename to standard format
                    pdb_files[0].rename(pdb_file)
                    self.logger.debug(f"Downloaded {pdb_id} using pdb_get")
                    return True

            except Exception as e:
                self.logger.debug(f"pdb_get failed for {pdb_id}: {e}, trying web download")

        # Method 2: Download from RCSB PDB (fallback or if pdb_get not available)
        try:
            url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
            self.logger.debug(f"Downloading {pdb_id} from RCSB: {url}")

            response = requests.get(url, timeout=30)
            response.raise_for_status()

            with open(pdb_file, 'w') as f:
                f.write(response.text)

            if pdb_file.exists() and pdb_file.stat().st_size > 0:
                self.logger.debug(f"Downloaded {pdb_id} from RCSB")
                return True

            return False

        except Exception as e:
            self.logger.error(f"Error downloading PDB {pdb_id}: {e}")
            return False

    def download_alphafold(self, uniprot_id: str, output_dir: Path) -> bool:
        """
        Download AlphaFold structure from EBI.

        Args:
            uniprot_id: UniProt accession ID
            output_dir: Directory to save AF model

        Returns:
            True if successful, False otherwise
        """
        try:
            url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v6.pdb"
            output_file = output_dir / f"AF-{uniprot_id}.pdb"

            result = subprocess.run(
                ['wget', '-O', str(output_file), url],
                capture_output=True,
                text=True,
                timeout=120
            )

            if output_file.exists() and output_file.stat().st_size > 0:
                return True

            return False

        except Exception as e:
            self.logger.error(f"Error downloading AlphaFold for {uniprot_id}: {e}")
            return False

    def parse_chain_from_pdb_info(self, pdb_info: PDBInfo) -> List[str]:
        """
        Extract chain IDs from PDB info chains string.

        Args:
            pdb_info: PDBInfo object

        Returns:
            List of chain IDs (e.g., ['A', 'B'])
        """
        chains = []
        if pdb_info.chains and pdb_info.chains != '-':
            # Format: "A=1-250,B=1-250" or "A/B=1-250"
            for part in pdb_info.chains.split(','):
                chain_part = part.split('=')[0].strip()
                # Handle multiple chains like "A/B"
                for ch in chain_part.split('/'):
                    if ch.strip():
                        chains.append(ch.strip())
        return chains

    def get_ligands_from_pdb(self, pdb_file: Path, target_chains: List[str]) -> List[LigandInfo]:
        """
        Extract biologically relevant ligands from PDB file.

        Args:
            pdb_file: Path to PDB file
            target_chains: List of target chain IDs

        Returns:
            List of LigandInfo objects with ligand details
        """
        ligands = {}  # Key: (ligand_name, chain), Value: list of coords
        protein_coords = {}  # Key: chain, Value: list of heavy atom coords

        try:
            with open(pdb_file, 'r') as f:
                for line in f:
                    if line.startswith('ATOM'):
                        chain = line[21:22].strip()
                        if chain in target_chains:
                            # Get heavy atom coordinates (not hydrogen)
                            atom_name = line[12:16].strip()
                            if not atom_name.startswith('H'):
                                x = float(line[30:38])
                                y = float(line[38:46])
                                z = float(line[46:54])

                                if chain not in protein_coords:
                                    protein_coords[chain] = []
                                protein_coords[chain].append(np.array([x, y, z]))

                    elif line.startswith('HETATM'):
                        ligand_name = line[17:20].strip()
                        chain = line[21:22].strip()

                        # Skip excluded ligands
                        if ligand_name in EXCLUDED_LIGANDS:
                            continue

                        # Get heavy atom coordinates
                        atom_name = line[12:16].strip()
                        if not atom_name.startswith('H'):
                            x = float(line[30:38])
                            y = float(line[38:46])
                            z = float(line[46:54])

                            key = (ligand_name, chain)
                            if key not in ligands:
                                ligands[key] = []
                            ligands[key].append(np.array([x, y, z]))

        except Exception as e:
            self.logger.error(f"Error parsing PDB file {pdb_file}: {e}")
            return []

        # Check which ligands have contact with target chains (within 4A)
        ligand_infos = []
        for (ligand_name, lig_chain), lig_coords in ligands.items():
            has_contact = False

            for target_chain in target_chains:
                if target_chain not in protein_coords:
                    continue

                prot_coords = protein_coords[target_chain]

                # Check distance between any ligand atom and any protein atom
                for lig_coord in lig_coords:
                    for prot_coord in prot_coords:
                        dist = np.linalg.norm(lig_coord - prot_coord)
                        if dist <= 4.0:
                            has_contact = True
                            break
                    if has_contact:
                        break
                if has_contact:
                    break

            if has_contact:
                # Calculate ligand center
                center = np.mean(lig_coords, axis=0)
                ligand_infos.append(LigandInfo(
                    ligand_name=ligand_name,
                    center=center,
                    chain=lig_chain
                ))

        return ligand_infos

    def get_chain_ligand_contacts(self, pdb_file: Path, target_chains: List[str]) -> Dict[str, List[LigandInfo]]:
        """
        Determine which chains contact which ligands.

        Args:
            pdb_file: Path to PDB file
            target_chains: List of target chain IDs

        Returns:
            Dictionary mapping chain ID to list of ligands that contact it
        """
        ligands = {}  # Key: (ligand_name, chain), Value: list of coords
        protein_coords = {}  # Key: chain, Value: list of heavy atom coords

        try:
            with open(pdb_file, 'r') as f:
                for line in f:
                    if line.startswith('ATOM'):
                        chain = line[21:22].strip()
                        if chain in target_chains:
                            # Get heavy atom coordinates (not hydrogen)
                            atom_name = line[12:16].strip()
                            if not atom_name.startswith('H'):
                                x = float(line[30:38])
                                y = float(line[38:46])
                                z = float(line[46:54])

                                if chain not in protein_coords:
                                    protein_coords[chain] = []
                                protein_coords[chain].append(np.array([x, y, z]))

                    elif line.startswith('HETATM'):
                        ligand_name = line[17:20].strip()
                        chain = line[21:22].strip()

                        # Skip excluded ligands
                        if ligand_name in EXCLUDED_LIGANDS:
                            continue

                        # Get heavy atom coordinates
                        atom_name = line[12:16].strip()
                        if not atom_name.startswith('H'):
                            x = float(line[30:38])
                            y = float(line[38:46])
                            z = float(line[46:54])

                            key = (ligand_name, chain)
                            if key not in ligands:
                                ligands[key] = []
                            ligands[key].append(np.array([x, y, z]))

        except Exception as e:
            self.logger.error(f"Error parsing PDB file {pdb_file}: {e}")
            return {}

        # Check which chains contact which ligands (within 4A)
        chain_ligands = {chain: [] for chain in target_chains}

        for (ligand_name, lig_chain), lig_coords in ligands.items():
            for target_chain in target_chains:
                if target_chain not in protein_coords:
                    continue

                prot_coords = protein_coords[target_chain]
                has_contact = False

                # Check distance between any ligand atom and any protein atom
                for lig_coord in lig_coords:
                    for prot_coord in prot_coords:
                        dist = np.linalg.norm(lig_coord - prot_coord)
                        if dist <= 4.0:
                            has_contact = True
                            break
                    if has_contact:
                        break

                if has_contact:
                    # Calculate ligand center
                    center = np.mean(lig_coords, axis=0)
                    chain_ligands[target_chain].append(LigandInfo(
                        ligand_name=ligand_name,
                        center=center,
                        chain=lig_chain
                    ))

        # Remove chains with no ligand contacts
        chain_ligands = {chain: ligs for chain, ligs in chain_ligands.items() if ligs}

        return chain_ligands

    def is_single_ligand_bound(self, ligands: List[LigandInfo], distance_threshold: float = 10.0) -> bool:
        """
        Check if structure has single ligand or clustered ligands.

        Args:
            ligands: List of LigandInfo objects
            distance_threshold: Maximum distance between ligand centers to be considered clustered

        Returns:
            True if single ligand or all ligands are clustered, False otherwise
        """
        if len(ligands) == 0:
            return False

        if len(ligands) == 1:
            return True

        # Check if all ligands are within distance threshold of each other
        for i in range(len(ligands)):
            for j in range(i + 1, len(ligands)):
                dist = np.linalg.norm(ligands[i].center - ligands[j].center)
                if dist > distance_threshold:
                    return False

        return True

    def align_pdb(self, query_pdb: Path, template_pdb: Path, output_pdb: Path, target_chain: str):
        """
        Align PDB structure to template using TMalign and save a single target protein chain and all ligands.

        Args:
            query_pdb: Query PDB file to align
            template_pdb: Template PDB file (reference)
            output_pdb: Output aligned PDB file
            target_chain: Single chain ID to save for protein atoms (all HETATMs are saved)
        """
        try:
            # Run TMalign
            cmd = ["./bin/TMalign", str(query_pdb), str(template_pdb)]
            result = subprocess.run(cmd, capture_output=True, text=True)
            stdout = result.stdout

            # Parse transformation matrix
            lines = stdout.splitlines()
            t = [0.0, 0.0, 0.0]
            u = [[0.0]*3 for _ in range(3)]

            idx = None
            for i, line in enumerate(lines):
                if "Rotation matrix to rotate Chain-1 to Chain-2" in line or \
                   "-------- rotation matrix to rotate Chain-1 to Chain-2" in line:
                    idx = i
                    break

            if idx is None:
                self.logger.error(f"Rotation matrix not found in TMalign output")
                return

            # Parse matrix (format: " t[i] u[i][0] u[i][1] u[i][2]")
            for i in range(3):
                parts = lines[idx+2+i].split()
                t[i] = float(parts[1])
                u[i][0] = float(parts[2])
                u[i][1] = float(parts[3])
                u[i][2] = float(parts[4])

            # Apply transformation and save target chain (ATOM) and all ligands (HETATM)
            with open(query_pdb, 'r') as fin, open(output_pdb, 'w') as fout:
                for line in fin:
                    if line.startswith("ATOM"):
                        chain = line[21:22].strip()
                        # Only save the specific target chain for protein atoms
                        if chain != target_chain:
                            continue
                    elif line.startswith("HETATM"):
                        # Save all HETATM lines regardless of chain
                        # The contact-checking logic will filter them later
                        pass
                    else:
                        # Skip non-ATOM/HETATM lines
                        continue

                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                    except:
                        fout.write(line)
                        continue

                    # Transform coordinates
                    X = t[0] + u[0][0]*x + u[0][1]*y + u[0][2]*z
                    Y = t[1] + u[1][0]*x + u[1][1]*y + u[1][2]*z
                    Z = t[2] + u[2][0]*x + u[2][1]*y + u[2][2]*z

                    new_line = (
                        line[:30]
                        + f"{X:8.3f}{Y:8.3f}{Z:8.3f}"
                        + line[54:]
                    )
                    fout.write(new_line)

        except Exception as e:
            self.logger.error(f"Error aligning PDB {query_pdb}: {e}")

    def check_same_pocket(self, ligand_centers: List[np.ndarray], distance_threshold: float = 10.0) -> bool:
        """
        Check if all ligands are in the same pocket.

        Args:
            ligand_centers: List of ligand center coordinates
            distance_threshold: Maximum distance to be considered same pocket

        Returns:
            True if all in same pocket, False otherwise
        """
        if len(ligand_centers) <= 1:
            return True

        for i in range(len(ligand_centers)):
            for j in range(i + 1, len(ligand_centers)):
                dist = np.linalg.norm(ligand_centers[i] - ligand_centers[j])
                if dist > distance_threshold:
                    return False

        return True

    def process_target(self, uniprot_id: str) -> bool:
        """
        Process a single target through the entire protein filtering pipeline.

        Steps:
        1. Fetch PDB list
        2. Download PDBs
        3. Download AlphaFold model
        4. Detect ligand-bound structures
        5. Align structures
        6. Check pocket clustering

        Args:
            uniprot_id: Target uniprot ID

        Returns:
            True if target passes all filters (single binding site), False otherwise
        """
        target_dir = self.curated_dir / uniprot_id

        self.logger.info(f"Processing target {uniprot_id}")

        # Step 1: Fetch PDB list and canonical sequence
        pdb_list, sequence = self.fetch_uniprot_info(uniprot_id)
        if not pdb_list:
            self.logger.warning(f"No PDB structures found for {uniprot_id}")
            return False

        self.save_pdb_list(target_dir, pdb_list)
        self.logger.info(f"Found {len(pdb_list)} PDB structures for {uniprot_id}")

        # Step 2 & 3: Create pdb directory and download structures
        pdb_dir = target_dir / "pdb"
        pdb_dir.mkdir(exist_ok=True)

        # Download AlphaFold model
        af_downloaded = self.download_alphafold(uniprot_id, pdb_dir)
        if not af_downloaded:
            self.logger.warning(f"No AlphaFold model available for {uniprot_id}")
            return False

        af_model = pdb_dir / f"AF-{uniprot_id}.pdb"

        # Download PDBs
        downloaded_pdbs = []
        for pdb_info in pdb_list:
            if self.download_pdb(pdb_info.pdb_id, pdb_dir):
                downloaded_pdbs.append(pdb_info)

        if not downloaded_pdbs:
            self.logger.warning(f"No PDB files downloaded for {uniprot_id}")
            return False

        self.logger.info(f"Downloaded {len(downloaded_pdbs)} PDB files")

        # Step 4: Detect chain-ligand contacts for each structure
        # Structure: list of (pdb_info, chain_ligands_dict)
        # where chain_ligands_dict: {chain_id: [LigandInfo, ...]}
        chain_bound_pdbs = []
        for pdb_info in downloaded_pdbs:
            pdb_file = pdb_dir / f"{pdb_info.pdb_id.lower()}.pdb"
            if not pdb_file.exists():
                continue

            # Get target chains
            target_chains = self.parse_chain_from_pdb_info(pdb_info)
            if not target_chains:
                continue

            # Get chain-ligand contacts
            chain_ligands = self.get_chain_ligand_contacts(pdb_file, target_chains)

            # Only keep structures where at least one chain has ligand contacts
            # and each chain with contacts has single/clustered ligands
            valid_chains = {}
            for chain, ligands in chain_ligands.items():
                if self.is_single_ligand_bound(ligands):
                    valid_chains[chain] = ligands

            if valid_chains:
                chain_bound_pdbs.append((pdb_info, valid_chains))

        if not chain_bound_pdbs:
            self.logger.info(f"No chain-ligand bound structures for {uniprot_id}")
            return False

        # Save ligand-bound PDB list
        ligand_bound_file = target_dir / "ligand_bound_pdbs.txt"
        with open(ligand_bound_file, 'w') as f:
            for pdb_info, chain_ligands in chain_bound_pdbs:
                for chain, ligands in chain_ligands.items():
                    ligand_names = ','.join([lig.ligand_name for lig in ligands])
                    f.write(f"{pdb_info.pdb_id}  Chain {chain}: {ligand_names}\n")

        total_chain_structures = sum(len(chain_ligs) for _, chain_ligs in chain_bound_pdbs)
        self.logger.info(f"Found {total_chain_structures} chain-ligand pairs across {len(chain_bound_pdbs)} PDB structures")

        # Step 5: Align structures - create separate file for each chain
        aligned_dir = target_dir / "aligned"
        aligned_dir.mkdir(exist_ok=True)

        # Structure: list of (pdb_info, chain, ligands, aligned_pdb_path)
        aligned_structures = []
        for pdb_info, chain_ligands in chain_bound_pdbs:
            pdb_file = pdb_dir / f"{pdb_info.pdb_id.lower()}.pdb"

            # Align each chain separately
            for chain, ligands in chain_ligands.items():
                aligned_pdb = aligned_dir / f"{pdb_info.pdb_id.lower()}_{chain}.pdb"
                self.align_pdb(pdb_file, af_model, aligned_pdb, chain)

                if aligned_pdb.exists():
                    aligned_structures.append((pdb_info, chain, ligands, aligned_pdb))

        if not aligned_structures:
            self.logger.warning(f"No structures aligned for {uniprot_id}")
            return False

        self.logger.info(f"Aligned {len(aligned_structures)} chain-specific structures")

        # Step 6: Check pocket clustering
        # Re-extract ligand centers from aligned structures
        all_ligand_centers = []
        pocket_info = []

        for pdb_info, chain, original_ligands, aligned_pdb in aligned_structures:
            # Extract ligands from aligned structure
            # Use only the specific chain since the file only contains that chain
            aligned_ligands = self.get_ligands_from_pdb(aligned_pdb, [chain])

            if aligned_ligands:
                for lig in aligned_ligands:
                    all_ligand_centers.append(lig.center)
                    pocket_info.append({
                        'pdb_id': pdb_info.pdb_id,
                        'chain': chain,
                        'aligned_file': aligned_pdb.name,
                        'ligand_name': lig.ligand_name,
                        'center': lig.center
                    })

        # Check if all ligands are in same pocket
        is_single_pocket = self.check_same_pocket(all_ligand_centers)

        # Save pocket info
        import csv
        pocket_csv = target_dir / "pocket_info.csv"
        with open(pocket_csv, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['PDB_ID', 'Chain', 'Aligned_File', 'Ligand_Name', 'Center_X', 'Center_Y', 'Center_Z'])
            for info in pocket_info:
                writer.writerow([
                    info['pdb_id'],
                    info['chain'],
                    info['aligned_file'],
                    info['ligand_name'],
                    f"{info['center'][0]:.3f}",
                    f"{info['center'][1]:.3f}",
                    f"{info['center'][2]:.3f}"
                ])

        if not is_single_pocket:
            self.logger.info(f"Multiple binding pockets detected for {uniprot_id}")
            return False

        # Write per-target sequence FASTA (combined into sequences.fasta at the end)
        if sequence:
            seq_fasta = target_dir / "sequence.fasta"
            with open(seq_fasta, 'w') as f:
                f.write(f">{uniprot_id}\n{sequence}\n")

        self.logger.info(f"Target {uniprot_id} passed all filters (single binding site)")
        return True

    def _get_best_structure(self, target_dir: Path) -> Optional[Tuple[str, float]]:
        """
        Return (pdbid_chain, resolution) for the best (lowest-resolution) structure.

        Reads pocket_info.csv for available aligned structures and pdbid.list for resolutions.
        pdbid_chain format matches the aligned filename stem, e.g. "3mms_A".
        """
        import csv as _csv
        pocket_csv = target_dir / "pocket_info.csv"
        if not pocket_csv.exists():
            return None

        pdbid_chains = set()
        with open(pocket_csv) as f:
            for row in _csv.DictReader(f):
                pdbid_chains.add(row['Aligned_File'].replace('.pdb', ''))

        if not pdbid_chains:
            return None

        # Parse resolutions from pdbid.list
        res_map: Dict[str, float] = {}
        pdbid_list = target_dir / "pdbid.list"
        if pdbid_list.exists():
            with open(pdbid_list) as f:
                for line in f:
                    if line.startswith('#') or not line.strip():
                        continue
                    parts = line.split()
                    if len(parts) >= 3:
                        pdbid = parts[0].upper()
                        try:
                            res_map[pdbid] = float(parts[2].rstrip('A'))
                        except ValueError:
                            res_map[pdbid] = 99.0

        best_chain, best_res = None, 99.0
        for pc in pdbid_chains:
            pdbid = pc.split('_')[0].upper()
            res = res_map.get(pdbid, 99.0)
            if res < best_res:
                best_res = res
                best_chain = pc

        return (best_chain, best_res) if best_chain else None

    def run_pipeline(self, n_processes: int = 1) -> List[str]:
        """
        Run protein filtering pipeline on all targets.

        Args:
            n_processes: Number of parallel processes (default: 1 for sequential)

        Returns:
            List of target IDs that passed filters
        """
        targets = self.get_target_list()
        self.logger.info(f"Processing {len(targets)} targets")

        passed_targets = []

        if n_processes == 1:
            # Sequential processing
            for target in targets:
                try:
                    if self.process_target(target):
                        passed_targets.append(target)
                except Exception as e:
                    self.logger.error(f"Error processing {target}: {e}")
        else:
            # Parallel processing
            import multiprocessing as mp
            with mp.Pool(n_processes) as pool:
                results = pool.map(self._process_target_wrapper, targets)
                passed_targets = [t for t, passed in zip(targets, results) if passed]

        # Save passed targets
        passed_file = self.curated_dir / "passed_targets.txt"
        with open(passed_file, 'w') as f:
            for target in passed_targets:
                f.write(f"{target}\n")

        self.logger.info(f"Pipeline complete: {len(passed_targets)}/{len(targets)} targets passed")

        # Combine per-target sequence.fasta into global sequences.fasta
        # and write best_structure.tsv for pocket RMSD computation
        self._write_global_outputs(passed_targets)

        # Optionally delete filtered-out targets
        for target in targets:
            if target not in passed_targets:
                import shutil
                target_dir = self.curated_dir / target
                if target_dir.exists():
                    shutil.rmtree(target_dir)

        return passed_targets

    def _write_global_outputs(self, passed_targets: List[str]) -> None:
        """Write sequences.fasta and best_structure.tsv for all passed targets."""
        seq_fasta = self.curated_dir / "sequences.fasta"
        best_str_tsv = self.curated_dir / "best_structure.tsv"

        n_seq = 0
        n_best = 0
        with open(seq_fasta, 'w') as sf, open(best_str_tsv, 'w') as bf:
            bf.write("uniprot\tpdbid_chain\tresolution\n")
            for target in passed_targets:
                target_dir = self.curated_dir / target

                # Append sequence
                per_seq = target_dir / "sequence.fasta"
                if per_seq.exists():
                    sf.write(per_seq.read_text())
                    n_seq += 1

                # Append best structure entry
                best = self._get_best_structure(target_dir)
                if best:
                    bf.write(f"{target}\t{best[0]}\t{best[1]:.2f}\n")
                    n_best += 1

        self.logger.info(f"Wrote sequences.fasta ({n_seq} sequences) and "
                         f"best_structure.tsv ({n_best} entries)")

    def _process_target_wrapper(self, target: str) -> bool:
        """Wrapper for multiprocessing."""
        try:
            return self.process_target(target)
        except Exception as e:
            self.logger.error(f"Error processing {target}: {e}")
            return False
