aaList = list("ACDEFGHIKLMNPQRSTVWY")
mapping = {
    'A': 'Ala',
    'H': 'His',
    'Y': 'Tyr',
    'R': 'Arg',
    'T': 'Thr',
    'K': 'Lys',
    'M': 'Met',
    'D': 'Asp',
    'N': 'Asn',
    'C': 'Cys',
    'Q': 'Gln',
    'E': 'Glu',
    'G': 'Gly',
    'I': 'Ile',
    'L': 'Leu',
    'F': 'Phe',
    'P': 'Pro',
    'S': 'Ser',
    'W': 'Trp',
    'V': 'Val'
    }

mapping_inv = {v.upper():k for k,v in mapping.items()}