import requests, urllib.parse

# Amino acid base names and their primary substitution sites
AMINO_ACIDS = {
    'lys': ('lysine', 'Nε'),
    'arg': ('arginine', 'Nω'),
    'his': ('histidine', 'N³'),
    'ser': ('serine', 'O'),
    'thr': ('threonine', 'O'),
    'tyr': ('tyrosine', 'O'),
    'cys': ('cysteine', 'S'),
    'asp': ('aspartic acid', 'β'),
    'glu': ('glutamic acid', 'γ'),
    'asn': ('asparagine', 'Nδ'),
    'gln': ('glutamine', 'Nε'),
    'trp': ('tryptophan', 'N¹'),
    'met': ('methionine', 'S'),
    'orn': ('ornithine', 'Nδ'),
    'hse': ('homoserine', 'O'),
}

# Protecting group definitions
PROTECTING_GROUPS = {
    'boc': '(tert-butoxycarbonyl)',
    'z': '(benzyloxycarbonyl)',
    'cbz': '(benzyloxycarbonyl)',
    'fmoc': '(9-fluorenylmethoxycarbonyl)',
    'ac': '(acetyl)',
    'alloc': '(allyloxycarbonyl)',
    'dde': '((1-(4,4-dimethyl-2,6-dioxocyclohexylidene)ethyl))',
    'ivdde': '((1-(4,4-dimethyl-2,6-dioxocyclohexylidene)-3-methylbutyl))',
    'mtt': '(4-methyltrityl)',
    'proc': '(propargyloxycarbonyl)',
    'pbf': '(2,2,4,6,7-pentamethyldihydrobenzofuran-5-sulfonyl)',
    'tos': '(tosyl)',
    'mtr': '(4-methoxy-2,3,6-trimethylbenzenesulfonyl)',
    'pmc': '(2,2,5,7,8-pentamethylchroman-6-sulfonyl)',
    'no2': '(nitro)',
    'trt': 'trityl',
    'bn': 'benzyl',
    'bzl': 'benzyl',
    'dnp': '2,4-dinitrophenyl',
    'tbu': 'tert-butyl',
    'acm': 'acetamidomethyl',
    'stbu': 'tert-butyl',
    'mob': '4-methoxybenzyl',
    'mmt': '4-methoxytrityl',
    'npys': '3-nitro-2-pyridinesulfenyl',
    'xan': '(xanthenyl)',
    'for': 'formyl',
    'hoc': '1-(1-adamantyl)-1-methylethoxycarbonyl',
    'o': 'oxide',
    'o2': ',S-dioxide',
    # Ester forms
    'otbu': 'tert-butyl ester',
    'obut': 'tert-butyl ester',  # typo variant
    'obn': 'benzyl ester',
    'ome': 'methyl ester',
    'oet': 'ethyl ester',
    'oall': 'allyl ester',
    'ofm': '9-fluorenylmethyl ester',
}

# Special cases that don't follow the standard pattern
SPECIAL_CASES = {
    'tyr(2-br-z)': ('tyrosine', 'O-2-bromobenzyloxycarbonyl'),
    'tyr(2,6-cl2-bn)': ('tyrosine', 'O-2,6-dichlorobenzyl'),
    'met(o2)': ('methionine', 'S,S-dioxide'),
}

def generate_side_chain_protections():
    """Generate the complete side chain protections dictionary"""
    protections = {}
    
    # Generate standard combinations
    for aa_code, (aa_name, site) in AMINO_ACIDS.items():
        for pg_code, pg_name in PROTECTING_GROUPS.items():
            key = f"{aa_code}({pg_code})"
            
            # Handle special formatting for different types
            if pg_code.startswith('o'):  # Ester forms
                value = f"{site}-{pg_name}"
            elif aa_code in ['asp', 'glu'] and not pg_code.startswith('o'):
                # For asp/glu without 'o' prefix, assume ester
                value = f"{site}-{pg_name} ester"
            elif aa_code == 'met' and pg_code == 'o':
                value = f"{site}-{pg_name}"
            else:
                value = f"{site}-{pg_name}"
            
            protections[key] = (aa_name, value)
    
    # Add special cases
    protections.update(SPECIAL_CASES)
    
    return protections

# Generate the complete dictionary
SIDE_CHAIN_PROTECTIONS = generate_side_chain_protections()
len(SIDE_CHAIN_PROTECTIONS)

AA_FULL = {'gly': 'glycine',
 'ala': 'alanine',
 'val': 'valine',
 'leu': 'leucine',
 'leuc': 'leucine',
 'ile': 'isoleucine',
 'ileu': 'isoleucine',
 'ser': 'serine',
 'thr': 'threonine',
 'cys': 'cysteine',
 'met': 'methionine',
 'phe': 'phenylalanine',
 'tyr': 'tyrosine',
 'trp': 'tryptophan',
 'asp': 'aspartic acid',
 'glu': 'glutamic acid',
 'asn': 'asparagine',
 'gln': 'glutamine',
 'lys': 'lysine',
 'arg': 'arginine',
 'his': 'histidine',
 'pro': 'proline',
 'hyp': '4-hydroxyproline',
 'orn': 'ornithine',
 'dab': '2,4-diaminobutyric acid',
 'dap': '2,3-diaminopropionic acid',
 'sar': 'sarcosine',
 'abu': '2-aminobutyric acid',
 'aib': '2-aminoisobutyric acid',
 'nva': 'norvaline',
 'nle': 'norleucine',
 'hse': 'homoserine',
 'hcy': 'homocysteine',
 'pen': 'penicillamine',
 'sec': 'selenocysteine',
 'pyl': 'pyrrolysine',
 'aad': '2-aminoadipic acid',
 'aan': 'α-asparagine',
 'aca': '2-aminocapric acid',
 'agn': 'α-glutamine',
 'ahx': '6-aminohexanoic acid',
 'apm': '2-aminopimelic acid',
 'app': 'γ-amino-β-hydroxybenzenepentanoic acid',
 'asu': '2-aminosuberic acid',
 'azagly': '2-azaglycine',
 'aze': '2-carboxyazetidine',
 'bal': 'β-alanine',
 'bas': 'β-aspartic acid',
 'bly': '3,6-diaminohexanoic acid',
 'bua': 'butanoic acid',
 'bux': '4-amino-3-hydroxybutanoic acid',
 'cap': 'γ-amino-β-hydroxycyclohexanepentanoic acid',
 'cha': '3-cyclohexylalanine',
 'cit': 'N5-aminocarbonylornithine',
 'cya': '3-sulfoalanine',
 'dpm': 'diaminopimelic acid',
 'dpr': '2,3-diaminopropanoic acid',
 'dsu': '2,7-diaminosuberic acid',
 'edc': 'S-ethylthiocysteine',
 'ggu': 'γ-glutamic acid',
 'gla': 'γ-carboxyglutamic acid',
 'glc': 'hydroxyacetic acid',
 'glp': 'pyroglutamic acid',
 'pglu': 'pyroglutamic acid',
 'har': 'homoarginine',
 'hhs': 'homohistidine',
 'hiv': '2-hydroxyisovaleric acid',
 'hva': '2-hydroxypentanoic acid',
 'hyl': '5-hydroxylysine',
 'igl': 'α-amino-2-indanacetic acid',
 'iva': 'isovaline',
 'lac': '2-hydroxypropanoic acid',
 'maa': 'mercaptoacetic acid',
 'mba': 'mercaptobutanoic acid',
 'mhp': '3-hydroxy-4-methylproline',
 'mpa': 'mercaptopropanoic acid',
 'nal': '3-naphthylalanine',
 'nphe': 'N-benzylglycine',
 'nty': 'nortyrosine',
 'oaa': 'ω-amino acid',
 'oic': '2-carboxyoctahydroindole',
 '3-pal': '3-(3-pyridyl)-alanine',
 'phg': '2-phenylglycine',
 'pip': '2-carboxypiperidine',
 'pmp': 'β,β-cyclopentamethylene-β-mercaptopropionic acid',
 'pser': 'phospho-L-serine',
 'pthr': 'phospho-L-threonine',
 'ptyr': 'phospho-L-tyrosine',
 'pya(4)': '3-(4-pyridinyl)-alanine',
 'spg': '1-amino-1-carboxycyclopentane',
 'sta': 'statine',
 'thi': '3-thienylalanine',
 'tac': '2-tolylaminophenylcabonyl',
 'tic': '3-carboxyisoquinoline',
 'tle': '3-methylvaline',
 'tml': 'ε-N-trimethyllysine',
 'tza': '3-thiazolylalanine',
 'wil': 'α-amino-2,4-dioxopyrimidinepropanoic acid'}

AA_YL = {'gly': 'glycyl',
 'ala': 'alanyl',
 'val': 'valyl',
 'leu': 'leucyl',
 'leuc': 'leucyl',
 'ile': 'isoleucyl',
 'ileu': 'isoleucyl',
 'ser': 'seryl',
 'thr': 'threonyl',
 'cys': 'cysteinyl',
 'met': 'methionyl',
 'phe': 'phenylalanyl',
 'tyr': 'tyrosyl',
 'trp': 'tryptophyl',
 'asp': 'aspartyl',
 'glu': 'glutamyl',
 'asn': 'asparaginyl',
 'gln': 'glutaminyl',
 'lys': 'lysyl',
 'arg': 'arginyl',
 'his': 'histidyl',
 'pro': 'prolyl',
 'hyp': '4-hydroxyprolyl',
 'orn': 'ornithyl',
 'dab': '2,4-diaminobutyryl',
 'dap': '2,3-diaminopropionyl',
 'sar': 'sarcosyl',
 'abu': '2-aminobutyryl',
 'aib': '2-aminoisobutyryl',
 'nva': 'norvalyl',
 'nle': 'norleucyl',
 'hse': 'homoserinyl',
 'hcy': 'homocysteinyl',
 'pen': 'penicillaminyl',
 'sec': 'selenocysteinyl',
 'pyl': 'pyrrolysinyl',
 'aad': '2-aminoadipyl',
 'aan': 'α-asparaginyl',
 'aca': '2-aminocapryl',
 'agn': 'α-glutaminyl',
 'ahx': '6-aminohexanoyl',
 'apm': '2-aminopimelyl',
 'app': 'γ-amino-β-hydroxybenzenepentanoyl',
 'asu': '2-aminosuberyl',
 'azagly': '2-azaglycyl',
 'aze': '2-carboxyazetidinyl',
 'bal': 'β-alanyl',
 'bas': 'β-aspartyl',
 'bly': '3,6-diaminohexanoyl',
 'bua': 'butanoyl',
 'bux': '4-amino-3-hydroxybutanoyl',
 'cap': 'γ-amino-β-hydroxycyclohexanepentanoyl',
 'cha': '3-cyclohexylalanyl',
 'cit': 'N5-aminocarbonylornithyl',
 'cya': '3-sulfoalanyl',
 'dpm': 'diaminopimelyl',
 'dpr': '2,3-diaminopropanoyl',
 'dsu': '2,7-diaminosuberyl',
 'edc': 'S-ethylthiocysteinyl',
 'ggu': 'γ-glutamyl',
 'gla': 'γ-carboxyglutamyl',
 'glc': 'hydroxyacetyl',
 'glp': 'pyroglutamyl',
 'pglu': 'pyroglutamyl',
 'har': 'homoarginyl',
 'hhs': 'homohistidinyl',
 'hiv': '2-hydroxyisovaleryl',
 'hva': '2-hydroxypentanoyl',
 'hyl': '5-hydroxylysinyl',
 'igl': 'α-amino-2-indanacetyl',
 'iva': 'isovalinyl',
 'lac': '2-hydroxypropanoyl',
 'maa': 'mercaptoacetyl',
 'mba': 'mercaptobutanoyl',
 'mhp': '3-hydroxy-4-methylprolyl',
 'mpa': 'mercaptopropanoyl',
 'nal': '3-naphthylalanyl',
 'nphe': 'N-benzylglycyl',
 'nty': 'nortyrosinyl',
 'oaa': 'ω-aminoacyl',
 'oic': '2-carboxyoctahydroindolyl',
 '3-pal': '3-(3-pyridyl)-alanyl',
 'phg': '2-phenylglycyl',
 'pip': '2-carboxypiperidyl',
 'pmp': 'β,β-cyclopentamethylene-β-mercaptopropionyl',
 'pser': 'phospho-L-serinyl',
 'pthr': 'phospho-L-threoninyl',
 'ptyr': 'phospho-L-tyrosinyl',
 'pya(4)': '3-(4-pyridinyl)-alanyl',
 'spg': '1-amino-1-carboxycyclopentyl',
 'sta': 'statinyl',
 'thi': '3-thienylalanyl',
 'tac': '2-tolylaminophenylcabonyl',
 'tic': '3-carboxyisoquinolinyl',
 'tle': '3-methylvalyl',
 'tml': 'ε-N-trimethyllysinyl',
 'tza': '3-thiazolylalanyl',
 'wil': 'α-amino-2,4-dioxopyrimidinepropanoyl'}

N_CAPS = {'h2n': '',
 'h': '',
 'nh2': '',
 'cbz': 'benzyloxycarbonyl',
 'z': 'benzyloxycarbonyl',
 'boc': 'tert-butoxycarbonyl',
 'tboc': 'tert-butoxycarbonyl',
 'fmoc': '9-fluorenylmethoxycarbonyl',
 'alloc': 'allyloxycarbonyl',
 'troc': '2,2,2-trichloroethoxycarbonyl',
 'teoc': '2-(trimethylsilyl)ethoxycarbonyl',
 'poc': 'propargyloxycarbonyl',
 'moc': 'methoxycarbonyl',
 'eoc': 'ethoxycarbonyl',
 'azoc': '2-azidobenzyloxycarbonyl',
 'nvoc': '6-nitroveratryloxycarbonyl',
 'ddz': '3,5-dimethoxybenzyloxycarbonyl',
 'ac': 'acetyl',
 'piv': 'pivaloyl',
 'bz': 'benzoyl',
 'tos': 'tosyl',
 'ts': 'tosyl',
 'ms': 'mesyl',
 'mes': 'mesyl',
 'tf': 'trifluoromethanesulfonyl',
 'ns': '4-nitrobenzenesulfonyl',
 'dts': 'dithiasuccinoyl',
 'for': 'formyl',
 'tfa': 'trifluoroacetyl',
 'mtr': '4-methoxy-2,3,6-trimethylbenzenesulfonyl',
 'pmc': '2,2,5,7,8-pentamethylchroman-6-sulfonyl',
 'pbf': '2,2,4,6,7-pentamethyldihydrobenzofuran-5-sulfonyl',
 'bn': 'benzyl',
 'pmb': '4-methoxybenzyl',
 'dmb': '2,4-dimethoxybenzyl',
 'tmb': '2,4,6-trimethoxybenzyl',
 'trt': 'trityl',
 'mmt': '4-methoxytrityl',
 'dmt': "4,4'-dimethoxytrityl",
 'myr': 'Myristoyl'}

C_CAPS = {'h': '', #typo?
 'oh': '',
 'ome': 'methyl ester',
 'och3': 'methyl ester',
 'oet': 'ethyl ester',
 'otbu': 'tert-butyl ester',
 'obn': 'benzyl ester',
 'obzl': 'benzyl ester',
 'oallyl': 'allyl ester',
 'oall': 'allyl ester',
 'oph': 'phenyl ester',
 'otmse': '2-(trimethylsilyl)ethyl ester',
 'otbdms': 'tert-butyldimethylsilyl ester',
 'otips': 'triisopropylsilyl ester',
 'ofm': '9-fluorenylmethyl ester',
 'odmb': '2,4-dimethoxybenzyl ester',
 'opmb': '4-methoxybenzyl ester',
 'onp': '4-nitrophenyl ester',
 'opfp': 'pentafluorophenyl ester',
 'osu': 'N-hydroxysuccinimidyl ester',
 'obt': 'benzotriazol-1-yl ester',
 'oat': '7-azabenzotriazol-1-yl ester',
 'nh2': 'amide',
 'nh': 'amine', #typo?
 'nhme': 'N-methylamide',
 'nme2': 'N,N-dimethylamide',
 'nhet': 'N-ethylamide',
 'net2': 'N,N-diethylamide',
 'nhbn': 'N-benzylamide',
 'nhph': 'N-phenylamide',
 'nhnh2': 'hydrazide',
 'nhnh': 'hydrazide', #typo?
 'conhnh2': 'hydrazide',
 'cho': 'aldehyde',
 'ch2oh': 'alcohol',
 'conh2': 'amide',
 'cf3': 'trifluoromethane'}

# SIDE_CHAIN_PROTECTIONS = {'lys(boc)': ('lysine', 'Nε-(tert-butoxycarbonyl)'),
#  'lys(z)': ('lysine', 'Nε-(benzyloxycarbonyl)'),
#  'lys(cbz)': ('lysine', 'Nε-(benzyloxycarbonyl)'),
#  'lys(fmoc)': ('lysine', 'Nε-(9-fluorenylmethoxycarbonyl)'),
#  'lys(ac)': ('lysine', 'Nε-(acetyl)'),
#  'lys(alloc)': ('lysine', 'Nε-(allyloxycarbonyl)'),
#  'lys(dde)': ('lysine',
#  'Nε-((1-(4,4-dimethyl-2,6-dioxocyclohexylidene)ethyl))'),
#  'lys(ivdde)': ('lysine',
#  'Nε-((1-(4,4-dimethyl-2,6-dioxocyclohexylidene)-3-methylbutyl))'),
#  'lys(mtt)': ('lysine', 'Nε-(4-methyltrityl)'),
#  'lys(proc)': ('lysine', 'Nε-(propargyloxycarbonyl)'),
#  'lys(cl-z)': ('lysine', 'Nε-(2-chlorobenzyloxycarbonyl)'),
#  'arg(pbf)': ('arginine',
#  'Nω-(2,2,4,6,7-pentamethyldihydrobenzofuran-5-sulfonyl)'),
#  'arg(boc)': ('arginine', 'Nω-(tert-butoxycarbonyl)'),
#  'arg(z)': ('arginine', 'Nω-(benzyloxycarbonyl)'),
#  'arg(tos)': ('arginine', 'Nω-(tosyl)'),
#  'arg(mtr)': ('arginine', 'Nω-(4-methoxy-2,3,6-trimethylbenzenesulfonyl)'),
#  'arg(pmc)': ('arginine', 'Nω-(2,2,5,7,8-pentamethylchroman-6-sulfonyl)'),
#  'arg(no2)': ('arginine', 'Nω-(nitro)'),
#  'his(trt)': ('histidine', 'N³-trityl'),
#  'his(boc)': ('histidine', 'N³-tert-butoxycarbonyl'),
#  'his(bn)': ('histidine', 'N³-benzyl'),
#  'his(tos)': ('histidine', 'N³-tosyl'),
#  'his(dnp)': ('histidine', 'N³-2,4-dinitrophenyl'),
#  'ser(tbu)': ('serine', 'O-tert-butyl'),
#  'ser(bn)': ('serine', 'O-benzyl'),
#  'ser(trt)': ('serine', 'O-trityl'),
#  'ser(ac)': ('serine', 'O-acetyl'),
#  'thr(tbu)': ('threonine', 'O-tert-butyl'),
#  'thr(bn)': ('threonine', 'O-benzyl'),
#  'thr(ac)': ('threonine', 'O-acetyl'),
#  'tyr(tbu)': ('tyrosine', 'O-tert-butyl'),
#  'tyr(bn)': ('tyrosine', 'O-benzyl'),
#  'tyr(bzl)': ('tyrosine', 'O-benzyl'),
#  'tyr(2-br-z)': ('tyrosine', 'O-2-bromobenzyloxycarbonyl'),
#  'tyr(2,6-cl2-bn)': ('tyrosine', 'O-2,6-dichlorobenzyl'),
#  'tyr(ac)': ('tyrosine', 'O-acetyl'),
#  'cys(trt)': ('cysteine', 'S-trityl'),
#  'cys(tbu)': ('cysteine', 'S-tert-butyl'),
#  'cys(bn)': ('cysteine', 'S-benzyl'),
#  'cys(acm)': ('cysteine', 'S-acetamidomethyl'),
#  'cys(stbu)': ('cysteine', 'S-tert-butyl'),
#  'cys(mob)': ('cysteine', 'S-4-methoxybenzyl'),
#  'cys(mmt)': ('cysteine', 'S-4-methoxytrityl'),
#  'cys(npys)': ('cysteine', 'S-3-nitro-2-pyridinesulfenyl'),
#  'asp(tbu)': ('aspartic acid', 'β-tert-butyl ester'),
#  'asp(otbu)': ('aspartic acid', 'β-tert-butyl ester'),
#  'asp(obut)': ('aspartic acid', 'β-tert-butyl ester'), #typo?
#  'asp(bn)': ('aspartic acid', 'β-benzyl ester'),
#  'asp(obn)': ('aspartic acid', 'β-benzyl ester'),
#  'asp(ome)': ('aspartic acid', 'β-methyl ester'),
#  'asp(oet)': ('aspartic acid', 'β-ethyl ester'),
#  'asp(oall)': ('aspartic acid', 'β-allyl ester'),
#  'asp(ofm)': ('aspartic acid', 'β-9-fluorenylmethyl ester'),
#  'glu(tbu)': ('glutamic acid', 'γ-tert-butyl ester'),
#  'glu(otbu)': ('glutamic acid', 'γ-tert-butyl ester'),
#  'glu(obut)': ('glutamic acid', 'γ-tert-butyl ester'), #typo?
#  'glu(bn)': ('glutamic acid', 'γ-benzyl ester'),
#  'glu(obn)': ('glutamic acid', 'γ-benzyl ester'),
#  'glu(ome)': ('glutamic acid', 'γ-methyl ester'),
#  'glu(oet)': ('glutamic acid', 'γ-ethyl ester'),
#  'glu(oall)': ('glutamic acid', 'γ-allyl ester'),
#  'glu(ofm)': ('glutamic acid', 'γ-9-fluorenylmethyl ester'),
#  'asn(trt)': ('asparagine', 'Nδ-(trityl)'),
#  'asn(xan)': ('asparagine', 'Nδ-(xanthenyl)'),
#  'gln(trt)': ('glutamine', 'Nε-(trityl)'),
#  'gln(xan)': ('glutamine', 'Nε-(xanthenyl)'),
#  'trp(boc)': ('tryptophan', 'N¹-tert-butoxycarbonyl'),
#  'trp(for)': ('tryptophan', 'N¹-formyl'),
#  'trp(hoc)': ('tryptophan', 'N¹-1-(1-adamantyl)-1-methylethoxycarbonyl'),
#  'met(o)': ('methionine', 'S-oxide'),
#  'met(o2)': ('methionine', 'S,S-dioxide'),
#  'orn(boc)': ('ornithine', 'Nδ-(tert-butoxycarbonyl)'),
#  'orn(z)': ('ornithine', 'Nδ-(benzyloxycarbonyl)'),
#  'orn(fmoc)': ('ornithine', 'Nδ-(9-fluorenylmethoxycarbonyl)'),
#  'hse(tbu)': ('homoserine', 'O-tert-butyl'),
#  'hse(bn)': ('homoserine', 'O-benzyl')}

COUNTER_ACIDS = {'hcl': 'hydrochloride',
                 '2hcl': 'dihydrochloride',
                 '3hcl': 'trihydrochloride',
                 '4hcl': 'tetrahydrochloride',
                 'hydrochloride': 'hydrochloride',
                 'hbr': 'hydrobromide',
                 '2hbr': 'dihydrobromide',
                 'hci': 'hydroiodide',
                 'tfa': 'trifluoroacetate',
                 'trifluoroacetate': 'trifluoroacetate',
                 'acetate': 'acetate',
                 'diacetate': 'diacetate',
                 'acoh': 'acetate',
                 '2acoh': 'diacetate'}

GREEK_LETTERS = ['α', 'β', 'γ', 'δ', 'ε', 'ζ', 'η', 'θ', 'ι', 'κ', 'λ', 'μ', 'ν', 'ξ', 'ο', 'π', 'ρ', 'σ', 'τ', 'υ', 'φ', 'χ', 'ψ', 'ω']

def split_peptide_shorthand(shorthand):
    tokens = []
    current_token = ""
    paren_depth = 0
    
    for char in shorthand.strip():
        if char == '(':
            paren_depth += 1
            current_token += char
        elif char == ')':
            paren_depth -= 1
            current_token += char
        elif char == '-' and paren_depth == 0:
            # Split here - we're not inside parentheses
            if current_token:
                tokens.append(current_token)
                current_token = ""
        else:
            current_token += char
    
    if current_token:
        tokens.append(current_token)
    
    return tokens

def parse_protected_residue(token):
    """Parse a token that might contain side chain protection"""
    if '(' in token and token.endswith(')'):
        if token.lower() in SIDE_CHAIN_PROTECTIONS:
            base_aa, protection = SIDE_CHAIN_PROTECTIONS[token.lower()]
            return base_aa, protection
        else:
            # Try to parse manually
            base = token.split('(')[0]
            protection_part = token[token.find('(')+1:-1]
            if base in AA_FULL:
                return AA_FULL[base], f"with {protection_part} protection"
    return None, None

def shorthand_to_iupac(shorthand: str) -> str:
    original_shorthand = shorthand.strip()
    original_shorthand = shorthand.strip('-')
    original_shorthand = original_shorthand.replace('--', '-')
    
    # Check if it's a cyclic peptide
    is_cyclic = False
    if original_shorthand.lower().startswith('cyclo(') and original_shorthand.endswith(')'):
        is_cyclic = True
        shorthand = original_shorthand[6:-1].strip()
    elif original_shorthand.lower().startswith('cyclo[') and original_shorthand.endswith(']'):
        is_cyclic = True
        shorthand = original_shorthand[6:-1].strip()
    elif original_shorthand.lower().startswith('cyclo-(') and original_shorthand.endswith(')'):
        is_cyclic = True
        shorthand = original_shorthand[7:-1].strip() 
    elif original_shorthand.lower().startswith('cyclo-[') and original_shorthand.endswith(']'):
        is_cyclic = True
        shorthand = original_shorthand[7:-1].strip()
    elif original_shorthand.lower().startswith('cyclo-{') and original_shorthand.endswith('}'):
        is_cyclic = True
        shorthand = original_shorthand[7:-1].strip()
    elif original_shorthand.lower().startswith('cyclo-{') and original_shorthand.endswith('}'):
        is_cyclic = True
        shorthand = original_shorthand[7:-1].strip() 
    elif original_shorthand.lower().startswith('cyclo'):
        is_cyclic = True
        shorthand = original_shorthand[5:].strip()

    shorthand = shorthand.strip()
    shorthand = shorthand.strip('-')
    shorthand = shorthand.replace('--', '-')

    items_to_remove = []
    counter_acid_suffix = None
    if len(shorthand.split('.')) == 2:
        for item in shorthand.split('.'):
            if item.lower().strip().replace(' ', '') in COUNTER_ACIDS:
                counter_acid_suffix = COUNTER_ACIDS[item.lower().strip().replace(' ', '')]
                items_to_remove.append(item)
        if len(items_to_remove) == 1:
            new_shorthand_split = shorthand.split('.')
            new_shorthand_split.remove(items_to_remove[0])
            shorthand = '.'.join(new_shorthand_split)
    
    tokens = split_peptide_shorthand(shorthand.strip())

    # N-cap
    prefix = ''
    if tokens and tokens[0].lower() in N_CAPS:
        prefix = N_CAPS[tokens.pop(0).lower()]

    # C-cap  
    suffix = ''
    if tokens and tokens[-1].lower() in C_CAPS:
        suffix = C_CAPS[tokens.pop().lower()]

    if not tokens:
        raise ValueError("No residues found")

    parts = []
    i = 0
    while i < len(tokens):
        t = tokens[i]
        t_lower_stripped = t.lower().strip().replace(' ', '')
        is_last = (i == len(tokens) - 1)
        
        # For cyclic peptides, all residues should be in -yl form except we need special handling
        # since there's no true "last" residue in a cycle
        if is_cyclic:
            is_last = False  # Treat all as intermediate residues

        # Check if this token is a stereochemistry indicator
        if t_lower_stripped in ['d', '(d)']:
            # Next token should be processed with D- prefix
            i += 1
            if i >= len(tokens):
                # 'd' was the last token, treat as unknown
                parts.append('d')
                break
            
            next_token = tokens[i]
            next_token_lower_stripped = next_token.lower().strip().replace(' ', '')
            greek_prefix=None
            if next_token[0] in GREEK_LETTERS:
                greek_prefix = next_token[0]
                next_token = next_token[1:]

            methyl_modifier=None
            if next_token != '' and next_token_lower_stripped.startswith('me'):
                if next_token_lower_stripped[2:] in AA_FULL:
                    methyl_modifier = 'methyl-'    
                    next_token = next_token[2:]
            if next_token != '' and next_token_lower_stripped.startswith('(me)'):
                if next_token_lower_stripped[4:] in AA_FULL:
                    methyl_modifier = 'methyl-'    
                    next_token = next_token[4:]

            is_last_updated = is_last and not is_cyclic
            
            # Try to parse the next token
            base_aa, protection = parse_protected_residue(next_token_lower_stripped)
            if not base_aa and next_token_lower_stripped in AA_FULL:
                base_aa = AA_FULL[next_token_lower_stripped]
            
            if base_aa:
                # Known amino acid
                if is_last_updated:
                    base = base_aa
                else:
                    # Convert to -yl form
                    base_key = None
                    for key, val in AA_FULL.items():
                        if val == base_aa:
                            base_key = key
                            break
                    if base_key and base_key.lower().strip().replace(' ', '') in AA_YL:
                        base = AA_YL[base_key.lower().strip().replace(' ', '')]
                    else:
                        base = base_aa.replace('ine', 'yl').replace('ic acid', 'yl')
                
                # Add D- prefix (not for glycine)
                if not base_aa.startswith('glycine'):
                    base = 'd-' + base

                if greek_prefix:
                    base = greek_prefix + '-' + base

                if methyl_modifier:
                    base = methyl_modifier + base
                    
                if protection:
                    base = f"{protection}-{base}"
                parts.append(base)
            else:
                # Unknown token after 'd'
                parts.append(f"d-{next_token}")
                
        elif t_lower_stripped in ['l', '(l)']:
            # Next token should be processed with L- prefix  
            i += 1
            if i >= len(tokens):
                # 'l' was the last token, treat as unknown
                parts.append('l')
                break
                
            next_token = tokens[i]
            next_token_lower_stripped = next_token.lower().strip().replace(' ', '')
            greek_prefix=None
            if next_token[0] in GREEK_LETTERS:
                greek_prefix = next_token[0]
                next_token = next_token[1:]

            methyl_modifier=None
            if next_token != '' and next_token_lower_stripped.startswith('me'):
                if next_token_lower_stripped[2:] in AA_FULL:
                    methyl_modifier = 'methyl-'    
                    next_token = next_token[2:]
            if next_token != '' and next_token_lower_stripped.startswith('(me)'):
                if next_token_lower_stripped[4:] in AA_FULL:
                    methyl_modifier = 'methyl-'    
                    next_token = next_token[4:]

            is_last_updated = is_last and not is_cyclic
            
            # Try to parse the next token
            base_aa, protection = parse_protected_residue(next_token_lower_stripped)
            if not base_aa and next_token_lower_stripped in AA_FULL:
                base_aa = AA_FULL[next_token_lower_stripped]
            
            if base_aa:
                # Known amino acid
                if is_last_updated:
                    base = base_aa
                else:
                    # Convert to -yl form
                    base_key = None
                    for key, val in AA_FULL.items():
                        if val == base_aa:
                            base_key = key
                            break
                    if base_key and base_key.lower().strip().replace(' ', '') in AA_YL:
                        base = AA_YL[base_key.lower().strip().replace(' ', '')]
                    else:
                        base = base_aa.replace('ine', 'yl').replace('ic acid', 'yl')
                
                # Add L- prefix (not for glycine)
                if not base_aa.startswith('glycine'):
                    base = 'l-' + base

                if greek_prefix:
                    base = greek_prefix + '-' + base

                if methyl_modifier:
                    base = methyl_modifier + base
                    
                if protection:
                    base = f"{protection}-{base}"
                parts.append(base)
            else:
                # Unknown token after 'l'
                parts.append(f"l-{next_token}")

        elif t_lower_stripped in ['dl', '(dl)', '(d/l)', 'd/l', 'd,l', '(d,l)']: 
            # Next token should be processed with L- prefix  
            i += 1
            if i >= len(tokens):
                # 'dl' was the last token, treat as unknown
                parts.append('d/l')
                break
                
            next_token = tokens[i]
            next_token_lower_stripped = next_token.lower().strip().replace(' ', '')
            greek_prefix=None
            if next_token[0] in GREEK_LETTERS:
                greek_prefix = next_token[0]
                next_token = next_token[1:]

            methyl_modifier=None
            if next_token != '' and next_token_lower_stripped.startswith('me'):
                if next_token_lower_stripped[2:] in AA_FULL:
                    methyl_modifier = 'methyl-'    
                    next_token = next_token[2:]
            if next_token != '' and next_token_lower_stripped.startswith('(me)'):
                if next_token_lower_stripped[4:] in AA_FULL:
                    methyl_modifier = 'methyl-'    
                    next_token = next_token[4:]

            is_last_updated = is_last and not is_cyclic
            
            # Try to parse the next token
            base_aa, protection = parse_protected_residue(next_token_lower_stripped)
            if not base_aa and next_token_lower_stripped in AA_FULL:
                base_aa = AA_FULL[next_token_lower_stripped]
            
            if base_aa:
                # Known amino acid
                if is_last_updated:
                    base = base_aa
                else:
                    # Convert to -yl form
                    base_key = None
                    for key, val in AA_FULL.items():
                        if val == base_aa:
                            base_key = key
                            break
                    if base_key and base_key.lower().strip().replace(' ', '') in AA_YL:
                        base = AA_YL[base_key.lower().strip().replace(' ', '')]
                    else:
                        base = base_aa.replace('ine', 'yl').replace('ic acid', 'yl')
                
                # Add L- prefix (not for glycine)
                if not base_aa.startswith('glycine'):
                    base = 'dl-' + base

                if greek_prefix:
                    base = greek_prefix + '-' + base

                if methyl_modifier:
                    base = methyl_modifier + base
                    
                if protection:
                    base = f"{protection}-{base}"
                parts.append(base)
            else:
                # Unknown token after 'l'
                parts.append(f"dl-{next_token}")
                
        else:
            # Regular token (not preceded by d or l)
            if t_lower_stripped.startswith('(d)'):
                base_prefix = 'd-'
                t = t[3:]
            elif t_lower_stripped.startswith('(l)'):
                base_prefix = 'l-'
                t = t[3:]
            elif t_lower_stripped.startswith('(d/l)'):
                base_prefix = '(d/l)-'
                t = t[5:]
            else:
               base_prefix = 'l-' 

            greek_prefix=None
            if t != '' and t[0] in GREEK_LETTERS:
                greek_prefix = t[0]
                t = t[1:]

            methyl_modifier=None
            if t != '' and t_lower_stripped.startswith('me'):
                if t_lower_stripped[2:] in AA_FULL:
                    methyl_modifier = 'methyl-'    
                    t = t[2:]
            if t != '' and t_lower_stripped.startswith('(me)'):
                if t_lower_stripped[4:] in AA_FULL:
                    methyl_modifier = 'methyl-'    
                    t = t[4:]
                
            
            base_aa, protection = parse_protected_residue(t_lower_stripped)
            if not base_aa and t_lower_stripped in AA_FULL:
                base_aa = AA_FULL[t_lower_stripped]

            if base_aa:
                # Known amino acid
                if is_last and not is_cyclic:
                    base = base_aa
                else:
                    # Convert to -yl form
                    base_key = None
                    for key, val in AA_FULL.items():
                        if val == base_aa:
                            base_key = key
                            break
                    if base_key and base_key.lower().strip().replace(' ', '') in AA_YL:
                        base = AA_YL[base_key.lower().strip().replace(' ', '')]
                    else:
                        base = base_aa.replace('ine', 'yl').replace('ic acid', 'yl')
                
                # Add L- prefix by default (not for glycine)
                if not base_aa.startswith('glycine'):
                    base = base_prefix + base
                    
                if greek_prefix:
                    base = greek_prefix + '-' + base

                if methyl_modifier:
                    base = methyl_modifier + base
                    
                if protection:
                    base = f"{protection}-{base}"
                parts.append(base)
            else:
                # Unknown token - keep as-is without adding stereochemistry
                parts.append(t)
        
        i += 1

    name = ''
    if is_cyclic:
        name += 'cyclo('
    
    if prefix:
        name += prefix + '-'
    name += '-'.join(parts)
    if suffix:
        name += ' ' + suffix
        
    if is_cyclic:
        name += ')'

    if counter_acid_suffix:
        name += ' ' + counter_acid_suffix
        
    return name


def opsin_smiles(iupac_name: str) -> str:
    url = f"https://opsin.ch.cam.ac.uk/opsin/{urllib.parse.quote(iupac_name)}.smi"
    r = requests.get(url, timeout=15)
    r.raise_for_status()
    return r.text.strip()

# Demo with your problematic examples
examples = [
    # "Fmoc-Cys-OAllyl",
    # "Fmoc-Lys(Proc)-OH",
    # "Cbz-Phe-Gly-OMe",
    # "H2N-Phe-Phe-Gly-Thr-Phe-Phe-Gly-OH",
    # "Arg(Pbf)-OH",
    # "Fmoc-Ser(tBu)-OH",
    # "Ac-Cys(Trt)-NH2",
    # "Boc-Phe-Gly-OMe",
    # "Fmoc-His(Trt)-OH",
    # "Fmoc-Lys(Proc)",
    # "(D)Leu-OBn", 
    # "boc-gly-l-ala-l-phe-l-ileu-gly-l-leu-l-met-nh",
    # "cyclo(arg-gly-asp-d-phe-lys)",
    # "cyclo-(l-meleu-d-vala-l-meleu-d-lac)",
    # "l-histidyl-d-βNal-l-alanine methyl ester",
    # "H-Tyr-D-Ala-Gly-Phe-NH-NH-D-Phe-D-Asp-D-Nle-Trp-H",
    # "H-Phe-His-Leuc -Val-Ile-His-NH2.2HCl",
    # "pglu-asn-trp",
    "Cbz-DL-Lys(Cbz)-DL-Val-DL-Pro-DL-Val-CF3",
    "Cbz-L-Phe-L-Lys(Cbz)-OMe",
    "H2N-Phe-Phe-Gly-Thr-Phe-Phe-Gly-OH",
    "t-Boc-Lys(CBZ)-Gln-Nle-Ala-Val-Lys(CBZ)-NH2",
]

for ex in examples:
    try:
        iupac = shorthand_to_iupac(ex)
        print(f"{ex} -> {iupac}")
        # Uncomment to test SMILES generation
        # smiles = opsin_smiles(iupac)
        # print(f"  -> {smiles}")
    except Exception as e:
        print(f"{ex} -> ERROR: {e}")
