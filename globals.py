import os

# GLOBAL VARIABLES
only_print_command = False
use_cluster = False #!
leave_one_out_xval = True #! False #
score_with_all_seeds = True #False #!
only_auc = False # In the analysis_xval if True auc.txt created only, other graphs are not drawn
navlakha_data = True #False #!

DEFAULT_TOP_SCORING_CUTOFF = "5%" #"10%" #"1%" #"5%" # If ends with "%" taken as percentage otherwise as score - At the time of analysis it was "10%"
N_LINKER_THRESHOLD = 2 # For Netlink method
DEFAULT_SEED_SCORE = 1.0 # Default score for seed nodes, used when no score given in assoc file
DEFAULT_NON_SEED_SCORE = 0.01 # Default score for non-seed nodes
ALLOWED_MAX_DEGREE = 100000 #175 #90 # Max degree allowed in the graph filtering
N_SAMPLE_GRAPH = 100 # Number of random graphs to be generated
N_X_VAL = None #! #5 #182 # Number of cross validation folds, readjusted if leave_one_out_xval is True
N_SEED = None #Will be set during run
N_RANDOM_NEGATIVE_FOLDS = None #10 # Number of non-seed scores to be averaged for negative score calculation, 
			    # If 0 all non seeds are included as they are, If None all non seeds are included averaged to scale with the size of test seeds
REPLICABLE = 123 #63826 #9871354 #123 #None # Assign a predefined seed in randomization for initial test folds creation and N_RANDOM_NEGATIVE_FOLD generation
ONLY_LARGEST_COMPONENT = True 
GO_ENRICHMENT_P_VALUE_CUTOFF = 0.05

# FOLLOWING LOCAL ONLY TO MAIN
ignore_experiment_failures = False
delay_experiment = True
tex_format = True #False 
functional_enrichment = False

MODE = "compare" # prepare, score, analyze, compare, summary, module
user_friendly_id = "navlakha" #"biana_no_tap_relevance-omim_alzheimer-ns_nz_nd_union-top5" # "biana_no_tap-omim" #biana_no_tap_no_reliability-omim_alzheimer-nd-top5 #"biana_no_tap-omim_perturbed_%d_10" % i_parameter #"biana_no_tap-all_seed20below" #"biana-all" #"biana_no_tap-omim" #"omim_alzheimer-diabetes" #"all_vs_all" # a.k.a. emre friendly id for compare and summary
summary_seed_cutoff = 1 #None #2 #20 # Seed cutoff considered for inclusion of an experiment in sum_up_experiments, if None seed.dat is not created. Also used in compare_experiments if analysis_type is user
prepare_mutated = None # Creates permutad/pruned networks 
analyze_network = False 
exclude_seeds_in_comparison = True


omim_phenotypes = ["alzheimer", "breast cancer", "diabetes", "insulin", "anemia", "myopathy", "neuropathy", "obesity", "parkinson disease", "prostate cancer", "hypertension", "leukemia", "lung cancer", "asthma", "ataxia", "epilepsy", "schizophrenia", "cardiomyopathy", "cataract", "spastic paraplegia", "lymphoma", "mental retardation", "systemic lupus erythematosus"] # "autism", "aneurysm",  
omim_phenotypes = [ "omim_" + "_".join(p.split()) for p in omim_phenotypes ]

goh_phenotypes = ["goh_developmental", "goh_connective_tissue", "goh_ear,nose,throat", "goh_endocrine", "goh_psychiatric", "goh_immunological", "goh_neurological", "goh_respiratory", "goh_multiple", "goh_renal", "goh_skeletal", "goh_bone", "goh_dermatological", "goh_cancer", "goh_ophthamological", "goh_metabolic", "goh_nutritional", "goh_muscular", "goh_hematological", "goh_gastrointestinal", "goh_cardiovascular"] #,"goh_unclassified"] 

chen_phenotypes = ["atherosclerosis",  "ischaemic_stroke",  "systemic_scleroderma",  "migraine",  "epilepsy",  "cirrhosis",  "ulcerative_colitis",  "cervical_carcinoma",  "osteoarthritis",  "inflammatory_bowel_disease",  "myocardial_ischemia",  "endometrial_carcinoma",  "pancreatitis",  "graves_disease",  "neural_tube_defects",  "lymphoma",  "endometriosis",  "autism",  "hypercholesterolaemia"]
chen_phenotypes = [ "chen_" + p for p in chen_phenotypes ]

hsdl_phenotypes = [ "INNER_CELL_MASS", "TROPHECTODERM", "Embryonic_stem_cell", "ECTODERM", "MESENDODERM", "ENDODERM", "MESODERM", "neural_progenitor_cell", "neural_stem_cell", "heart", "cardiomyocyte", "digestive_tube", "pancreas", "liver" ]
hsdl_phenotypes = [ "hsdl_" + p.replace(" ", "_").lower() for p in hsdl_phenotypes ]


ppis = []
ppis += ["hprd"] #, "ophid"]
#ppis += ["goh", "entrez", "biana_no_tap_no_reliability", "biana_no_tap_relevance", "biana_no_reliability"] 
#ppis += ["rivasi"]
#ppis += ["biana_no_tap_no_reliability", "biana_no_tap_relevance", "biana_no_reliability"] 
#ppis += ["goh"]
#ppis += ["entrez"]
#ppis += ["biana_no_reliability"]
#ppis += ["biana_no_tap_no_reliability"] 
#ppis += ["biana_no_tap_relevance"]
#ppis += ["biana_no_tap_no_reliability_permuted_p%s_%s" % (p, i) for p in xrange(10,110,10) for i in xrange(1,101) ] 
#ppis += ["biana_no_tap_no_reliability_permuted_p%s_%s" % (p, i) for p in xrange(10,110,10) for i in xrange(1,11) ] 
#ppis += ["biana_no_tap_no_reliability_pruned_p%s_%s" % (p, i) for p in xrange(i_parameter,i_parameter+10,10) for i in xrange(1,11) ] 
#ppis += ["biana_no_tap_no_reliability_pruned_non_seed_interactions_p%s_%s" % (p, i) for p in xrange(10,100,10) for i in xrange(1,101) ] 
#ppis += ["biana_no_tap_no_reliability_pruned_p%s_%s" % (p, i) for p in xrange(10,100,10) for i in xrange(1,11) ] 
#ppis += ["david"]
#ppis += ["javi"] #["goh"] #["piana_joan_exp", "piana_joan_all"] #["david"] #["goh", "biana_no_tap_no_reliability", "biana_no_reliability", "biana_no_tap_relevance"]
#ppis = ["ori_coexpression_1e-2", "ori_network", "ori_coexpression", "ori_coexpression_colocalization", "ori_colocalization", "ori_coexpression_colocalization_1e-2"]
#ppis += ["ori_no_tap_coexpression_1e-2", "ori_no_tap_network", "ori_no_tap_coexpression", "ori_no_tap_coexpression_colocalization", "ori_no_tap_colocalization", "ori_no_tap_coexpression_colocalization_1e-2"]
#ppi += ["goh_1e5", "biana_coexpression"]

phenotypes = []
#phenotypes += navlakha_phenotypes # Now located at the bottom of the page
#phenotypes += chen_phenotypes + omim_phenotypes + goh_phenotypes 
#phenotypes += hsdl_phenotypes
#phenotypes += omim_phenotypes 
#phenotypes += [ "perturbed_%s_p%i_%i" % (d, p, i) for d in omim_phenotypes for p in xrange(10,100,10) for i in xrange(1,101) ]
#phenotypes += [ "perturbed_%s_p%i_%i" % (d, p, i) for d in omim_phenotypes for p in xrange(i_parameter,i_parameter+10,10) for i in xrange(1,11) ]  
#phenotypes += goh_phenotypes 
#phenotypes += chen_phenotypes 
#phenotypes += ["omim_prostate_cancer"]
#phenotypes += ["omim_breast_cancer", "omim_lung_cancer"]
#phenotypes += ["omim_leukemia"]
#phenotypes += ["omim_alzheimer"] 
#phenotypes += ["omim_insulin"] 
#phenotypes += ["omim_diabetes"]
#phenotypes += ["omim_parkinson_disease"]
#phenotypes += ["apoptosis_joan"]
#phenotypes += ["custom"] #["aneurysm"] #["apoptosis_joan"] #["alzheimer_david_CpOGU", "alzheimer_david_CpOIN", "alzheimer_david_RpOGU", "alzheimer_david_RpOIN"] #["aneurysm", "breast_cancer"]

scoring_parameters = []
scoring_parameters += [("nr", 1, 1), ("ff", 1, 5)]
scoring_parameters += [("nz", 1, 5), ("ns", 3, 2)] 
#scoring_parameters += [("nr", 1, 1)]
scoring_parameters += [("nd", 1, 1)]
scoring_parameters += [("mcl", 1, 1)]
#scoring_parameters += [("nz", 1, 5)]
#scoring_parameters += [("ns", 3, 2)]
#scoring_parameters += [("ff", 1, 5)]
#scoring_parameters += [("ns", 2, 3), ("ns", 2, 2)]
#scoring_parameters += [("nw",1, 1)]
#scoring_parameters += [("nx", 1, 1)]
#scoring_parameters += [("ns", 2, 2), ("ns", 2, 3), ("ns", 2, 4), ("ns", 3, 3)]
#scoring_parameters += [("ff", 1, i) for i in xrange(1,9)]
#scoring_parameters += [("nz", 1, i) for i in xrange(4,6)]
##scoring_parameters += [("nz", 1, i) for i in xrange(1,9)]
##scoring_parameters += [("ns", r, i) for r in xrange(1,9) for i in xrange(1,5)]
##scoring_parameters += [("ns", r, i) for r in xrange(4,9) for i in xrange(1,3)]
##scoring_parameters += [("nh", r, i) for r in (1,2,3) for i in xrange(1,5)]
##scoring_parameters += [("n1", r, i) for r in (1,2,3) for i in xrange(1,5)]


# Directory of the project
base_dir = ".."
base_dir = os.path.abspath(base_dir) + os.sep
data_dir = base_dir + "data" + os.sep
data_dir = os.path.abspath(data_dir) + os.sep
src_dir = base_dir + "src" + os.sep

# BIANA node & network files
biana_node_file_prefix = data_dir + "human_interactome_biana" + os.sep + "human_nodes"
biana_network_file_prefix = data_dir + "human_interactome_biana" + os.sep + "human_network"

# PPIs from existing studies
goh_network_file = data_dir + "goh_human_ppi" + os.sep + "ppi.sif"
rhodes_network_file = data_dir + "rhodes_human_probabilistic_ppi" + os.sep + "ppi.sif"
goh_network_file_filtered_by_degree = goh_network_file[:-4] + "_degree_filtered.sif"
rhodes_network_file_filtered_by_degree = rhodes_network_file[:-4] + "_degree_filtered.sif"
rhodes_interaction_relevance_file = rhodes_network_file[:-4] + ".eda"

# Gene info file 
gene_info_file = data_dir + "gene_info" + os.sep + "genes.tsv"

COMPARISON_GOLD_STANDARD_FILE = data_dir + "alzheimer_gold" + os.sep + "gene_list.txt" 
#COMPARISON_GOLD_STANDARD_FILE = data_dir + "netage" + os.sep + "AD_genes.txt" 
#COMPARISON_GOLD_STANDARD_FILE = data_dir + "uwaging" + os.sep + "aging.txt" # 22, 9 from genage
#COMPARISON_GOLD_STANDARD_FILE = data_dir + "genage" + os.sep + "aging_candidates.txt" # 261, 9 from uwaging, 91 from netage
#COMPARISON_GOLD_STANDARD_FILE = data_dir + "uwaging" + os.sep + "uwaging_genage_netage.txt" # 5 intersection of uwaging - genage - netage 
#COMPARISON_GOLD_STANDARD_FILE = data_dir + "uwaging" + os.sep + "uwaging_genage.txt" # 9 intersection of uwaging - genage 
#COMPARISON_GOLD_STANDARD_FILE = data_dir + "uwaging" + os.sep + "uwaging_mutex_genage_netage.txt" # 12 intersection of uwaging - genage & uwagin - netage
#COMPARISON_GOLD_STANDARD_FILE = data_dir + "uwaging" + os.sep + "mutex_uwaging_genage_netage.txt" # 99 intersection of uwaging - genage & uwaging - netage & genage - netage
#COMPARISON_GOLD_STANDARD_FILE = data_dir + "netage" + os.sep + "longetivity.txt" # 456, 8 from uwaging, 91 from genage

scoring_methods = ["nd", "nz", "ns", "ff", "nr", "nw", "nl", "nx", "nh", "n1", "nb"]

THRESHOLDS = { "nr": [ 4e-6, 2e-5, 5e-5, 1e-4, 2e-4, 3e-4, 4e-4, 5e-4, 1e-3, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 ],
		"ff": [ 1e-3, 1e-2, 2e-2, 5e-2, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3 ], 
		"nd": [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9 ],  
		"nz": [ 0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.018, 0.02, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9 ],
		"ns": [ 0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.018, 0.02, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9 ] }

THRESHOLDS = { "nr": [ j*10**-i for i in xrange(2,7) for j in xrange(1,10) ] + [ 0.05*i for i in xrange(1,21) ],
		"ff": [ j*10**-i for i in xrange(2,5) for j in xrange(1,10) ] + [ 0.05*i for i in xrange(1,31) ], 
		"nd": [ 0.01*i for i in xrange(1,101) ],  
		"nz": [ 0.01*i for i in xrange(1,101) ],
		"ns": [ 0.01*i for i in xrange(1,101) ],
		"mcl": [ 0.01*i for i in xrange(1,101) ] } 


# Scoring related parameters 
#PPI = "biana" # biana output network as it is (do not forget to revert _degree_filtered.sif.original to _degree_filtered.sif, this one is _degree_filtere_disconnected_only.sif)
#PPI = "biana_no_reliability" # only largest connected component (lcc)
#PPI = "biana_reliability" # reliability filtered lcc
#PPI = "biana_no_tap_no_reliability" # tap filtered lcc
#PPI = "biana_no_tap_no_reliability_1e-5" # tap filtered lcc with non seed scores of 1e-5
#PPI = "biana_no_tap_reliability" # tap & reliability filtered lcc
#PPI = "biana_no_tap_relevance" # tap filtered & string edge score assigned lcc # manually assigned +1 to all edge scores to reduce max/min edge score ratio
#PPI = "biana_no_tap_exp_db_relevance" tap filtered & string exp & db edge score assigned lcc
#PPI = "biana_no_tap_corelevance" # tap filtered & string co-exp score assigned lcc
#PPI = "biana_no_tap_reliability_relevance" # tap & reliability filtered & string edge score assigned lcc
#PPI = "goh" 
#PPI = "rhodes"

#SCORING = "ns" #"netscore"
#SCORING = "nz" #"netzcore"
#SCORING = "nh" #"netzscore"
#SCORING = "n1" #"netz1score"
#SCORING = "nd" #"netshort"
#SCORING = "nw" #"netween"
#SCORING = "nl" #"netlink"
#SCORING = "nr" #"netrank"
#SCORING = "nx" #"netrandom"
#SCORING = "nb" #"netZscore" (cortesy of baldo)
#SCORING = "ff" #"FunctionalFlow" 


navlakha_phenotypes = ['abacavir', 'abdominal', 'acampomelic', 'acromesomelic', 'acth', 'adenocarcinoma', 'adenomas', 'adrenal', 'adrenocortical', 'adrenoleukodystrophy', 'adrenomyeloneuropathy', 'afibrinogenemia', 'agammaglobulinemia', 'aids', 'alagille', 'albinism', 'alcohol', 'aldosterone', 'alexander', 'alport', 'alzheimer', 'amelogenesis', 'amyloidosis', 'amyotrophic', 'anemia', 'angelman', 'angioedema', 'anorexia', 'anterior', 'antley-bixler', 'aortic', 'aplastic', 'apolipoprotein', 'arrhythmogenic', 'arthrogryposis', 'asthma', 'ataxia', 'atelosteogenesis', 'atherosclerosis', 'atopy', 'atrial', 'atrioventricular', 'attention', 'autism', 'autoimmune', 'azoospermia', 'bamforth-lazarus', 'bardet-biedl', 'bare', 'bartter', 'basal', 'bcg', 'becker', 'beckwith-wiedemann', 'bernard-soulier', 'bethlem', 'bladder', 'bleeding', 'blepharophimosis', 'blood', 'blue-cone', 'boomerang', 'brachydactyly', 'bradyopsia', 'brain', 'branchiootorenal', 'breast', 'brugada', 'budd-chiari', 'butterfly', 'c1q', 'c1r,c1s', 'c4', 'c8', 'campomelic', 'cardiomyopathy', 'carnitine', 'cataract', 'caudal', 'celiac', 'central', 'cerebellar', 'cerebral', 'cerebrooculofacioskeletal', 'charcot-marie-tooth', 'cholestasis', 'chondrodysplasia', 'chondrosarcoma', 'chorea', 'choreoathetosis', 'chromosome', 'chronic', 'cirrhosis', 'cleft', 'cockayne', 'coenzyme', 'cold-induced', 'coloboma', 'colon', 'colorblindness', 'colorectal', 'combined', 'complement', 'complex', 'cone', 'cone-rod', 'congenital', 'corneal', 'cornelia', 'coronary', 'corpus', 'craniofacial', 'craniosynostosis', 'creatine', 'creutzfeldt-jakob', 'crohn', 'crouzon', 'cutis', 'cystic', 'deafness', 'dejerine-sottas', 'dementia', 'dent', 'dermatitis', 'diabetes', 'dna', 'dysfibrinogenemia', 'dyskeratosis', 'dyslexia', 'dystonia', 'ectodermal', 'ectopia', 'ehlers-danlos', 'elliptocytosis', 'emery-dreifuss', 'emphysema', 'encephalopathy', 'endometrial', 'enolase', 'epidermolysis', 'epidermolytic', 'epilepsy', 'epileptic', 'epiphyseal', 'episodic', 'erythremia', 'erythrocytosis', 'esophageal', 'exostoses', 'exudative', 'factor', 'fanconi', 'fetal', 'fibromatosis', 'fletcher', 'focal', 'foveomacular', 'fructose', 'fundus', 'gastric', 'gastrointestinal', 'gaucher', 'generalized', 'germ', 'giant', 'glanzmann', 'glaucoma', 'glioblastoma', 'glomerulosclerosis', 'glycine', 'glycogen', 'glycogenosis', 'gm2-gangliosidosis', 'goiter', 'gonadal', 'graves', 'griscelli', 'growth', 'h.', 'hdl', 'heinz', 'hemangioma', 'hematuria', 'hemochromatosis', 'hemolytic', 'hemolytic-uremic', 'hemophagocytic', 'hemophilia', 'hemorrhagic', 'hepatic', 'hepatitis', 'hepatoblastoma', 'hepatocellular', 'hereditary', 'hermansky', 'heterotaxy', 'high', 'von_hippel-lindau', 'hirschsprung', 'histiocytoma', 'hiv', 'holoprosencephaly', 'homocystinuria', 'huntington', 'hypercholanemia', 'hypercholesterolemia', 'hyperekplexia', 'hyperinsulin', 'hyperlipoproteinemia', 'hyperoxaluria', 'hyperparathyroidism', 'hypertension', 'hyperthyroidism', 'hypertriglyceridemia', 'hypodontia', 'hypogonadotropic', 'hypokalemic', 'hypomagnesemia', 'hypoparathyroidism', 'hypophosphat', 'hypothyroidism', 'hypotrichosis', 'ichthyosiform', 'ichthyosis', 'iga', 'immunodeficiency', 'inclusion', 'inflammatory', 'insomnia', 'insulin', 'intervertebral', 'intracranial', 'intrauterine', 'invasive', 'iridogoniodysgenesis', 'iron', 'jackson-weiss', 'jervell', 'joubert', 'juvenile', 'kallmann', 'kaposi', 'keratitis', 'keratosis', 'ladd', 'larsen', 'leber', 'leigh', 'leiomyomatosis', 'leopard', 'leprosy', 'lethal', 'leukemia', 'leuko', 'li', 'liddle', 'lipodystrophy', 'lipoid', 'lipoma', 'lipoprotein', 'lissencephaly', 'loeys-dietz', 'long', 'longevity', 'lumbar', 'lung', 'lupus', 'lymphangioleiomyomatosis', 'lymphoma', 'lymphoproliferative', 'macrothrombocytopenia', 'macular', 'major', 'malaria', 'male', 'malignant', 'marfan', 'mast', 'maturity-onset', 'meckel', 'medullary', 'medulloblastoma', 'megakaryoblastic', 'megaloblastic', 'melanoma', 'memory', 'meningioma', 'mental', 'metachromatic', 'metaphyseal', 'methemoglobinemia', 'microcephaly', 'microphthalmia', 'migraine', 'mismatch', 'mitochondrial', 'mody', 'mowat-wilson', 'mucoepidermoid', 'mucopolysaccharidosis', 'muir-torre', 'multiple', 'muscle', 'muscular', 'myasthenic', 'mycobacterial', 'mycobacterium', 'myelodysplastic', 'myelogenous', 'myeloid', 'myeloproliferative', 'myocardial', 'myoclonic', 'myopathy', 'myotonia', 'nasu-hakola', 'nephrolithiasis', 'nephronophthisis', 'nephrotic', 'neural', 'neuroblastoma', 'neurofibromatosis', 'neuropathy', 'neutral', 'neutropenia', 'nevus', 'nicotine', 'niemann', 'night', 'nonsmall', 'noonan', 'obesity', 'obsessive-compulsive', 'oguchi', 'oligodontia', 'omenn', 'opitz', 'optic', 'orofacial', 'orthostatic', 'ossification', 'osteoarthritis', 'osteogenesis', 'osteolysis', 'osteopetrosis', 'osteoporosis', 'osteosarcoma', 'ovarian', 'ovarioleukodystrophy', 'pachyonychia', 'paget', 'pancreatic', 'pancreatitis', 'parathyroid', 'parietal', 'parkinson', 'paroxysmal', 'peroxisomal', 'persistent', 'pfeiffer', 'phenylketonuria', 'pheochromocytoma', 'phosphoglycerate', 'pick', 'pituitary', 'placental', 'platelet', 'polycystic', 'polycythemia', 'polyposis', 'pontocerebellar', 'porphyria', 'prader-willi', 'premature', 'progressive', 'prostate', 'pseudohermaphroditism', 'pseudohypoaldosteronism', 'pseudohypoparathyroidism', 'psoriasis', 'psoriatic', 'pulmonary', 'pyogenic', 'pyruvate', 'qt', 'refsum', 'renal', 'retinal', 'retinitis', 'rett', 'rhabdomyosarcoma', 'rheumatoid', 'rhizomelic', 'rieger', 'roussy-levy', 'rubenstein-taybi', 'sarcoma', 'sars', 'scapuloperoneal', 'schizophrenia', 'scid', 'senior-loken', 'sensory', 'severe', 'short', 'simpson-golabi-behmel', 'skin,hair,eye', 'sleep', 'smith-magenis', 'spastic', 'spherocytosis', 'spinal', 'spinocerebellar', 'spondylocarpotarsal', 'spondyloepimetaphyseal', 'spondyloepiphyseal', 'squamous', 'stature', 'stevens-johnson', 'stickler', 'stroke', 'subcortical', 'symphalangism', 'syndactyly', 'synpolydactyly', 'systemic', 't-cell', 'tetralogy', 'thalassemia', 'thrombocythemia', 'thrombocytopenia', 'thrombophilia', 'thyroid', 'thyrotropin-releasing', 'transient', 'trichothiodystrophy', 'tuberculosis', 'tuberous', 'ullrich', 'usher', 'uv', 'venous', 'ventricular', 'vohwinkel', 'waardenburg', 'von_willebrand', 'williams-beuren', 'wilms', 'xeroderma', 'zellweger']

#navlakha_phenotypes = ["afibrinogenemia"] 
#navlakha_phenotypes = ["myocardial"] 

navlakha_phenotypes = [ "navlakha_" + p.replace(" ", "_").lower() for p in navlakha_phenotypes ] 

phenotypes += navlakha_phenotypes #!

