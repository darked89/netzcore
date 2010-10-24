
# GLOBAL VARIABLES
only_print_command = False
use_cluster = True
#use_cluster = False
leave_one_out_xval = False #True 
score_with_all_seeds = False
only_auc = True # In the analysis_xval if True auc.txt created only, other graphs are not drawn

DEFAULT_TOP_SCORING_PERCENTAGE = 10 #1 #5 #10 # At the time of analysis it was 10
N_LINKER_THRESHOLD = 2 # For Netlink method
DEFAULT_SEED_SCORE = 1.0 # Default score for seed nodes, used when no score given in assoc file
DEFAULT_NON_SEED_SCORE = 0.01 # Default score for non-seed nodes
ALLOWED_MAX_DEGREE = 100000 #175 #90 # Max degree allowed in the graph filtering
N_SAMPLE_GRAPH = 100 # Number of random graphs to be generated
N_X_VAL = 5 #None #5 #182 # Number of cross validation folds, readjusted if leave_one_out_xval is True
N_SEED = None #Will be set during run
N_RANDOM_NEGATIVE_FOLDS = None #10 # Number of non-seed scores to be averaged for negative score calculation, 
			    # If 0 all non seeds are included as they are, If None all non seeds are included averaged to scale with the size of test seeds
REPLICABLE = 123 #63826 #9871354 #123 #None # Assign a predefined seed in randomization for initial test folds creation and N_RANDOM_NEGATIVE_FOLD generation

# FOLLOWING LOCAL ONLY TO MAIN
ignore_experiment_failures = False
delay_experiment = True
tex_format = True #False 
functional_enrichment = False

MODE = "summary" # prepare, score, analyze, compare, summary
user_friendly_id = "biana_no_tap_permuted_p10-omim" #"biana_no_tap-omim_perturbed_%d_10" % i_parameter #"biana_no_tap-all_seed20below" #"biana-all" #"biana_no_tap-omim" #"omim_alzheimer-diabetes" #"all_vs_all" # a.k.a. emre friendly id for compare and summary
summary_seed_cutoff = 1 #None #2 #20 # Seed cutoff considered for inclusion of an experiment in sum_up_experiments, if None seed.dat is not created
prepare_mutated = False #True # Creates permutad/pruned networks 
analyze_network = False #True

ppis = []
#ppis += ["hprd"] #, "ophid"]
#ppis += ["goh", "entrez", "biana_no_tap_no_reliability", "biana_no_tap_relevance", "biana_no_reliability"] 
#ppis += ["rivasi"]
#ppis += ["biana_no_tap_no_reliability", "biana_no_tap_relevance", "biana_no_reliability"] 
#ppis += ["goh"]
#ppis += ["entrez"]
#ppis += ["biana_no_tap_no_reliability"] 
#ppis += ["biana_no_tap_no_reliability_permuted_p%s_%s" % (p, i) for p in xrange(10,110,10) for i in xrange(1,101) ] 
#ppis += ["biana_no_tap_no_reliability_permuted_p%s_%s" % (p, i) for p in xrange(10,110,10) for i in xrange(1,11) ] 
#ppis += ["biana_no_tap_no_reliability_pruned_p%s_%s" % (p, i) for p in xrange(i_parameter,i_parameter+10,10) for i in xrange(1,11) ] 
#ppis += ["biana_no_tap_no_reliability_pruned_p%s_%s" % (p, i) for p in xrange(10,100,10) for i in xrange(1,101) ] 
#ppis += ["biana_no_tap_no_reliability_pruned_p%s_%s" % (p, i) for p in xrange(10,100,10) for i in xrange(1,11) ] 
ppis += ["biana_no_tap_relevance"]
#ppis += ["biana_no_reliability"]
#ppis += ["david"]
#ppis += ["javi"] #["goh"] #["piana_joan_exp", "piana_joan_all"] #["david"] #["goh", "biana_no_tap_no_reliability", "biana_no_reliability", "biana_no_tap_relevance"]
#ppis = ["ori_coexpression_1e-2", "ori_network", "ori_coexpression", "ori_coexpression_colocalization", "ori_colocalization", "ori_coexpression_colocalization_1e-2"]
#ppis += ["ori_no_tap_coexpression_1e-2", "ori_no_tap_network", "ori_no_tap_coexpression", "ori_no_tap_coexpression_colocalization", "ori_no_tap_colocalization", "ori_no_tap_coexpression_colocalization_1e-2"]
#ppi += ["goh_1e5", "biana_coexpression"]

phenotypes = []
#phenotypes += ["navlakha_abdominal"]
#phenotypes += perturbed_omim_phenotypes 
#phenotypes += navlakha_phenotypes
phenotypes += chen_phenotypes + omim_phenotypes + goh_phenotypes 
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
#scoring_parameters += [("nr", 1, 1), ("ff", 1, 5)]
#scoring_parameters += [("nz", 1, 5), ("ns", 3, 2)] 
scoring_parameters += [("nr", 1, 1)]
#scoring_parameters += [("nd", 1, 1)]
#scoring_parameters += [("nz", 1, 5)]
#scoring_parameters += [("ns", 3, 2)]
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
		"ns": [ 0.01*i for i in xrange(1,101) ] }


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


goh_phenotypes = ["goh_developmental", "goh_connective_tissue", "goh_ear,nose,throat", "goh_endocrine", "goh_psychiatric", "goh_immunological", "goh_neurological", "goh_respiratory", "goh_multiple", "goh_renal", "goh_skeletal", "goh_bone", "goh_dermatological", "goh_cancer", "goh_ophthamological", "goh_metabolic", "goh_nutritional", "goh_muscular", "goh_hematological", "goh_gastrointestinal", "goh_cardiovascular"] #,"goh_unclassified"] 

omim_phenotypes = ["alzheimer", "breast cancer", "diabetes", "insulin", "anemia", "myopathy", "neuropathy", "obesity", "parkinson disease", "prostate cancer", "hypertension", "leukemia", "lung cancer", "asthma", "ataxia", "epilepsy", "schizophrenia", "cardiomyopathy", "cataract", "spastic paraplegia", "lymphoma", "mental retardation", "systemic lupus erythematosus"] # "autism", "aneurysm",  
omim_phenotypes = [ "omim_" + "_".join(p.split()) for p in omim_phenotypes ]

chen_phenotypes = ["atherosclerosis",  "ischaemic_stroke",  "systemic_scleroderma",  "migraine",  "epilepsy",  "cirrhosis",  "ulcerative_colitis",  "cervical_carcinoma",  "osteoarthritis",  "inflammatory_bowel_disease",  "myocardial_ischemia",  "endometrial_carcinoma",  "pancreatitis",  "graves_disease",  "neural_tube_defects",  "lymphoma",  "endometriosis",  "autism",  "hypercholesterolaemia"]
chen_phenotypes = [ "chen_" + p for p in chen_phenotypes ]

#navlakha_phenotypes = ["Macular","Pick","Pontocerebellar","Larsen","Opitz","CPT","Pituitary","Sensory","Leber","Thrombophilia","Corneal","Orofacial","Hemophilia","Chondrodysplasia","Noonan","Hypoglycemia","Meckel","Brachydactyly","Hypophosphat","Lethal","Endometrial","Nasopharyngeal","Long","Generalized","Aldosterone","Phenylketonuria","Bethlem","Cystinuria","Cerebellar","Atelosteogenesis","IgA","Myasthenic","Carnitine","Myoclonic","Short","Hyperthyroidism","Osteopetrosis","Nonsmall","Esophageal","Nephrotic","Colon","Urolithiasis","Juvenile","Nicotine","Pfeiffer","Becker","Keratosis","Cone","Megakaryoblastic","von_Hippel-Lindau","Insulin","Thrombocytopenia","Melanoma","Hereditary","Chorea","Trifunctional","Sick","Brugada","Pheochromocytoma","Myocardial","Creatine","Polyposis","Porphyria","Neutral","Ullrich","Pachyonychia","Orthostatic","Prader-Willi","Diabetes","Oguchi","Neuropathy","Osteolysis","Mast","Joubert","Hepatic","Paragangliomas","Ichthyosiform","Angelman","Agammaglobulinemia","Epilepsy","Crohn","Lung","Polycystic","SARS","Skin,hair,eye","Anorexia","Glomerulosclerosis","Jervell","Mucolipidosis","Li","Glucocorticoid","Longevity","Osteogenesis","HIV","Enolase","Attention","Scapuloperoneal","H.","Carcinoid","Focal","Lupus","Hepatocellular","Simpson-Golabi-Behmel","Dyskeratosis","Pseudohypoaldosteronism","Myeloproliferative","Jackson-Weiss","Pulmonary","Homocystinuria","Fetal","Arrhythmogenic","Autoimmune","Optic","Rieger","Muscular","Syndactyly","Histiocytoma","Hemolytic","Dandy-Walker","Coloboma","Giant","Psoriasis","Insensitivity","Multiple","Bladder","Hypodontia","Blood","LADD","Pancreatic","Hyperinsulin","Retinitis","Leprosy","Dystonia","Obesity","Keratitis","Nephrolithiasis","ACTH","Cerebrooculofacioskeletal","Choreoathetosis","Arthrogryposis","Epiphyseal","Cone-rod","Hepatitis","Neutropenia","Factor","Corpus","HDL","Malaria","Leopard","Breast","Ceroid-lipofuscinosis","Congenital","Renal","Paraganglioma","Adrenoleukodystrophy","Trichothiodystrophy","Male","Heterotaxy","Spastic","Bleeding","Symphalangism","Gonadal","Huntington","Hemolytic-uremic","Abdominal","Migraine","Cold-induced","Refsum","AIDS","Myotonia","Atherosclerosis","Osteoarthritis","Lipoprotein","Iridogoniodysgenesis","Beckwith-Wiedemann","Alexander","Tetralogy","Abacavir","Premature","Parietal","Immunodeficiency","Metachromatic","Stature","Boomerang","Mucoepidermoid","Medullary","Anterior","Propionicacidemia","Zellweger","Pancreatitis","Amyotrophic","Hypercholesterolemia","Usher","Fletcher","Pelvic","Tuberous","Alcohol","C1r,C1s","Germ","Emphysema","Campomelic","Adenocarcinoma","Griscelli","Vohwinkel","Pseudohypoparathyroidism","Intracranial","Sitosterolemia","Mitochondrial","Intervertebral","Myelogenous","Exudative","Nasu-Hakola","Amyloidosis","Pigmented","QT","Glycine","Lumbar","High","Azoospermia","Coenzyme","Growth","SCID","Hemangioma","Waardenburg","Dementia","Chromosome","Psoriatic","Hyperlipoproteinemia","Glycogen","Cerebral","Cholestasis","Cutis","Coronary","Amelogenesis","Bernard-Soulier","Colorectal","Myopathy","Pseudohermaphroditism","Fanconi","Phosphoglycerate","Bardet-Biedl","Carney","Blue-cone","Progressive","Pseudoxanthoma","Thalassemia","Squamous","Muscle","Glutaricaciduria","Protoporphyria","Mowat-Wilson","Lymphoproliferative","Mycobacterium","Glycogenosis","Invasive","Sarcoma","Adrenocortical","Craniofacial","Bamforth-Lazarus","Graves","Aplastic","Erythremia","T-cell","Elliptocytosis","Branchiootorenal","Xeroderma","Smith-Magenis","Liddle","Rhizomelic","Insomnia","Hypertriglyceridemia","Tuberculosis","Kallmann","Hypermethioninemia","Central","Wilms","Emery-Dreifuss","Myeloid","Glanzmann","Hypomagnesemia","Omenn","Rheumatoid","Stevens-Johnson","Albinism","Methylmalonic","Intrauterine","Inclusion","Oligodontia","Charcot-Marie-Tooth","Parkinson","Chronic","Hypertension","Atopy","Adrenal","Senior-Loken","Autism","Ceroid","BCG","Epileptic","Paroxysmal","Platelet","Thyrotropin-releasing","Fibromatosis","Dyslexia","Hepatoblastoma","Carbohydrate-deficient","Complement","Hemochromatosis","Leukemia","Hyperparathyroidism","Transient","Fructose","Prostate","Placental","C4","Subcortical","Ehlers-Danlos","Episodic","Glaucoma","Dent","Antley-Bixler","Night","Atrioventricular","Ciliary","Creutzfeldt-Jakob","Heinz","Colorblindness","Ectodermal","Hyperoxaluria","Hypogonadotropic","Schizophrenia","Muir-Torre","Ataxia","Spondyloepimetaphyseal","Neurofibromatosis","Combined","Afibrinogenemia","Asthma","Atrial","Stickler","Loeys-Dietz","Ovarian","Peters","Cystic","Hypotrichosis","Cleft","Convulsions","Hyperekplexia","MODY","Adenomas","Pyogenic","Synpolydactyly","Obsessive-compulsive","Roussy-Levy","Medulloblastoma","Myelodysplastic","Caudal","Celiac","Iron","Tyrosinemia","C1q","Methemoglobinemia","Lymphangioleiomyomatosis","Lipoid","Peroxisomal","Hirschsprung","Leuko","Stroke","Hypercholanemia","Exostoses","Gaucher","Walker-Warburg","Cockayne","Thyroid","Kaposi","Pyruvate","Ectopia","Thrombocythemia","Spondylocarpotarsal","Acromesomelic","Angioedema","Amyotrophy","Blepharophimosis","Dysfibrinogenemia","Epidermolytic","Hyperprolinemia","Hemorrhagic","Nephronophthisis","Williams-Beuren","Epidermolysis","Microphthalmia","Inflammatory","Spinal","Systemic","Hypothyroidism","Hermansky","Osteoporosis","Aicardi-Goutieres","Bradyopsia","Meningioma","Niemann","Rett","Butterfly","Leiomyomatosis","Hypokalemic","Megaloblastic","Adrenomyeloneuropathy","Alzheimer","Hemophagocytic","Rubenstein-Taybi","Microcephaly","Ventricular","Mucopolysaccharidosis","Ossification","Lipoma","Acyl-CoA","Mental","Leigh","Apolipoprotein","UV","Ovarioleukodystrophy","Bare","Dejerine-Sottas","Hypoparathyroidism","Bartter","Retinal","Mandibuloacral","Venous","Nevus","Spondyloepiphyseal","von_Willebrand","Ichthyosis","Basal","C8","Acampomelic","Cardiomyopathy","Maturity-onset","DNA","Holoprosencephaly","Mycobacterial","Persistent","Severe","Fraser","Cataract","Aortic","Foveomacular","Deafness","Marfan","Complex","Krabbe","Erythrocytosis","Gastrointestinal","Spherocytosis","Fundus","Dermatitis","Polycythemia","Rhabdomyosarcoma","Parathyroid","Memory","Hematuria","Alport","GM2-gangliosidosis","Goiter","Brain","Neural","Mismatch","Cirrhosis","Budd-Chiari","Metaphyseal","Crouzon","Lipodystrophy","Paget","Lymphoma","Anemia","Gastric","Macrothrombocytopenia","Cornelia","Alagille","Major","Craniosynostosis","Osteosarcoma","Neuroblastoma","Glioblastoma","Chondrosarcoma","Lissencephaly","Sleep","Spinocerebellar","Malignant","Encephalopathy","Maple"]

navlakha_phenotypes = ["glycogenosis", "keratosis", "roussy-levy", "ectodermal", "hypogonadotropic", "marfan", "systemic", "hypothyroidism", "pseudohypoaldosteronism", "osteoporosis", "emphysema", "graves", "oligodontia", "campomelic", "rhabdomyosarcoma", "neurofibromatosis", "nasu-hakola", "antley-bixler", "jackson-weiss", "anemia", "intervertebral", "sleep", "obesity", "ichthyosiform", "heinz", "bleeding", "hyperekplexia", "syndactyly", "cornelia", "hyperlipoproteinemia", "cockayne", "ovarian", "hereditary", "amelogenesis", "spondyloepiphyseal", "malaria", "ichthyosis", "scapuloperoneal", "osteolysis", "abacavir", "exudative", "insulin", "pick", "atrioventricular", "microphthalmia", "cataract", "angelman", "pancreatitis", "leuko", "spinocerebellar", "bradyopsia", "von_hippel-lindau", "germ", "t-cell", "elliptocytosis", "psoriasis", "liddle", "fanconi", "thrombophilia", "mismatch", "transient", "lipoprotein", "tetralogy", "myasthenic", "stroke", "celiac", "stature", "epidermolysis", "male", "glaucoma", "osteosarcoma", "keratitis", "high", "coenzyme", "endometrial", "kaposi", "niemann", "acth", "cold-induced", "chorea", "larsen", "venous", "thyrotropin-releasing", "leukemia", "megaloblastic", "blood", "carnitine", "prostate", "short", "colorblindness", "leigh", "epilepsy", "brain", "polycystic", "gonadal", "exostoses", "gastrointestinal", "psoriatic", "hemolytic-uremic", "aplastic", "erythremia", "cholestasis", "epileptic", "cutis", "li", "thrombocythemia", "alexander", "ehlers-danlos", "dna", "atherosclerosis", "insomnia", "acampomelic", "microcephaly", "crohn", "spinal", "mody", "nonsmall", "retinitis", "metachromatic", "mitochondrial", "bardet-biedl", "lipoma", "albinism", "nicotine", "mental", "apolipoprotein", "rieger", "hirschsprung", "nephrolithiasis", "choreoathetosis", "hermansky", "aids", "huntington", "alcohol", "cerebral", "paroxysmal", "griscelli", "xeroderma", "hemorrhagic", "phenylketonuria", "bamforth-lazarus", "blepharophimosis", "cirrhosis", "sensory", "alzheimer", "pachyonychia", "lipodystrophy", "becker", "pancreatic", "simpson-golabi-behmel", "nephrotic", "senior-loken", "major", "neuropathy", "parkinson", "hypophosphat", "azoospermia", "symphalangism", "arthrogryposis", "myelogenous", "meningioma", "acromesomelic", "von_willebrand", "corpus", "blue-cone", "leiomyomatosis", "goiter", "lipoid", "pseudohypoparathyroidism", "retinal", "hypertriglyceridemia", "basal", "hemophagocytic", "kallmann", "gm2-gangliosidosis", "congenital", "histiocytoma", "mycobacterial", "chondrodysplasia", "alagille", "iron", "bare", "medullary", "intrauterine", "encephalopathy", "glioblastoma", "hemangioma", "dementia", "squamous", "neutropenia", "generalized", "cleft", "aldosterone", "butterfly", "multiple", "severe", "iridogoniodysgenesis", "paget", "hyperthyroidism", "nephronophthisis", "dent", "dystonia", "agammaglobulinemia", "foveomacular", "longevity", "lissencephaly", "fletcher", "malignant", "lethal", "hypomagnesemia", "tuberous", "adrenocortical", "giant", "epiphyseal", "fundus", "dermatitis", "hypercholesterolemia", "opitz", "factor", "branchiootorenal", "leber", "hiv", "hyperparathyroidism", "cerebellar", "breast", "tuberculosis", "orofacial", "creutzfeldt-jakob", "qt", "progressive", "cystic", "esophageal", "pituitary", "craniosynostosis", "pfeiffer", "brachydactyly", "chondrosarcoma", "schizophrenia", "angioedema", "spondyloepimetaphyseal", "thrombocytopenia", "coloboma", "spherocytosis", "platelet", "myopathy", "complement", "pheochromocytoma", "vohwinkel", "creatine", "enolase", "attention", "phosphoglycerate", "hypodontia", "macular", "myoclonic", "noonan", "macrothrombocytopenia", "ventricular", "meckel", "nevus", "inflammatory", "invasive", "peroxisomal", "zellweger", "colon", "atopy", "ataxia", "pulmonary", "smith-magenis", "abdominal", "dyslexia", "mycobacterium", "immunodeficiency", "budd-chiari", "combined", "rett", "sars", "alport", "adenomas", "osteogenesis", "rhizomelic", "polyposis", "porphyria", "c8", "emery-dreifuss", "lung", "placental", "c4", "hepatocellular", "obsessive-compulsive", "thalassemia", "premature", "adrenoleukodystrophy", "leprosy", "hdl", "arrhythmogenic", "lymphoproliferative", "cerebrooculofacioskeletal", "cone", "sarcoma", "hypoparathyroidism", "afibrinogenemia", "glycogen", "migraine", "cone-rod", "coronary", "complex", "colorectal", "hypotrichosis", "hemochromatosis", "myocardial", "optic", "pyruvate", "leopard", "osteoarthritis", "metaphyseal", "focal", "diabetes", "lupus", "oguchi", "episodic", "neural", "hyperinsulin", "craniofacial", "glanzmann", "trichothiodystrophy", "muscle", "mucoepidermoid", "fetal", "spastic", "methemoglobinemia", "glycine", "charcot-marie-tooth", "loeys-dietz", "hypercholanemia", "hemolytic", "adrenal", "thyroid", "asthma", "atrial", "ectopia", "pontocerebellar", "brugada", "pseudohermaphroditism", "ladd", "parathyroid", "dysfibrinogenemia", "atelosteogenesis", "lymphoma", "epidermolytic", "holoprosencephaly", "parietal", "caudal", "hepatic", "myeloproliferative", "c1q", "rheumatoid", "homocystinuria", "amyotrophic", "hyperoxaluria", "williams-beuren", "growth", "usher", "megakaryoblastic", "hypokalemic", "hepatitis", "refsum", "hepatoblastoma", "myotonia", "bethlem", "night", "crouzon", "cardiomyopathy", "central", "mowat-wilson", "mucopolysaccharidosis", "medulloblastoma", "myeloid", "anterior", "aortic", "inclusion", "deafness", "chronic", "h.", "waardenburg", "erythrocytosis", "gaucher", "anorexia", "glomerulosclerosis", "stickler", "long", "spondylocarpotarsal", "polycythemia", "adenocarcinoma", "rubenstein-taybi", "prader-willi", "memory", "corneal", "intracranial", "neutral", "ullrich", "pyogenic", "orthostatic", "maturity-onset", "renal", "subcortical", "uv", "joubert", "muir-torre", "boomerang", "dejerine-sottas", "ossification", "juvenile", "heterotaxy", "autoimmune", "synpolydactyly", "skin,hair,eye", "omenn", "muscular", "lumbar", "ovarioleukodystrophy", "hypertension", "lymphangioleiomyomatosis", "autism", "chromosome", "melanoma", "jervell", "c1r,c1s", "mast", "scid", "iga", "hematuria", "fibromatosis", "neuroblastoma", "bladder", "adrenomyeloneuropathy", "bernard-soulier", "fructose", "stevens-johnson", "gastric", "hemophilia", "wilms", "osteopetrosis", "bcg", "persistent", "myelodysplastic", "dyskeratosis", "bartter", "beckwith-wiedemann", "amyloidosis"]
navlakha_phenotypes = [ "navlakha_" + p.replace(" ", "_").lower() for p in navlakha_phenotypes ]


hsdl_phenotypes = ["lymphoblast", "prolymphocyte", "bronchi",
			"INTERMEDIATE MESODERM", "LATERAL MESODERM", "PARAXIAL MESODERM",
			"lymphoid progenitor", "hemoctyoblast", "T lymphocyte", "B lymphocyte",
			"PHARYNX", "digestive tube", "small lymphocyte", 
			"myeloid progenitor", "endothelium of blood vessels",
			"ECTODERM", "INNER CELL MASS", "PRIMITIVE GUT", "BLASTOCYST", "ZYGOTE",
			"MORULA", "hemangioblast tissue", "mesenchyme", "heart", 
			"SPLANCHNIC MESODERM", "MESENDODERM", "ENDODERM", "MESODERM"]
#  "lungs", "OUTER EPITHELIUM OF BODY", 
hsdl_phenotypes = [ "hsdl_" + p.replace(" ", "_").lower() for p in hsdl_phenotypes ]


perturbed_omim_phenotypes = [ "perturbed_%s_p%i_%i" % (d, p, i) for d in omim_phenotypes for p in xrange(10,100,10) for i in xrange(1,11) ] #101

