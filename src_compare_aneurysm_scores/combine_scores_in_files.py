
def concatanate_scores_from_file(file_in, file_scores, file_out):
    """
	Concanates 2 files where first columns are exactly the same (common and sorted without any missing element)
	in both and second column of second file is added to first one and written in a new file
    """
    fIn = open(file_in)
    fScore = open(file_scores)
    fOut = open(file_out, "w")

    for line, line2 in zip(fIn.readlines(), fScore.readlines()):
	id, score = line[:-1].split()
	id2, score2 = line2[:-1].split()
	if id != id2:
	    print "Format error"
	    return
	else:
	    fOut.write("%s\t%s\t%s\n" % (id, score, score2))
    fIn.close()
    fScore.close()
    fOut.close()
    return

if __name__ == "__main__":
    concatanate_scores_from_file("/home/emre/arastirma/data/aneurist/aneurist_Jul_05/aneursym_scores.txt", "/home/emre/arastirma/data/HEFalMp_predicted_gene_disease_associations/aneurysm_seed_pvalues.txt", "aneurysm_seed_scores.txt")

