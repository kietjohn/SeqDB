import sys
import time
from Bio import Entrez

Entrez.email = "your.email@domain.tld"
if not Entrez.email:
    print "you must add your email address"
    sys.exit(2)

# create an empty list we will fill with the gene names
genes = []
# open the file and begin to iterate over the lines
for line in open('genes.txt', 'r'):
    # at each line, get the value, strip the line ending
    # (\n) from the string, and add the gene name 
    # to our list
    genes.append(line.strip('\n'))
print "genes: ", genes

# create an empty list we will fill with the gene names
species = []
# open the file and begin to iterate over the lines
for line in open('species.txt', 'r'):
    # at each line, get the value, strip the line ending
    # (\n), split the string wherever there are '|' characters
    # and add the fourth element (the GI) of the split string
    # to our species list
    name = line.strip('\n').split('|')[3]
    species.append(name)
print "species: ", species

# for each gene in the gene list
for gene in genes:
    # Print out the gene name
    print "gene: ", gene
    # Open a gene-specific output file, to which we will write our 
    # results. The string substitution makes the output file name 
    # start with the name of the each gene, as we iterate 
    # `for gene in genes`
    output = open('{0}.fasta'.format(gene), 'w')
    # for each taxon in the species list
    for taxon in species:
        try:
            # get rid of the taxon version number
            taxon = taxon.split('.')[0]
            # print the taxon name.  The "\t" means print a tab
            # preceding the taxon name
            print "\ttaxon: ", taxon
            # We know the *name* of the gene that we want, and the
            # accession number of the species in which it is 
            # found, but we do not know the GI number of the 
            # gene (and protein product) we want. So we need to 
            # get that first, since it is the best reference to 
            # other items in NCBI, like the protein sequence.  
            # So, we first build a search  string for an Entrez  
            # search query. The result will contain the GI of  
            # the given gene in the given species.
            #
            # Substitute the taxon and gene name into the string
            terms = "{0}[accn] AND {1}[Gene]".format(taxon, gene)
            # Use the biopython Entrez class and esearch method to
            # search the Gene db using the terms we've defined 
            # above. Entrez Esearch's function is to return 
            # primary identifier (GIs) of records. The results 
            # are returned as XML.
            result = Entrez.esearch(db = 'Gene', term=terms)
            # Parse the XML using `read()` method of the Entrez 
            # class
            record = Entrez.read(result)
            # Get the GI from the first (and only, hopefully) 
            # record that was returned.
            gi = record["IdList"][0]
            # Now that we have the GI number, let's get the actual
            # record for the gene which is similar to a page like
            # this (but in XML versus HTML):
            #
            # http://www.ncbi.nlm.nih.gov/gene/4025118
            gene_record = Entrez.efetch(db="gene", id=gi, retmode="xml")
            # Parse the XML using `read()` method of the Entrez 
            # class
            xml = Entrez.read(gene_record)
            # The `read()` method parses the returned XML to give
            # us many, many nested dictionaries.  This is slightly 
            # confusing because the structure returned actually 
            # contains lists nested within dictionaries (for a 
            # number of reasons). Something like this:
            #
            # xml = {
            #            [
            #                {'record':
            #                    [
            #                        {'name1':value}
            #                    ]
            #                }
            #            ]
            #        }
            #
            # So, what we are doing below is getting the first 
            # element of every list then using the "key" that we
            # know (e.g. "Entrezgene_locus") to return the value 
            # it is associated with, which is usually another 
            # list. Then, we get the first (0th) element of that 
            # list, which is a dictionary. Then, we get use the 
            # "key" ('Gene-commentary_products') to get the next
            # "value", and so forth.
            products = xml[0]['Entrezgene_locus'][0]['Gene-commentary_products'][0]
            # Ensure we've grabbed the correct thing by making 
            # sure that the "type" entry of the item we are on is
            # "peptide".
            assert products['Gene-commentary_type'].attributes['value']== 'peptide'
            # Get the accession number for the protein
            prot_accession = products['Gene-commentary_accession']
            # Print the accession number, preceded by two tab 
            # characters
            print "\t\taccession: ", prot_accession
            # Get the GI of the protein - this is in a different 
            # place that the Accession number - so we traverse
            # the hierarchy again
            prot_gi = xml[0]['Entrezgene_locus'][0]['Gene-commentary_products'][0]['Gene-commentary_seqs'][0]['Seq-loc_whole']['Seq-id']['Seq-id_gi']
            print "\t\tprotein GI: ", prot_gi
            # Finally, use the `efetch()` method of Entrez to grab 
            # the fasta sequence of the protein we're after.
            seq = Entrez.efetch(db="protein", id = prot_gi, retmod='text', rettype='fasta')
            # In one fell swoop, read in the result, and write it 
            # out to our gene-specific file
            output.write(seq.read().replace('\n\n', '\n'))
            # Sleep for 0.3 seconds, to keep NCBI servers happy
            time.sleep(0.3)
        except IndexError:
            error = "Unable to locate record.\n"
            print error
            output.write(error)
            # Sleep for 0.3 seconds, to keep NCBI servers happy
            time.sleep(0.3)
    # Flush the buffer (push the results into the file and close 
    # the file. This essentially "writes" the results to the file.
    output.close()