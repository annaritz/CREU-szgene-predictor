## 11/11/17 by Anna Ritz
## Computes statistics about the brain_top network
## from GIANT (http://hb.flatironinstitute.org/download)

## IMPORT STATEMENTS
import matplotlib.pyplot as plt

## GLOBAL VARIABLES
FILENAME = 'brain_top' ## entrez_1, entrez_2, prob
NUMBER_OF_LINES = 41668864 ## number of lines in file

## Trims the network at multiple places, reporting the size.
def main():
	## list of probabilities to split into files.
	probs = [0.125,0.15,0.175,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]

	## list of files to write.  (names have probability in them)
	outfiles = [open('../source/infiles/networkfiles/brain_top_geq_%.3f.txt' % (p),'w') for p in probs]

	## counts for the number of lines for each file.
	counts = [0]*len(probs)

	## read_file
	i = 0 ## use a counter to print to console occaisonally
	with open(FILENAME) as fin: ## open the file.
		for line in fin: ## for every line...
			## print something to screene every million lines.
			if i % 1000000 == 0:
				print('reading line %d (%.4f)' % (i,i/NUMBER_OF_LINES))
			i+=1

			## get the probability (strip whitespace & split on whitespace first)
			prob = float(line.strip().split()[2]) 

			## this is a list of indices of outfiles to write this line to.
			to_write = [p for p in range(len(probs)) if prob >= probs[p]]

			## write the line for each of these outfiles.
			for j in to_write:
				outfiles[j].write(line)	
				counts[j] += 1	## increment counts 

	## print the final counts to screen
	## and close outfiles.
	print('Prob\tCount\tProportion')
	for i in range(len(probs)):
		print('%.3f\t%d\t%.2e' % (probs[i],counts[i],counts[i]/NUMBER_OF_LINES))
		outfiles[i].close()

	print('done')
	return

if __name__ == '__main__':
	main()