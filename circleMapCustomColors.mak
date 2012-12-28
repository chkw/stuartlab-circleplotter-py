# 08JUN12	chrisw
# This makefile draws CircleMap images using UCSC Stuart lab's circlePlot.py.
# circlePlot.py was written/modified by:
# Steve Benz, Zack Sanborn, Sam Ng, Chris Wong
#
# circlePlot.py plots the image files using  the matplotlib python library

# where to look for data files
DATA_DIR = Data

# columns -> samples
# rows -> features
# This is the file where samples and features are specified.
MATRIX_FILE = $(DATA_DIR)/kegg-slea.expr-level.cdm

# each additional matrix file in this list will add a ring
# file format is same as for MATRIX_FILE
OTHER_DATA_FILES = $(DATA_DIR)/kegg-slea.gbm-subtype-codes.cdm

# A file containing the color codings specified via the -k option.
# When -k is specified, the file must contain color coding information 
# for each ring
COLOR_CODING = $(DATA_DIR)/kegg-slea.colorcoding.tsv

# CircleMap images are output to this directory.
OUTPUT_DIR = Circle_output

# little script meant to transpose a single row into a single column
PERL_TRANSPOSE = perl -ne 'chomp;@l=split(/\t/);foreach $$word (@l){print "$$word\n" if $$word ne ""};'
test:

#: draw CircleMap images
plot:
	head -n 1 $(MATRIX_FILE) \
	| $(PERL_TRANSPOSE) \
	| tail -n +2 \
	> sampleIDs.tmp ;
	\
	cut -f 1 $(MATRIX_FILE) \
	| tail -n +2 \
	> features.tmp ;
	\
	mkdir -p $(OUTPUT_DIR) ;
	\
	python circlePlot.py \
		-s sampleIDs.tmp \
		-f features.tmp \
		-l \
		-k $(COLOR_CODING) \
		$(OUTPUT_DIR) \
		$(MATRIX_FILE) \
		$(OTHER_DATA_FILES) ;
	\
	rm -rf sampleIDs.tmp features.tmp ;
	\

#: delete output directory
clean:
	rm -rf $(OUTPUT_DIR)
	rm *tmp
	rm *pyc
