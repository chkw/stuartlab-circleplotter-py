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
MATRIX_FILE = $(DATA_DIR)/sample_matrix.tab

# feature to use for determining sample ordering
ORDERING_FEATURE = APC

# each additional matrix file in this list will add a ring
# file format is same as for MATRIX_FILE
OTHER_DATA_FILES = $(DATA_DIR)/sample_matrix_2.tab

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
	circlePlot.py \
		-s sampleIDs.tmp \
		-f features.tmp \
		-o $(ORDERING_FEATURE) \
		-l \
		$(OUTPUT_DIR) \
		$(MATRIX_FILE) \
		$(OTHER_DATA_FILES) ;
	\
	rm -rf sampleIDs.tmp features.tmp ;
	\

#: delete output directory
clean:
	rm -rf $(OUTPUT_DIR)
