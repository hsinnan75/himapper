.KEEP_STAT:

all: main

CXX             = g++
FLAGS		= -Wall -D NDEBUG -O3 -m64
LIB		= -lz -lm -lpthread -lbz2 -lstdc++
SOURCE		= main.cpp GetData.cpp Mapping.cpp AlignmentCandidates.cpp tools.cpp bwt_index.cpp bwt_search.cpp Digest.cpp Normalization.cpp
HEADER		= structure.h
OBJECT		= $(SOURCE:%.cpp=%.o)

main:		$(OBJECT)
			$(CXX) $(FLAGS) $(OBJECT) -o himapper BWT_Index/libbwa.a $(LIB) 
			
%.o:		%.cpp $(HEADER)
			$(CXX) $(FLAGS) -c $<

clean:
		rm -f *.o *~
