.KEEP_STAT:

all:		himapper

bwt_index:
		$(MAKE) -C src/BWT_Index libbwa.a

himapper: bwt_index
		$(MAKE) -C src
		mkdir -p bin/ && cp -f src/$@ bin/

clean:
		rm -f bin/himapper bin/bwt_index
		$(MAKE) clean -C src
		$(MAKE) clean -C src/BWT_Index

