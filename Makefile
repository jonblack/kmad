OBJECTS = MSA.o Sequences.o ScoringMatrix.o FeaturesProfile.o Profile.o Residue.o substitutionMatrix.o misc.o txtProc.o vecUtil.o findVal.o
GPP = g++ -O3 -std=c++11
CPPFLAGS =  -c 

KMAN: $(OBJECTS)
	$(GPP) -o $@ $(OBJECTS) /usr/local/lib/libboost_program_options.a
MSA.o: MSA.cpp
	$(GPP) $(CPPFLAGS) MSA.cpp
Sequences.o: Sequences.h  Sequences.cpp
	$(GPP) $(CPPFLAGS) Sequences.cpp
ScoringMatrix.o: ScoringMatrix.h ScoringMatrix.cpp
	$(GPP) $(CPPFLAGS) ScoringMatrix.cpp
FeaturesProfile.o: FeaturesProfile.h FeaturesProfile.cpp
	$(GPP) $(CPPFLAGS) FeaturesProfile.cpp
Profile.o: Profile.h Profile.cpp
	$(GPP) $(CPPFLAGS) Profile.cpp
Residue.o: Residue.h Residue.cpp
	$(GPP) $(CPPFLAGS) Residue.cpp
substitutionMatrix.o: substitutionMatrix.h substitutionMatrix.cpp
	$(GPP) $(CPPFLAGS) substitutionMatrix.cpp
misc.o: misc.h misc.cpp
	$(GPP) $(CPPFLAGS) misc.cpp
txtProc.o: txtProc.h txtProc.cpp
	$(GPP) $(CPPFLAGS) txtProc.cpp
vecUtil.o: vecUtil.h vecUtil.cpp
	$(GPP) $(CPPFLAGS) vecUtil.cpp
findVal.o: findVal.h findVal.cpp
	$(GPP) $(CPPFLAGS) findVal.cpp
clean: 
	rm -f *.o KMAN *.h.gch
install:
	cp KMAN /usr/local/bin/KMAN

