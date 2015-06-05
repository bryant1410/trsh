include ../../environment

INSTANCIA = ../../ProblemInstances/instancias/A_DU_RM_CL_01.txt

all: MainLan MainSeq FuerzaBruta

clean:
	rm -f MainLan MainSeq MainWan *.o *% *~

MainLan: newGA.req.o newGA.pro.o MainLan.o
	$(CXX) $(LDFLAGS) $^ $(LOADLIBES) $(CPPFLAGS) -o $@

MainWan: newGA.req.o newGA.pro.o MainWan.o
	$(CXX) $(LDFLAGS) $^ $(LOADLIBES) $(CPPFLAGS) -o $@

MainSeq: newGA.req.o newGA.pro.o MainSeq.o
	$(CXX) $(LDFLAGS) $^ $(LOADLIBES) $(CPPFLAGS) -o $@

FuerzaBruta: newGA.req.o newGA.pro.o fuerzabruta.o
	$(CXX) $(LDFLAGS) $^ $(LOADLIBES) $(CPPFLAGS) -o $@

Greedy: newGA.req.o newGA.pro.o greedy.o
	$(CXX) $(LDFLAGS) $^ $(LOADLIBES) $(CPPFLAGS) -o $@

LAN: MainLan
	$(RUN) -np 10 ./MainLan Config.cfg

WAN: MainWan
	$(RUN) -v  -p4pg pgfileWan MainWan

SEQ: MainSeq
	./MainSeq newGA.cfg $(INSTANCIA) res/sol.txt

EjecFuerzaBruta: FuerzaBruta
	./FuerzaBruta $(INSTANCIA)

EjecGreedy: Greedy
	./Greedy $(INSTANCIA)

.PHONY: all clean
