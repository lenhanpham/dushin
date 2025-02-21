module load intel/2013_sp1.2.144
echo "compiling dushin.f" > report_compile.txt
ifort -c -g dushin.f
echo "compiling bmatred.f" >> report_compile.txt
ifort -c -O bmatred.f
echo "compiling project0freq.f" >> report_compile.txt
ifort -c -O proj0freq.f
echo "compiling recalc-freq.f" >> report_compile.txt
ifort -c -O recalc-freq.f
echo "compiling subs.f" >> report_compile.txt
ifort -c -g subs.f
echo "compiling subs-dftb.f" >> report_compile.txt
ifort -c -g subs-dftb.f
echo "compiling ddiag.f dmpower.f dmatinv.f" >> report_compile.txt
ifort -c -O ddiag.f dmpower.f dmatinv.f
echo './dushin line' >> report_compile.txt
ifort -g dushin.o subs.o recalc-freq.o proj0freq.o bmatred.o ddiag.o dmpower.o dmatinv.o
-o ./dushin
echo './displace line' >> report_compile.txt
ifort -g displace.f subs.o recalc-freq.o proj0freq.o bmatred.o ddiag.o dmpower.o dmatinv.o
-o ./displace
#echo './plot-modes line' >> report_compile.txt
#ifort -g plot-modes.f subs.o recalc-freq.o proj0freq.o bmatred.o ddiag.o dmpower.o dmatinv.o
-o ./plot-modes
echo './compare-geom line' >> report_compile.txt
ifort -g compare-geom.f subs.o recalc-freq.o proj0freq.o bmatred.o ddiag.o dmpower.o dmatinv.o
-o ./compare-geom