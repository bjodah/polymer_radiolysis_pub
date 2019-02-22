# Example usage:
#
#  $ CXXFLAGS=-Og make -B DEBUG=1 CXX=g++-7
#  $ CXXFLAGS="-O3 -DENSEMBLE_VERBOSE -funroll-loops -ffast-math" make demo
#  $ CMAKE_PREFIX_PATH=/opt/symengine-262478e/ CPLUS_INCLUDE_PATH=~/vc/kinetgine:~/.local/include:~/vc/anyode/include:~/vc/pycvodes/pycvodes/include/:/opt/symengine-262478e/include CXX=g++-7 make WITH_KINETGINE=1 demo -B


varied-W1.json: mk_varied.py
	./mk_varied.py --outfile $@
prodW1.mk: prodW1.mk.mako
	mako-render --var N0=100 --var TOKEN=W1 --var N_LEVELS="[1,10,100,1000,10000]" $< >$@
-include prodW1.mk

varied-W2.json: mk_varied.py
	./mk_varied.py --outfile $@
prodW2.mk: prodW2.mk.mako
	mako-render --var N0=100 --var TOKEN=W2 --var N_LEVELS="[1,10,100,1000,10000]" $< >$@
-include prodW2.mk

varied-W3.json: mk_varied.py
	./mk_varied.py --outfile $@
prodW3.mk: prodW2.mk.mako
	mako-render --var N0=100 --var TOKEN=W3 --var N_LEVELS="[1,10,100,1000,10000]" $< >$@
-include prodW3.mk

varied-S1.json: mk_varied.py
	./mk_varied.py --outfile $@  --dr-pulse-Gy-s 16e6 --t-pulse-s 50e-9 --f-pulse-Hz 25 --tot-dose-Gy 4.8e3
prodS1.mk: prodW2.mk.mako
	mako-render --var N0=100 --var TOKEN=S1 --var N_LEVELS="[1,10,100,1000,10000]" $< >$@
-include prodS1.mk

varied-W4.json: mk_varied.py
	./mk_varied.py --outfile $@ --count-wait "[(2,0.090)]"
prodW4.mk: prodW2.mk.mako
	mako-render --var N0=100 --var TOKEN=W4 --var N_LEVELS="[1,10,100,1000,10000]" $< >$@
-include prodW4.mk


varied-Q1.json: mk_varied.py
	./mk_varied.py --outfile $@  --tot-dose-Gy 320

prodQ1.mk: prodQ1.mk.mako
	mako-render --var N0=10 --var TOKEN=Q1 --var N_LEVELS="[1,10,100]" $< >$@
-include prodQ1.mk

include prodX.mk

prodP.mk: prodP.mk.mako
	mako-render --var DOSERATE=1.0e7 --var TOKEN=P --var MW_kDa=410.0 --var prim=', "loop_scaling_inter": -1, "chain_rad_chain_rad_M_s": 2.5e7' $< >$@
-include prodP.mk

prodP2.mk: prodP.mk.mako
	mako-render --var DOSERATE=1.0e7 --var TOKEN=P2 --var MW_kDa=41.0 --var prim=', "loop_scaling_inter": -1, "chain_rad_chain_rad_M_s": 9.95e7' $< >$@
-include prodP2.mk

prodP3.mk: prodP.mk.mako
	mako-render --var DOSERATE=1.0e7 --var TOKEN=P3 --var MW_kDa=4.1 --var prim=', "loop_scaling_inter": -1, "chain_rad_chain_rad_M_s": 39.6e7' $< >$@
-include prodP3.mk


prodPa.mk: prodP.mk.mako
	mako-render --var DOSERATE=1.0e7 --var TOKEN=Pa --var MW_kDa=410.0 --var prim=', "loop_scaling_inter": -0.5, "chain_rad_chain_rad_M_s": 2.5e7' $< >$@
-include prodPa.mk

prodPa2.mk: prodP.mk.mako
	mako-render --var DOSERATE=1.0e7 --var TOKEN=Pa2 --var MW_kDa=41.0 --var prim=', "loop_scaling_inter": -0.5, "chain_rad_chain_rad_M_s": 9.95e7' $< >$@
-include prodPa2.mk

prodPa3.mk: prodP.mk.mako
	mako-render --var DOSERATE=1.0e7 --var TOKEN=Pa3 --var MW_kDa=4.1 --var prim=', "loop_scaling_inter": -0.5, "chain_rad_chain_rad_M_s": 39.6e7' $< >$@
-include prodPa3.mk


prodPb.mk: prodP.mk.mako
	mako-render --var DOSERATE=1.0e7 --var TOKEN=Pb --var MW_kDa=410.0 --var prim=', "loop_scaling_inter": -2, "chain_rad_chain_rad_M_s": 2.5e7' $< >$@
-include prodPb.mk

prodPb2.mk: prodP.mk.mako
	mako-render --var DOSERATE=1.0e7 --var TOKEN=Pb2 --var MW_kDa=41.0 --var prim=', "loop_scaling_inter": -2, "chain_rad_chain_rad_M_s": 9.95e7' $< >$@
-include prodPb2.mk

prodPb3.mk: prodP.mk.mako
	mako-render --var DOSERATE=1.0e7 --var TOKEN=Pb3 --var MW_kDa=4.1 --var prim=', "loop_scaling_inter": -2, "chain_rad_chain_rad_M_s": 39.6e7' $< >$@
-include prodPb3.mk

prodPnointer.mk: prodP.mk.mako
	mako-render --var DOSERATE=1.0e7 --var TOKEN=Pnointer --var MW_kDa=410.0 --var prim=', "loop_scaling_inter": -1, "chain_rad_chain_rad_M_s": 0.0' $< >$@
-include prodPnointer.mk

prodPnointra.mk: prodP.mk.mako
	mako-render --var DOSERATE=1.0e7 --var TOKEN=Pnointra --var MW_kDa=410.0 --var prim=', "loop_scaling_inter": -1, "t_half_loop_closure_intra_s_94kDa": 1e10' $< >$@
-include prodPnointra.mk

prodPnointernoloopscaling.mk: prodP.mk.mako
	mako-render --var DOSERATE=1.0e7 --var TOKEN=Pnointernoloopscaling --var MW_kDa=410.0 --var prim=', "loop_scaling_inter": 0.0, "chain_rad_chain_rad_M_s": 0.0, "loop_scaling_intra": 0.0' $< >$@
-include prodPnointernoloopscaling.mk

prodPnoloopscaling.mk: prodP.mk.mako
	mako-render --var DOSERATE=1.0e7 --var TOKEN=Pnoloopscaling --var MW_kDa=410.0 --var prim=', "loop_scaling_inter": 0.0, "loop_scaling_intra": 0.0' $< >$@
-include prodPnoloopscaling.mk

prodR.mk: prodP.mk.mako
	mako-render --var DOSERATE=2.0e7 --var TOKEN=R --var MW_kDa=410.0 --var prim=', "loop_scaling_inter": -1, "chain_rad_chain_rad_M_s": 2.5e7' $< >$@
-include prodR.mk

prodR2.mk: prodP.mk.mako
	mako-render --var DOSERATE=2.0e7 --var TOKEN=R2 --var MW_kDa=41.0 --var prim=', "loop_scaling_inter": -1, "chain_rad_chain_rad_M_s": 9.95e7' $< >$@
-include prodR2.mk

prodR3.mk: prodP.mk.mako
	mako-render --var DOSERATE=2.0e7 --var TOKEN=R3 --var MW_kDa=4.1 --var prim=', "loop_scaling_inter": -1, "chain_rad_chain_rad_M_s": 39.6e7' $< >$@
-include prodR3.mk


prodEP.mk: prodEP.mk.mako
	mako-render --var DOSERATE=1.0e7 --var TOKEN=EP --var MW_kDa=410.0 --var prim=', "loop_scaling_inter": -1, "chain_rad_chain_rad_M_s": 2.5e7' $< >$@
-include prodEP.mk

prodEP2.mk: prodEP.mk.mako
	mako-render --var DOSERATE=1.0e7 --var TOKEN=EP2 --var MW_kDa=41.0 --var prim=', "loop_scaling_inter": -1, "chain_rad_chain_rad_M_s": 9.95e7' $< >$@
-include prodEP2.mk

prodEP3.mk: prodEP.mk.mako
	mako-render --var DOSERATE=1.0e7 --var TOKEN=EP3 --var MW_kDa=4.1 --var prim=', "loop_scaling_inter": -1, "chain_rad_chain_rad_M_s": 39.6e7' $< >$@
-include prodEP3.mk

prodER.mk: prodEP.mk.mako
	mako-render --var DOSERATE=2.0e7 --var TOKEN=ER --var MW_kDa=410.0 --var prim=', "loop_scaling_inter": -1, "chain_rad_chain_rad_M_s": 2.5e7' $< >$@
-include prodER.mk

prodER2.mk: prodEP.mk.mako
	mako-render --var DOSERATE=2.0e7 --var TOKEN=ER2 --var MW_kDa=41.0 --var prim=', "loop_scaling_inter": -1, "chain_rad_chain_rad_M_s": 9.95e7' $< >$@
-include prodER2.mk

prodER3.mk: prodEP.mk.mako
	mako-render --var DOSERATE=2.0e7 --var TOKEN=ER3 --var MW_kDa=4.1 --var prim=', "loop_scaling_inter": -1, "chain_rad_chain_rad_M_s": 39.6e7' $< >$@
-include prodER3.mk

prodES.mk: prodEP.mk.mako
	mako-render --var DOSERATE=2.0e7 --var TOKEN=ES --var MW_kDa=410.0 --var prim=', "loop_scaling_inter": -1, "chain_rad_chain_rad_M_s": 2.5e7, "loop_scaling_intra": -2' $< >$@
-include prodES.mk


varied-MPPL1a.json: mk_varied.py
	./mk_varied.py --outfile $@  --dr-pulse-Gy-s 32e6 --t-pulse-s 5e-6 --f-pulse-Hz 400 --count-wait "[(2,0.090),(125,60)]" --tot-dose-Gy 80e3 
prodMPPL1a.mk: prodMP.mk.mako
	mako-render --var TOKEN=MPPL1a --var MW_kDa=410.0 $< >$@
-include prodMPPL1a.mk

varied-MPPL1aZ.json: mk_varied.py
	./mk_varied.py --outfile $@  --dr-pulse-Gy-s 32e6 --t-pulse-s 5e-6 --f-pulse-Hz 400 --count-wait "[(2,0.090),(125,60)]" --tot-dose-Gy 320 
prodMPPL1aZ.mk: prodMP.mk.mako
	mako-render --var TOKEN=MPPL1aZ --var MW_kDa=410.0 $< >$@
-include prodMPPL1aZ.mk

varied-MP3.json: mk_varied.py
	./mk_varied.py --outfile $@  --dr-pulse-Gy-s 32e6 --t-pulse-s 5e-6 --f-pulse-Hz 400 --count-wait "[(2,0.090),(125,60)]" --tot-dose-Gy 3.2e3 
prodMP3.mk: prodMP3.mk.mako
	mako-render --var TOKEN=MP3 --var MW_kDa=410.0 $< >$@
-include prodMP3.mk

varied-MP4.json: mk_varied.py
	./mk_varied.py --outfile $@  --dr-pulse-Gy-s 32e6 --t-pulse-s 5e-6 --f-pulse-Hz 400 --count-wait "[(2,0.090),(125,60)]" --tot-dose-Gy 3.2e3 
prodMP4.mk: prodMP3.mk.mako
	mako-render --var TOKEN=MP4 --var MW_kDa=410.0 $< >$@
-include prodMP4.mk


varied-MP5a.json: mk_varied.py
	./mk_varied.py --outfile $@  --dr-pulse-Gy-s 32e6 --t-pulse-s 5e-6 --f-pulse-Hz 400 --count-wait "[(2,0.090),(125,60)]" --tot-dose-Gy 3.2e3 
prodMP5a.mk: prodMP5.mk.mako
	mako-render --var TOKEN=MP5a --var MW_kDa=410.0 --var prim=', "t_half_loop_closure_intra_s_94kDa": 8e-3' $< >$@
-include prodMP5a.mk

varied-MP5b.json: mk_varied.py
	./mk_varied.py --outfile $@  --dr-pulse-Gy-s 32e6 --t-pulse-s 5e-6 --f-pulse-Hz 400 --count-wait "[(2,0.090),(125,60)]" --tot-dose-Gy 3.2e3 
prodMP5b.mk: prodMP5.mk.mako
	mako-render --var TOKEN=MP5b --var MW_kDa=410.0 --var prim=', "t_half_loop_closure_intra_s_94kDa": 8e-4' $< >$@
-include prodMP5b.mk

varied-MP7a.json: mk_varied.py
	./mk_varied.py --outfile $@  --dr-pulse-Gy-s 32e6 --t-pulse-s 5e-6 --f-pulse-Hz 400 --count-wait "[(2,0.090),(125,60)]" --tot-dose-Gy 3.2e3 
prodMP7a.mk: prodMP7.mk.mako
	mako-render --var TOKEN=MP7a --var MW_kDa=410.0 --var t_half_loop_closure_intra_s_94kDa='8e-3' $< >$@
-include prodMP7a.mk

varied-MP7b.json: mk_varied.py
	./mk_varied.py --outfile $@  --dr-pulse-Gy-s 32e6 --t-pulse-s 5e-6 --f-pulse-Hz 400 --count-wait "[(2,0.090),(125,60)]" --tot-dose-Gy 3.2e3 
prodMP7b.mk: prodMP7.mk.mako
	mako-render --var TOKEN=MP7b --var MW_kDa=410.0 --var t_half_loop_closure_intra_s_94kDa='8e-4' $< >$@
-include prodMP7b.mk

varied-MP7c.json: mk_varied.py
	./mk_varied.py --outfile $@  --dr-pulse-Gy-s 32e6 --t-pulse-s 5e-6 --f-pulse-Hz 400 --count-wait "[(2,0.090),(125,60)]" --tot-dose-Gy 3.2e3 
prodMP7c.mk: prodMP7.mk.mako
	mako-render --var TOKEN=MP7c --var MW_kDa=410.0 --var t_half_loop_closure_intra_s_94kDa='8e-5' $< >$@
-include prodMP7c.mk

varied-MP7b2.json: mk_varied.py
	./mk_varied.py --outfile $@  --dr-pulse-Gy-s 32e6 --t-pulse-s 5e-6 --f-pulse-Hz 400 --count-wait "[(2,0.090),(125,60)]" --tot-dose-Gy 3.2e3 
prodMP7b2.mk: prodMP7.mk.mako
	mako-render --var TOKEN=MP7b2 --var MW_kDa=410.0 --var t_half_loop_closure_intra_s_94kDa='8e-4' $< >$@
-include prodMP7b2.mk

varied-MP7b3.json: mk_varied.py
	./mk_varied.py --outfile $@  --dr-pulse-Gy-s 32e6 --t-pulse-s 5e-6 --f-pulse-Hz 400 --count-wait "[(2,0.090),(125,60)]" --tot-dose-Gy 3.2e3 
prodMP7b3.mk: prodMP7.mk.mako
	mako-render --var TOKEN=MP7b3 --var MW_kDa=410.0 --var t_half_loop_closure_intra_s_94kDa='8e-4' $< >$@
-include prodMP7b3.mk

varied-MP7d.json: mk_varied.py
	./mk_varied.py --outfile $@  --dr-pulse-Gy-s 32e6 --t-pulse-s 5e-6 --f-pulse-Hz 400 --count-wait "[(2,0.090),(125,60)]" --tot-dose-Gy 3.2e3 
prodMP7d.mk: prodMP7.mk.mako
	mako-render --var TOKEN=MP7d --var MW_kDa=410.0 --var t_half_loop_closure_intra_s_94kDa='8e-5' $< >$@
-include prodMP7d.mk

varied-MP8d.json: mk_varied.py
	./mk_varied.py --outfile $@  --dr-pulse-Gy-s 32e6 --t-pulse-s 5e-6 --f-pulse-Hz 400 --count-wait "[(2,0.090),(125,60)]" --tot-dose-Gy 3.2e3 
prodMP8d.mk: prodMP8.mk.mako
	mako-render --var TOKEN=MP8d --var MW_kDa=410.0 --var t_half_loop_closure_intra_s_94kDa='8e-5' $< >$@
-include prodMP8d.mk

varied-MP8e.json: mk_varied.py
	./mk_varied.py --outfile $@  --dr-pulse-Gy-s 32e6 --t-pulse-s 5e-6 --f-pulse-Hz 400 --count-wait "[(2,0.090),(125,60)]" --tot-dose-Gy 3.2e3 
prodMP8e.mk: prodMP8.mk.mako
	mako-render --var TOKEN=MP8e --var MW_kDa=410.0 --var t_half_loop_closure_intra_s_94kDa='8e-6' $< >$@
-include prodMP8e.mk

varied-MP9e.json: mk_varied.py
	./mk_varied.py --outfile $@  --dr-pulse-Gy-s 32e6 --t-pulse-s 5e-6 --f-pulse-Hz 400 --count-wait "[(2,0.090),(125,60)]" --tot-dose-Gy 3.2e3 
prodMP9e.mk: prodMP9.mk.mako
	mako-render --var TOKEN=MP9e --var MW_kDa=410.0 --var t_half_loop_closure_intra_s_94kDa='8e-6' --var alkyl_inter_M_s_410kDa='2.5e7' --var peroxy_factor='1' $< >$@
-include prodMP9e.mk

varied-MP9f.json: mk_varied.py
	./mk_varied.py --outfile $@  --dr-pulse-Gy-s 32e6 --t-pulse-s 5e-6 --f-pulse-Hz 400 --count-wait "[(2,0.090),(125,60)]" --tot-dose-Gy 3.2e3 
prodMP9f.mk: prodMP9.mk.mako
	mako-render --var TOKEN=MP9f --var MW_kDa=410.0 --var t_half_loop_closure_intra_s_94kDa='8e-6' --var alkyl_inter_M_s_410kDa='2.5e6' --var peroxy_factor='10' $< >$@
-include prodMP9f.mk

varied-MP10d.json: mk_varied.py
	./mk_varied.py --outfile $@  --dr-pulse-Gy-s 32e6 --t-pulse-s 5e-6 --f-pulse-Hz 400 --count-wait "[(2,0.090),(125,60)]" --tot-dose-Gy 20e3 
prodMP10d.mk: prodMP10.mk.mako
	mako-render --var TOKEN=MP10d --var MW_kDa=410.0 --var t_half_loop_closure_intra_s_94kDa='8e-5' --var alkyl_inter_M_s_410kDa='2.5e6'  $< >$@
-include prodMP10d.mk

varied-MP10d40.json: mk_varied.py
	./mk_varied.py --outfile $@  --dr-pulse-Gy-s 32e6 --t-pulse-s 5e-6 --f-pulse-Hz 400 --count-wait "[(2,0.090),(125,60)]" --tot-dose-Gy 40e3 
prodMP10d40.mk: prodMP10.mk.mako
	mako-render --var TOKEN=MP10d40 --var MW_kDa=410.0 --var t_half_loop_closure_intra_s_94kDa='8e-5' --var alkyl_inter_M_s_410kDa='2.5e6'  $< >$@
-include prodMP10d40.mk

varied-MP10d40mf.json: mk_varied.py
	./mk_varied.py --outfile $@  --dr-pulse-Gy-s 32e6 --t-pulse-s 5e-6 --f-pulse-Hz 400 --count-wait "[(2,0.090),(125,60)]" --tot-dose-Gy 40e3 
prodMP10d40mf.mk: prodMP10.mk.mako
	mako-render --var TOKEN=MP10d40mf --var MW_kDa=410.0 --var t_half_loop_closure_intra_s_94kDa='8e-5' --var alkyl_inter_M_s_410kDa='2.5e6'  $< >$@
-include prodMP10d40mf.mk

varied-MP10d032mf.json: mk_varied.py
	./mk_varied.py --outfile $@  --dr-pulse-Gy-s 32e6 --t-pulse-s 5e-6 --f-pulse-Hz 0320 --count-wait "[(2,0.090),(125,60)]" --tot-dose-Gy 0.32e3 
prodMP10d032mf.mk: prodMP10d032.mk.mako
	mako-render --var TOKEN=MP10d032mf --var MW_kDa=410.0 --var t_half_loop_closure_intra_s_94kDa='8e-5' --var alkyl_inter_M_s_410kDa='2.5e6'  $< >$@
-include prodMP10d032mf.mk

include prodXB.mk
