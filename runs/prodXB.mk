prodXB1.mk: prodXB1.mk.mako
	mako-render --var TOKEN=XB1  --var MW_kDa=410.0 --var t_half_loop_closure_intra_s_94kDa='8e-5' --var alkyl_inter_M_s_410kDa='2.5e6' $< >$@  || rm $@
-include prodXB1.mk

prodXB2.mk: prodXB2.mk.mako
	mako-render --var TOKEN=XB2  --var MW_kDa=410.0 --var t_half_loop_closure_intra_s_94kDa='8e-5' --var alkyl_inter_M_s_410kDa='2.5e6' $< >$@ || rm $@
-include prodXB2.mk

prodXB3.mk: prodXB3.mk.mako
	mako-render --var TOKEN=XB3  --var MW_kDa=410.0 --var t_half_loop_closure_intra_s_94kDa='8e-5' --var alkyl_inter_M_s_410kDa='2.5e6' $< >$@ || rm $@
-include prodXB3.mk

.PHONY: all-XB
all-XB: all-XB1 all-XB2 all-XB3
