MCC=mcc

gr_ode_ln: GR_ODE_LN.m GR_ODE_LN_Model.m GR_ODE_LN_parameter_list.m
	$(MCC) -m $^ -o $@

clean:
	rm -f *.log *.c *.txt *.sh *.prj gr_ode_ln

