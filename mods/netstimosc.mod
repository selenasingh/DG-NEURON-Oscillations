: Modified netstim125 to include oscillating interval between spikes, currently testing
: ISSUES: use of net_move here (commented out)

NEURON	{ 
  ARTIFICIAL_CELL NetStimOsc
  RANGE number, start, forcestop, freq, status, nspk
  THREADSAFE : only true if every instance has its own distinct Random
  POINTER donotuse
}

PARAMETER {
	:interval	= 10 (ms) <1e-9,1e9>: time between spikes (msec)
	number		= 50 <0,1e9>	: number of spikes (independent of noise)
	start		= 0 (ms)	: start of first spike
	forcestop 	= 1200 (ms)	: stop firing spikes
	PI 			= 3.14159265358979323846
	freq 		= 3
	status 		= 1 : unused
	nspk 		= 10 : unused
}

ASSIGNED {
	event (ms)
	on
	ispike
	donotuse
}

PROCEDURE seed(x) {
	set_seed(x)
}

INITIAL {
	on = 0 : off
	ispike = 0
	if (start >= 0 && number > 0) {
		on = 1
		: randomize the first spike so on average it occurs at
		: start + noise*interval
		event = start + invl(interval(start)) - interval(start)
		: but not earlier than 0
		if (event < 0) {
			event = 0
		}
		: no event after time "forcestop"...
		if (event < forcestop) {
			net_send(event, 3)
		}
	}
}	

PROCEDURE init_sequence(t(ms)) {
	if (number > 0) {
		on = 1
		event = 0
		ispike = 0
	}
}

FUNCTION interval (t (ms)) {
	interval = 500*sin(2*PI*freq*(t)/1000 + 1.6)+500
}

FUNCTION invl(mean (ms)) (ms) {
	if (mean <= 0.) {
		mean = .01 (ms) : I would worry if it were 0.
	}
	invl = mean + mean*erand()
}
VERBATIM
double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);
ENDVERBATIM

FUNCTION erand() {
VERBATIM
	if (_p_donotuse) {
		/*
		:Supports separate independent but reproducible streams for
		: each instance. However, the corresponding hoc Random
		: distribution MUST be set to Random.negexp(1)
		*/
		_lerand = nrn_random_pick(_p_donotuse);
	}else{
		/* only can be used in main thread */
		if (_nt != nrn_threads) {
hoc_execerror("multithread random in NetStim"," only via hoc Random");
		}
ENDVERBATIM
		: the old standby. Cannot use if reproducible parallel sim
		: independent of nhost or which host this instance is on
		: is desired, since each instance on this cpu draws from
		: the same stream
		erand = exprand(1)
VERBATIM
	}
ENDVERBATIM
}

PROCEDURE noiseFromRandom() {
VERBATIM
 {
	void** pv = (void**)(&_p_donotuse);
	if (ifarg(1)) {
		*pv = nrn_random_arg(1);
	}else{
		*pv = (void*)0;
	}
 }
ENDVERBATIM
}

PROCEDURE next_invl() {
	if (number > 0) {
		event = invl(interval(t))
	}
	if (ispike >= number) {
		on = 0
	}
}

NET_RECEIVE (w) {
	if (flag == 0) { : external event
		if (w > 0 && on == 0) { : turn on spike sequence
			: but not if a netsend is on the queue
			init_sequence(t)
			: randomize the first spike so on average it occurs at
			: noise*interval (most likely interval is always 0)
			next_invl()
			event = event + interval(t)
			net_send(event, 1)
		:}else if (w > 0 && on == 1) {			: DOESN'T WORK
			: assume interval has changed, recalculate time of next event
			:next_invl()
			:net_move(t + event)
		}else if (w < 0) { : turn off spiking definitively
			on = 0
		}
	}
	if (flag == 3) { : from INITIAL
		if (on == 1) { : but ignore if turned off by external event
			if (t< forcestop) {
				init_sequence(t)
				net_send(0, 1)
			}
		}
	}
	if (flag == 1 && on == 1 && t< forcestop) {
		ispike = ispike + 1
		net_event(t)
		next_invl()
		if (on == 1) {
			net_send(event, 1)
		}
	}
}

