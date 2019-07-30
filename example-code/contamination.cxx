// -*- outline-regexp: "/\\*-+";  -*-
/*-  Identification and Changes  */

/*
  contamination.cxx -- Written by Randall Gray 
  Initial coding: 
  Date: 2007.10.22
  Location: floyd.hba.marine.csiro.au:/home/gray/Projects/InVitro/contamination.cxx


  2019-07-30-07:21:00 -- minor explanatory note

  There is a "member cube" referred to below.  This is a mechanism
  which allows me to implement mortality on a population exposed to
  multiple causes of mortality without the problem of double counting.
  It makes a necessary assumption of linear independence in mortality
  (that is it implicitly rejects "potentiating" or "ameliorating"
  interactions), but it is a lot better than iterating through
  (potenitally) hundreds of thousands of individual interactions
  within a population.  The cube is actually a hypercube where each
  axis corresponds to the influence of one source of mortality.

  The code used to calculate susceptibility (the parameterisation
  which includes the O.D.E. stuff) and the other Evaluate() and
  Calculate() calls which may appear are part of a C library I wrote
  many years ago to allow me to move some of the traditionally "hard
  coded" aspects of a model into its parameterisation.  This means
  that it is easy to explore alternative functional responses to
  environmental conditions. In practice, the p-coded functions are
  slower, but usually not enough slower to warrant change; when it
  does need change, it is a simple matter to replace a call to the
  evaluator with a hardcoded function.


  Compile with: g++ -ggdb -DDEBUGGING -Wall -o contamination contamination.cxx
  History:

  $Log: contamination.cxx,v $
  Revision 1.22  2009/04/17 01:25:25  roger
  Major patch bomb which might break stuff

  Revision 1.21  2009/03/18 03:54:18  gray
  Code writers committed for psychiatric assessment w.r.t. stasis/migration

  Revision 1.20  2009/02/20 01:27:40  gray
  Whopping big patch bomb

  Revision 1.19  2008/10/17 05:15:15  gray
  Wow, what can I say?  Animal, AnimalP, dynamics are the culprits, changed UERBOSE to VERBOSE

  Revision 1.18  2008/10/06 03:47:14  gray
  Not this one.

  Revision 1.17  2008/08/29 04:11:21  gray
  More messing with indenting!

  Revision 1.16  2008/08/26 08:21:39  gray
  Massive mods.  Lotsa space fixing and such. Merged Beths stuff in

  Revision 1.15  2008/08/01 01:04:54  gray
  Adding a set of celestial bodies

  Revision 1.14  2008/07/10 23:32:08  gray
  debuggering

  Revision 1.13  2008/04/30 00:50:16  gray
  Moved repository to Njal, converted from "se-model" to "invitro"

  Revision 1.12  2008/04/01 01:06:19  roger
  shuffle PrmExpr stuff so it is a bit cleaner and easier to add stuff to

  Revision 1.11  2008/03/31 00:10:55  gray
  Incorfporated script files

  Revision 1.10  2008/03/07 00:56:00  gray
  enable state stuff

  Revision 1.9  2008/02/27 06:09:20  gray
  This commits all the agents files that have been wobbled in the
  great PARAM_NOR adventure.  The attributions will go away.

  Revision 1.8  2008/02/26 22:44:23  gray
  Eat toxic scum!

  Revision 1.7  2008/02/21 04:21:26  gray
  more contaminant propagating

  Revision 1.6  2008/02/21 01:12:11  gray
  contaminant propagation

  Revision 1.5  2008/02/20 02:42:10  gray
  Dunno what we were doing, but it needs to go in

  Revision 1.4  2007/11/29 05:02:56  gray
  Doodled the *ERBOSE, fatal and warning macros.  Yahoo.

  Revision 1.3  2007/10/31 03:58:43  gray
  working on testing the contamination influences on fecundity

  Revision 1.2  2007/10/24 00:21:03  gray
  Minor emacs tweaks -- no code changes yet

  Revision 1.1  2007/10/22 03:16:40  gray
  initial coding


*/

/*-  Copyright  */

/*
  (C) 2007 CSIRO Australia
  All rights reserved
*/

/*-  Discussion  */
/*
   This class is inherited by entities which are susceptible to contaminants.
	
	
*/
/*

  THIS USES THE FOLLOWING

  cinfo[i].DT = cc->GetVarRef2("dt");
  cinfo[i].concv = cc->GetVarRef2("conc");
  cinfo[i].imassv = cc->GetVarRef2("imass");
  cinfo[i].atev = cc->GetVarRef2("ate");
  cinfo[i].currentloadv = cc->GetVarRef2("current_load");

*/





/*-  Configuration stuff  */
#pragma implementation

#ifndef __contamination_cxx
#define __contamination_cxx
#endif

/*-  Included files  */

#include "contamination.hxx"
#include "contsrc.hxx"

#include "memchk.h"

/*-  Local variables, constants, and defines  */

/*-  Code  */

/*-- serialisation code for the individual contaminants */

/*--- void *Contamination::Get_cinfo_State(int i, int *len) -- package things */
void *Contamination::Get_cinfo_State(int i, int *len) {
	assert(len);
	void *v[2];
	int l[2];

	v[0] = &cinfo[i].current_load;
	l[0] = sizeof(cinfo[i].current_load);
	assert(cinfo[i].name);
	v[1] = cinfo[i].name;
	l[1] = strlen(cinfo[i].name)+1;
	return pack_mem(v, l, 2, len);
}

/*--- Contamination::Set_cinfo_State(void *data, int len, int i) -- restore state */
void Contamination::Set_cinfo_State(void *data, int len, int i) {
	void **v;
	int *l;

	int n = unpack_mem(data, len, &v, &l);
	assert(n == 2);

	assert(v[0]);
	assert(l[0] == sizeof(double));
	cinfo[i].current_load = *(double *)v[0];
	if (cinfo[i].name) Free(cinfo[i].name);

	assert(v[1]);
	assert(l[1] > 1);
	cinfo[i].name = Strdup((char *)v[1]);
	Free(v);
	Free(l);
}

/*-- double Contamination::Level(char *name) -- Return the level of indicated contaminant */
double Contamination::Level(char *name) {
	int i;

	if (!member_cube) return DNaN;
	assert(cinfo);
	for (i = 0; i < n_cinfo; i++) {
		if (!strcmp(name, cinfo[i].name)) return cinfo[i].current_load;
	}
	return DNaN;
}


/*-- serialisation code for the whole set of  contaminants */

/*--- void *Contamination::GetState(int* len) --  package things */
void *Contamination::GetState(int* len) {
	assert(len);
	void **v = 0, *d = 0;
	int i, *l = 0;

	v = (void **)Calloc(n_cinfo + 4, sizeof(void *));
	if (!v) abort();
	l = (int *)Calloc(n_cinfo + 4, sizeof(int));
	if (!l) abort();


	v[0] = ctaxon;
	l[0] = ctaxon?strlen(ctaxon)+1:0;
	v[1] = &n_cinfo;
	l[1] = sizeof(n_cinfo);

	if (member_cube) {
		v[2] = member_cube->GetState(&l[2]);
	}
	else v[2] = 0, l[2] = 0;
	v[3] = cname;
	l[3] = cname?strlen(cname)+1:0;

	for (i = 0; i < n_cinfo; i++) {
		v[i+4] = Get_cinfo_State(i, &l[i+4]);
	}

	d = pack_mem(v, l, 4 + n_cinfo, len);
	if (v[2]) Free(v[2]);
	for (i = 0; i < n_cinfo; i++) {
		if (v[i+4]) Free(v[i+4]);
	}

	Free(v);
	Free(l);

	return d;
}

/*--- Contamination::SetState(void *data, int len) -- unpackage things */
void Contamination::SetState(void *data, int len) {
	assert(data);
	void **v;
	int *l, n;
	n = unpack_mem(data, len, &v, &l);
	assert(n >= 4);

	assert(v[1]);
	assert(l[1] == sizeof(int));
			 
	n_cinfo = *(int*)v[1];
	assert(n_cinfo == n-4);

	if (cinfo) free_cinfo();
	if (member_cube) delete member_cube;

	if (ctaxon) {
		Free(ctaxon);
		ctaxon = 0;
	}
	if (v[0]) ctaxon = Strdup((char *)v[0]);

	if (cname) {
		Free(cname);
		cname = 0;
	}
	if (v[3]) cname = Strdup((char *)v[3]);

	member_cube = 0;
	if (v[2]) {
		assert(l[2] > 0);
		member_cube = new Cube(n_cinfo+1);
		member_cube->SetState(v[2], l[2]);
	}
	if (n_cinfo > 0) {

		cinfo = (_cinfo*)Calloc(n_cinfo, sizeof(_cinfo));
		if (!cinfo) abort();
		
		for (int i = 0; i < n_cinfo; i++) {
			cinfo[i].vbid = -1;
			cinfo[i].name = 0;
			cinfo[i].current_load = 0;
		}

		for (int i = 0; i < n_cinfo; i++) {
			Set_cinfo_State(v[i+4], l[i+4], i);
		}
	}

	Free(v); Free(l);
	// Most of the state setting is actually read only parameter data
	// and should be done in ReInit or PostInit
}

/*-- Contamination::Reset()  -- reset (zero) contamination data */
void 
	PrmEnvExpr::Reset();
	ContaminantSink::Reset();
	DeathLogger::Reset();
	cname = 0;
	ctaxon = 0;
	cinfo = 0;
	n_cinfo = 0;
	member_cube = 0;
}

/*-- int Contamination::ReInit(int attach) -- reinitialise after moving between kernels */
int Contamination::ReInit(int attach) {
	if (!ContaminantSink::Init(ctaxon)) return 0;;
	if (!DeathLogger::Init(ctaxon,cname)) return 0;;
	if (!init_contaminant_stuff()) return 0;

	return 1;
}


/*-- Contamination::zero() -- zero out data without freeing */
void Contamination::zero()
{
	cname = 0;
	ctaxon = 0;
	cinfo = 0;
	n_cinfo = 0;
	member_cube = 0;
}

/*-- Constructors / destructors  for Contamination */

/*--- Contamination constructor */
Contamination::Contamination()
{
	zero();
}

/*--- Contamination destructor */
Contamination::~Contamination()
{
	if (ctaxon) Free(ctaxon);
	if (cname) Free(cname);
	if (cinfo) free_cinfo();
	if (member_cube) {
		delete member_cube;
		member_cube = 0;
	}
}


/*
  The accumulation and depuration are handled by differential equations specified like so:


  ode format:	ode(dy/dt = expression, y(lb) = y0, stepsize, value)			 
  in y and t

  In the parameter files of sink agents there will be blocks which 
  tell the agent what kinds of agent to query for pertinent information;
  thus: 

ContaminantSink { # Juvenile sharks
	contaminants {
		AdministriviumIncapacitate {
			# concentration in the water ("conc") and interval ("t") are passed into the system on contact
			# ingested mass ("ate") is passed into the system on ingestion
			# "imass" is initialised
			# "current_load" is initialised to the current mass of contaminant associated with the organism
	
			#acute_lethal = "50% 60[ug/l] @ 200[hours], 50% 122[ug/l] @ 96[hours]"	
			#chronic_lethal = "50% 60[ug/l] @ 400[hours], 50% 122[ug/l] @ 192[hours]"	

			acute_lethal = "none"	(source for data)
			chronic_lethal = "none"	(source for value)

			const { 
				#decay_rate = 0.003	
				decay_rate = 0.03	
				#decay_rate = 0.1	
				#decay_rate = 0.3	
	
				resp_filter_rate = 28[ml/(sec*kg)]	
				resp_uptake = 0.05
				organic_uptake = 0.95	
				stepsize = 1[ug/l]	
			}

			variables {
				parameters {
					volume = "imass * l/kg"
				}
				exposure_rate = "resp_uptake * conc * resp_filter_rate * imass"
			}
	
			load_update = "ate * organic_uptake + volume * ode(dC/dt = exposure_rate  / volume - decay_rate * C(t), C(0) = current_load/volume, t/20, t)"	(Source goes here)
			reproductive_impairment = "current_load>0?0.9:0.0"	

			contaminant_tick = 20[min]
		}
	}
}  

*/

/*-- close up shop, dump data appropriately, free memory */

/*--- Contamination::Shutdown() -- make deathlogger spit out any unsaved data and close files */
int Contamination::Shutdown()
{
	if (!DeathLogger::Shutdown()) return 0;;
	for (int i=0;i<n_cinfo;i++) {
		VERBOSE("Contamination::Shutdown", "ctaxon %s contaminant %s load %g",
			ctaxon, cinfo[i].name, cinfo[i].current_load);
	}
	return 1;
}

/*--- Contamination::free_cinfo() -- free data */
void Contamination::free_cinfo()
{
	assert(cinfo);
	assert(n_cinfo > 0);
	assert(member_cube);

	for (int i=0;i<n_cinfo;i++) {
		assert(cinfo[i].name);
		Free(cinfo[i].name);
		if (cinfo[i].DT) CCalc::FreeCalcVar(cinfo[i].DT);
		if (cinfo[i].imassv) CCalc::FreeCalcVar(cinfo[i].imassv);
		if (cinfo[i].concv) CCalc::FreeCalcVar(cinfo[i].concv);
		if (cinfo[i].atev) CCalc::FreeCalcVar(cinfo[i].atev);
		if (cinfo[i].currentloadv) CCalc::FreeCalcVar(cinfo[i].currentloadv);
	}
	Free(cinfo);
	cinfo = 0;
	n_cinfo = 0;
}


/*-- Contamination::OverrideLocalMembers() -- boolean for presence of a member_cube  */
int Contamination::OverrideLocalMembers()
{
	return member_cube?1:0;
}

/*-- double Contamination::cgetMembers()  -- get the number of members represented */
double Contamination::cgetMembers() {
	if (!member_cube) abort();

	return member_cube->Value();
}

/*-- Contamination::csetMembers()  -- decrease the number of members represented */
void Contamination::csetMembers(double m) {
	if (!member_cube) abort();
	member_cube->AdjustN(cgetMembers() - m, 0);
}


/*-- Contamination::load_LC(char *s, char *tag, EndpointSurf *ES) -- load the contaminant profile for the indicated vulnerable taxon */

int Contamination::load_LC(char *s, char *tag, EndpointSurf *ES) {
	char *points = PGetS(PARAM_OPT, ctaxon, GetCName(CLASS_CONTSINK), "contaminants", s, tag, (char *)0);

	if (!points) {
		// It's ok -- we don't have anything to load.
		VERBOSE("Poisoning", "No contaminant stuff for %s", ctaxon);		
		return 0;
	}
	else {
		char conc0[30] = "", conc1[30] = "";
		char time0[30] = "", time1[30] = "";
		double p0= 0, p1 = 0, c0 = 0, c1 = 0, t0 = 0, t1 = 0;

		int n = sscanf(points, " %lf %% %[^@ ] @ %[^, ], %lf %% %[^@ ] @ %s", &p0, conc0, time0, &p1, conc1, time1);
		if (n != 6) {
			if (!strcasecmp(points,"none")) return 1;
	 
			fatal(1,"The format of a LC string needs to be 'none' or like\n\t40%% 120[mg/l] @ 48[hours], 75%% 150[mg/l] @ 92[hours]\n"
				"I parsed %d things in the string '%s'\n	#1 %f '%s' '%s'\n	#2 %f '%s' '%s'\n", n, points,
				p0, conc0, time0, p1, conc1, time1);
		}
		p0 = p0/100.0;
		p1 = p1/100.0;

		if (strlen(conc0) >= 30) fatal(1,"String too long in contaminant spec");
		if (strlen(conc1) >= 30) fatal(1,"String too long in contaminant spec");
		if (strlen(time0) >= 30) fatal(1,"String too long in contaminant spec");
		if (strlen(time1) >= 30) fatal(1,"String too long in contaminant spec");

		CCalc cc;
		c0 = cc.UnitEvaluate(conc0);
		t0 = cc.UnitEvaluate(time0);
		c1 = cc.UnitEvaluate(conc1);
		t1 = cc.UnitEvaluate(time1);
		if (!ES->SetSurface(p0,c0,t0,p1,c1,t1)) abort();
	}
	return 1;
}


/*-- Contamination::ZapContaminantSetup(int i) -- get rid of all the machinery */

void Contamination::ZapContaminantSetup(int i) {
	assert(i >= 0 && i < n_cinfo);
	
//	cinfo[i].vbid = -1;

	cinfo[i].update.string = 0;
	cinfo[i].update.id = -1;

	cinfo[i].reproduce.string = 0;
	cinfo[i].reproduce.id = -1;

	cinfo[i].forage.string = 0;
	cinfo[i].forage.id = -1;

	cinfo[i].move.string = 0;
	cinfo[i].move.id = -1;

	if (cinfo[i].DT) cinfo[i].DT = CCalc::FreeCalcVar(cinfo[i].DT);
	if (cinfo[i].concv) cinfo[i].concv = CCalc::FreeCalcVar(cinfo[i].concv);
	if (cinfo[i].imassv) cinfo[i].imassv = CCalc::FreeCalcVar(cinfo[i].imassv);
	if (cinfo[i].atev) cinfo[i].atev = CCalc::FreeCalcVar(cinfo[i].atev);
	if (cinfo[i].currentloadv) cinfo[i].currentloadv = CCalc::FreeCalcVar(cinfo[i].currentloadv);

	cinfo[i].tick = DNaN;
	cinfo[i].conc = DNaN;
	cinfo[i].ate = 0;
}


/*-- Contamination::ContaminantSetup(char *s, int i) -- set up data for vulnerable taxa using Load_LC */
int Contamination::ContaminantSetup(char *s, int i) {
	cinfo[i].vbid = PrmEnvExpr::LoadBigBlock(cinfo[i].vbid, ctaxon, GetCName(CLASS_CONTSINK), "contaminants", s, (char*)0);

	if (!PrmEnvExpr::EnvBlockOk(cinfo[i].vbid)) {
		VERBOSE("ContaminationContaminantInit", "Missing environment agent");
//	 free_cinfo();
		return 0;
	}
	
	RCCalc *cc = PrmEnvExpr::GetCCalc(cinfo[i].vbid);
	assert(cc);
	
	if (!cinfo[i].name) {
		cinfo[i].name = Strdup(s);
		cinfo[i].current_load = 0;
	}
	else if (cinfo[i].name != s) { 
		abort();
	}
	
	// Get parameters from the parameterisation corpus
	cinfo[i].cont_tick = PGetN(PARAM_NOR|PARAM_REQ, ctaxon, GetCName(CLASS_CONTSINK), "contaminants", s, "contaminant_tick", (char *)0);
	cinfo[i].update.string = PGetS(PARAM_REQ, ctaxon, GetCName(CLASS_CONTSINK), "contaminants", s, "load_update", (char *)0);

	cinfo[i].reproduce.string = PGetS(PARAM_OPT, ctaxon, GetCName(CLASS_CONTSINK), "contaminants", s, "reproductive_impairment", (char *)0);
	//if (!cinfo[i].reproduce.string) cinfo[i].reproduce.string = "0";

	cinfo[i].forage.string = PGetS(PARAM_OPT, ctaxon, GetCName(CLASS_CONTSINK), "contaminants", s, "foraging_impairment", (char *)0);
	//if (!cinfo[i].forage.string) cinfo[i].forage.string = "0";

	cinfo[i].move.string = PGetS(PARAM_OPT, ctaxon, GetCName(CLASS_CONTSINK), "contaminants", s, "movement_impairment", (char *)0);
	//if (!cinfo[i].move.string) cinfo[i].move.string = "0";
	
	// Load the LC% data here
	int caught = 0;
	caught += load_LC(s, "acute_lethal", &cinfo[i].acute_lethal);
	caught += load_LC(s, "chronic_lethal", &cinfo[i].chronic_lethal);
	caught += load_LC(s, "reproduction", &cinfo[i].reproduction);
	caught += load_LC(s, "movement", &cinfo[i].movement);
	caught += load_LC(s, "foraging", &cinfo[i].foraging);
	if (!caught) fatal(1,"You didn't specify a response for %s to contaminant %s", ctaxon, s);
	else VERBOSE("Poisoning", "%s is sensitive to %d different pathologies for %s", ctaxon, caught, s);

	
	// Specific impairments
	cinfo[i].update.id = cc->AddProgram(cinfo[i].update.string);
	if (cinfo[i].forage.string) cinfo[i].forage.id = cc->AddProgram(cinfo[i].forage.string);
	if (cinfo[i].move.string) cinfo[i].move.id = cc->AddProgram(cinfo[i].move.string);
	if (cinfo[i].reproduce.string) cinfo[i].reproduce.id = cc->AddProgram(cinfo[i].reproduce.string);

	assert(cinfo[i].update.id >= 0);

	// and initialise the rest
	cinfo[i].DT = cc->GetVarRef2("dt");
	cinfo[i].concv = cc->GetVarRef2("conc");
	cinfo[i].imassv = cc->GetVarRef2("imass");
	cinfo[i].atev = cc->GetVarRef2("ate");
	cinfo[i].currentloadv = cc->GetVarRef2("current_load");
	cinfo[i].tick = 0;
	cinfo[i].conc = 0;
	cinfo[i].ate = 0;

	return 1;
}


/*-- Contamination::init_contaminant_stuff() -- gets the machinery set for the more detailed routines */
int Contamination::init_contaminant_stuff() {
	int Nnew;
	int i, j, k;
	char *s;

	if (!contaminants) abort();
	Nnew = contaminants->NumInterest();
	if (!Nnew) {
		warning("%s is not interested in any contaminants", ctaxon);
		return 1;
	}

	if (n_cinfo > 0 && n_cinfo != Nnew) fatal(1,"Contaminant lists must remain constant through ctaxon changes");
	else n_cinfo = Nnew;

	if (!member_cube) {// populations & schools, y'know
		member_cube = new Cube(n_cinfo+1, PgetMembers());
		if (!member_cube) abort();
	}

	for (j = 0; cinfo && j < n_cinfo; j++) {
		k = 1;
		s = contaminants->GetInterest(j);
		for (i = 0; k && i < n_cinfo; i++) {
			if (!strcmp(cinfo[i].name, s)) k = 0;
		}	 
		if (k) fatal(1,"Contaminant lists must really remain constant through ctaxon changes");
	}

	if (!cinfo) {
		cinfo = (_cinfo*)Calloc(n_cinfo, sizeof(_cinfo));
		if (!cinfo) abort();
		
		for (i = 0; i < n_cinfo; i++) {
			cinfo[i].vbid = -1;
			cinfo[i].name = 0;
			cinfo[i].current_load = 0;
		}
	}

	for (i = 0; i < n_cinfo; i++) {
		ZapContaminantSetup(i);

		if (cinfo[i].name) s = cinfo[i].name;
		else s = contaminants->GetInterest(j);

		ContaminantSetup(s, i);
	}

	PsetMembers(DNaN);

	return 1;
}

/* Contamination::Init(char *taxon, char *name) -- initialise the contaminant data for a given taxon (the agent name is copied to the cname) */
int Contamination::Init(char *taxon, char *name)
{
	assert(taxon && *taxon);
	if (ctaxon) Free(ctaxon);
	ctaxon = Strdup(taxon);
	assert(name);
	if (cname) Free(cname);
	cname = Strdup(name);
	if (!DeathLogger::Init(ctaxon,cname)) return 0;;
	if (!ContaminantSink::Init(ctaxon)) return 0;
	if (!init_contaminant_stuff()) return 0;

	// Do we really need a get big block thing here? ... Worth thinking about

	return 1;
}

/*-- Contamination::Isa(int mask) -- queries role of agent */
int Contamination::Isa(int mask)
{
	return ContaminantSink::Isa(mask);
}

/*-- indicates compatibility with kernel settings and capabilities */
int Contamination::Compatibility()
{
	return ContaminantSink::Compatibility() & (ATTR|RESET|REINIT|STATE);
}

// Usually no need to touch Intoxicate
//double Contamination::Intoxicate(double t, double dt)
//{
//	// enable the intoxication code
//	return ContaminantSink::Intoxicate(t, dt);
//}


/*-- CommitIntoxicate(double t, double dt, double actual_dt) --  Commit any intoxication post behaviour */
int Contamination::CommitIntoxicate(double t, double dt, double actual_dt)
{
	if (!n_cinfo) return 1; // nothing to commit

	if (isnan(actual_dt)) { // Oops, we've popped our cogs.
		return 1;
	}

	int i;
	double new_load = 0;
	double k, *K = 0;
	double old_members = member_cube->Value();

	K = (double *)Malloc(n_cinfo * sizeof(*K));
	if (!K) abort();

	// Collect K
	for (i = 0, k = 0; i < n_cinfo; i++) {
		RCCalc *cc = PrmEnvExpr::GetCCalc(cinfo[i].vbid);
		assert(cc);

		K[i] = cinfo[i].acute_lethal.value(cinfo[i].conc, actual_dt);
		if (K[i] > 0) {
			VERBOSE("Poisoning", "%s conc = %f, load = %f K = %f", cinfo[i].name, cinfo[i].conc, cinfo[i].current_load, K[i]);
		}
		k += K[i];

		cc->SetVarRef2(cinfo[i].imassv, getIMass());
		cc->SetVarRef2(cinfo[i].DT, actual_dt);

		cc->SetVarRef2(cinfo[i].concv, cinfo[i].conc);
		cc->SetVarRef2(cinfo[i].atev, cinfo[i].ate); // We aren't eating t, stuff
		cc->SetVarRef2(cinfo[i].currentloadv, cinfo[i].current_load);

		/* get environment info */
		PrmEnvExpr::Configure(t, cinfo[i].vbid);
		PrmEnvExpr::ValidateVariables(cinfo[i].vbid);

		new_load = cc->Calculate(cinfo[i].update.id); // update load level
		VERBOSE("CommitIntoxicate", "%s %f -> %f  conc = %f dt = %f imass = %f ate = %f", 
			cinfo[i].name, cinfo[i].current_load, new_load, 
			cinfo[i].conc, actual_dt, 
			getIMass(), cinfo[i].ate);
		cinfo[i].current_load = new_load; // change load for contaminant
	
		for (int iq = 0; profile && iq < profile->N; iq++) {
			if (!strcmp(profile->c_list[iq].name, cinfo[i].name)) {
				profile->c_list[iq].mass = new_load;
				break;
			}
		}
	}
	 
	if (k > 0) { // This is mostly pertinent for populations and schools
		// Adjust acute mortality based on water concentration here

		double dk = cgetMembers(), ddk;
		member_cube->AdjustLevels(K, 1, n_cinfo);
		ddk = cgetMembers();

		LogDeath(t, dk - ddk, ddk, getIMass(), "AcutePoisoning");
	}
	
	// Adjust chronic mortality based on tissue load here
	for (k = 0, i = 0; i < n_cinfo; i++) {
		K[i] = cinfo[i].chronic_lethal.value(cinfo[i].current_load, actual_dt);
		k += K[i];
	}

	// Change the levels in the member_cube to reflect any mortality we've inflicted
	if (k > 0) {
		double dk = cgetMembers(), ddk;
		member_cube->AdjustLevels(K, 1, n_cinfo);
		ddk = cgetMembers();

		LogDeath(t, dk - ddk, ddk, getIMass(), "ChronicPoisoning");
	}

#if defined(MAINTAIN_THINGS_MEMBERS)
	PsetMembers(cgetMembers());
#endif

	
	// and don't carry *this* lot of contaminant across to the next iteration
	for (i = 0; i < n_cinfo; i++) {
		cinfo[i].conc = 0;
		cinfo[i].ate = 0;
	}

	if (cgetMembers() - old_members < 0) {
		VERBOSE("Poisoning", "%f %s died due to contaminants", old_members - cgetMembers(), ctaxon);
	}

	Free(K);

	return 1; // For now we'll say it worked
}

/*-- Contaminantion::LocalIntoxicate(agent, t, dt, contaminant) -- An individual has been hit */
// We're about to get nuked by something
// Note that dt is an estimate and the intoxication may need to be adjusted
// in CommitIntoxicate if it is used
double Contamination::LocalIntoxicate(int agent, double t, double dt, char *contaminant)
{
	VERBOSE("LocalIntoxicate", "Intoxicating %s", contaminant);
	int cx = -1;
	
	for (int i = 0; i < n_cinfo; i++) {
		if (!strcmp(cinfo[i].name, contaminant)) {
			cx = i;
			break;
		}
	}

	assert(cx >= 0);
	assert(KISA(agent, CLASS_CONTSRC));

	// Get the contaminant index
	int cidx = ContaminantSource::GetCSNum(KID(agent), contaminant);
	assert(cidx >= 0);
	// Now get the value
	double d = ContaminantSource::GetCSValue(KID(agent), 
		t, getLocation(), cidx);
	if (isnan(d)) return dt; // not applicable

	cinfo[cx].tick = Min(cinfo[cx].cont_tick, dt);
	cinfo[cx].conc = Max(cinfo[cx].conc, d);
//	cinfo[cx].ate = 0;

	return cinfo[cx].tick;
}


/*-- Contamination::Ingest(char *contaminant, double mass, double t) -- contamination by eating something a bit funny */
// We're about to get nuked by something
// Note that dt is an estimate and the intoxication may need to be adjusted
// in CommitIntoxicate if it is used
int Contamination::Ingest(char *contaminant, double mass, double t)
{
	VERBOSE("Ingest", "Ingesting %s", contaminant);
	int cx = -1;
	
	assert(mass > 0);

	for (int i = 0; i < n_cinfo; i++) {
		if (!strcmp(cinfo[i].name, contaminant)) {
			cx = i;
			break;
		}
	}
	if (cx <0) return 0;

	assert(cx >= 0);

	cinfo[cx].ate += mass;
	return 1;
}

/*-- Contamination::getReproductiveImpairment(double t) -- service routine */
double Contamination::getReproductiveImpairment(double t) {
	double d = 1.0;
	double v = 0;

	for (int i = 0; i < n_cinfo; i++) {
		if (!cinfo[i].reproduce.string) continue;

		RCCalc *cc = PrmEnvExpr::GetCCalc(cinfo[i].vbid);
		assert(cc);

		cc->SetVarRef2(cinfo[i].imassv, getIMass());
		cc->SetVarRef2(cinfo[i].currentloadv, cinfo[i].current_load);

		PrmEnvExpr::Configure(t, cinfo[i].vbid);
		PrmEnvExpr::ValidateVariables(cinfo[i].vbid);

		v = cc->Calculate(cinfo[i].reproduce.id);

		d *= (1.0 - v);
	}

	return 1.0 - d;
}

/*-- Contamination::getForagingImpairment(double t) -- service routine */
double Contamination::getForagingImpairment(double t) {
	double d = 1.0;
	double v = 0;

	for (int i = 0; i < n_cinfo; i++) {
		if (!cinfo[i].forage.string) continue;

		RCCalc *cc = PrmEnvExpr::GetCCalc(cinfo[i].vbid);
		assert(cc);

		cc->SetVarRef2(cinfo[i].imassv, getIMass());
		cc->SetVarRef2(cinfo[i].currentloadv, cinfo[i].current_load);

		PrmEnvExpr::Configure(t, cinfo[i].vbid);
		PrmEnvExpr::ValidateVariables(cinfo[i].vbid);

		v = cc->Calculate(cinfo[i].forage.id);

		d *= (1.0 - v);
	}

	return 1.0 - d;
}

/*-- Contamination::getMovementImpairment(double t) -- service routine */
double Contamination::getMovementImpairment(double t) {
	double d = 1.0;
	double v = 0;

	for (int i = 0; i < n_cinfo; i++) {
		if (!cinfo[i].move.string) continue;

		RCCalc *cc = PrmEnvExpr::GetCCalc(cinfo[i].vbid);
		assert(cc);

		cc->SetVarRef2(cinfo[i].imassv, getIMass());
		cc->SetVarRef2(cinfo[i].currentloadv, cinfo[i].current_load);

		PrmEnvExpr::Configure(t, cinfo[i].vbid);
		PrmEnvExpr::ValidateVariables(cinfo[i].vbid);

		v = cc->Calculate(cinfo[i].move.id);

		d *= (1.0 - v);
	}

	return 1.0 - d;
}

/*-  The End  */
