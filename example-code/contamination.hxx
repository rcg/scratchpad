// -*- outline-regexp: "/\\*-+";  -*-
/*-  Identification and Changes  */

/*
  contamination.hxx -- Written by Randall Gray 
  Initial coding: 
  Date: 2007.10.22
  Location: floyd.hba.marine.csiro.au:/home/gray/Projects/InVitro/contamination.hxx

  History:

  $Log: contamination.hxx,v $
  Revision 1.10  2008/10/06 03:47:14  gray
  Not this one.

  Revision 1.9  2008/08/29 04:11:21  gray
  More messing with indenting!

  Revision 1.8  2008/08/28 05:55:05  gray
  More reformatting.  Closer ... I think.	More testing is needed

  Revision 1.7  2008/08/26 08:21:39  gray
  Massive mods.  Lotsa space fixing and such. Merged Beths stuff in

  Revision 1.6  2008/04/30 00:50:16  gray
  Moved repository to Njal, converted from "se-model" to "invitro"

  Revision 1.5  2008/02/21 04:21:26  gray
  more contaminant propagating

  Revision 1.4  2008/02/20 02:42:10  gray
  Dunno what we were doing, but it needs to go in

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

/*-  Configuration stuff  */
#pragma interface

#ifndef __contamination_hxx
#define in_contamination_hxx
#define __contamination_hxx

#include "prmenvexpr.hxx"
#include "contsink.hxx"
#include "cube.hxx"
#include "endpointsurf.hxx"
#include "deathlogger.hxx"


class Contamination: virtual public PrmEnvExpr, virtual public ContaminantSink,
	virtual public DeathLogger
{
public:
	Contamination();
	virtual ~Contamination();
	virtual int Init(char *tax, char *nm);
	virtual int Isa(int mask);
	virtual int Compatibility();
	virtual int Shutdown();

	virtual void *GetState(int*);
	virtual void SetState(void*, int);
	virtual void Reset();
	virtual int ReInit(int);

	double Level(char *name);
	

protected:
	virtual int load_LC(char*, char*, EndpointSurf*);
	virtual void ZapContaminantSetup(int i);
	virtual int ContaminantSetup(char *name, int i);

	virtual int CommitIntoxicate(double t, double dt, double actual_dt);
	//virtual double Intoxicate(double t, double dt);
	virtual double LocalIntoxicate(int agentid, double t, double dt, char* contaminant);
	virtual int Ingest(char* contaminant, double mass, double t);

	virtual void PsetMembers(double m)=0;
	virtual double PgetMembers()=0;
	void csetMembers(double m);
	double cgetMembers();
	int OverrideLocalMembers();

	virtual int init_contaminant_stuff();

	virtual double getReproductiveImpairment(double t);
	virtual double getForagingImpairment(double t);
	virtual double getMovementImpairment(double t);

	Cube *member_cube;
	char *ctaxon, *cname;

	typedef struct {
		EndpointSurf acute_lethal, chronic_lethal, foraging, reproduction, movement;

		double tick, conc, ate;
		double current_load;
		double cont_tick;

		struct {
			int id;
			char *string;
		} update, forage, reproduce, move;

		int vbid;
		CCalc::CalcVar *concv, *imassv, *atev, *currentloadv, *DT;
		char *name;
	} _cinfo;

	_cinfo *cinfo;
	int n_cinfo;

private:
	void zero();
	void free_cinfo();
	void *Get_cinfo_State(int, int*);
	void Set_cinfo_State(void*, int, int);

Attribute:
	virtual R3 getLocation()=0;
	virtual double getIMass()=0;

};

#endif
/*-  The End  */

