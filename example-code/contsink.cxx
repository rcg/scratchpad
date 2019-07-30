#include <stdlib.h>
#include <string.h>
#include "contsink.hxx"
#include "contsrc.hxx"
#include "memchk.h"

/* 
  2019-07-30-07:21:00 -- minor explanatory note

  This class provides an agent which inherits it the ability to be
  recognised as a source of contaminants.  Most of the file consists
  of the serialisation code it needs to run with the kernels.  The
  KID2() call you see in GetProfile() supports the interface to the
  multitasking kernels and the KID() call in Intoxicate() is a call to
  whichever kernel is in control at the time.

  It does provide the Intoxicate() call, which is how the contaminants
  are actually introduced into an agent.

*/
/*-- ContaminantSink() --  */
ContaminantSink::ContaminantSink() {
	taxname = 0;
	contaminants = 0;
	profile = 0;
}
/*-- ~ContaminantSink() --  */
ContaminantSink::~ContaminantSink() {
	if (taxname) Free(taxname);
	if (contaminants) delete contaminants;
	if (profile) delete profile;
}
/*-- Reset() --  */
void ContaminantSink::Reset() {
	PrmAgent::Reset();
	contaminants = 0;
	profile = 0;
	taxname = 0;
}
/*-- Init(char *taxon) --  */
int ContaminantSink::Init(char *taxon) {
	if (taxname) Free(taxname);
	taxname = 0;
	if (!taxon) return 0;
	taxname = Strdup(taxon);

//#warning also need to load up parameter based stuff here
	if (contaminants) contaminants->ClearInterests();
	else contaminants = new ContaminantList();

	if (profile) delete profile;
	profile = 0;

	if (PGetI(0, taxon, GetCName(CLASS_CONTSINK),	"contaminants", "disable", (char *)0)) return 1;

	char **clist = PGetNodes(taxon, GetCName(CLASS_CONTSINK),"contaminants", (char *)0);
	if (!clist) return 1; // nothing to do

	profile = new ContaminantProfile();
	if (!profile) abort();
	
	for (int i=0;clist[i];i++) {
		VERBOSE("ContaminantSink", "Init %s", clist[i]);
		contaminants->RegisterInterest(clist[i]);
		profile->AddContaminant(clist[i], DNaN);
	}
	Free(clist);
	assert(profile->N > 0);
	return 1;
}
/*-- Isa(int mask) --  */
int ContaminantSink::Isa(int mask) {
	if ((mask&CLASS_MASK) == CLASS_CONTSINK) return 1;
	return PrmAgent::Isa(mask);
}
/*-- Compatibility() --  */
int ContaminantSink::Compatibility() {
	return (ATTR|RESET|STATE|REINIT);
}
/*-- GetState(int *sz) --  */
void *ContaminantSink::GetState(int *sz) {
	void *v[3];
	int l[3];

	if (taxname) {
		v[0] = taxname;
		l[0] = strlen(taxname)+1;
	}
	else v[0] = 0, l[0] = 0;

	if (profile) v[1] = profile->GetState(&l[1]);
	else v[1] = 0, l[1] = 0;

	if (contaminants) v[2] = contaminants->GetState(&l[2]);
	else v[2] = 0, l[2] = 0;

	void *d = pack_mem(v, l, 3, sz);

	if (v[1]) Free(v[1]);
	if (v[2]) Free(v[2]);
	return d;
}
/*-- SetState(void *d, int sz) --  */
void ContaminantSink::SetState(void *d, int sz) {
	if (taxname) Free(taxname);
	if (profile) delete profile;
	if (contaminants) delete contaminants;
	taxname = 0;
	profile = 0;
	contaminants = 0;

	void **v;
	int *l;

	int n = unpack_mem(d, sz, &v, &l);
	assert(n == 3);
	if (v[0]) taxname = Strdup((char*)v[0]);
	if (v[1]) profile = new ContaminantProfile(v[1], l[1]);
	if (v[2]) {
		contaminants = new ContaminantList();
		contaminants->SetState(v[2], l[2]);
	}
	Free(v); Free(l);
}
/*-- Get(int attribute, void *args, int args_size,	void *data, int *size) */
void *ContaminantSink::Get(int attribute, void *args, int args_size,	void *data, int *size) {
	if (PrmAgent::Isa(attribute)) 
		return PrmAgent::Get(attribute, args, args_size, data, size);
	switch(attribute) {
	case ATTR_CONTSINK_PROFILE:
		if (!profile) {
			assert(size);
			*size = 0;
			return 0;
		}
		assert(!args);
		assert(!args_size);
		return profile->GetState(size);
	}
	abort();
}

double ContaminantSink::Intoxicate(double t, double dt) {
	// figure out which agents will intoxicate us and
	// call LocalIntoxicate for each one

	if (!contaminants) return dt;
	int num;
	int *ia = FindAgentsByClass(CLASS_CONTSRC, &num);
	if (!ia) return dt;

	for (int i=0;i<num;i++) {
		int nc = contaminants->NumInterest();
		for (int j=0;j<nc;j++) {
			char *c = contaminants->GetInterest(j);
			if (ContaminantSource::GetCSNum(KID(ia[i]), c) >= 0)
				dt = LocalIntoxicate(ia[i], t, dt, c);
		}
	}
	Free(ia);
	return dt;
}

int ContaminantSink::SetProfile(ContaminantProfile *p) {
	assert(p);
	if (profile) delete profile;
	profile = new ContaminantProfile(p);
	if (!profile) abort();
	return 1;
}

ContaminantProfile *ContaminantSink::GetProfile(KID2(xid)) {
#ifdef PRODUCTION_KERNEL
	return PKDACCESS(ContaminantSink,xid)getProfile();
#else
	int sz = 0;
	void *v = KGET(xid, ATTR_CONTSINK_PROFILE, 0, 0, 0, &sz);
	if (!v) abort();
	ContaminantProfile *p = new ContaminantProfile(v, sz);
	Free(v);
	return p;
#endif
}

ContaminantProfile *ContaminantSink::getProfile()
{
	return new ContaminantProfile(profile);
}

