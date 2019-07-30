#include <stdlib.h>
#include <string.h>
#include "contsrc.hxx"
#include "memchk.h"

/* 
  2019-07-30-07:21:00 -- minor explanatory note

  This class really only provides an agent which inherits it the
  ability to be recognised as a source of contaminants.  Most of the
  file is devoted to the serialisation code and the mechanisms that 
  allow it to run with the kernels.

*/

ContaminantSource::ContaminantSource()
{
	taxname = 0;
	contaminants = 0;
}

ContaminantSource::~ContaminantSource()
{
	if (taxname) Free(taxname);
	if (contaminants) delete contaminants;
}

void ContaminantSource::Reset()
{
	PrmAgent::Reset();
	taxname = 0;
	contaminants = 0;
}

int ContaminantSource::Init(char *taxon)
{
	if (taxname) Free(taxname);
	taxname = 0;
	if (!taxon) return 0;
	taxname = Strdup(taxon);

//#warning also need to load up parameter based stuff here
	if (contaminants) contaminants->ClearSources();
	else contaminants = new ContaminantList();

	char **clist = PGetNodes(taxon, GetCName(CLASS_CONTSRC), "contaminants", (char *)0);
	if (!clist) return 1; // nothing to do
	for (int i=0;clist[i];i++) {
		VERBOSE("ContaminantSource", "Init %s", clist[i]);
		contaminants->RegisterAsContaminantSource(clist[i]);
		
	}
	Free(clist);
	return 1;
}


int ContaminantSource::Isa(int mask)
{
	if ((mask&CLASS_MASK) == CLASS_CONTSRC) return 1;
	return PrmAgent::Isa(mask);
}

void *ContaminantSource::Get(int attribute, void *args, int args_size,
	void *data, int *size)
{
	double d;
	if ((attribute&CLASS_MASK) != CLASS_CONTSRC)
		return PrmAgent::Get(attribute, args, args_size, data, size);
	switch(attribute) {
	case ATTR_CONTSRC_CSVALUE: {
		assert(args);
		assert(args_size > 0);
		void **v;
		int *l;
		int n = unpack_mem(args, args_size, &v, &l);
		assert(n == 3);
		assert(v[0]);
		assert(v[1]);
		assert(v[2]);
		assert(l[0] == sizeof(double));
		assert(l[1] == sizeof(R3));
		assert(l[2] == sizeof(int));
		d = getCSValue(*(double*)v[0], *(R3*)v[1], *(int*)v[2]);
		Free(v); Free(l);
		return GetReturn(data, size, &d, sizeof(d));
	}
	case ATTR_CONTSRC_CSNUM: {
		assert(args);
		assert(args_size > 0);
		int i = getCSNum((char*)args);
		return GetReturn(data, size, &i, sizeof(i));
	}
	case ATTR_CONTSRC_CSID:
		abort();	// how did we get here?

	}
	fatal(1, "Unknown attribute 0x%08x", attribute);
}

int ContaminantSource::IsSource(KID2(xid), char *contaminant)
{
	return (ContaminantSource::GetCSNum(KID(xid), contaminant) >= 0);
}

int ContaminantSource::GetCSNum(KID2(xid), char *contaminant)
{
	assert(contaminant);
#ifdef PRODUCTION_KERNEL
	return PKDACCESS(ContaminantSource,xid)getCSNum(contaminant);
#else
	int i, sz = sizeof(int);
	if (!KGET(xid, ATTR_CONTSRC_CSNUM, 
			(void*)contaminant, strlen(contaminant)+1, &i, &sz)) abort();
	return i;
#endif
}

double ContaminantSource::GetCSValue(KID2(xid), double t, R3 loc, int cid)
{
#ifdef PRODUCTION_KERNEL
	return PKDACCESS(ContaminantSource,xid)getCSValue(t, loc, cid);
#else
	double d;
	// pack up arguments and send through to agent
	void *v[3], *data;
	int l[3], sz, dsz;
	v[0] = &t; l[0] = sizeof(t);
	v[1] = &loc; l[1] = sizeof(loc);
	v[2] = &cid; l[2] = sizeof(cid);
	data = pack_mem (v, l, 3, &sz);
	dsz = sizeof(double);
	if (!KGET(xid, ATTR_CONTSRC_CSVALUE, data, sz, &d, &dsz)) abort();
	Free(data);
	return d;
#endif
}

int ContaminantSource::getCSNum(char* contaminant)
{
	assert(contaminants);
	assert(contaminant);
	for (int i=0;i<contaminants->NumSource();i++) {
		char *str = contaminants->GetSource(i);
		if (!strcmp(str, contaminant)) {
			return i;
		}
	}
	return -1;
}

