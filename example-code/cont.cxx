// -*- outline-regexp: "/\\*-+";  -*-
/*-  Identification and Changes  */

/*
  cont.cxx -- Written by Randall Gray 
  Initial coding: 
  Date: 2007.10.22
  Location: floyd.hba.marine.csiro.au:/home/gray/Projects/InVitro/contamination.cxx

  This file is responsible for much of the communication between the
  components that deal with contaminants.



  2019-07-30-07:21:00 -- minor explanatory note

  The "pack" and "unpack" routines in this file were there so we could
  pass agents between different operating environments.  These
  routines were implemented for *every* agent.

  One of the side benefits was that we could save snapshots of the
  model during execution and resume them.  This was useful in bug
  tracking, comparing divergent scenarios (we could start from a
  common point), and static analysis of the state of the simulated
  system.

  You may notice the comment lines which precede most functions; these
  are structured so the "outline" package in emacs can give me just
  the bits I am working on, rather than the whole file.

*/


/*-- Includes, definitions and variables */
/*--- includes */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "cont.hxx"
#include "packmem.h"
#include "memchk.h"


/*--- definitions */
#define NO_STRDUP  

/*--- variables */


/*-- Constructors/destructors */

/*--- ContaminantList() */
ContaminantList::ContaminantList() {
	source = new StringTable();
	interest = new StringTable();
	if (!source || !interest) abort();
}

/*--- ~ContaminantList() */
ContaminantList::~ContaminantList() {
	delete source;
	delete interest;
}

/*-- RegisterAsContaminantSource(char *contaminant) */
int ContaminantList::RegisterAsContaminantSource(char *contaminant) {
	return source->Insert(contaminant, "0");
}

/*-- RegisterAsContaminantSource(char *contaminant) -- a potential sink*/
int ContaminantList::RegisterInterest(char *contaminant) {
	return interest->Insert(contaminant, "0");
}

/*-- IsSource(char *contaminant) -- predicate */
int ContaminantList::IsSource(char *contaminant) {
	if (source->GetValue(contaminant)) return 1;
	return 0;
}

/*-- IsInterested(char *contaminant) -- predicate, potential sink/logger */
int ContaminantList::IsInterested(char *contaminant) {
	if (interest->GetValue(contaminant)) return 1;
	return 0;
}

/*-- NumSource() --  return the number of sources for contact*/
int ContaminantList::NumSource() {
	return source->Num();
}

/*-- NumInterest() -- return the number of entities interested in contaminants */
int ContaminantList::NumInterest() {
	return interest->Num();
}

// we can use NO_STRDUP to speed things a little (comparing pointers
// rather than strings) at the cost of reduced flexibility and potentially
// "missing" contaminants
//
// TESTED AND OK FOR PRODUCTION

ar *ContaminantList::GetInterest(int rec) {  // return key for contaminant number "rec"
#ifdef NO_STRDUP
	return interest->GetKey(rec);
#else
	return Strdup(interest->GetKey(rec));
#endif
}

/*-- GetSource(int rec) --  return the actual source  */
char *ContaminantList::GetSource(int rec) {
#ifdef NO_STRDUP
	return source->GetKey(rec);
#else
	return Strdup(source->GetKey(rec));
#endif
}


/*-- GetState(int *sz) -- The following two routines serialise / unserialise for migration/rollback */
void *ContaminantList::GetState(int *sz) {
	void *v[2], *d;
	int l[2];
	assert(sz);
	v[0] = source->GetState(&l[0]);
	v[1] = interest->GetState(&l[1]);
	d =  pack_mem(v, l, 2, sz);
	if (v[0]) Free(v[0]);
	if (v[1]) Free(v[1]);
	return d;
}
/*-- SetState(void *d, int sz) --  */
void ContaminantList::SetState(void *d, int sz) {
	void **v;
	int *l, n;
	n = unpack_mem(d, sz, &v, &l);
	assert(n == 2);
	source->SetState(v[0], l[0]);
	interest->SetState(v[1], l[1]);
	Free(v);
	Free(l);
}


/*-- ClearInterests() --  */
void ContaminantList::ClearInterests() {
	assert(interest);
	delete interest;
	interest = new StringTable();
}
/*-- ClearSources() --  */
void ContaminantList::ClearSources() {
	assert(source);
	delete source;
	source = new StringTable();
}
/*-- ContaminantProfile() -- (Re-)Initialise */
ContaminantProfile::ContaminantProfile() {
	N = 0;
	c_list = 0;
}


/*-- ContaminantProfile(ContaminantProfile *c) --  Copy profile  */
ContaminantProfile::ContaminantProfile(ContaminantProfile *c) {
	assert(c);
	N = c->N;
	c_list = (Contaminant*)Calloc(N, sizeof(Contaminant));
	if (!c_list) abort();
	for(int i=0;i<N;i++) {
		assert(c->c_list[i].name);
		c_list[i].name = Strdup(c->c_list[i].name);
		if (!c_list[i].name) abort();
		c_list[i].mass = c->c_list[i].mass;
	}
}


/*-- ~ContaminantProfile() -- dispose of contaminants (OHS is watching!) */
ContaminantProfile::~ContaminantProfile() {
	for(int i=0;i<N;i++) {
		assert(c_list);
		assert(c_list[i].name);
		Free(c_list[i].name);
	}
	if (c_list) Free(c_list);
}


/*-- unpack_struct(int i, void *d, int sz) -- Used in serialisation/extraction  */
void ContaminantProfile::unpack_struct(int i, void *d, int sz) {
	void **v;
	int *l;
	int n = unpack_mem(d, sz, &v, &l);
	assert(n == 2);
	assert(v[0]);
	assert(l[0] == sizeof(double));
	assert(v[1]);
	assert(l[1] > 1);

	assert(c_list);
	c_list[i].mass = *(double*)v[0];
	c_list[i].name = Strdup((char*)v[1]);
	if (!c_list[i].name) abort();
	Free(v); Free(l);
}

/*-- unpack_list(void *d, int sz) -- Used in serialisation/extraction */
void ContaminantProfile::unpack_list(void *d, int sz) {
	void **v;
	int *l;

	assert(N > 0);
	int n = unpack_mem(d, sz, &v, &l);
	assert(n == N);
	c_list = (Contaminant*)Calloc(N, sizeof(Contaminant));
	if (!c_list) abort();
	for(int i=0;i<N;i++) {
		unpack_struct(i, v[i], l[i]);
	}
	Free(v); Free(l);
}



/*-- ContaminantProfile(void *d, int sz) -- Used in serialisation/extraction */
ContaminantProfile::ContaminantProfile(void *d, int sz) {
	void **v;
	int *l;

	int n = unpack_mem(d, sz, &v, &l);
	assert(n == 2);
	assert(v[0]);
	assert(l[0] == sizeof(int));
	N = *(int*)v[0];
	if (v[1]) {
		assert(N > 0);
		unpack_list(v[1], l[1]);
	}
	else {
		assert(N == 0);
		c_list = 0;
	}
	Free(v); Free(l);
}


/*-- pack_struct(int i, int *sz) -- Used in serialisation/extraction */
void *ContaminantProfile::pack_struct(int i, int *sz) {
	void *v[2];
	int l[2];
	v[0] = &c_list[i].mass; l[0] = sizeof(c_list[i].mass);
	assert(c_list[i].name && (strlen(c_list[i].name) > 0));
	v[1] = c_list[i].name; l[1] = strlen(c_list[i].name)+1;
	return pack_mem(v, l, 2, sz);
}


/*-- pack_list(int *sz) -- Used in serialisation/extraction */
void *ContaminantProfile::pack_list(int *sz) {
	if (!N) {
		*sz = 0;
		return 0;
	}
	void **v = (void**)Calloc(N, sizeof(void*));
	int *l = (int*)Calloc(N, sizeof(int));
	if (!v || !l) abort();
	for(int i=0;i<N;i++) {
		v[i] = pack_struct(i, &l[i]);
	}
	void *d = pack_mem(v, l, N, sz);
	for(int i=0;i<N;i++) if (v[i]) Free(v[i]);
	Free(v); Free(l);
	return d;
}


/*-- GetState(int *sz) -- Used in serialisation/extraction */
void *ContaminantProfile::GetState(int *sz) {
	void *v[2];
	int l[2];
	v[0] = &N; l[0] = sizeof(N);
	v[1] = pack_list(&l[1]);

	void *d = pack_mem(v, l, 2, sz);
	if (v[1]) Free(v[1]);
	return d;
}

/*-- AddContaminant(char *name, double mass) -- */
int ContaminantProfile::AddContaminant(char *name, double mass) {
	assert(N >= 0);
	c_list = (Contaminant*)Realloc(c_list, (N+1)*sizeof(Contaminant));
	if (!c_list) abort();
	c_list[N].mass = mass;
	if (!name || (strlen(name)<1)) abort();
	c_list[N].name = Strdup(name);
	if (!c_list[N].name) abort();
	N++;
	return 1;
}
