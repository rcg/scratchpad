#ifndef _CONT_HXX_INCLUDED_
#define _CONT_HXX_INCLUDED_

#include "stringtable.hxx"

// Probably need to add a heap of stuff to this later
// so the agent can store the information on
// how to deal with each contaminant

class ContaminantList
{
public:
	ContaminantList();
	~ContaminantList();
	int RegisterAsContaminantSource( char *contaminant ); // This makes the agent ... well, a source
	int RegisterInterest( char *contaminant );            // This is used to indicate which contaminants are pertinent
	int IsSource( char *contaminant );                    // boolean test for the given contaminant
	int IsInterested( char *contaminant );                // boolean test to see if an agent is interested in a contaminant
	int NumSource();                                      // returns the number of sources in the system
	int NumInterest();                                    // returns the number of entities in the system that are interested in contaminants
	char *GetInterest( int rec );
	char *GetSource( int rec );  
	void *GetState( int *sz );
	void SetState( void *d, int sz );
	void ClearInterests();
	void ClearSources();
private:
	StringTable *source, *interest;
};

class ContaminantProfile
{
public:
	ContaminantProfile();
	ContaminantProfile( ContaminantProfile* );
	ContaminantProfile( void*, int );
	~ContaminantProfile();
	void *GetState( int *sz );

	typedef struct _Contaminant {
		char *name;
		double mass;
	} Contaminant;
	int N;
	Contaminant *c_list;

	int AddContaminant( char *name, double mass );

private:
	void *pack_list( int* );
	void *pack_struct( int, int* );
	void unpack_list( void*, int );
	void unpack_struct( int, void*, int );
};

#endif
