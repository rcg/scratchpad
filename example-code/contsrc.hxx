#ifndef _CONTSRC_HXX_INCLUDED_
#define _CONTSRC_HXX_INCLUDED_

#include "prmagent.hxx"
#include "r3.hxx"
#include "cont.hxx"

class ContaminantSource : virtual public PrmAgent 
{
public:
	ContaminantSource();
	virtual ~ContaminantSource();
	int Isa(int);
	int Init(char*);
	virtual void Reset();
	virtual void *Get(int attribute,
		void *args, int args_size,
		void *data, int *size);

protected:
	ContaminantList *GetContaminantList() { return contaminants; };

public:
	static int IsSource(KID2(xid), char *contaminant);
	static int GetCSNum(KID2(xid), char *contaminant);
	static double GetCSValue(KID2(xid), double, R3, int);
Attribute:
	virtual double getCSValue(double t, R3 location, int cid)=0;
	virtual int getCSNum(char* contaminant);
private:
	char *taxname;
	ContaminantList *contaminants;
};

// Number of source contaminants
#define ATTR_CONTSRC_CSNUM		(CLASS_CONTSRC|0x0001) // int
// Contaminant source name to id
#define ATTR_CONTSRC_CSID		(CLASS_CONTSRC|0x0002) // int
// Contaminant value at location
#define ATTR_CONTSRC_CSVALUE	(CLASS_CONTSRC|0x0003) // double

#endif


