#ifndef _CONTSINK_HXX_INCLUDED_
#define _CONTSINK_HXX_INCLUDED_

#include "prmagent.hxx"
#include "searchagent.hxx"
#include "cont.hxx"

class ContaminantSink : virtual public PrmAgent, virtual public SearchAgent
{
public:
	ContaminantSink();
	virtual ~ContaminantSink();
	int Init(char *taxon);
	virtual int Isa(int);
	virtual int Compatibility();
	virtual void *Get(int attribute, void *args, int args_size,
		void *data, int *size);
	virtual void Reset();
	virtual void *GetState(int*);
	virtual void SetState(void*, int);

	virtual double Intoxicate(double t, double dt);
	virtual int CommitIntoxicate(double t, double dt, double dt2)=0;
	
	virtual int SetProfile(ContaminantProfile*);
	static ContaminantProfile *GetProfile(KID2(xid));
protected:
	virtual double LocalIntoxicate(int agent, double t, double dt, char *contaminant)=0;

	ContaminantList *contaminants;
	ContaminantProfile *profile;

private:
	int *get_contaminant_agent_list(int *num);
	char *taxname;
Attribute:
	virtual ContaminantProfile *getProfile();
};

#define ATTR_CONTSINK_PROFILE		(CLASS_CONTSINK|0x0001)

#endif
