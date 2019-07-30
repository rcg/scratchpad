// -*- outline-regexp: "/\\*-+";  -*-
/*-  Identification and Changes  */

/*
  cube.hxx -- Written by Randall Gray 
  Initial coding: 
  Date: 2007.09.06
  Location: odin.valhalla.asgard:/usr/home/gray/Projects/Ningaloo/cube.hxx

  History:

  $Log: cube.hxx,v $
  Revision 1.7  2008/08/29 04:11:21  gray
  More buggering with indenting!

  Revision 1.6  2008/08/26 08:21:39  gray
  Massive mods.  Lotsa space fixing and such. Merged Beths stuff in

  Revision 1.5  2007/10/24 00:21:03  gray
  Minor emacs tweaks -- no code changes yet

  Revision 1.4  2007/10/12 02:58:12  gray
  muchous state save/restore stuff

  Revision 1.3  2007/10/09 01:29:04  gray
  Reworked init_contaminants to do it incrementally.  Will need to test with a taxon change or some such thing

  Revision 1.2  2007/10/08 03:12:34  gray
  Why dont you come with me while we poison the fishes in the park, la da de da da de da

  Revision 1.1  2007/09/21 00:44:12  gray
  poisoning support


*/

/*-  Copyright  */

/*
  (C) 2007 CSIRO Australia
  All rights reserved
*/

/*-  Discussion  */

/*-  Configuration stuff  */

#ifndef __cube_hxx
#define in_cube_hxx
#define __cube_hxx

/*-  Types, defines, includes, externs and code  */

class Cube {
private:
	int n;
	double value, *v;

	double LValue();

public:
	Cube(int N);
	Cube(int N, double val);
	virtual ~Cube();
	
	double Value();
	double level(int i);
	void setMembers(double d);
	int add_dimension();
	double AdjustN(double K, int i); // Adjusts by a number removed
	double AdjustLevels(double *level, int base, int n); // Adjusts the levels by the proportion of the remaining range	
	double proportion_of_box(double *v, int dim);
#if 0
	double AdjustLevel(double level, int i); // Adjust according to a level in one of the other indices
#endif
  
	virtual void *GetState(int *sz);
	virtual void SetState(void *v, int sz);

};




/*-  The End  */

#endif
