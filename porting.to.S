Porting this code to S-PLUS should be tedious but not too difficult, at
least if a GPL license is sufficient so that functions missing in S-PLUS
can be copied from R.  I do not have the time or interest to maintain an
S-PLUS port, but would be happy to answer questions. 

It would be sufficient to port the functions for survey.design2 and
svyrep.design objects, as survey.design objects are deprecated.


The main differences will be
   -- class() needs to be replaced by oldClass()
   -- eval.parent() needs to be replaced by eval(   ,sys.parent())
   -- The names of some of the components of a glm object may need to
       be changed in svyglm and its methods
   -- In the recursive variance estimation you need special handling of 
	empty strata because NULL[,1] is an error in S-PLUS
   -- svytable() and postStratify() use the xtabs() function. The R xtabs() 
	function could be added to the package (which would change the licence 
	from LGPL to  GPL), or a substitute found.
   -- svyglm and svycoxph use with(), which needs to be changed to eval. 
      Some care may be needed to get the right scoping behaviour
   -- hadamard() uses lexical scope to hide its data objects. These would 
	need to be put in the global search path.
   -- bquote() needs to be ported from R or replaced by substitute()
   -- "split<-" and "colnames<-" need to be ported from R or rewritten
   -- The syntax of UseMethod() when dispatch is not on the first 
	argument appears to be different in S-PLUS
   -- explicit argument passing needs to be added to many *apply calls, to 
      quantile functions in svyquantile and svrepquantile, and to the 
      objective and gradient functions in svymle()

