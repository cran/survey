
With the exception of svymle, porting this code to S-PLUS should not  be 
too difficult. I do not have the time or interest to maintain an S-PLUS 
port, but would be happy to answer questions.

The main differences will be
   --  eval.parent() needs to be replaced with eval(   ,sys.parent())
   --  The names of some of the components of a glm object may need to
       be changed in svyglm and its methods
   --  svytable() and postStratify() use the xtabs() function. The R xtabs() 
	function could be added to the package (which would change the licence 
	from LGPL to  GPL), or a substitute found.
   -- svyglm and svycoxph use with(), which needs to be changed to eval. 
      Some care may be needed to get the right scoping behaviour
   -- hadamard() uses lexical scope to hide its data objects. These would 
	need to be put in the global search path.
   -- svrepquantile() and svyquantile() use lexical scope to avoid 
	passing data to nested functions.  The data would have to be 
	passed explicitly.
   -- svychisq() and summary.svytable() uses bquote(), which can be 
	replaced by substitute()

svymle() will take more effort but should be feasible. It makes fairly
heavy use of lexical scope in defining the objective function and
gradient.  Lots of explicit argument passing will have to be added. 

