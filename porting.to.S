
With the exception of svymle, porting this code to S-PLUS should not
be difficult if anyone is interested in doing so. The main differences
will be
   --  eval.parent() needs to be replaced with eval(   ,sys.parent())
   --  The names of some of the components of a glm object may need to
       be changed in svyglm and its methods
   --  svytable() uses the xtabs() function. The R xtabs() function could 
       be added to the package (which would change the licence from LGPL to
       GPL), or a substitute found.
   --  svyglm and svycoxph use with(), which needs to be changed to eval.

svymle() will take more effort but should be feasible
   --  optim(), used in svymle(), needs to be replaced by nlmin or nlminb, 
       which have different calling formats.
   --  svymle() also makes fairly heavy use of lexical scope in defining
       the objective function and gradient.  Lots of explicit argument
       passing will have to be added.
