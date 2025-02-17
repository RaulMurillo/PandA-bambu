<a name="library"></a> 
<div class="divName"> 
<h2 class="name">Name:</h2> <?php linkTo("command","library","library");?> 
<span class="smallDescription">binds an external mathematical function to a variable in Sollya 
</span> 
</div> 
<div class="divLibraryName"> 
<h2 class="libraryname">Library names:</h2> 
<span class="commandline type">sollya_obj_t sollya_lib_libraryfunction(sollya_obj_t, char *,</span> 
<span class="commandline type">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;int (*)(mpfi_t, mpfi_t, int))</span> 
<span class="commandline type">sollya_obj_t sollya_lib_build_function_libraryfunction(sollya_obj_t, char *,</span> 
<span class="commandline type">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;int (*)(mpfi_t,</span> 
<span class="commandline type">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;mpfi_t, int))</span> 
<span class="commandline type">sollya_obj_t sollya_lib_libraryfunction_with_data(</span> 
<span class="commandline type">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;sollya_obj_t, char *,</span> 
<span class="commandline type">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;int (*)(mpfi_t, mpfi_t, int, void *),</span> 
<span class="commandline type">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;void *, void (*)(void *))</span> 
<span class="commandline type">sollya_obj_t sollya_lib_build_function_libraryfunction_with_data(</span> 
<span class="commandline type">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;sollya_obj_t, char *,</span> 
<span class="commandline type">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;int (*)(mpfi_t,</span> 
<span class="commandline type">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;mpfi_t, int, void *),</span> 
<span class="commandline type">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;void *, void (*)(void *))</span> 
</div> 
<div class="divUsage"> 
<h2 class="category">Usage: </h2> 
<span class="commandline"><?php linkTo("command","library","library");?>(<span class="arg">path</span>) : <span class="type">string</span> -&gt; <span class="type">function</span></span> 
 
</div> 
<div class="divDescription"> 
<h2 class="category">Description: </h2><ul> 
<li>The command <?php linkTo("command","library","library");?> lets you extend the set of mathematical 
functions known to Sollya. 
By default, Sollya knows the most common mathematical functions such 
as <?php linkTo("command","exp","exp");?>, <?php linkTo("command","sin","sin");?>, <?php linkTo("command","erf","erf");?>, etc. Within Sollya, these functions may be 
composed. This way, Sollya should satisfy the needs of a lot of 
users. However, for particular applications, one may want to 
manipulate other functions such as Bessel functions, or functions 
defined by an integral or even a particular solution of an ODE. 
</li><li><?php linkTo("command","library","library");?> makes it possible to let Sollya know about new functions. In 
order to let it know, you have to provide an implementation of the 
function you are interested in. This implementation is a C file containing 
a function of the form: 
              int my_ident(sollya_mpfi_t result, sollya_mpfi_t op, int n)  
The semantic of this function is the following: it is an implementation of 
the function and its derivatives in interval arithmetic. 
my_ident(result, I, n) shall store in result an enclosure  
of the image set of the n-th derivative 
of the function f over I: f^(n)(I) C result. 
</li><li>The integer value returned by the function implementation currently has no 
meaning. 
</li><li>You do not need to provide a working implementation for any n. Most 
functions of Sollya requires a relevant implementation only for f, 
f' and f''. For higher derivatives, its is not so critical and the 
implementation may just store [-Inf, +Inf] in result whenever n&gt;2. 
</li><li>Note that you should respect somehow interval-arithmetic standards in your 
implementation: result has its own precision and you should perform the  
intermediate computations so that result is as tight as possible. 
</li><li>You can include sollya.h in your implementation and use library  
functionnalities of Sollya for your implementation. However, this requires to 
have compiled Sollya with -fPIC in order to make the Sollya executable code 
position independent and to use a system on with programs, using dlopen to 
open dynamic routines can dynamically open themselves. IMPORTANT NOTICE: as 
the code will be run in a context where a sollya session is already opened, 
the library functions must be used directly, without calling sollya_lib_init 
and sollya_lib_close (calling these functions would conflict with the 
current session, leading to weird and hard to debug behaviors). 
</li><li>To bind your function into Sollya, you must use the same identifier as the 
function name used in your implementation file (my_ident in the previous 
example). Once the function code has been bound to an identifier, you can use 
a simple assignment to assign the bound identifier to yet another identifier. 
This way, you may use convenient names inside Sollya even if your 
implementation environment requires you to use a less convenient name. 
</li><li>The dynamic object file whose name is given to <?php linkTo("command","library","library");?> for binding of an 
external library function may also define a destructor function 
  int sollya_external_lib_close(void). 
If Sollya finds such a destructor function in the dynamic object file, it 
will call that function when closing the dynamic object file again. 
This happens when Sollya is terminated or when the current Sollya session 
is restarted using <?php linkTo("command","restart","restart");?>. 
The purpose of the destructor function is to allow the dynamically bound code 
to free any memory that it might have allocated before Sollya is terminated  
or restarted.  
The dynamic object file is not necessarily needed to define a destructor 
function. This ensure backward compatibility with older Sollya external  
library function object files. 
When defined, the destructor function is supposed to return an integer 
value indicating if an error has happened. Upon success, the destructor 
functions is to return a zero value, upon error a non-zero value. 
</ul> 
</div> 
<div class="divExamples"> 
<div class="divExample"> 
<h2 class="category">Example 1: </h2> 
&nbsp;&nbsp;&nbsp;&gt; bashexecute("gcc -fPIC -Wall -c libraryexample.c -I$HOME/.local/include");<br> 
&nbsp;&nbsp;&nbsp;&gt; bashexecute("gcc -shared -o libraryexample libraryexample.o -lgmp -lmpfr");<br> 
&nbsp;&nbsp;&nbsp;&gt; myownlog = library("./libraryexample");<br> 
&nbsp;&nbsp;&nbsp;&gt; evaluate(log(x), 2);<br> 
&nbsp;&nbsp;&nbsp;0.69314718055994530941723212145817656807550013436025<br> 
&nbsp;&nbsp;&nbsp;&gt; evaluate(myownlog(x), 2);<br> 
&nbsp;&nbsp;&nbsp;0.69314718055994530941723212145817656807550013436025<br> 
</div> 
</div> 
<div class="divSeeAlso"> 
<span class="category">See also: </span><?php linkTo("command","function","function");?>, <?php linkTo("command","bashexecute","bashexecute");?>, <?php linkTo("command","externalproc","externalproc");?>, <?php linkTo("command","externalplot","externalplot");?>, <?php linkTo("command","diff","diff");?>, <?php linkTo("command","evaluate","evaluate");?>, <?php linkTo("command","libraryconstant","libraryconstant");?> 
</div> 
