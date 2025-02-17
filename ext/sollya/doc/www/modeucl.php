<a name="modeucl"></a> 
<div class="divName"> 
<h2 class="name">Name:</h2> <?php linkTo("command","modeucl","mod");?> 
<span class="smallDescription">Computes the euclidian division of polynomials or numbers and returns the rest 
</span> 
</div> 
<div class="divLibraryName"> 
<h2 class="libraryname">Library name:</h2> 
<span class="commandline type">sollya_obj_t sollya_lib_euclidian_mod(sollya_obj_t, sollya_obj_t)</span> 
</div> 
<div class="divUsage"> 
<h2 class="category">Usage: </h2> 
<span class="commandline"><?php linkTo("command","modeucl","mod");?>(<span class="arg">a</span>, <span class="arg">b</span>) : (<span class="type">function</span>, <span class="type">function</span>) -&gt; <span class="type">function</span></span> 
 
</div> 
<div class="divParameters"> 
<h2 class="category">Parameters: </h2> 
<ul> 
<li><span class="arg">a</span> is a constant or a polynomial.</li> 
<li><span class="arg">b</span> is a constant or a polynomial.</li> 
</ul> 
</div> 
<div class="divDescription"> 
<h2 class="category">Description: </h2><ul> 
<li><?php linkTo("command","modeucl","mod");?>(<span class="arg">a</span>,<span class="arg">b</span>) computes <span class="arg">a</span> <?php linkTo("command","minus","-");?> (<span class="arg">b</span> <?php linkTo("command","mult","*");?> <?php linkTo("command","diveucl","div");?>(<span class="arg">a</span>,<span class="arg">b</span>)). 
In other words, it returns the remainder of the Euclidian division 
of <span class="arg">a</span> by <span class="arg">b</span>. 
</li><li>See <?php linkTo("command","diveucl","div");?> for subtle cases involving polynomials whose degree can 
not easily be computed by the tool as their leading coefficient is 
given as a constant expression that is mathematically zero but for 
which the tool is unable to detect this fact. 
</ul> 
</div> 
<div class="divExamples"> 
<div class="divExample"> 
<h2 class="category">Example 1: </h2> 
&nbsp;&nbsp;&nbsp;&gt; mod(1001, 231);<br> 
&nbsp;&nbsp;&nbsp;77<br> 
&nbsp;&nbsp;&nbsp;&gt; mod(13, 17);<br> 
&nbsp;&nbsp;&nbsp;13<br> 
&nbsp;&nbsp;&nbsp;&gt; mod(-14, 15);<br> 
&nbsp;&nbsp;&nbsp;1<br> 
&nbsp;&nbsp;&nbsp;&gt; mod(-213, -5);<br> 
&nbsp;&nbsp;&nbsp;-3<br> 
&nbsp;&nbsp;&nbsp;&gt; print(mod(23/13, 11/17));<br> 
&nbsp;&nbsp;&nbsp;105 / 221<br> 
&nbsp;&nbsp;&nbsp;&gt; print(mod(exp(13),-sin(17)));<br> 
&nbsp;&nbsp;&nbsp;exp(13) + 460177 * sin(17)<br> 
</div> 
<div class="divExample"> 
<h2 class="category">Example 2: </h2> 
&nbsp;&nbsp;&nbsp;&gt; mod(24 + 68 * x + 74 * x^2 + 39 * x^3 + 10 * x^4 + x^5, 4 + 4 * x + x^2);<br> 
&nbsp;&nbsp;&nbsp;0<br> 
&nbsp;&nbsp;&nbsp;&gt; mod(24 + 68 * x + 74 * x^2 + 39 * x^3 + 10 * x^4 + x^5, 2 * x^3);<br> 
&nbsp;&nbsp;&nbsp;24 + x * (68 + x * 74)<br> 
&nbsp;&nbsp;&nbsp;&gt; mod(x^2, x^3);<br> 
&nbsp;&nbsp;&nbsp;x^2<br> 
</div> 
<div class="divExample"> 
<h2 class="category">Example 3: </h2> 
&nbsp;&nbsp;&nbsp;&gt; mod(exp(x), x^2);<br> 
&nbsp;&nbsp;&nbsp;exp(x)<br> 
&nbsp;&nbsp;&nbsp;&gt; mod(x^3, sin(x));<br> 
&nbsp;&nbsp;&nbsp;x^3<br> 
</div> 
</div> 
<div class="divSeeAlso"> 
<span class="category">See also: </span><?php linkTo("command","gcd","gcd");?>, <?php linkTo("command","diveucl","div");?>, <?php linkTo("command","numberroots","numberroots");?> 
</div> 
