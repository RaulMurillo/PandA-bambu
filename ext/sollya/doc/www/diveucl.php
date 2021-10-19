<a name="diveucl"></a> 
<div class="divName"> 
<h2 class="name">Name:</h2> <?php linkTo("command","diveucl","div");?> 
<span class="smallDescription">Computes the euclidian division of polynomials or numbers and returns the quotient 
</span> 
</div> 
<div class="divLibraryName"> 
<h2 class="libraryname">Library name:</h2> 
<span class="commandline type">sollya_obj_t sollya_lib_euclidian_div(sollya_obj_t, sollya_obj_t)</span> 
</div> 
<div class="divUsage"> 
<h2 class="category">Usage: </h2> 
<span class="commandline"><?php linkTo("command","diveucl","div");?>(<span class="arg">a</span>, <span class="arg">b</span>) : (<span class="type">function</span>, <span class="type">function</span>) -&gt; <span class="type">function</span></span> 
 
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
<li>When both <span class="arg">a</span> and <span class="arg">b</span> are constants, <?php linkTo("command","diveucl","div");?>(<span class="arg">a</span>,<span class="arg">b</span>) computes 
<?php linkTo("command","floor","floor");?>(<span class="arg">a</span> <?php linkTo("command","divide","/");?> <span class="arg">b</span>). In other words, it returns the quotient of the Euclidian 
division of <span class="arg">a</span> by <span class="arg">b</span>. 
</li><li>When both <span class="arg">a</span> and <span class="arg">b</span> are polynomials with at least one being non-constant, 
<?php linkTo("command","diveucl","div");?>(<span class="arg">a</span>,<span class="arg">b</span>) computes a polynomial <span class="arg">q</span> such that the polynomial <span class="arg">r</span> equal to 
<span class="arg">a</span> - <span class="arg">q</span> * <span class="arg">b</span> is of degree strictly smaller than the degree of <span class="arg">b</span> (see 
exception below). In order to recover <span class="arg">r</span>, use the <?php linkTo("command","modeucl","mod");?> command. 
</li><li><?php linkTo("command","diveucl","div");?> works on polynomials whose coefficients are constant 
expressions that cannot be simplified (by the tool) to rational 
numbers. In most cases, the tool is able to perform the Euclidian 
polynomial division for such polynomials and stop the Euclidian 
division algorithm only when <span class="arg">r</span> is of degree strictly smaller than 
the degree of <span class="arg">b</span>. In certain cases, when the polynomials involve 
coefficients given as constant expressions that are mathematically 
zero but for which the tool is unable to detect this fact, the tool 
may be unable to correctly determine that <span class="arg">r</span> is actually of degree 
stricly smaller than the degree of <span class="arg">b</span>. The issue arises in particular 
for polynomials whose leading coefficient is a constant expression 
which is zero without the tool being able to detect this. In these 
cases, <?php linkTo("command","diveucl","div");?>, together with <?php linkTo("command","modeucl","mod");?>, just guarantee that <span class="arg">q</span> and 
<span class="arg">r</span>, as returned by the two commands, satisfy that <span class="arg">r</span> added to the 
product of <span class="arg">q</span> and <span class="arg">b</span> yields <span class="arg">a</span>, and that <span class="arg">r</span> is of the smallest 
degree the tool can admit. However, there might exist another pair of 
a quotient and remainder polynomial for which the remainder polynomial 
is of a degree less than the one of <span class="arg">r</span>. 
</li><li>When at least one of <span class="arg">a</span> or <span class="arg">b</span> is a function that is no polynomial, 
<?php linkTo("command","diveucl","div");?>(<span class="arg">a</span>,<span class="arg">b</span>) returns 0. 
</ul> 
</div> 
<div class="divExamples"> 
<div class="divExample"> 
<h2 class="category">Example 1: </h2> 
&nbsp;&nbsp;&nbsp;&gt; div(1001, 231);<br> 
&nbsp;&nbsp;&nbsp;4<br> 
&nbsp;&nbsp;&nbsp;&gt; div(13, 17);<br> 
&nbsp;&nbsp;&nbsp;0<br> 
&nbsp;&nbsp;&nbsp;&gt; div(-14, 15);<br> 
&nbsp;&nbsp;&nbsp;-1<br> 
&nbsp;&nbsp;&nbsp;&gt; div(-213, -5);<br> 
&nbsp;&nbsp;&nbsp;42<br> 
&nbsp;&nbsp;&nbsp;&gt; div(23/13, 11/17);<br> 
&nbsp;&nbsp;&nbsp;2<br> 
&nbsp;&nbsp;&nbsp;&gt; div(exp(13),-sin(17));<br> 
&nbsp;&nbsp;&nbsp;460177<br> 
</div> 
<div class="divExample"> 
<h2 class="category">Example 2: </h2> 
&nbsp;&nbsp;&nbsp;&gt; div(24 + 68 * x + 74 * x^2 + 39 * x^3 + 10 * x^4 + x^5, 4 + 4 * x + x^2);<br> 
&nbsp;&nbsp;&nbsp;6 + x * (11 + x * (6 + x))<br> 
&nbsp;&nbsp;&nbsp;&gt; div(24 + 68 * x + 74 * x^2 + 39 * x^3 + 10 * x^4 + x^5, 2 * x^3);<br> 
&nbsp;&nbsp;&nbsp;19.5 + x * (5 + x * 0.5)<br> 
&nbsp;&nbsp;&nbsp;&gt; div(x^2, x^3);<br> 
&nbsp;&nbsp;&nbsp;0<br> 
</div> 
<div class="divExample"> 
<h2 class="category">Example 3: </h2> 
&nbsp;&nbsp;&nbsp;&gt; div(exp(x), x^2);<br> 
&nbsp;&nbsp;&nbsp;0<br> 
</div> 
</div> 
<div class="divSeeAlso"> 
<span class="category">See also: </span><?php linkTo("command","gcd","gcd");?>, <?php linkTo("command","modeucl","mod");?>, <?php linkTo("command","numberroots","numberroots");?> 
</div> 
